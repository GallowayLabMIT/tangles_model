import scipy.integrate
import numpy as np
from numba import njit, vectorize, float32, float64, boolean
from numba.experimental import jitclass
import multiprocessing
import pandas as pd
import tempfile
import os
import enum
import gc
import subprocess
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

DEFAULT_PARAMS = {
    'mRNA_drag': 1/20, # pN nm^(alpha / 1)
    'mRNA_exponent': 1, # the value of alpha
    'DNA_twist_mobility': 10, # s pN nm
    'RNAP_radius': 15, # nm
    'RNAP_velocity': 20, # nm / s
    'RNAP_torque_cutoff': 12, # pN nm
    'RNAP_stall_torque_width': 3, #pN
    'DNA_force': 1, # pN
    'DNA_bend_plength': 50, # pN
    'DNA_twist_plength': 95, # pN
    'DNA_plectonome_twist_plength': 24, # pN
    'temperature': 298, # K
    'base_promoter_initiation_rate': 1 / 120, # 1 / sec
    'topo_rate': 1 / 1200, # 1 / sec
    'mRNA_deg_rate': 1 / 1200 # 1 / sec
}

# Most member variables are float64's except for left_bc_free/right_bc_free
supercoiling_varspec = [(name,float64) for name in
                       ['kb', 'w0', 'tau_s', 'tau_0', 'tau_p', 'sc_s', 'sc_p', 'rnac_r',
                        'v_0', 't_c', 'stall_torque_width', 'chi', 'eta', 'alpha',
                       'left_bc_position', 'right_bc_position', 'left_bc_val', 'right_bc_val']]
supercoiling_varspec.append(('left_bc_free', boolean))
supercoiling_varspec.append(('right_bc_free', boolean))

@jitclass(supercoiling_varspec)
class SupercoilingPolymeraseModel(object):
    """
    Class that encapsulates all physical constants and
    calculations needed to compute an ODE model for
    supercoiling induced by transcription.
    
    It exposes the system_derivatives member function,
    which can be used in an ODE solver to represent the motion
    of each of the polymerases and the induced supercoiling.
    """
    def __init__(self, A, C, P, T, f,
                 v_0, critical_torque, stall_torque_width, rnac_r,
                 chi, eta, alpha, left_bc_free, right_bc_free,
                 left_bc_loc, right_bc_loc, left_bc_val, right_bc_val):
        """
        Constructor that takes all relevant system constants and initializes
        internal parameters used by other functions.
    
        Args:
        -----
        A: DNA bend persistence length (nm)
        C: DNA twist persistance length (nm)
        P: DNA plectonome twist persistence length (nm)
        T: Temperature of the system (K)
        f: Constant force held on DNA (pN)
        v_0: The maximum speed of the polymerase (nm / s)
        critical_torque: Torque (in pN nm) at which the
            polymerases begin to stall.
        stall_torque_width: How wide the stall torque distribution is (pN)
        rnac_r: The radius of a RNAP (used for hard-sphere repulsion)
        chi: The DNA twisting mobility (s pN nm)
        eta: The The mRNA drag coefficient (pN (nm^(alpha - 1)))
        alpha: The mRNA power-law scaling exponent.
        left_bc_free: A boolean representing if the left BC is free.
        right_bc_free: A boolean representing if the right BC is free.
        left_bc_loc: If the left BC is not free, encodes the location of the BC.
        right_bc_loc: If the right BC is not free, encodes the location of the BC.
        left_bc_val: If the left BC is not free, encodes the value of phi at the BC.
        right_bc_val: If the right BC is not free, encodes the value of phi at the BC.
        """
        self.kb =  1380649 / 100000000 # pN nm / K
        self.w0 =  1.85 # 1/nm

        # Init torque model constants
        c = self.kb * T * C * self.w0 ** 2 # pN
        p = self.kb * T * P * self.w0 ** 2 # pN
        g = f - np.sqrt(self.kb * T * f / A) # pN
        cs = c * (1 - C / (4 * A) * np.sqrt(self.kb * T / (A * f))) # pN

        # Compute critical torque model values
        self.tau_s = cs / self.w0 # pN nm
        self.tau_0 = np.sqrt(2 * p * g / (self.w0 ** 2 * (1 - p / cs))) # pN nm
        self.tau_p = p / self.w0 # pN nm
        self.sc_s = 1 / cs * np.sqrt(2 * p * g  / (1 - p / cs)) # dimless
        self.sc_p = 1 / p * np.sqrt(2 * p * g / (1 - p / cs)) # dimless
        
        self.rnac_r = rnac_r
        self.v_0 = v_0
        self.t_c = critical_torque
        self.stall_torque_width = stall_torque_width
        self.chi = chi
        self.eta = eta
        self.alpha = alpha
        self.left_bc_free = left_bc_free
        self.right_bc_free = right_bc_free
        self.left_bc_position = left_bc_loc
        self.right_bc_position = right_bc_loc
        self.left_bc_val = left_bc_val
        self.right_bc_val = right_bc_val
    
    def evaluate_twist(self, x, states):
        """
        Given the current polymerase states,
        computes the excess DNA twist phi at each
        of the points in x.
        
        Args:
        -----
        x: An ndarray representing the points to evaluate phi at.
        states: A (3N,) shape ndarray encoding (location, mRNA, phi) for each polymerase.
            This vector is assumed to be sorted in increasing order of location.
        
        Returns:
        --------
        An ndarray matching the shape of x, where twist has been evaluated.
        """
        augmented_loc = np.concatenate((np.array([self.left_bc_position]),
                                        states[::3],
                                        np.array([self.right_bc_position])))
        augmented_phi = np.concatenate((np.array([self.left_bc_val]),
                                        states[2::3],
                                        np.array([self.right_bc_val])))
        # Replace these with free BCs if needed
        if self.left_bc_free:
            # Step w0 away from BC, to avoid divide by zero errors.
            augmented_loc[0] = augmented_loc[1] - 1
            augmented_phi[0] = augmented_phi[1]
        if self.right_bc_free:
            augmented_loc[-1] = augmented_loc[-2] + 1
            augmented_phi[-1] = augmented_phi[-2]
        
        return np.interp(x, augmented_loc, augmented_phi)
    
    def torque_response(self, sc):
        """
        Given a ndarray of supercoiling values, calculates
        the torque response on each using vectorized Numpy functions.
        
        The torque response is a function of the provided
        DNA persistence lengths, the temperature and force,
    
        Args:
        -----
        sc: A ndarray containing supercoiling density values.
        
        Returns:
        --------
        The calculated torque, using a two-phase constant-force
        torque equation.
        """
        abs_sc = np.abs(sc)
        result = self.tau_0 * np.sign(sc)
        result[abs_sc < self.sc_s] = self.tau_s * sc[abs_sc < self.sc_s]
        result[abs_sc >= self.sc_p] = self.tau_p * sc[abs_sc >= self.sc_p]
        return result
    
    def polymerase_velocity(self, sc_behind, sc_ahead):
        """
        Given supercoiling densities behind and ahead of
        each polymerase, returns the velocity of each
        polymerase.
        
        Args:
        -----
        sc_behind: An ndarray containing the supercoiling
            density behind each polymerase.
        sc_ahead: An ndarray containing the supercoiling
            density ahead of each polymerase.
        
        Returns:
        --------
        The velocity of each polymerase (in nm/s).
        """
        # Restrict exponential argument to 10, to ensure we don't overflow
        return self.v_0 / (
            (1 + np.exp(np.minimum(10.0,
                           (np.abs(self.torque_response(sc_behind)) - self.t_c) / self.stall_torque_width))) *
            (1 + np.exp(np.minimum(10.0,
                           (np.abs(self.torque_response(sc_ahead)) - self.t_c) / self.stall_torque_width))))
    
    def system_derivatives(self, states, polymerase_directions):
        """
        Given the state of the system, in terms of a location, nascant mRNA length,
        and local excess DNA rotation, computes the derivatives of all of these states.
    
        Args:
        -----
        states: A (3N,) shape ndarray encoding (location, mRNA, phi) for each polymerase.
            This vector is assumed to be sorted in increasing order of location
        polymerase_directions: A (N,) shape ndarray encoding +-1, depending on the
            direction that the polymerase is moving.

        Returns:
        --------
        A 3N-length vector encoding the time derivatives of each of the states.
        """
        augmented_loc = np.concatenate((np.array([self.left_bc_position]),
                                        states[::3],
                                        np.array([self.right_bc_position])))
        augmented_phi = np.concatenate((np.array([self.left_bc_val]),
                                        states[2::3],
                                        np.array([self.right_bc_val])))
        # Replace these with free BCs if needed
        if self.left_bc_free:
            # Step w0 away from BC, to avoid divide by zero errors.
            augmented_loc[0] = augmented_loc[1] - 1
            augmented_phi[0] = augmented_phi[1]
        if self.right_bc_free:
            augmented_loc[-1] = augmented_loc[-2] + 1
            augmented_phi[-1] = augmented_phi[-2]
        # Compute supercoiling in each N + 1 region
        supercoiling = np.diff(augmented_phi) / (np.diff(augmented_loc) * -self.w0)
        rnac_torque = self.torque_response(supercoiling[1:]) - self.torque_response(supercoiling[:-1])
        velocities = polymerase_directions * self.polymerase_velocity(supercoiling[1:], supercoiling[:-1])
        # Calculate inter-polymerase distances (augmented with large distances)
        inter_distances = np.concatenate((np.array([np.inf]),
                                          np.diff(augmented_loc[1:-1]),
                                          np.array([np.inf])))
        # Halt if the polymerase is within the 2 RNAC radii in the direction it is going
        velocities[(polymerase_directions > 0) & (inter_distances[1:] < 2 * self.rnac_r)] = 0
        velocities[(polymerase_directions < 0) & (inter_distances[:-1] < 2 * self.rnac_r)] = 0
        
        drag = self.eta * (states[1::3] ** self.alpha)
        angular_changes =  drag * velocities * self.w0 / (self.chi + drag) - rnac_torque / (self.chi + drag)

        derivatives = np.zeros(states.shape)
        derivatives[::3] = velocities
        derivatives[1::3] = polymerase_directions * velocities # Ensure that transcript size is strictly positive
        derivatives[2::3] = angular_changes
        return derivatives

class BoundaryCondition(object):
    """
    Simple tuple-like class that stores boundary conditions,
    with simple helper static methods for creating BCs.
    """

    def __init__(self, location, twist_value, is_free):
        """
        Sets internal state for this BoundaryCondition

        Args:
        -----
        location: Location in nm of this boundary condition.
        twist_values: Excess DNA twist, in radians, at the given location.
        is_free: If the BC is "free", e.g. has unrestricted twist.
        """
        self.loc = location
        self.excess_twist = twist_value
        self.is_free = is_free
    
    @staticmethod
    def free(location):
        """
        Returns a free BC at the given location.

        Args:
        -----
        location: Location in nm of this free BC.
        """
        return BoundaryCondition(location, 0, True)
    
    @staticmethod
    def fixed(location, excess_twist=0):
        """
        Returns a fixed-twist BC at the given location.

        Args:
        -----
        location: Location in nm of this fixed-twist BC.
        excess_twist: The fixed excess twist at this location. Default is zero excess twist.
        """
        return BoundaryCondition(location, excess_twist, False)

class PromoterType(enum.Enum):
    STATIC = enum.auto()
    SC_DEPENDENT = enum.auto()

class PAS_Type(enum.Enum):
    NO_RT = enum.auto()
    LENGTH_RT = enum.auto()
    TIME_RT = enum.auto()

class Promoter(object):
    """
    Stores the information required to represent a promoter
    """
    def __init__(self, type, rate, location, direction):
        """
        Args:
        -----
        type: a PromoterType specifying the specific type of promoter.
        rate: A base promoter on-rate, specified in units of 1/s
        location: The location in nm of the promoter.
        direction: Either +/-1, representing the direction polymerases start
        """
        self.type = type
        self.base_rate = rate
        self.location = location
        self.direction = direction

class PAS(object):
    """
    Stores information required to represent a polyA signal
    """
    def __init__(self, type, location, direction, readthrough_duration=0.0):
        """
        Args:
        -----
        type: a PAS_Type specifying the type of polyA-signal.
        location: A location in nm of the PAS.
        direction: Either +/-1, representing the direction in which we terminate.
        readthrough_duration: Either a distance in nm or time in seconds to readthrough, on average.
        """
        self.type = type
        self.location = location
        self.direction = direction
        self.readthrough_duration = readthrough_duration


class SimulationState(object):
    """
    Class that stores all relevant information about the system state needed for a simulation.

    This class includes the base SupercoilngPolymeraseModel, which stores the dynamic information
    necessary for a simulation, as well as the normal (3N,) vector that describes the current
    dynamic motion of the relevant LNCs.

    This class also includes dynamic information that is _not_ directly used in the ODE solver,
    such as the location of genes, direction of each polymerase, and termination conditions.

    Finally, this class contains mutation functions that describe common operations, such
    as adding a polymerase, removing a polymerase, and topoisomerase action.

    Member variables:
    -----------------
    bcs: A 2-tuple containing the bounding BoundaryCondition's
    params: A dictionary containing all relevant physical constants.
    physical_model: A SupercoilingPolymeraseModel used to evaluate all ODE equations.
    time: the current simulation time.
    state: a (3N,) ndarray representing triples of (location, mRNA_length, excess_twist)
        for each active polymerase.
    genes: A list of the form [(gene_name, (gene_start_in_nm, gene_end_in_nm))].
        Genes are sorted by location!
    polymerase_state: a N x 3 ndarray that stores [direction, start_loc, uid] for each active polymerase.
        This is used to identify which transcripts to create when termination occurs.
    transcripts: a G x 1 ndarray that stores the number of transcripts for each gene. This can be used
        to simulate higher-order interactions between genes.
    """
    def __init__(self, start_time, bcs, genes, params, physical_model):
        """
        Given the system parameters and boundary conditions,
        initializes a SimulationState.
        
        Args:
        -----
        start_time: The initial time to start the simulation at.
        bcs: A 2-tuple encoding (Left_BoundaryCondition, Right_BoundaryCondition). 
            Each BC should be an element of type BoundaryCondition.
        genes: A dictionary, where keys are gene names and values are (start, end) tuples
            encoding the start/end of genes.
        params: A dictionary encoding relevant keyvalue pairs representing physical constants.
            See the docstring for `create_physical_model`.
        physical_model: A SupercoilingPolymeraseModel representing the underlying physics.
            Can be created with SimulationState.create_physical_model
        
        Returns:
        --------
        A SimulationState with all relevant values set.
        """
        self.bcs = bcs
        self.params = params
        
        # Init model
        self.physical_model = physical_model
        
        # Init helper things to track.
        self.time = start_time
        self.state = np.zeros(0)
        self.genes = sorted(genes.items(), key=lambda x:min(x[1]))
        self.polymerase_state = np.zeros((0,3))
        self.transcripts = np.zeros((len(self.genes),))
        self.next_uid = 0
        pass

    @staticmethod
    def create_physical_model(bcs, params=DEFAULT_PARAMS):
        """
        Given a dictionary of parameters, creates a SupercoilingPolymeraseModel.
        Use this static function to create the physical model, so that we don't
        have to recompute all of the time.

        Args:
        -----
        bcs: A 2-tuple encoding (Left_BoundaryCondition, Right_BoundaryCondition). 
            Each BC should be an element of type BoundaryCondition.
        params: A dictionary containing the following keys:
            mRNA_drag:            pN nm^(alpha / 1)
            mRNA_exponent:        the value of alpha
            DNA_twist_mobility:   s pN nm
            RNAP_radius:           nm
            RNAP_velocity:        nm / s
            RNAP_torque_cutoff:   pN nm
            RNAP_stall_torque_width: pN
            DNA_force:            pN
            DNA_bend_plength:     pN
            DNA_twist_plength:    pN
            DNA_plectonome_twist_plength: pN
            temperature:          K
        """
        return SupercoilingPolymeraseModel(
            params['DNA_bend_plength'], params['DNA_twist_plength'],
            params['DNA_plectonome_twist_plength'], params['temperature'],
            params['DNA_force'], params['RNAP_velocity'],
            params['RNAP_torque_cutoff'], params['RNAP_stall_torque_width'],
            params['RNAP_radius'],
            params['DNA_twist_mobility'], params['mRNA_drag'],
            params['mRNA_exponent'], bcs[0].is_free, bcs[1].is_free,
            bcs[0].loc, bcs[1].loc, bcs[0].excess_twist, bcs[1].excess_twist)
        
    def calculate_torque(self, location):
        augmented_loc = np.concatenate((np.array([self.physical_model.left_bc_position]),
                                        self.state[::3],
                                        np.array([self.physical_model.right_bc_position])))
        augmented_phi = np.concatenate((np.array([self.physical_model.left_bc_val]),
                                        self.state[2::3],
                                        np.array([self.physical_model.right_bc_val])))
        # Replace these with free BCs if needed
        if self.physical_model.left_bc_free:
            # Step w0 away from BC, to avoid divide by zero errors.
            augmented_loc[0] = augmented_loc[1] - 1
            augmented_phi[0] = augmented_phi[1]
        if self.physical_model.right_bc_free:
            augmented_loc[-1] = augmented_loc[-2] + 1
            augmented_phi[-1] = augmented_phi[-2]
        # Compute supercoiling in each N + 1 region
        supercoiling = np.diff(augmented_phi) / (np.diff(augmented_loc) * -self.physical_model.w0)
        rnac_torque = self.physical_model.torque_response(supercoiling[1:]) - self.physical_model.torque_response(supercoiling[:-1])

        relative_loc = augmented_loc - location
        region_idx = np.where(augmented_loc - location < 0)[0][-1]
        return rnac_torque[region_idx] if len(rnac_torque) > 0 else 0

    def remove_supercoiling(self):
        """
        Randomly removes all supercoiling in an intergenic region.
        """
        if len(self.state) == 0:
            return
        # Genes are already sorted
        # We now rewrite all polymerases that are sitting on the range
        # (idx).start to (idx + 1).end
        # We do this by removing all relevant polymerases from state, then
        # using the model.compute_twist function to recalculate, at each
        # of the polymerases to rewrite
        if len(self.genes) == 1:
            lower_bound = -np.inf
            upper_bound = np.inf
        else:
            relax_after_idx = np.random.choice(range(len(self.genes) - 1))
            lower_bound = min(self.genes[relax_after_idx][1])
            upper_bound = max(self.genes[relax_after_idx + 1][1])
        polymerase_rewrite_idx = np.where(
            (self.state[::3] > lower_bound) & (self.state[::3] < upper_bound))
        polymerase_positions = self.state[::3][polymerase_rewrite_idx]
        if len(polymerase_positions) == 0:
            # No polymerases in this range: just return
            return
        # Exploit polymerases always being sorted
        effective_state = np.delete(self.state, range(np.min(polymerase_rewrite_idx) * 3,
                                            (np.max(polymerase_rewrite_idx) + 1) * 3))
        self.state[2::3][polymerase_rewrite_idx] = \
            self.physical_model.evaluate_twist(polymerase_positions, effective_state)
    
    def add_polymerase(self, location, direction):
        """
        Adds a new polymerase, given a location and direction.

        If the promoter site is already occupied, nothing happens.

        Args:
        -----
        location: double, the location in nm to start the polymerase.
        direction: +/-1, an integer representing which direction the polymerase should move.
        """
        # Skip if promoter site already occupied
        if np.min(
            np.abs(np.concatenate((np.array([np.inf]), self.state[::3])) - location)) < self.physical_model.rnac_r:
            return

        new_twist = self.physical_model.evaluate_twist(location, self.state)
        # Find the correct insertion index
        possible_inserts = np.append(self.state[::3], np.inf)
        insert_idx = np.where(possible_inserts > location)[0][0]
        # Insert new polymerase state into the correct location
        self.state = np.insert(self.state, insert_idx * 3, [location, 0, new_twist])
        self.polymerase_state = np.insert(self.polymerase_state, insert_idx,
            [direction, location, self.next_uid], axis=0)
        self.next_uid += 1

    def remove_polymerase(self, polymerase_id):
        """
        Removes a polymerase by index. Given the polymerase status, this acts as termination (e.g.
        creates any necessary transcripts)
        
        Args:
        -----
        polymerase_id: An integer representing which polymerase to remove.
        """
        self.cleave_transcript(polymerase_id)
        self.state = np.delete(self.state, np.s_[(polymerase_id * 3):((polymerase_id + 1) * 3)])
        self.polymerase_state = np.delete(self.polymerase_state, polymerase_id, axis=0)

    def cleave_transcript(self, polymerase_id):
        """
        Simulates a polyA-site cleavage of of a transcript. This means that we update
        the transcripts counter for any genes included, and set the mRNA length back to zero.

        Args:
        -----
        polymerase_id: The polymerase ID to cleave.
        """
        # Compute transcript region
        region = sorted([self.polymerase_state[polymerase_id,1], self.state[3 * polymerase_id]])
        direction = self.polymerase_state[polymerase_id,0]
        # Update transcripts for any gene that falls within this region that is in the right direction
        for i, gene in enumerate(self.genes):
            if min(gene[1]) >= region[0] and max(gene[1]) <= region[1] and np.sign(gene[1][1] - gene[1][0]) == direction:
                # Transcript encompasses this gene
                self.transcripts[i] += 1
        # Update polymerase state, by cleaving mRNA length and bumping slightly forward
        self.state[(polymerase_id * 3) + 1] = 0
        self.state[polymerase_id * 3] += .0001 * direction
    
    def degrade_transcript(self, gene_idx):
        """
        Removes a transcript from a specific gene, reducing its count by one.

        Args:
        -----
        gene_idx: The index of the gene transcript to degrade.
        """
        if self.transcripts[gene_idx] > 0:
            self.transcripts[gene_idx] -= 1


@njit
def sample_discrete(probs):
    """Randomly sample an index with probability given by probs."""
    # Generate random number
    q = np.random.rand()
    
    # Find index
    i = 0
    p_sum = 0.0
    while p_sum < q:
        p_sum += probs[i]
        i += 1
    return i - 1

class SimulationResult(object):
    """
    An all-in-one container for the results of a stochastic simulation. Multiple helper functions
    are available that summarize this data, or give other outputs (such as a created movie)

    Member variables:
    -----------------
    gene_map: A list of gene names, where the i'th entry is the name of the i'th gene.
    transcript_snapshots: A (t x (g + 1)) ndarray storing the number of each gene transcript
        at a given simulation time. Each row is a snapshot, with the first entry being the time and the 
        next entries being the number of transcripts.
    polymerase_history: A list of 2-tuples, encoding (time_array, polymerase_state), e.g. a concatnation of 
        ODE solution results.
    events: A dictionary of the form {event_name_str : list}, where a list of times
        at which point the events occurred are listed.
    """
    def __init__(self, genes, bcs, tspan):
        """
        Create a SimulationResult container, with the necessary variables describing
        the simulation.

        Args:
        -----
        genes: A dictionary, where keys are gene names and values are (start, end) tuples
            encoding the start/end of genes.
        bcs: A 2-tuple encoding (Left_BoundaryCondition, Right_BoundaryCondition). 
            Each BC should be an element of type BoundaryCondition.
        tspan: A 2-tuple encoding (start_time, end_time)
        """
        self.genes = genes
        self.bcs = bcs
        self.tspan = tspan

        self.gene_map = [x[0] for x in sorted(genes.items(), key=lambda y: min(y[1]))]

        self.transcript_snapshots = np.array([[0] * (1 + len(self.gene_map))])
        self.polymerase_history = []
        self.events = {}
    
    def add_event(self, event_name, time):
        """
        Adds an event to the event record.

        Args:
        -----
        event_name: str, the type of event to record.
        time: double, the time at which the event occurred.
        """
        if event_name not in self.events:
            self.events[event_name] = [time]
        else:
            self.events[event_name].append(time)
    
    def add_transcript_snapshot(self, transcripts, time):
        """
        Adds an updated transcript count at the given time.

        Args:
        -----
        transcripts: A (g x 1) ndarray specifying the number of transcripts.
        time: double, the time at which this snapshot is valid.
        """
        self.transcript_snapshots = np.append(self.transcript_snapshots,
            np.concatenate(([time], transcripts))[np.newaxis,:], axis=0)
    
    def add_polymerase_history(self, time_arr, state_arr):
        """
        Update the polymerase activity history.

        Args:
        -----
        time_arr: A ndarray representing the times at which the state was simulated
        state:arr: A (t x 3N) ndarray representing the simulated state of the system.
        """
        self.polymerase_history.append((time_arr, state_arr))



class SimulationRunner(object):
    """
    This class is responsible for setting up and running simulations. It handles multiprocessing-based runs as well.

    After initialization, call the `simulate` function to generate SimulationResults
    """
    def __init__(self, genes, promoters, PASes, bcs, tspan, params=DEFAULT_PARAMS):
        """
        Create a new SimulationRunner, with all internal state necessary to run simulations.

        Args:
        -----
        genes: A dictionary, where keys are gene names and values are (start, end) tuples
            encoding the start/end of genes.
        promoter: A list of Promoter objects.
        PASes: A list of PAS objects.
        bcs: A 2-tuple encoding (Left_BoundaryCondition, Right_BoundaryCondition). 
            Each BC should be an element of type BoundaryCondition.
        tspan: A 2-tuple encoding (start_time, end_time)
        params: A dictionary encoding relevant keyvalue pairs representing physical constants.
            See the docstring for `create_physical_model`.
        """

        self.genes = genes
        self.promoters = promoters
        self.pas = PASes
        self.bcs = bcs
        self.tspan = tspan
        self.params = params
        self.physical_model = SimulationState.create_physical_model(bcs, params)
    
    def multi_simulation(self, n):
        """
        Calculates multiple runs in parallel
        """
        return [self.single_simulation() for _ in range(n)]

    def single_simulation(self ,idx=0):
        """
        Runs a simulation, returning a SimulationResult.
        """
        while True:
            try:
                state = SimulationState(self.tspan[0], self.bcs, self.genes, self.params, self.physical_model)
                result = SimulationResult(self.genes, self.bcs, self.tspan)

                timed_terminations  = [] # stores (polymerase_uid, stop_time)
                length_terminations = [] # stores (polymerase_uid, stop_distance)
                
                while state.time < self.tspan[1]:
                    # Perform a Gillespie draw to determine the next event to occur.
                    # For N promoters and G genes, we have N + G + 1 stochastic events to pull from.
                    # The first N are promoter initiation events, the next G are mRNA degradation events,
                    # and the last is topoisomerase activity.
                    #
                    # For this draw, we use base rates to determine when to check a condition.
                    stochastic_propensities = np.concatenate((np.array([x.base_rate for x in self.promoters])
                                                                * self.params["base_promoter_initiation_rate"],
                                                             state.transcripts * self.params["mRNA_deg_rate"],
                                                             np.array([self.params["topo_rate"]])))
                    for i, promoter in enumerate(self.promoters):
                        if promoter.type == PromoterType.SC_DEPENDENT:
                            # Multiply base rate by 3*kb*T
                            stochastic_propensities[i] *= np.exp(2 * state.physical_model.kb * state.params['temperature'])
                    gillespie_mean_time = 1 / np.sum(stochastic_propensities)
                    attempt_offset = np.random.exponential(gillespie_mean_time)


                    # Try to simulate until the next offset
                    if state.state.shape[0] > 0:
                        # Run an ODE model. We need stop conditions.
                        # The first stop condition is if any polymerase reached a PAS in the correct direction.
                        def reached_pas(t, y, i):
                            # Checks if polymerase i is close to a PAS.
                            return min([y[3*i] - pas.location for pas in self.pas if pas.direction == state.polymerase_state[i,0]])
                        def reached_time_termination(t, y, j):
                            return t - timed_terminations[j][1]
                        def reached_length_termination(t,y,j):
                            return y[3 * uid_mapping[length_terminations[j][0]]] - length_terminations[j][1]

                        def rebind_stop_func(func, index):
                            def bound(t, x):
                                return func(t, x, index)
                            return bound
                        

                        num_polymerases = state.polymerase_state.shape[0]
                        # Calculate UID to index mapping
                        uid_mapping = {}
                        for i in range(num_polymerases):
                            uid_mapping[state.polymerase_state[i,2]] = i
                        # Remove any unnecessary termination commands
                        timed_terminations = [x for x in timed_terminations if x[0] in uid_mapping]
                        length_terminations = [x for x in length_terminations if x[0] in uid_mapping]



                        num_time_terminations = len(timed_terminations)
                        num_length_terminations = len(length_terminations)
                        stop_funcs = ([rebind_stop_func(reached_pas, i) for i in range(num_polymerases)]
                                    + [rebind_stop_func(reached_time_termination,i) for i in range(num_time_terminations)]
                                    + [rebind_stop_func(reached_length_termination, i) for i in range(num_length_terminations)])

                        for i in range(len(stop_funcs)):
                            stop_funcs[i].terminal = True

                        ode_result = scipy.integrate.solve_ivp(
                            lambda t,y:self.physical_model.system_derivatives(y,
                                                    state.polymerase_state[:,0].flatten()),
                            (state.time, state.time + attempt_offset),
                            state.state,
                            events=stop_funcs,
                            method='RK45')
                        
                        if ode_result.status == -1:
                            print('-', end='', flush=True)
                            raise RuntimeError('Integration error')
                        
                        state.time = ode_result.t[-1]
                        state.state = ode_result.y[:,-1]
                        result.add_polymerase_history(ode_result.t, ode_result.y)

                        if ode_result.status == 1:
                            # Handle an event and continue
                            event_idx = np.where(np.array(
                                [len(x) for x in ode_result.t_events]) > 0)[0][0]
                            
                            if event_idx < num_polymerases:
                                # Identify which PAS was closest:
                                pas_idx = np.argmin(np.abs(np.array([pas.location - state.state[event_idx * 3] for pas in self.pas])))
                                pas = self.pas[pas_idx]
                                # Handle PAS
                                if pas.type == PAS_Type.NO_RT:
                                    state.remove_polymerase(event_idx)
                                elif pas.type == PAS_Type.TIME_RT:
                                    state.cleave_transcript(event_idx)
                                    timed_terminations.append((state.polymerase_state[event_idx,2],
                                        state.time + np.random.exponential(pas.readthrough_duration)))
                                elif pas.type == PAS_Type.LENGTH_RT:
                                    state.cleave_transcript(event_idx)
                                    length_terminations.append((state.polymerase_state[event_idx,2],
                                        state.state[event_idx * 3] + np.random.exponential(
                                        pas.readthrough_duration)))
                                else:
                                    raise RuntimeError('Invalid PAS type')

                                pass
                            elif event_idx < num_polymerases + num_time_terminations:
                                # Time termination
                                idx = event_idx - num_polymerases
                                state.remove_polymerase(uid_mapping[timed_terminations[idx][0]])
                                del timed_terminations[idx]
                            else:
                                # Space termination
                                idx = event_idx - num_polymerases - num_time_terminations
                                state.remove_polymerase(uid_mapping[length_terminations[idx][0]])
                                del length_terminations[idx]
                            
                            result.add_transcript_snapshot(state.transcripts, state.time)
                            # We're done, skip the Gillespie draw.
                            continue

                    # If we are here, we reached the final offset. This means our initial Gillespie draw
                    # was valid.
                    event_idx = sample_discrete(stochastic_propensities / np.sum(stochastic_propensities))

                    assert(event_idx < len(stochastic_propensities))
                    if event_idx == len(stochastic_propensities) - 1:
                        # Topo rate
                        result.add_event('topo', state.time)
                        state.remove_supercoiling()
                    elif event_idx < len(self.promoters):
                        # Try initating a promoter
                        promoter = self.promoters[event_idx]
                        if promoter.type == PromoterType.STATIC:
                            pass
                        elif promoter.type == PromoterType.SC_DEPENDENT:
                            # Check for supercoiling dependence, that is 1.2 * 2 Pi * torque at that location. Do this by finding a probability
                            # and multiplying by the exp[-E / (kb * t)]
                            energy = state.calculate_torque(promoter.location) * 1.2 * 2.0 * np.pi
                            # Adjust for 2kB T shift
                            if np.random.random() > np.exp((-energy / (state.physical_model.kb * state.params['temperature'])) - 2):
                                # Skip; the initiation doesn't happen
                                continue
                        else:
                            raise RuntimeError('Invalid promoter type')

                        state.add_polymerase(promoter.location, promoter.direction)
                    else:
                        # Degrade a mRNA
                        transcript_idx = event_idx - len(self.promoters)
                        state.degrade_transcript(transcript_idx)
                        result.add_transcript_snapshot(state.transcripts, state.time)
                print('+', end='', flush=True)
                return result
            except (RuntimeError, FloatingPointError):
                continue # Retry

def write_single_frame(run_result, interp_expression, frame_times, frame_idx, genes, x_range_override, png_name):
    """
    Writes a single animation frame to disk as a PNG. This
    is intended as a helper function to be run in parallel.
    
    Args:
    -----
    run_result: A postprocessed simulation dataset.
    interp_expression: A interpolated expression dataset to use for plotting
    frame_times: A list of available frame times to plot at.
    frame_idx: The current frame being generated
    genes: A list of genes to draw on the plot.
    x_range_override: A 2-tuple encoding an overridden x-range, if specified
    png_name: The filename to save the figure at.
    
    Returns:
    ---------
    The saved figure size, in pixels.
    """
    time = frame_times[frame_idx]
    i = frame_idx
    # Inspired, but not identical to https://stackoverflow.com/a/43775224
    def multi_interp(x, xp, fp):
        j = np.searchsorted(xp, x) - 1
        d = (x - xp[j]) / (xp[j + 1] - xp[j])
        return (1 - d) * fp[j, :] + fp[j + 1, :] * d
    
    fig, axs = plt.subplots(1, 2, figsize=(10, 3), gridspec_kw={'width_ratios': [2, 1]})
       
    # Locate which raw segment this frame is in
    lookup_results =  np.where([np.any(r[0] <= time) & np.any(r[0] >= time) for r in run_result['raw']])
       
    excess_twist = multi_interp(time, run_result['time'], run_result['excess_twist'])
    sc_density = multi_interp(time, run_result['time'], run_result['sc_density'])

    x_domain = run_result['x_domain']
        
    # Plot genes
    for gene in genes:
        axs[0].plot([gene[0], gene[1]], [0, 0],
                linewidth=20, alpha=.5, zorder=1)

    # Plot DNA segments
    domain_points = np.array([x_domain, np.zeros(x_domain.shape)]).T.reshape(-1,1,2)
    segments = np.concatenate([domain_points[:-1], domain_points[1:]], axis=1)
    max_excursion = np.max(np.abs([np.min(run_result['sc_density']), np.max(run_result['sc_density'])]))
    cmap_norm = plt.Normalize(-.125, .125)
    dna_segments = LineCollection(segments, cmap='Spectral', norm=cmap_norm)
    dna_segments.set_array(sc_density)
    dna_segments.set_linewidth(8)
    line = axs[0].add_collection(dna_segments)
    fig.colorbar(line, ax=axs[0])
    if len(lookup_results[0]) > 0:
        in_range = lookup_results[0][0]
        polymerase_state = multi_interp(time, run_result['raw'][in_range][0], run_result['raw'][in_range][1].T)

        # Plot polymerases
        polymerase_locations = polymerase_state[::3]
        polymerase_transcripts = polymerase_state[1::3]
        axs[0].plot(polymerase_locations, np.zeros(polymerase_locations.shape), 'k.')

        n_polymerases = len(polymerase_locations)
        transcript_segments = np.zeros((n_polymerases, 2, 2))
        transcript_segments[:,:,0] = polymerase_locations[:,np.newaxis]
        transcript_segments[:,1,1] = polymerase_transcripts

        transcript_line = axs[0].add_collection(LineCollection(transcript_segments))
    axs[0].set_xlim(x_domain[0], x_domain[-1])
    if x_range_override is not None:
        axs[0].set_xlim(*x_range_override)
    axs[0].set_ylim(-0.5 * max([abs(g[1] - g[0]) for g in genes]), max([abs(g[1] - g[0]) for g in genes]))
    axs[0].get_yaxis().set_visible(False)
    axs[0].set_xlabel('Distances (nm)')
    # Plot transcript levels
    axs[1].plot(frame_times[:i], interp_expression[:,:i].T)
    axs[1].set_xlim(frame_times[0], frame_times[-1])
    axs[1].set_ylim(0, np.max(interp_expression))
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Expression level')

    #plt.show()
    plt.tight_layout()
    plt.savefig(png_name)
    pixel_size = plt.gcf().get_size_inches()*plt.gcf().dpi 
    plt.close()
    return pixel_size

def gen_movie(run_result, genes, n_frames, name, outfolder, make_mp4=False, x_range_override=None):
    """
    Given a postprocessed run result, makes a movie of the state of the system.
    The output movie is either a series of PNGs or a generated .mp4
    
    Args:
    -----
    run_result: A postprocessed simulation run result.
    genes: Gene definitions, for plotting the gene boundaries.
    n_frames: The number of frames to generate images at.
    name: The name of the output movie. If outputting PNGs, filenames are appended with the frame number.
        If outputing an mp4, the filename is directly output.
    outfolder: The folder to save the movie to.
    make_mp4: If mp4 output is preferred. This requires ffmpeg installed and accessible on PATH!
    x_range_override: Optionally gives the x axis range for the DNA view.
    """
    def expression_interp(t, t_known, known_expression):
        time_idx = np.array([np.where(t_known <= tp)[0][-1] for tp in t])
        return known_expression[time_idx, :].T

    if not make_mp4:
        png_folder = outfolder
    else:
        tempdir = tempfile.TemporaryDirectory()
        png_folder = tempdir.name
    n_digits = int(np.ceil(np.log10(n_frames)))
    
    frame_times = np.linspace(run_result['time'][0] + .1, run_result['time'][-1] - .1, n_frames)
    interp_expression = expression_interp(frame_times, run_result['time'], run_result['gene_expression'])
    
    for i, time in enumerate(frame_times):
        if i == 0:
            continue
        pixel_size = write_single_frame(run_result, interp_expression,
                                        frame_times, i, genes,
                                        x_range_override,
                                        os.path.join(png_folder, '{}_{:0{}}.png'.format(name, i, n_digits)))
    if make_mp4:
        gc.collect()
        output = subprocess.run(['ffmpeg', '-y', '-r', '60', '-f', 'image2', '-s',
                                      '{}x{}'.format(int(pixel_size[0]), int(pixel_size[1])),
                                      '-i',
                                      os.path.join(png_folder, '{}_%0{}d.png'.format(name, n_digits)),
                                      '-vcodec', 'libx264', '-crf',
                                      '25', '-pix_fmt', 'yuv420p',
                                      os.path.join(outfolder, '{}.mp4'.format(name))],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        tempdir.cleanup()
        if output.returncode != 0:
            print(output)
            raise RuntimeError('ffmpeg call failed!')