import scipy.integrate
import numpy as np
from numba import njit, vectorize, float32, float64, boolean
from numba.experimental import jitclass
import multiprocessing
import pandas as pd
import tempfile
import os
import gc
import subprocess
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

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
        Constructor that takes all relevant system constants and initalizes
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
        rnac_r: The radius of a RNAP (used for hard-sphere replusion)
        chi: The DNA twisting mobility (s pN nm)
        eta: The The mRNA drag coefficent (pN (nm^(alpha - 1)))
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

def bind_supercoiling_model(params, bcs):
    """
    Given the system parameters and boundary conditions,
    returns a PolymeraseSupercoilingModel representing the system.
    
    Args:
    -----
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
    bcs: A 2-tuple encoding (left_bc, right_bc). A BC is either
        the string "free", or a tuple of (location, value).
    
    Returns:
    --------
    A SupercoilingPolymeraseModel with all relevant values set.
    """
    # Calculate BCs
    left_bc_free = bcs[0] == 'free'
    right_bc_free = bcs[1] == 'free'
    left_bc = (0, 0)
    right_bc = (0, 0)
    if not left_bc_free:
        left_bc = bcs[0]
    if not right_bc_free:
        right_bc = bcs[1]
    
    return SupercoilingPolymeraseModel(
            params['DNA_bend_plength'], params['DNA_twist_plength'],
            params['DNA_plectonome_twist_plength'], params['temperature'],
            params['DNA_force'], params['RNAP_velocity'],
            params['RNAP_torque_cutoff'], params['RNAP_stall_torque_width'],
            params['RNAP_radius'],
            params['DNA_twist_mobility'], params['mRNA_drag'],
            params['mRNA_exponent'], left_bc_free, right_bc_free,
            left_bc[0], right_bc[0], left_bc[1], right_bc[1])

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

def attach_on_call(func):
    """Decorator for adding a stochastic event."""
    def wrap(self):
        self.stochastic_events.append(func(self))
    return wrap

class SupercoilingSimulation(object):
    """
    A base class that encodes a Gillepse-algorithm method
    for stochastically simulating transcript initiation, alongside
    a supercoiling model.
    
    This code tracks the location of polymerases as they move over
    the genome. Events such as polymerase addition and removal are
    natively tracked using built-in functions. Other genomic activity
    such as topoisomerase activity can be extended with this model.
    
    A decorator @attach_on_call is provided that allows arbitrary hooking
    into the event system. Events can be predicated on an ODE
    integrated state or as a Gillepsie-style stochastic event.
    """
    
    def __init__(self, params, bcs, genes):
        """
        Takes the problem formulation and initalizes the
        simulation internals.
        
        Additionally initalizes an stochastic_events list, that contains
        events in the form of tuples
        (base_rate, state_multiplier_func,state_mutate_function)
        
        Args:
        -----
        params: A dictionary containing all relevant physical paramteres. See
            the docstring of `bind_supercoiling_model` for details.
            The extra parameters expected by this method are:
                base_promoter_initiation_rate: Rate that promoters add polymerases (1/sec)
        bcs: A 2-tuple encoding (left_bc, right_bc). A BC is either
            the string "free", or a tuple of (location, value).
        genes: An optional list of tuples encoding the arguments to add_stochastic_gene_init_event
            (gene_start, gene_end, promoter_strength, readthrough_func, terminate_at_tes).
            The gene direction is inferred from these values.
        """
        self.model = bind_supercoiling_model(params, bcs)
        self.genes = [(g[0], g[1]) for g in genes]
        self.base_promoter_rate = params['base_promoter_initiation_rate']
        self.topo_rate = params['topo_rate']
        self.mrna_deg_rate = params['mRNA_deg_rate']
        
        self.stochastic_events = []
        
        for g in genes:
            self.add_stochastic_gene_init_event(*g)
    
    @attach_on_call
    def enable_topo_relaxation(self):
        """
        Returns the stochastic event information
        required to relax DNA.
        
        Returns:
        --------
        A tuple (base_rate, state_dependent_rate, mutate_func).
        Here, the base rate is set to the relaxation rate constant
        and the state_dependent_rate is always set to 1. Mutate_func is
        the key function that removes supercoiling from a gene pair.
        """
        state_dependent_func = lambda x: 1.0
        def remove_supercoiling(self, t, state):
            if len(state) == 0:
                return
            # Sort genes by physical location
            sorted_genes = sorted(self.genes, key=lambda x: min(x))
            # We now rewrite all polymerases that are sitting on the range
            # (idx).start to (idx + 1).end
            # We do this by removing all relevant polymerases from state, then
            # using the model.compute_twist function to recalculate, at each
            # of the polymerases to rewrite
            if len(sorted_genes) == 1:
                lower_bound = -np.inf
                upper_bound = np.inf
            else:
                relax_after_idx = np.random.choice(range(len(sorted_genes) - 1))
                lower_bound = min(sorted_genes[relax_after_idx])
                upper_bound = max(sorted_genes[relax_after_idx + 1])
            polymerase_rewrite_idx = np.where(
                (state[::3] > lower_bound) & (state[::3] < upper_bound))
            polymerase_positions = state[::3][polymerase_rewrite_idx]
            if len(polymerase_positions) == 0:
                # No polymerases in this range: just return
                return
            # Exploit polymerases always being sorted
            new_state = np.delete(state, range(np.min(polymerase_rewrite_idx) * 3,
                                               (np.max(polymerase_rewrite_idx) + 1) * 3))
            self.simstate[2::3][polymerase_rewrite_idx] = \
                self.model.evaluate_twist(polymerase_positions, new_state)
            
            if 'topo' not in self.event_times:
                self.event_times['topo'] = []
            self.event_times['topo'].append(t)
        
        return (self.topo_rate, state_dependent_func, remove_supercoiling)
    
    def add_polymerase(self, location, direction,
                       termination_func, terminate_gene_names=[], event_functions=[]):
        """
        Adds a polymerase to the current state. A termination function is passed
        which specifies when the polymerase should be removed from the simulation.
        
        Additional event functions can be passed which encode when extra ODE events
        should be triggered, in addition to a function
        specifying what happens when that event occurs. This can be used to simulate
        things such as transcript cleavage, stalling, etc.
        
        Args:
        -----
        location: A float encoding the start location
        direction: An integer, either +1 or -1, to indicate polymerase motion direction
        termination_func: A function that takes two arguments, returns 0 when termination should occur
            X: the 3N number of states encoding (position, transcript_length, excess_twist)
            i: The index of the polymerase you belong to
        terminate_gene_names: A list of strings specifying gene names. When RNAP termination
            occurs, expression of these gene names is increased by one.
        event_functions: A list of tuples (event_func, mutate_func)
            event_func: A function with the same signature as above.
            mutate_func: A function that takes the time, state, and 
               polymerase index and mutates them somehow, returning nothing.
        """
        new_twist = self.model.evaluate_twist(location, self.simstate)
        # Find the correct insertion index
        possible_inserts = np.append(self.simstate[::3], np.inf)
        insert_idx = np.where(possible_inserts > location)[0][0]
        # Insert new polymerase state into the correct location
        self.simstate = np.insert(self.simstate, insert_idx * 3, [location, 0, new_twist])
        self.tes_time = np.insert(self.tes_time, insert_idx, np.inf)
        self.polymerase_directions = np.insert(self.polymerase_directions, insert_idx, direction)
        
        # Add ODE functions
        def terminate_transcription(self, t, state, i):
            self.simstate = np.delete(self.simstate, np.s_[(i * 3):((i + 1) * 3)])
            self.polymerase_directions = np.delete(self.polymerase_directions, i)
            self.tes_time = np.delete(self.tes_time, i)
            del self.simstate_ode_funcs[i]
            # Add event time
            if 'genes' not in self.event_times:
                self.event_times['genes'] = {}
            for name in terminate_gene_names:
                if name not in self.event_times['genes']:
                    self.event_times['genes'][name] = []
                self.event_times['genes'][name].append(t)
                
        self.simstate_ode_funcs.insert(insert_idx,
                                       [(termination_func, terminate_transcription)] + event_functions)

    def add_stochastic_gene_init_event(self, gene_start, gene_end, strength,
                                       gene_names=None,
                                       readthrough_mean=0.0,
                                       readthrough_in_nm=True,
                                       state_func=lambda x: 1.0):
        """
        Initalizes a simple gene which has a state-independent initiation
        rate, and whose spawned polymerases debind exactly at the gene end.
        
        Polymerases are not allowed to add if there is a polymerase already
        occupying the starting site.
        
        Args:
        -----
        gene_start: The location of the gene starting location
        gene_end: The location of the gene ending location
        strength: The state-independent multiplier on the base polymerase addition rate
        gene_names: A list of gene names that this gene spans.
        readthrough_mean: The mean of an exponential distribution used to calcualte
            readthrough distance/time.
        readthrough_in_nm: If true, the readthrough func is a distance from the TES.
            If false, the readthrough func is a time in seconds past the TES.
        state_func: An optional function that takes the state of the system (3N vector)
            and returns an additional state-dependent multiplier
        
        Side effects:
        -------------
        Adds this gene to the stochastic_events struct.
        """
        if gene_names is None:
            gene_names = [str(gene_start)]
        direction = np.sign(gene_end - gene_start)
        def mutate_func(self, t, state):
            # Our mutate function simply calls add_polymerase.
            # We terminate at gene_end + readthrough_func() if in nm
            # and at tes_time + readthrough_func() if in seconds
            # If sthere is > 0nm of readthrough,
            # we add a state lambda that increases expression levels.
            readthrough = np.random.exponential(readthrough_mean)
            if readthrough == 0:
                self.add_polymerase(gene_start, direction,
                                    lambda s, t, x, i: x[3*i] - (gene_end + readthrough),
                                    gene_names)
            else:
                # Otherwise, we need to add a mutate function
                # that updates transcript levels at the TES site.
                def tes_update_helper(self, t, state, i):
                    # Set transcript length to zero, bumping polymerase slightly forward.
                    self.simstate[(i * 3) + 1] = 0
                    self.simstate[i * 3] += .0001 * direction
                    self.tes_time[i] = t # store transcript release time.
                    # Update genes expressed from released transcript
                    if 'genes' not in self.event_times:
                        self.event_times['genes'] = {}
                    for name in gene_names:
                        if name not in self.event_times['genes']:
                            self.event_times['genes'][name] = []
                        self.event_times['genes'][name].append(t)

                if readthrough_in_nm:
                    termination_func = lambda s, t, x, i: x[3 * i] - (gene_end + direction * readthrough)
                else:
                    termination_func = lambda s, t, x, i: t - (s.tes_time[i] + readthrough)
                self.add_polymerase(gene_start, direction,
                                    termination_func,
                                    event_functions=[
                                        (lambda s, t, x, i: x[3 * i] - gene_end,
                                        tes_update_helper)
                                    ])
        # Define a function that stops polymerase addition if there is already a polymerase
        # occupying the site
        reject_occupied_site = lambda x: np.min(
            np.abs(np.concatenate((np.array([np.inf]), x[::3])) - gene_start)
        ) > self.model.rnac_r
        final_state_func = lambda x: reject_occupied_site(x) * state_func(x)
        self.stochastic_events.append((self.base_promoter_rate * strength,
                                       final_state_func,
                                       mutate_func))
        
    def simulate(self, tspan):
        """
        Given a timespan, integrates using a Gillepse-inspired ODE
        solution.
        
        Args:
        -----
        tspan: A tuple containing (t_start, t_end) [s]
        
        Returns:
        --------
        A list of tuples, containing (times, states) for each integration interval.
        """
        # Reset sim state
        
        success = False
        while not success:
            try:
                self.simstate = np.zeros(0)
                self.simstate_ode_funcs = []
                self.polymerase_directions = np.zeros(0)
                self.event_times = {'genes': {}}
                self.tes_time = np.zeros(0)
                result = []

                last_t = tspan[0]

                while last_t < tspan[1]:
                    # Make a Gillepse draw to find the next event to occur
                    gillepse_mean_time = 1 / np.sum([x[0] for x in self.stochastic_events])
                    next_attempt_time = last_t + np.random.exponential(gillepse_mean_time)

                    while last_t < next_attempt_time:
                        if len(self.simstate) > 0:
                            # Do an ODE simulation
                            stop_events = []
                            mutate_funcs = []
                            def rebind_event_func(event, p_idx):
                                def bound(t, x):
                                    return event[0](self, t, x, p_idx)
                                return bound
                            def rebind_mutate_func(event, p_idx):
                                def bound(s, t, x):
                                    return event[1](s, t, x, p_idx)
                                return bound
                            
                            for polymerase_idx, events in enumerate(self.simstate_ode_funcs):
                                # Weird lambda currying. MWE that explains the problem:
                                # [(lambda curried_i:(lambda x: x * curried_i))(i) for i in range(10)]
                                # vs 
                                # [lambda x: x * i for i in range(10)]
                                # https://stackoverflow.com/a/34021333
                                stop_events += [rebind_event_func(event, polymerase_idx)
                                                for event in events]
                                mutate_funcs += [rebind_mutate_func(event, polymerase_idx)
                                                 for event in events]
                            for i in range(len(stop_events)):
                                stop_events[i].terminal = True

                            ode_result = scipy.integrate.solve_ivp(
                                lambda t,y:self.model.system_derivatives(y,
                                                        self.polymerase_directions),
                                (last_t, next_attempt_time),
                                self.simstate,
                                events=stop_events,
                                method='RK45')

                            if ode_result.status == -1:
                                print(f'INTEGRATION FAILURE: {ode_result.message}, restarting')
                                raise RuntimeError('Integration error')
                            result.append((ode_result.t, ode_result.y))
                            last_t = ode_result.t[-1]
                            self.simstate = ode_result.y[:,-1]

                            if ode_result.status == 1:
                                # A polymerase hit a stop event site. See which one it is:
                                event_idx = np.where(np.array(
                                    [len(x) for x in ode_result.t_events]) > 0)[0][0]
                                # Mutate our state
                                mutate_funcs[event_idx](self, last_t, self.simstate)
                        else:
                            # Just advance the time directly; no polymerases right now
                            last_t = next_attempt_time

                    # Perform a Gillepsie draw
                    event_probs = gillepse_mean_time * np.array([x[0] for x in self.stochastic_events])

                    #self.init_mean_time = 1 / np.sum(self.promoter_strength * self.base_promoter_rate)
                    #self.init_probability = self.promoter_strength * self.base_promoter_rate * self.init_mean_time
                    next_stochastic_idx = sample_discrete(event_probs)
                    # Apply the state-specific state value to see if we perform this event
                    if np.random.random() < self.stochastic_events[next_stochastic_idx][1](self.simstate):
                        # Perform the stochastic state mutation
                        self.stochastic_events[next_stochastic_idx][2](self, last_t, self.simstate)
            except (RuntimeError, FloatingPointError):
                continue # Retry
            success = True
        return (result, self.event_times)

    def postprocess_run(self, sim_run, domain_points=1000, domain_endpoints=None):
        """
        Given the raw run data, calculates parameters
        of interest. Namely, this is the supercoiling density,
        the number of mRNAs over time (calculated using the degradation rate)
        
        Args:
        -----
        sim_run: A tuple of the form (raw, event_times):
            raw: A list of tuples, of the form (time_ndarray, 3N-state_ndarray)
            event_times: A dictionary containing details on when events occured
        domain_ponts: The number of points used to discritize space
        domain_endpoints: (optional) A tuple containing (domain_start, domain_end)
            to explicitly discritize over.
        
        Returns:
        --------
        A dictionary containing the following key-value pairs:
            raw: The input sim_run
            time: A single ndarray containing all timepoints
            x_domain: A vector of the locations that were sampled
            excess_twist: A sample of the excess twist, phi, across the domain for each timepoint
            sc_density: The supercoiling density across the domain for each timepoint
            gene_expression: The expression of each gene over time (calculated with a secondary
                Gillepsie run)
        """
        raw, events = sim_run
        result = {'raw': raw}
        
        result['time'] = np.concatenate(list(x[0] for x in raw))
        min_x = min([np.min(x[1][::3]) for x in raw])
        max_x = max([np.max(x[1][::3]) for x in raw])
        if domain_endpoints is None:
            domain_endpoints = (min_x, max_x)
        result['x_domain'] = np.linspace(*domain_endpoints, domain_points)
        result['excess_twist'] = np.zeros((len(result['time']), domain_points))
        compute_idx = 0
        for run in raw:
            for i in range(len(run[0])):
                result['excess_twist'][compute_idx, :] = self.model.evaluate_twist(
                    result['x_domain'], run[1][:,i])
                compute_idx += 1
        dx = (max_x - min_x) / domain_points
        result['sc_density'] = np.diff(result['excess_twist']) / (dx * -self.model.w0)
        
        # Now compute mRNA at each timepoint
        # The propensity of mRNA degrading is mrna_deg_rate * num_of_mrna
        gene_names = sorted(events['genes'].keys())
        gene_expression = np.zeros((len(result['time']),len(gene_names)))
        gene_idx = 0
        for gene_name in gene_names:
            time_accum = np.array([0])
            mRNA_accum = np.array([0])
            
            creation_events = events['genes'][gene_name]
            time_accum = np.append(time_accum, np.array([creation_events[0]]))
            mRNA_accum = np.append(mRNA_accum, np.array([1]))
            creation_idx = 1
            while time_accum[-1] < result['time'][-1]:
                # Draw the time to the next degradation
                if mRNA_accum[-1] > 0:
                    next_deg = time_accum[-1] + \
                        np.random.exponential(1.0 / (mRNA_accum[-1] * self.mrna_deg_rate))
                else:
                    next_deg = result['time'][-1]
                # Check if a creation event occurs
                if creation_idx < len(creation_events):
                    if next_deg < creation_events[creation_idx]:
                        # Degradation occured!
                        time_accum = np.append(time_accum, np.array([next_deg]))
                        mRNA_accum = np.append(mRNA_accum, np.array([mRNA_accum[-1] - 1]))
                        continue
                    # Otherwise, a creation event occurs
                    time_accum = np.append(time_accum, np.array([creation_events[creation_idx]]))
                    mRNA_accum = np.append(mRNA_accum, np.array(mRNA_accum[-1] + 1))
                    creation_idx += 1
                else:
                    time_accum = np.append(time_accum, np.array([next_deg]))
                    mRNA_accum = np.append(mRNA_accum, np.array([mRNA_accum[-1] - 1]))
            
            def expression_interp(t, t_known, known_expression):
                time_idx = np.array([np.where(t_known <= tp)[0][-1] for tp in t])
                return known_expression[time_idx]
            gene_expression[:,gene_idx] = expression_interp(result['time'], time_accum, mRNA_accum)
            gene_idx += 1
        result['gene_expression'] = gene_expression
        result['gene_names'] = gene_names
        return result

def expression_single_run(params, bcs, genes, gene_names, sim_time, i=0):
    sim = SupercoilingSimulation(params, bcs, genes)
    sim.enable_topo_relaxation()
    result = sim.postprocess_run(sim.simulate((sim_time[0], sim_time[1])))
    times = np.linspace(*sim_time)
    id_val = np.ones(times.shape, dtype=int) * i
    df_result = pd.DataFrame(data={'id': id_val, 'time': times})
    for name in gene_names:
        if name not in result['gene_names']:
            df_result[name + '_expression'] = np.zeros(times.shape)
        else:
            df_result[name + '_expression'] = np.interp(times,
                                                        result['time'],
                                                        result['gene_expression'][:,result['gene_names'].index(name)])
    return df_result
    
def bulk_simulation(params, bcs, genes, gene_names, sim_time, n_runs):
    """
    Given the relevant parameters to a supercoiling simulation,
    runs several rounds of simulations in order to reach aggregate averages.
    
    Args:
    -----
    params, bcs, genes: Parameters required by the SupercoilingSimulation constructor
    sim_time: A tuple containing (start_time, end_time, n_points) over which the simulation should be run.
    n_runs: The number of runs to aggregate
    
    Returns:
    --------
    A single dataframe containing the columns:
        run_id: Integer containing which run it came from
        time: The time value of the datapoint
        gene_names: A series of columns with name = each of the names in gene_names
    """
    with multiprocessing.Pool() as p:
        runs = [p.apply_async(expression_single_run, args=(params, bcs, genes, gene_names, sim_time, i))
                for i in range(n_runs)]
        p.close()
        p.join()
        return pd.concat([r.get() for r in runs])

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