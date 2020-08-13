import numpy as np
import pandas as pd
import scipy
import sys

import stochastic_sc_model as ssc_model

def bulk_topology_simulation(geometry, dx_in_bp, expression_levels, simtime, n_simulations):
    """
    For a given gene geometry, with a interleaving space, runs a series of bulk simulations and
    returns the output dataframe.
    
    Args:
    -----
    geometry: A string containing:
        'tandem': Both gene A and B face from left to right: A-> |dx| B->
        'convergent': Gene A and B face each other: A-> |dx| <-B
        'divergent': Gene A and B face away from each other: <-A |dx| B->
    dx_in_bp: The inter-gene spacing in base pairs.
    expression_levels: A tuple containing expression levels for genes A and B.
    simtime: The amount of time to run the simulation for.
    n_simulations: The number of simulations to run

    Returns:
    --------
    A list of summary values, which are [geometry, dx_in_bp, A_expression, A_mean, A_std, B_mean, B_std]
    """
    endcap_distances = 3000 * .34
    gene_distance = 2000 * .34
    dx = dx_in_bp * .34
    
    gene_A_boundaries = (endcap_distances, endcap_distances + gene_distance)
    gene_B_boundaries = (gene_A_boundaries[1] + dx, gene_A_boundaries[1] + dx + gene_distance)
    end_barrier = gene_B_boundaries[1] + endcap_distances
    
    barriers = ((0,0), (end_barrier, 0))
    genes = []
    # Gene A faces to the right if tandem/convergent, otherwise left for divergent
    genes.append((gene_A_boundaries[1], gene_A_boundaries[0], expression_levels[0], 'A')
                 if geometry == 'divergent' else
                 (gene_A_boundaries[0], gene_A_boundaries[1], expression_levels[0], 'A'))
    # Gene B faces to the right except in the convergent case
    genes.append((gene_B_boundaries[1], gene_B_boundaries[0], expression_levels[1], 'B')
                 if geometry == 'convergent' else
                 (gene_B_boundaries[0], gene_B_boundaries[1], expression_levels[1], 'B'))
    
    raw_data = ssc_model.bulk_simulation(params, barriers, genes, ['A', 'B'], (0, simtime, 1000), n_simulations)
    mean_data = raw_data.groupby('time').mean().reset_index()
    mean_data = mean_data[mean_data['time']==simtime]
    std_data = raw_data.groupby('time').std().reset_index()
    std_data = std_data[std_data['time']==simtime]

    return [geometry, dx_in_bp, expression_levels[0], expression_levels[1],
            mean_data['A_expression'].iloc[0], std_data['A_expression'].iloc[0],
            mean_data['B_expression'].iloc[0], std_data['B_expression'].iloc[0]]

if __name__ == '__main__':
    # Ensure we don't hit floating point errors
    np.seterr('raise')
    np.seterr(under='ignore')

    params = {
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

    # We are given an arg with inter-gene spacing. We will run at 10 different upstream expression levels
    accum = []
    for geometry in ['tandem', 'convergent', 'divergent']:
        for expression in np.logspace(-2,0,10):
            accum.append(bulk_topology_simulation(geometry, int(sys.argv[1]), (expression, 1), 12000, 5000))
            print(f'Finished {geometry}/expression:{expression}')
    df = pd.DataFrame(accum, columns=['geometry', 'spacing', 'A_promoter_strength', 'B_promoter_strength', 'A_mean', 'A_std', 'B_mean', 'B_std'])
    df.to_feather(f'output/dataframes/dx{expression}.feather')