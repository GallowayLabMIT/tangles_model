import numpy as np
import pandas as pd
import scipy
import sys
import argparse

import stochastic_sc_model as ssc_model

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


def bulk_topology_simulation(geometry, dx_in_bp, topo_rate_multiplier, expression_levels, simtime, n_simulations):
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
    topo_rate_multiplier: The topoisomerase base rate to use. This value is multiplied by 1/1200 1/s.
    expression_levels: A tuple containing expression levels for genes A and B.
    simtime: The amount of time to run the simulation for.
    n_simulations: The number of simulations to run

    Returns:
    --------
    A list of summary values, which are
    [geometry, dx_in_bp, topo_rate, A_expression, B_expression, A_mean, A_std, B_mean, B_std]
    """
    endcap_distances = 3000 * .34
    gene_distance = 2000 * .34
    dx = dx_in_bp * .34
    
    gene_A_boundaries = (endcap_distances, endcap_distances + gene_distance)
    gene_B_boundaries = (gene_A_boundaries[1] + dx, gene_A_boundaries[1] + dx + gene_distance)
    end_barrier = gene_B_boundaries[1] + endcap_distances
    
    barriers = ((0,0), (end_barrier, 0))
    genes = []

    # Update the state function (last gene argument) to use supercoiling dependent.
    # Gene A faces to the right if tandem/convergent, otherwise left for divergent
    genes.append((gene_A_boundaries[1], gene_A_boundaries[0], expression_levels[0], 'A')
                 if geometry == 'divergent' else
                 (gene_A_boundaries[0], gene_A_boundaries[1], expression_levels[0], 'A'))
    # Gene B faces to the right except in the convergent case
    genes.append((gene_B_boundaries[1], gene_B_boundaries[0], expression_levels[1], 'B')
                 if geometry == 'convergent' else
                 (gene_B_boundaries[0], gene_B_boundaries[1], expression_levels[1], 'B'))
    
    params_copy = dict(params)
    params_copy['topo_rate'] *= topo_rate_multiplier

    raw_data = ssc_model.bulk_simulation(params_copy, barriers, genes, ['A', 'B'], (0, simtime, 1000), n_simulations)
    mean_data = raw_data.groupby('time').mean().reset_index()
    mean_data = mean_data[mean_data['time']==simtime]
    std_data = raw_data.groupby('time').std().reset_index()
    std_data = std_data[std_data['time']==simtime]

    return [geometry, dx_in_bp, params_copy['topo_rate'],
            expression_levels[0], expression_levels[1],
            mean_data['A_expression'].iloc[0], std_data['A_expression'].iloc[0],
            mean_data['B_expression'].iloc[0], std_data['B_expression'].iloc[0]]
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Python helper code to run a parameter sweep")
    parser.add_argument('--runner_id', type=int, required=True, help='The runner index of this instance')
    parser.add_argument('--num_runners', type=int, required=True, help='The total number of runners that work is divided over')
    parser.add_argument('--n_simulations', type=int, default=5000, required=True, help='The number of simulations to run for each parameter set')
    args = parser.parse_args()

    # Ensure we don't hit floating point errors
    np.seterr('raise')
    np.seterr(under='ignore')
    # We are given an arg with inter-gene spacing. We will run at 10 different upstream expression levels
    accum = []
    run_id = 0
    for geometry in ['tandem', 'convergent', 'divergent']:
        for expression in np.logspace(-2,0,10):
            for dx in [400, 500, 750, 1000, 1500, 2000, 2500, 3000]:
                for topo in [1, 2, 5, 10]:
                    if (run_id % args.num_runners) == args.runner_id:
                        accum.append(bulk_topology_simulation(geometry, dx, topo, (1, expression), 12000, args.n_simulations))
                        print(f'Finished {geometry}/expression:{expression}/dx:{dx}/topo_multiplier:{topo}')
                    run_id += 1
    df = pd.DataFrame(accum, columns=['geometry', 'spacing', 'topo_rate', 'A_promoter_strength', 'B_promoter_strength', 'A_mean', 'A_std', 'B_mean', 'B_std'])
    df.to_feather(f'output/dataframes/worker_{args.runner_id}.feather')