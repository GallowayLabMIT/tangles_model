import itertools
import stochastic_sc_model as ssc
from pathos.pools import ProcessPool

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
    'mRNA_deg_rate': 1 / 4000 # 1 / sec
}

def run_sweep(dx_in_bp, topo_multiplier, n_runs, ptype, id):
    params_copy = dict(params)
    params_copy['topo_rate'] *= topo_multiplier
    endcap_distances = 3000 * .34
    gene_distance = 2000 * .34
    dx = dx_in_bp * .34
    
    gene_A_boundaries = (endcap_distances, endcap_distances + gene_distance)
    gene_B_boundaries = (gene_A_boundaries[1] + dx, gene_A_boundaries[1] + dx + gene_distance)
    end_barrier = gene_B_boundaries[1] + endcap_distances

    genes = {
        'A': [gene_A_boundaries[0], gene_A_boundaries[1]],
        'B': [gene_B_boundaries[1], gene_B_boundaries[0]]
    }
    promoters = [
        ssc.Promoter(ptype, .8, gene_A_boundaries[0], 1),
        ssc.Promoter(ptype, .8, gene_B_boundaries[1], -1),
    ]
    pases = [
        ssc.PAS(ssc.PAS_Type.NO_RT, gene_A_boundaries[1], 1),
        ssc.PAS(ssc.PAS_Type.NO_RT, gene_B_boundaries[0], -1)
    ]

    bcs = [
        ssc.BoundaryCondition.fixed(0),
        ssc.BoundaryCondition.fixed(end_barrier)
    ]

    runner = ssc.SimulationRunner(genes, promoters, pases, bcs, (0, 12000))
    print(f'Starting run dx={dx_in_bp}/topo={topo_multiplier}/ptype={ptype}/id={id}')
    result_pd = runner.multi_simulation(n_runs)
    result_pd['dx'] = dx_in_bp
    result_pd['topo_rate'] = topo_multiplier
    result_pd['sc_dependent'] = ptype == ssc.PromoterType.SC_DEPENDENT
    out_file = f'output/dataframes/dx-{dx_in_bp}_topo-{topo_multiplier}_id-{id}_sc={ptype!=ssc.PromoterType.STATIC}.feather'
    print(f'Finished run dx={dx_in_bp}/topo={topo_multiplier}/ptype={ptype}/id={id}, saving to file{out_file}')
    result_pd.to_feather(out_file)

if __name__ == '__main__':
    pool = ProcessPool()
    pool.map(lambda args: run_sweep(args[0], args[1], 50, args[2], args[3]), itertools.product(
        [400, 500, 750, 1000, 1500, 2000, 2500, 3000],
        [1, 2, 3, 5, 7, 9, 10],
        [ssc.PromoterType.STATIC,ssc.PromoterType.SC_DEPENDENT],
        range(100)))
