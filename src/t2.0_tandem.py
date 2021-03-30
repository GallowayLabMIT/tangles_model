import itertools
from pathos.pools import ProcessPool
import sys
from modules import stochastic_sc_model as ssc

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

def run_sweep(topo_multiplier, induction, n_runs, orientation, rTTa_present, ptype, id):
    params_copy = dict(params)
    params_copy['topo_rate'] *= topo_multiplier
    endcap_distances = 3000 * .34
    fp_distance = 720 * .34
    spacer_distance = 649 * .34
    rTTa_distance = 904 * .34
    rTTa_spacer = 205 * .34

    EGFP_boundaries = (0, fp_distance)
    mRuby2_boundaries = (fp_distance + spacer_distance, fp_distance * 2 + spacer_distance)
    rTTa_boundaries = (-(rTTa_distance + rTTa_spacer) * .34, -(rTTa_spacer) * .34)

    genes = {}
    promoters = []
    pases = []

    if orientation == 'divergent':
        genes['EGFP'] = [EGFP_boundaries[1], EGFP_boundaries[0]]
        promoters.append(ssc.Promoter(ptype, induction, EGFP_boundaries[1], -1))
        pases.append(ssc.PAS(ssc.PAS_Type.NO_RT, EGFP_boundaries[0], -1))
    else:
        genes['EGFP'] = [EGFP_boundaries[0], EGFP_boundaries[1]]
        promoters.append(ssc.Promoter(ptype, induction, EGFP_boundaries[0], 1))
        pases.append(ssc.PAS(ssc.PAS_Type.NO_RT, EGFP_boundaries[1], 1))

    if orientation == 'convergent':
        genes['mRuby2'] = [mRuby2_boundaries[1], mRuby2_boundaries[0]]
        promoters.append(ssc.Promoter(ptype, .8, mRuby2_boundaries[1], -1))
        pases.append(ssc.PAS(ssc.PAS_Type.NO_RT, mRuby2_boundaries[0], -1))
    else:
        genes['mRuby2'] = [mRuby2_boundaries[0], mRuby2_boundaries[1]]
        promoters.append(ssc.Promoter(ptype, .8, mRuby2_boundaries[0], 1))
        pases.append(ssc.PAS(ssc.PAS_Type.NO_RT, mRuby2_boundaries[1], 1))

    if rTTa_present == True:
        if orientation == 'divergent':
            genes['rTTa'] = [rTTa_boundaries[1], rTTa_boundaries[0]]
            promoters.append(ssc.Promoter(ptype, .5, rTTa_boundaries[1], -1))
            pases.append(ssc.PAS(ssc.PAS_Type.NO_RT, rTTa_boundaries[0], -1))
        else:
            genes['rTTa'] = [rTTa_boundaries[0], rTTa_boundaries[1]]
            promoters.append(ssc.Promoter(ptype, .5, rTTa_boundaries[0], 1))
            pases.append(ssc.PAS(ssc.PAS_Type.NO_RT, rTTa_boundaries[1], 1))
        bcs = [
                ssc.BoundaryCondition.fixed(rTTa_boundaries[0] - endcap_distances),
                ssc.BoundaryCondition.fixed(mRuby2_boundaries[1] + endcap_distances)
            ]
    else:
        bcs = [
                ssc.BoundaryCondition.fixed(-endcap_distances),
                ssc.BoundaryCondition.fixed(mRuby2_boundaries[1] + endcap_distances)
            ]
    
    runner = ssc.SimulationRunner(genes, promoters, pases, bcs, (0, 12000))
    print(f'Starting run orientation={orientation}/rTTa={rTTa_present}/induction={induction}/topo={topo_multiplier}/ptype={ptype}/id={id}')
    result_pd = runner.multi_simulation(n_runs)
    result_pd['orientation'] = orientation
    result_pd['rTTa_present'] = rTTa_present
    result_pd['EGFP_induction'] = induction
    result_pd['topo_rate'] = topo_multiplier
    result_pd['sc_dependent'] = ptype == ssc.PromoterType.SC_DEPENDENT
    out_file = f'output/dataframes/t2.0_rTTa-{rTTa_present}_orientation-{orientation}_topo-{topo_multiplier}_id-{id}_sc={ptype!=ssc.PromoterType.STATIC}.feather'
    print(f'Finished run orientation={orientation}/rTTa={rTTa_present}/induction={induction}/topo={topo_multiplier}/ptype={ptype}/id={id}, saving to file{out_file}')
    result_pd.to_feather(out_file)

if __name__ == '__main__':
    pool = ProcessPool()
    pool.map(lambda args: run_sweep(args[0], args[1], 50, args[2], args[3], args[4], args[5]), itertools.product(
        [1],
        [0, .1, .2, .3, .4, .5, .6, .7, .8],
        ['tandem'],
        [False],
        [ssc.PromoterType.STATIC,ssc.PromoterType.SC_DEPENDENT],
        range(100)))
