import stochastic_sc_model as ssc
from pathos.pools import ProcessPool

if __name__ == '__main__':
    pool = ProcessPool(nodes=4)

    genes = {"RFP": [2000, 2600]}
    promoter = [ssc.Promoter(ssc.PromoterType.STATIC, 1, 2000, 1)]
    pas = [ssc.PAS(ssc.PAS_Type.NO_RT, 2600, 1, 10)]
    bcs = (ssc.BoundaryCondition.fixed(0), ssc.BoundaryCondition.fixed(4000))

    runner = ssc.SimulationRunner(genes, promoter, pas, bcs, (0, 10000))
    print(runner.multi_simulation(50))