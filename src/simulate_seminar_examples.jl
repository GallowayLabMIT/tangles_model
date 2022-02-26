using TanglesModel
filename = "output/seminar_examples.h5"

overall_start = time()

base_rate = 1.0 / 120.0
n_repeats = 3
i = 0

function gen_sim_params(;
    topo_rate_factor::Float64=1.0,
    sc_dependent::Bool=DEFAULT_SIM_PARAMS.sc_dependent,
    σ2_coeff::Float64=0.0)
    return SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate * topo_rate_factor,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate,
        sc_dependent,
        σ2_coeff
    )
end

# simulate small number of examples of turning on genes in a single run
σ2 = 0.025
induction = 1.0
params = gen_sim_params(sc_dependent=true, σ2_coeff = σ2)

bcs = LinearBoundaryParameters(10000.0 * 0.34, false, false)

step_time = 10000.0
coupling_func(_mRNA,t) = (t > step_time) ? induction : 0.0
tandem_reporter_up = DiscreteConfig([CoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34, (_,_)->1.0), CoupledGene(base_rate, 1, 6000 * 0.34, 7000 * 0.34, coupling_func)])
tandem_reporter_down =   DiscreteConfig([CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)])
convergent =  DiscreteConfig([CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34, (_,_)->1.0)])
divergent =   DiscreteConfig([CoupledGene(base_rate, 1, 4000 * 0.34, 3000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)])

start_time = time()
simulate_full_examples(filename, n_repeats, "fig.bm_examples.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 40000.0, 5000)
simulate_full_examples(filename, n_repeats, "fig.bm_examples.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 40000.0, 5000)
simulate_full_examples(filename, n_repeats, "fig.bm_examples.convergent", params, bcs, convergent, 40000.0, 5000)
simulate_full_examples(filename, n_repeats, "fig.bm_examples.divergent", params, bcs, divergent, 40000.0, 5000)
println("Done with base-model examples with params:\n\tinduction: ", induction, "\n\t σ2: ", σ2)
println("Ran round in ", time() - start_time, " seconds")

println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
