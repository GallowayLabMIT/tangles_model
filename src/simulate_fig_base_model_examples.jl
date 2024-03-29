using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig_bm_examples_sims-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_full_examples_per_node = 50
n_repeats = 5
i = 0

function gen_sim_params(;
    topo_rate_factor::Float64=0.0,
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
for _ in 1:n_repeats,
    induction in exp10.(range(-2,0.5,length=31)),
    σ2 in [0.02, 0.025, 0.03]

    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = gen_sim_params(sc_dependent=true, σ2_coeff = σ2)

    bcs = LinearBoundaryParameters(10000.0 * 0.34, false, false)

    step_time = 10000.0
    coupling_func(_mRNA,t) = (t > step_time) ? induction : 0.0
    tandem_reporter_up = DiscreteConfig([CoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34, (_,_)->1.0), CoupledGene(base_rate, 1, 6000 * 0.34, 7000 * 0.34, coupling_func)])
    tandem_reporter_down =   DiscreteConfig([CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)])
    convergent =  DiscreteConfig([CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34, (_,_)->1.0)])
    divergent =   DiscreteConfig([CoupledGene(base_rate, 1, 4000 * 0.34, 3000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)])

    start_time = time()
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig.bm_examples.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 40000.0, 4000, Dict("step_time" => step_time, "step_induction" => induction))
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig.bm_examples.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 40000.0, 4000, Dict("step_time" => step_time, "step_induction" => induction))
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig.bm_examples.convergent", params, bcs, convergent, 40000.0, 4000, Dict("step_time" => step_time, "step_induction" => induction))
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig.bm_examples.divergent", params, bcs, divergent, 40000.0, 4000, Dict("step_time" => step_time, "step_induction" => induction))
    println("Done with base-model examples with params:\n\tinduction: ", induction, "\n\t σ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
