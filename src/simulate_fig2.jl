using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig2_sims-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_examples_per_node = 100
n_full_examples_per_node = 10
n_repeats = 1
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


# Fig 2: simulate small number of examples of turning on genes in a single run 
for _ in 1:n_repeats,
    induction in exp10.(range(-2,0.5,length=30)),
    σ2 in [0.0, 0.02]

    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = gen_sim_params(sc_dependent=true, σ2_coeff = σ2)

    bcs = LinearBoundaryParameters(10000.0 * 0.34, false, false)

    step_time = 10000.0
    coupling_func(_mRNA,t) = (t > step_time) ? induction : 0.0
    tandem_down = DiscreteConfig([CoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34, (_,_)->1.0), CoupledGene(base_rate, 1, 6000 * 0.34, 7000 * 0.34, coupling_func)])
    tandem_up =   DiscreteConfig([CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)])
    convergent =  DiscreteConfig([CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34, (_,_)->1.0)])
    divergent =   DiscreteConfig([CoupledGene(base_rate, 1, 4000 * 0.34, 3000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)])

    start_time = time()
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig2.tandem_down", params, bcs, tandem_down, 20000.0, 1000, {"step_time" => step_time})
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig2.tandem_up", params, bcs, tandem_up, 20000.0, 1000, {"step_time" => step_time})
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig2.convergent", params, bcs, convergent, 20000.0, 1000, {"step_time" => step_time})
    simulate_discrete_runs(filename, n_full_examples_per_node, "fig2.divergent", params, bcs, divergent, 20000.0, 1000, {"step_time" => step_time})
    println("Done with fig2 with params:\n\tinduction: ", induction, "\n\t σ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
