using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig_hyperparams_sims-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_examples_per_node = 40
n_full_examples_per_node = 40
n_repeats = 50
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

# Simulate changing hyperparameter sensitivity
for _ in 1:n_repeats,
    drag_coeff in [1/50, 1/20, 1/10, 1/5],
    drag_exp in [0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 3.0],
    stall_torque in [8, 10, 12, 14, 16],
    stall_width in [1, 3, 5],
    σ2 in [0.02, 0.025, 0.03],
    induction in exp10.(range(-2,0.5,length=10))

    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = SimulationParameters(
        TanglesModel.mRNA_Parameters(drag_coeff, drag_exp),
        TanglesModel.RNAP_Parameters(15, 20, stall_torque, stall_width),
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate,
        true,
        σ2
    )
    bcs = LinearBoundaryParameters(10000.0 * 0.34, false, false)
    tandem_reporter_up   = DiscreteConfig([UncoupledGene(base_rate * induction, 1, 6000 * 0.34, 7000 * 0.34), UncoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34)])
    tandem_reporter_down =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])
    convergent =  DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34)])
    divergent =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])

    start_time = time()
    simulate_summarized_runs(filename, n_examples_per_node, "hyperparams.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "hyperparams.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "hyperparams.convergent", params, bcs, convergent, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "hyperparams.divergent", params, bcs, divergent, 15000.0)
    println("Done with hyperparam sweep with params:\n\tdrag_coeff: ", drag_coeff,
        "\n\tdrag_exp: ", drag_exp, "\n\tstall_torque: ", stall_torque, "\n\tstall_width: ", stall_width,
        "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
