using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig_base_model_spacing_sims-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_examples_per_node = 100
n_full_examples_per_node = 100
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


# Simulate 1f: sweeping the intergenic spacing
for _ in 1:n_repeats,
    is_plasmid in [false, true],
    induction in exp10.(range(-2,0.5,length=30)),
    σ2 in [0.0, 0.02, 0.025, 0.03],
    spacing in [500, 1000, 1500, 2000, 3000, 4000, 5000, 10000]

    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = gen_sim_params(sc_dependent=true, σ2_coeff = σ2)

    bcs = is_plasmid ? CircularBoundaryParameters((8000.0 + spacing) * 0.34) : LinearBoundaryParameters((8000.0 + spacing) * 0.34, false, false)
    tandem_reporter_up   = DiscreteConfig([UncoupledGene(base_rate * induction, 1, (4000 + spacing) * 0.34, (5000 + spacing) * 0.34), UncoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34)])
    tandem_reporter_down =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, (4000 + spacing) * 0.34, (5000 + spacing) * 0.34)])
    convergent =  DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, (5000 + spacing) * 0.34, (4000 + spacing) * 0.34)])
    divergent =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, (4000 + spacing) * 0.34, (5000 + spacing) * 0.34)])

    start_time = time()
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.spacing.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.spacing.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.spacing.convergent", params, bcs, convergent, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.spacing.divergent", params, bcs, divergent, 15000.0)
    println("Done with base-model spacing sweep with params:\n\tis_plasmid: ", is_plasmid, "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end
println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
