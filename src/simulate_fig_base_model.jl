using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig_base_model_sims-node" * lpad(node_idx, 5, "0") * ".h5"

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


# steady-state response of three constructs with four different boundary conditions.
for _ in 1:n_repeats,
    is_plasmid in [false, true],
    sc_initiation in [false, true],
    induction in exp10.(range(-2,0.5,length=30)),
    σ2 in [0.0, 0.01, 0.02, 0.025, 0.029, 0.03, 0.031, 0.0316, 0.035, 0.05, 0.1]


    # If we aren't doing sc-dependent initiation, skip if σ2 is not zero
    if ~sc_initiation
        if σ2 != 0.0
            continue
        end
    end



    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    # Perform simulations with default hyperparameters
    params = gen_sim_params(sc_dependent=sc_initiation, σ2_coeff = σ2)

    bcs = is_plasmid ? CircularBoundaryParameters(10000.0 * 0.34) : LinearBoundaryParameters(10000.0 * 0.34, false, false)
    tandem_reporter_up   = DiscreteConfig([UncoupledGene(base_rate * induction, 1, 6000 * 0.34, 7000 * 0.34), UncoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34)])
    tandem_reporter_down =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])
    convergent =  DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34)])
    divergent =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])

    start_time = time()
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.convergent", params, bcs, convergent, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig.bm.divergent", params, bcs, divergent, 15000.0)
    println("Done with base model sims with params:\n\tis_plasmid: ", is_plasmid, "\n\tsc_dependent: ", sc_initiation, "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
