using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig1_fig2_sims-node" * lpad(node_idx, 5, "0") * ".h5"

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
n_repeats = 50
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


# Fig 1c-e (+ supplemental fig): steady-state response of three constructs with four different boundary conditions.
for _ in 1:n_repeats,
    is_plasmid in [false, true],
    sc_initiation in [false, true],
    induction in exp10.(range(-2,0.5,length=30)),
    σ2 in vcat([0.0], exp10.(range(-3.0,-1.0,length=9)))


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
    tandem_down = [UncoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate * induction, 1, 6000 * 0.34, 7000 * 0.34)]
    tandem_up =   [UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)]
    convergent =  [UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34)]
    divergent =   [UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)]

    start_time = time()
    simulate_summarized_runs(filename, n_examples_per_node, "fig1.tandem_downstream", params, bcs, tandem_down, 2, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig1.tandem_upstream", params, bcs, tandem_up, 2, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig1.convergent", params, bcs, convergent, 2, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig1.divergent", params, bcs, divergent, 2, 15000.0)
    println("Done with fig 1c-d-e with params:\n\tis_plasmid: ", is_plasmid, "\n\tsc_dependent: ", sc_initiation, "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

# Simulate 1f: sweeping the intergenic spacing
for _ in 1:n_repeats,
    is_plasmid in [false, true],
    induction in exp10.(range(-2,0.5,length=30)),
    σ2 in [0.0, 0.02],
    spacing in [500, 1000, 1500, 2000, 3000, 4000, 5000, 10000]

    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = gen_sim_params(sc_dependent=true, σ2_coeff = σ2)

    bcs = is_plasmid ? CircularBoundaryParameters((8000.0 + spacing) * 0.34) : LinearBoundaryParameters((8000.0 + spacing) * 0.34, false, false)
    tandem_down = [UncoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate * induction, 1, (4000 + spacing) * 0.34, (5000 + spacing) * 0.34)]
    tandem_up =   [UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, (4000 + spacing) * 0.34, (5000 + spacing) * 0.34)]
    convergent =  [UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, (5000 + spacing) * 0.34, (4000 + spacing) * 0.34)]
    divergent =   [UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, (4000 + spacing) * 0.34, (5000 + spacing) * 0.34)]

    start_time = time()
    simulate_summarized_runs(filename, n_examples_per_node, "fig1f.spacing.tandem_downstream", params, bcs, tandem_down, 2, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig1f.spacing.tandem_upstream", params, bcs, tandem_up, 2, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig1f.spacing.convergent", params, bcs, convergent, 2, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "fig1f.spacing.divergent", params, bcs, divergent, 2, 15000.0)
    println("Done with fig 1f with params:\n\tis_plasmid: ", is_plasmid, "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

# Simulate changing hyperparameter sensitivity
for _ in 1:n_repeats,
    drag_coeff in [1/100, 1/50, 1/20, 1/10, 1/5, 1/2],
    drag_exp in [0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 2.5, 3.0],
    stall_torque in [8, 10, 12, 14, 16],
    stall_width in [1, 3, 5],
    σ2 in [0.0, 0.02],
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
    convergent =  [UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34)]
    divergent =   [UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)]

    start_time = time()
    simulate_summarized_runs(filename, n_examples_per_node, "hyperparams.convergent", params, bcs, convergent, 2, 15000.0)
    simulate_summarized_runs(filename, n_examples_per_node, "hyperparams.divergent", params, bcs, divergent, 2, 15000.0)
    println("Done with hyperparam sweep with params:\n\tdrag_coeff: ", drag_coeff,
        "\n\tdrag_exp: ", drag_exp, "\n\tstall_torque: ", stall_torque, "\n\tstall_width: ", stall_width,
        "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

# Fig 2: simulate small number of examples of turning on genes in a single run 
for _ in 1:n_repeats,
    induction in exp10.(range(-2,0.5,length=30)),
    σ2 in [0.0, 0.02]

    params = gen_sim_params(sc_dependent=true, σ2_coeff = σ2)

    bcs = LinearBoundaryParameters(10000.0 * 0.34, false, false)

    coupling_func(_mRNA,t) = (t > 10000.0) ? induction : 0.0
    tandem_down = [CoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34, (_,_)->1.0), CoupledGene(base_rate, 1, 6000 * 0.34, 7000 * 0.34, coupling_func)]
    tandem_up =   [CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)]
    convergent =  [CoupledGene(base_rate, 1, 3000 * 0.34, 4000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34, (_,_)->1.0)]
    divergent =   [CoupledGene(base_rate, 1, 4000 * 0.34, 3000 * 0.34, coupling_func), CoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34, (_,_)->1.0)]

    start_time = time()
    simulate_full_examples(filename, n_full_examples_per_node, "fig2.tandem_down", params, bcs, tandem_down, 2, 20000.0)
    simulate_full_examples(filename, n_full_examples_per_node, "fig2.tandem_up", params, bcs, tandem_up, 2, 20000.0)
    simulate_full_examples(filename, n_full_examples_per_node, "fig2.convergent", params, bcs, convergent, 2, 20000.0)
    simulate_full_examples(filename, n_full_examples_per_node, "fig2.divergent", params, bcs, divergent, 2, 20000.0)
    println("Done with fig2 with params:\n\tinduction: ", induction, "\n\t σ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")