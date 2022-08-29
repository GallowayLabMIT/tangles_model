using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig_topo_dependence-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_summaries_per_node = 20
n_burst_datasets_per_node = 5
n_repeats = 10
i = 0

function gen_sim_params(;
    topo_rate_factor::Float64=0.0,
    mRNA_deg_rate_factor::Float64=1.0,
    sc_dependent::Bool=DEFAULT_SIM_PARAMS.sc_dependent,
    σ2_coeff::Float64=0.0,
    topo_type::TanglesModel.TopoisomeraseType=NoTopoisomerase(),
    torque_perturb::TanglesModel.TorqueFunctionPerturbation=NoTorqueFunctionPerturbation(),
    rnap_perturb::TanglesModel.RNAPInitPerturbation=NoRNAPInitPerturbation())
    return SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate * topo_rate_factor,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate * mRNA_deg_rate_factor,
        sc_dependent,
        σ2_coeff,
        topo_type::TanglesModel.TopoisomeraseType,
        torque_perturb,
        rnap_perturb
    )
end

# steady-state response of three constructs with four different boundary conditions.
for _ in 1:n_repeats,
    sc_initiation in [true],
    is_plasmid in [false, true],
    induction in [0.0, 0.1, 0.33, 0.5, 1.0, 1.5, 2.0],
    σ2 in [0.025],
    topo in [NoTopoisomerase(), IntragenicTopoisomerase(), IntergenicTopoisomerase()],
    topo_multiplier = [0.0, 0.1, 1.0, 10.0]

    # Skip topo multiplier if we aren't using a topoisomerase
    if typeof(topo) == NoTopoisomerase
        if topo_multiplier != 0.0
            continue
        end
    end

    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    # Perform simulations with default hyperparameters
    params = gen_sim_params(sc_dependent=sc_initiation, σ2_coeff = σ2, topo_rate_factor=topo_multiplier, topo_type=topo)

    bcs = is_plasmid ? CircularBoundaryParameters(10000.0 * 0.34) : LinearBoundaryParameters(10000.0 * 0.34, false, false)
    tandem_reporter_up   = DiscreteConfig([UncoupledGene(base_rate * induction, 1, 6000 * 0.34, 7000 * 0.34), UncoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34)])
    tandem_reporter_down =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])
    convergent =  DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34)])
    divergent =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])

    start_time = time()
    simulate_sc_rnap_dynamics(filename, n_burst_datasets_per_node, "fig.topo.sc_bursts.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 200, 10000.0, 20, Dict{String,Float64}())
    simulate_sc_rnap_dynamics(filename, n_burst_datasets_per_node, "fig.topo.sc_bursts.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 200, 10000.0, 20, Dict{String,Float64}())
    simulate_sc_rnap_dynamics(filename, n_burst_datasets_per_node, "fig.topo.sc_bursts.convergent", params, bcs, convergent, 200, 10000.0, 20, Dict{String,Float64}())
    simulate_sc_rnap_dynamics(filename, n_burst_datasets_per_node, "fig.topo.sc_bursts.divergent", params, bcs, divergent, 200, 10000.0, 20, Dict{String,Float64}())
    simulate_summarized_runs(filename, n_summaries_per_node, "fig.topo.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 15000.0)
    simulate_summarized_runs(filename, n_summaries_per_node, "fig.topo.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 15000.0)
    simulate_summarized_runs(filename, n_summaries_per_node, "fig.topo.convergent", params, bcs, convergent, 15000.0)
    simulate_summarized_runs(filename, n_summaries_per_node, "fig.topo.divergent", params, bcs, divergent, 15000.0)
    println("Done with base model sc bursting sims with params:\n\tis_plasmid: ", is_plasmid, "\n\tsc_dependent: ", sc_initiation, "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end
