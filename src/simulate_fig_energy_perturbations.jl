using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/energy_torque_perturbations-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_summarized_runs_per_repeat = 15
n_sc_burst_runs_per_repeat = 5
n_repeats = 50
i = 0

function gen_sim_params(;
    σ2_coeff::Float64=0.0,
    torque_perturb::TanglesModel.TorqueFunctionPerturbation=NoTorqueFunctionPerturbation(),
    rnap_perturb::TanglesModel.RNAPInitPerturbation=NoRNAPInitPerturbation())
    return SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate * 0.0,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate,
        true,
        σ2_coeff,
        NoTopoisomerase(),
        torque_perturb,
        rnap_perturb
    )
end


# steady-state response and SC/burst dynamics of four constructs under two perturbative scenarios
for _ in 1:n_repeats,
    is_plasmid in [false, true],
    sc_initiation in [true],
    induction in vcat([0.0,1.0],exp10.(range(-2,0.5,length=30))),
    σ2 in [0.025],
    scenario in ["none", "buffering", "energy_well"]


    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    # Perform simulations with default hyperparameters, other than perturbations
    if scenario == "buffering"
        params = gen_sim_params(σ2_coeff = σ2, torque_perturb=PositiveSupercoilingBuffering(0.031))
    elseif scenario == "energy_well"
        params = gen_sim_params(σ2_coeff = σ2, rnap_perturb=RNAPInitEnergyWell(-0.06, 0.125))
    else
        params = gen_sim_params(σ2_coeff = σ2)
    end

    bcs = is_plasmid ? CircularBoundaryParameters(10000.0 * 0.34) : LinearBoundaryParameters(10000.0 * 0.34, false, false)
    tandem_reporter_up   = DiscreteConfig([UncoupledGene(base_rate * induction, 1, 6000 * 0.34, 7000 * 0.34), UncoupledGene(base_rate, 2, 3000 * 0.34, 4000 * 0.34)])
    tandem_reporter_down =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])
    convergent =  DiscreteConfig([UncoupledGene(base_rate * induction, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(base_rate, 2, 7000 * 0.34, 6000 * 0.34)])
    divergent =   DiscreteConfig([UncoupledGene(base_rate * induction, 1, 4000 * 0.34, 3000 * 0.34), UncoupledGene(base_rate, 2, 6000 * 0.34, 7000 * 0.34)])

    start_time = time()
    simulate_summarized_runs(filename, n_summarized_runs_per_repeat, "fig.perturb.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 15000.0)
    simulate_summarized_runs(filename, n_summarized_runs_per_repeat, "fig.perturb.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 15000.0)
    simulate_summarized_runs(filename, n_summarized_runs_per_repeat, "fig.perturb.convergent", params, bcs, convergent, 15000.0)
    simulate_summarized_runs(filename, n_summarized_runs_per_repeat, "fig.perturb.divergent", params, bcs, divergent, 15000.0)
    if (induction == 0.0) || (induction == 1.0)
        simulate_sc_rnap_dynamics(filename, n_sc_burst_runs_per_repeat, "fig.perturb.sc_bursts.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 200, 15000.0, 20, Dict{String,Float64}())
        simulate_sc_rnap_dynamics(filename, n_sc_burst_runs_per_repeat, "fig.perturb.sc_bursts.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 200, 15000.0, 20, Dict{String,Float64}())
        simulate_sc_rnap_dynamics(filename, n_sc_burst_runs_per_repeat, "fig.perturb.sc_bursts.convergent", params, bcs, convergent, 200, 15000.0, 20, Dict{String,Float64}())
        simulate_sc_rnap_dynamics(filename, n_sc_burst_runs_per_repeat, "fig.perturb.sc_bursts.divergent", params, bcs, divergent, 200, 15000.0, 20, Dict{String,Float64}())
    end

    println("Done with perturbation sims with params:\n\tis_plasmid: ", is_plasmid, "\n\tinduction: ", induction, "\n\tscenario: ", scenario)
    println("Ran round in ", time() - start_time, " seconds")
end

println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
