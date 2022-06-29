using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig_sc_bursting-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_repeats_per_node = 10
n_repeats = 10
i = 0

function gen_sim_params(;
    topo_rate_factor::Float64=0.0,
    mRNA_deg_rate_factor::Float64=1.0,
    sc_dependent::Bool=DEFAULT_SIM_PARAMS.sc_dependent,
    σ2_coeff::Float64=0.0)
    return SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate * topo_rate_factor,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate * mRNA_deg_rate_factor,
        sc_dependent,
        σ2_coeff
    )
end

# steady-state response of three constructs with four different boundary conditions.
for _ in 1:n_repeats,
    is_plasmid in [false, true],
    sc_initiation in [false, true],
    induction in [0.0, 0.1, 0.5, 1.0, 1.5, 2.0],
    σ2 in [0.0, 0.02, 0.025, 0.03]


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
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm.sc_bursts.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 200, 15000.0, 20, Dict{String,Float64}())
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm.sc_bursts.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 200, 15000.0, 20, Dict{String,Float64}())
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm.sc_bursts.convergent", params, bcs, convergent, 200, 15000.0, 20, Dict{String,Float64}())
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm.sc_bursts.divergent", params, bcs, divergent, 200, 15000.0, 20, Dict{String,Float64}())
    println("Done with base model sc bursting sims with params:\n\tis_plasmid: ", is_plasmid, "\n\tsc_dependent: ", sc_initiation, "\n\tinduction: ", induction, "\n\tσ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end

# Simulate turn-on dynamics
for _ in 1:n_repeats,
    induction in [0.0, 0.1, 0.5, 1.0, 1.5, 2.0],
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
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm_examples.sc_bursts.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 200, 40000.0, 20, Dict("step_time" => step_time, "step_induction" => induction))
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm_examples.sc_bursts.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 200, 40000.0, 20, Dict("step_time" => step_time, "step_induction" => induction))
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm_examples.sc_bursts.convergent", params, bcs, convergent, 200, 40000.0, 20, Dict("step_time" => step_time, "step_induction" => induction))
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.bm_examples.sc_bursts.divergent", params, bcs, divergent, 200, 40000.0, 20, Dict("step_time" => step_time, "step_induction" => induction))
    println("Done with base-model examples (sc bursting) with params:\n\tinduction: ", induction, "\n\t σ2: ", σ2)
    println("Ran round in ", time() - start_time, " seconds")
end
# Simulate toggles
for _ in 1:n_repeats,
    hill_coeff in 1.0:0.25:5.0,
    k_factor in exp10.(range(-1.0,1.0,length=3)),
    mRNA_deg_rate_factor in exp10.(range(0.0,2.0,length=5)),
    σ2 in [0.02, 0.025, 0.03]

    # Distribute work between nodes
    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = gen_sim_params(
        sc_dependent = true,
        σ2_coeff = σ2,
        mRNA_deg_rate_factor = mRNA_deg_rate_factor)


    bcs = LinearBoundaryParameters(10000.0 * 0.34, false, false)

    k_val = 200.0 * k_factor / mRNA_deg_rate_factor
    basal_k = round(Int32, 15.0 / mRNA_deg_rate_factor)
    mRNA_ic = convert(Array{Int32,1}, [basal_k, 0])
    hill_func(inhibitor) = (k_val) / (k_val + inhibitor^hill_coeff)
    extra_attrs = Dict("hill_coeff"=>hill_coeff, "K_val"=>k_val, "K_factor"=>k_factor, "mRNA_deg_factor"=>mRNA_deg_rate_factor)

    # Generate genes with a 10000-second "burn-in" period
    tandem_reporter_up =      DiscreteConfig([
        CoupledGene(base_rate, 1, 3500 * 0.34, 4500 * 0.34, (mRNA,_)->hill_func(mRNA[2])),
        CoupledGene(base_rate, 2, 5500 * 0.34, 6500 * 0.34, (mRNA,t)->(t > 10000) * hill_func(mRNA[1]))
    ])
    tandem_reporter_down =      DiscreteConfig([
        CoupledGene(base_rate, 1, 5500 * 0.34, 6500 * 0.34, (mRNA,_)->hill_func(mRNA[2])),
        CoupledGene(base_rate, 2, 3500 * 0.34, 4500 * 0.34, (mRNA,t)->(t > 10000) * hill_func(mRNA[1]))
    ])
    convergent =  DiscreteConfig([
        CoupledGene(base_rate, 1, 3500 * 0.34, 4500 * 0.34, (mRNA,_)->hill_func(mRNA[2])),
        CoupledGene(base_rate, 2, 6500 * 0.34, 5500 * 0.34, (mRNA,t)->(t > 10000) * hill_func(mRNA[1]))
    ])
    divergent =   DiscreteConfig([
        CoupledGene(base_rate, 1, 4500 * 0.34, 3500 * 0.34, (mRNA,_)->hill_func(mRNA[2])),
        CoupledGene(base_rate, 2, 5500 * 0.34, 6500 * 0.34, (mRNA,t)->(t > 10000) * hill_func(mRNA[1]))
    ])


    start_time = time()
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.toggles.sc_bursts.tandem_reporter_upstream", params, bcs, tandem_reporter_up, 200, 20000.0, 10, mRNA_ic, extra_attrs)
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.toggles.sc_bursts.tandem_reporter_downstream", params, bcs, tandem_reporter_down, 200, 20000.0, 10, mRNA_ic, extra_attrs)
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.toggles.sc_bursts.convergent", params, bcs, convergent, 200, 20000.0, 10, mRNA_ic, extra_attrs)
    simulate_sc_rnap_dynamics(filename, n_repeats_per_node, "fig.toggles.sc_bursts.divergent", params, bcs, divergent, 200, 20000.0, 10, mRNA_ic, extra_attrs)
    println("Done with toggles (sc-bursts) with params:\n\thill_coeff: ", hill_coeff, "\n\tK_factor: ", k_factor, "\n\tmRNA_deg_fac: ", mRNA_deg_rate_factor)
    println("Ran round in ", time() - start_time, " seconds")
end
println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
