using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig5_sims-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 160.0
n_examples_per_node = 100
n_repeats = 100
i = 0

function gen_sim_params(;
    topo_rate_factor::Float64=1.0,
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

# Figure 5: simulate the toggle switch in the three orientations, while changing
# both the input hill coefficient and the size of the reservoir.
for _ in 1:n_repeats,
    hill_coeff in 1.0:0.1:5.0,
    k_factor in exp10.(range(-1.0,1.0,length=3)),
    mRNA_deg_rate_factor in exp10.(range(0.0,2.0,length=5)),
    topo_factor in exp10.(range(-0.5,2.5,length=21))

    # Distribute work between nodes
    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = gen_sim_params(
        sc_dependent = true,
        σ2_coeff = 0.02,
        topo_rate_factor = topo_factor,
        mRNA_deg_rate_factor = mRNA_deg_rate_factor)


    bcs = LinearBoundaryParameters(10000.0 * 0.34, false, false)

    k_val = 200.0 * k_factor / mRNA_deg_rate_factor
    basal_k = round(Int32, 15.0 / mRNA_deg_rate_factor)
    mRNA_ic = convert(Array{Int32,1}, [basal_k, 0])
    hill_func(inhibitor) = (k_val) / (k_val + inhibitor^hill_coeff)
    extra_attrs = Dict("hill_coeff"=>hill_coeff, "K_val"=>k_val, "K_factor"=>k_factor, "mRNA_deg_factor"=>mRNA_deg_rate_factor, "topo_factor"=>topo_factor)

    # Generate genes with a 10000-second "burn-in" period
    tandem =      DiscreteConfig([
        CoupledGene(base_rate, 1, 3500 * 0.34, 4500 * 0.34, (mRNA,_)->hill_func(mRNA[2])),
        CoupledGene(base_rate, 2, 5500 * 0.34, 6500 * 0.34, (mRNA,t)->(t > 10000) * hill_func(mRNA[1]))
    ])
    tandem_up =      DiscreteConfig([
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
    simulate_discrete_runs(filename, n_examples_per_node, "fig5.tandem", params, bcs, tandem, 25000.0, 500, mRNA_ic, extra_attrs)
    simulate_discrete_runs(filename, n_examples_per_node, "fig5.tandem_up", params, bcs, tandem_up, 25000.0, 500, mRNA_ic, extra_attrs)
    simulate_discrete_runs(filename, n_examples_per_node, "fig5.convergent", params, bcs, convergent, 25000.0, 500, mRNA_ic, extra_attrs)
    simulate_discrete_runs(filename, n_examples_per_node, "fig5.divergent", params, bcs, divergent, 25000.0, 500, mRNA_ic, extra_attrs)
    println("Done with fig 5 with params:\n\thill_coeff: ", hill_coeff, "\n\tK_factor: ", k_factor, "\n\tmRNA_deg_fac: ", mRNA_deg_rate_factor, "\n\ttopo_fac: ", topo_factor)
    println("Ran round in ", time() - start_time, " seconds")
end

println("Done with ALL fig5 simulations; node shutting down after ", time() - overall_start, " seconds!")
