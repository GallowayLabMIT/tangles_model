using TanglesModel
node_idx = parse(Int64, ARGS[1])
n_nodes = parse(Int64, ARGS[2])
filename = "output/modeling_paper/fig_zinani_topo_sims-node" * lpad(node_idx, 5, "0") * ".h5"

println("""
======================================================
Run summary:
""", "\tprocess_idx: ", node_idx, "\n\ttotal number of processes: ", n_nodes,"\n",
"""
======================================================""")

overall_start = time()

base_rate = 1.0 / 120.0
n_examples_per_node = 10
n_repeats = 100
i = 0

function gen_sim_params(;
    temperature::Float64=273.15+37.0,
    topo_rate_factor::Float64=1.0,
    mRNA_deg_rate::Float64=1.0,
    sc_dependent::Bool=DEFAULT_SIM_PARAMS.sc_dependent,
    σ2_coeff::Float64=0.0,
    topo_type::TanglesModel.TopoisomeraseType)
    return SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate * topo_rate_factor,
        mRNA_deg_rate,
        sc_dependent,
        σ2_coeff,
        topo_type
    )
end

# Figure 4:

# 6017 bp before next gene to the left, then her1
# Her1 is 6405 bp long
# They are separated by 12107 bp
# Her7 is 1315 bp long
# Her7 is 8549 bp before the next gene to the right
her1_space = 6017 * 0.34
her1_len   = 6405 * 0.34
her1_her7_space = 12107 * 0.34
her7_len   = 1315 * 0.34
her7_space = 8549 * 0.34

# Figure 3: simulate the her1/her7 system, in the three following conditions:
#   1. Uncoupled: separate her1/her7 to be very distant from each other
#   2. Fully coupled: replicate their results by tying her1/her7 mRNA production precisely together (just use 1).
#   3. TANGLES coupled: use the real topology.
base_zinani_rates = Dict(
    "mRNA_synthesis" => 5, "protein_synthesis" => 5,
    "mRNA_degradation" => 0.2, "protein_degradation" =>0.3,
    "dimer_association" => 0.01, "dimer_dissociation" => 400)
# Rescale Zinani rates to per-second instead of per-minute
time_rescaled_rates = Dict(key => val / 60.0 for (key, val) in base_zinani_rates)
print("Time rescaled rates:\n\t")
println(time_rescaled_rates)
# We have the following discrete components:
# 1:  her1 mRNA
# 2:  her7 mRNA
# 3:  her1 protein
# 4:  her7 protein
# 5:  hes6 protein (fixed at value of 100 in Zinani)
# 6:  unoccupied_her1_promoter (binary!)
# 7:  unoccupied_her7_promoter (binary!)
# 8:  Her1/Her1_occupied_her1_promoter (binary!)
# 9:  Her7/Hes6_occupied_her1_promoter (binary!)
# 10: Her1/Her1_occupied_her7_promoter (binary!)
# 11: Her7/Hes6_occupied_her7_promoter (binary!)
dmap = Dict(
    "her1_mRNA" => 1,
    "her7_mRNA" => 2,
    "her1_protein" => 3,
    "her7_protein" => 4,
    "hes6_protein" => 5,
    "her1_promoter_empty" => 6,
    "her7_promoter_empty" => 7,
    "her1_promoter_with_11" => 8,
    "her1_promoter_with_76" => 9,
    "her7_promoter_with_11" => 10,
    "her7_promoter_with_76" => 11
)
discrete_ic = convert(Array{Int32,1}, [0, 0, 0, 0, 100, 1, 1, 0, 0, 0, 0])


for _ in 1:n_repeats,
    temperature in [273.15 + 21.5],
    σ2 in [0.02, 0.025, 0.03],
    state in ["uncoupled", "fully-coupled", "tangles-coupled"],
    (topo, topo_str) in zip(
        [NoTopoisomerase(), OriginalTopoisomerase(), IntragenicTopoisomerase(), IntergenicTopoisomerase()],
        ["none", "original", "intragenic", "intergenic"]
    )

    # Distribute work between nodes
    if i % n_nodes != node_idx
        global i += 1
        continue
    end
    global i += 1

    params = gen_sim_params(
        temperature = temperature,
        sc_dependent = true,
        σ2_coeff = σ2,
        mRNA_deg_rate = time_rescaled_rates["mRNA_degradation"],
        topo_type = topo)

    if state == "uncoupled"
        # Free-end BCs, with 10 million basepairs total
        endpoint = 1000000000.0 * 0.34
        bcs = LinearBoundaryParameters(endpoint, true, true)
        config = DiscreteConfig([
            # Genes are active when they have an unoccupied promoter site
            CoupledGene(time_rescaled_rates["mRNA_synthesis"], 1, her1_space + her1_len, her1_space,
                (discrete,_)->Float64(discrete[dmap["her1_promoter_empty"]])),
            CoupledGene(time_rescaled_rates["mRNA_synthesis"], 2, endpoint - (her7_space + her7_len), endpoint - her7_space,
                (discrete,_)->Float64(discrete[dmap["her7_promoter_empty"]])),
        ], 9, [
            ((discrete,_)->time_rescaled_rates["protein_synthesis"] * discrete[dmap["her1_mRNA"]]) => [dmap["her1_protein"] => 1],
            ((discrete,_)->time_rescaled_rates["protein_synthesis"] * discrete[dmap["her7_mRNA"]]) => [dmap["her7_protein"] => 1],
            ((discrete,_)->time_rescaled_rates["protein_degradation"] * discrete[dmap["her1_protein"]]) => [dmap["her1_protein"] => -1],
            ((discrete,_)->time_rescaled_rates["protein_degradation"] * discrete[dmap["her7_protein"]]) => [dmap["her7_protein"] => -1],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her1_promoter_empty"]] * discrete[dmap["her1_protein"]] * (discrete[dmap["her1_protein"]] - 1) / 2) => [
                dmap["her1_promoter_empty"] => -1, dmap["her1_protein"] => -2, dmap["her1_promoter_with_11"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her7_promoter_empty"]] * discrete[dmap["her1_protein"]] * (discrete[dmap["her1_protein"]] - 1) / 2) => [
                dmap["her7_promoter_empty"] => -1, dmap["her1_protein"] => -2, dmap["her7_promoter_with_11"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her1_promoter_empty"]] * discrete[dmap["hes6_protein"]] * discrete[dmap["her7_protein"]]) => [
                dmap["her1_promoter_empty"] => -1, dmap["hes6_protein"] => -1, dmap["her7_protein"] => -1, dmap["her1_promoter_with_76"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her7_promoter_empty"]] * discrete[dmap["hes6_protein"]] * discrete[dmap["her7_protein"]]) => [
                dmap["her7_promoter_empty"] => -1, dmap["hes6_protein"] => -1, dmap["her7_protein"] => -1, dmap["her7_promoter_with_76"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her1_promoter_with_11"]]) => [
                dmap["her1_promoter_with_11"] => -1, dmap["her1_promoter_empty"] => 1, dmap["her1_protein"] => 2
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her7_promoter_with_11"]]) => [
                dmap["her7_promoter_with_11"] => -1, dmap["her7_promoter_empty"] => 1, dmap["her1_protein"] => 2
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her1_promoter_with_76"]]) => [
                dmap["her1_promoter_with_76"] => -1, dmap["her1_promoter_empty"] => 1, dmap["hes6_protein"] => 1, dmap["her7_protein"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her7_promoter_with_76"]]) => [
                dmap["her7_promoter_with_76"] => -1, dmap["her7_promoter_empty"] => 1, dmap["hes6_protein"] => 1, dmap["her7_protein"] => 1
            ]
        ])
        start_time = time()
        simulate_discrete_runs(filename, n_examples_per_node, "fig.zinani.uncoupled", params, bcs, config, 15000.0, 1000, discrete_ic, Dict{String,Float64}("temperature"=>temperature), Dict{String,String}("topo"=>topo_str))
        println("Done with Zinani fig with params:\n\ttemperature: ", temperature, "\n\ttype: ", state)
        println("Ran round in ", time() - start_time, " seconds")
    elseif state == "fully-coupled"
        # Free-end BCs, with 10 million basepairs total
        endpoint = 1000000000.0 * 0.34
        bcs = LinearBoundaryParameters(endpoint, true, true)
        config = DiscreteConfig([
            # Both genes are active, but only one of the mRNAs makes both proteins
            MultiCoupledGene(time_rescaled_rates["mRNA_synthesis"], Array{UInt32,1}([1,2]), her1_space + her1_len, her1_space,
                (discrete,_)->Float64(discrete[dmap["her1_promoter_empty"]])),
            MultiCoupledGene(time_rescaled_rates["mRNA_synthesis"], Array{UInt32,1}(), endpoint - (her7_space + her7_len), endpoint - her7_space,
                (discrete,_)->Float64(discrete[dmap["her7_promoter_empty"]])),
        ], 9, [
            ((discrete,_)->time_rescaled_rates["protein_synthesis"] * discrete[dmap["her1_mRNA"]]) => [dmap["her1_protein"] => 1],
            ((discrete,_)->time_rescaled_rates["protein_synthesis"] * discrete[dmap["her7_mRNA"]]) => [dmap["her7_protein"] => 1],
            ((discrete,_)->time_rescaled_rates["protein_degradation"] * discrete[dmap["her1_protein"]]) => [dmap["her1_protein"] => -1],
            ((discrete,_)->time_rescaled_rates["protein_degradation"] * discrete[dmap["her7_protein"]]) => [dmap["her7_protein"] => -1],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her1_promoter_empty"]] * discrete[dmap["her1_protein"]] * (discrete[dmap["her1_protein"]] - 1) / 2) => [
                dmap["her1_promoter_empty"] => -1, dmap["her1_protein"] => -2, dmap["her1_promoter_with_11"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her7_promoter_empty"]] * discrete[dmap["her1_protein"]] * (discrete[dmap["her1_protein"]] - 1) / 2) => [
                dmap["her7_promoter_empty"] => -1, dmap["her1_protein"] => -2, dmap["her7_promoter_with_11"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her1_promoter_empty"]] * discrete[dmap["hes6_protein"]] * discrete[dmap["her7_protein"]]) => [
                dmap["her1_promoter_empty"] => -1, dmap["hes6_protein"] => -1, dmap["her7_protein"] => -1, dmap["her1_promoter_with_76"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her7_promoter_empty"]] * discrete[dmap["hes6_protein"]] * discrete[dmap["her7_protein"]]) => [
                dmap["her7_promoter_empty"] => -1, dmap["hes6_protein"] => -1, dmap["her7_protein"] => -1, dmap["her7_promoter_with_76"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her1_promoter_with_11"]]) => [
                dmap["her1_promoter_with_11"] => -1, dmap["her1_promoter_empty"] => 1, dmap["her1_protein"] => 2
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her7_promoter_with_11"]]) => [
                dmap["her7_promoter_with_11"] => -1, dmap["her7_promoter_empty"] => 1, dmap["her1_protein"] => 2
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her1_promoter_with_76"]]) => [
                dmap["her1_promoter_with_76"] => -1, dmap["her1_promoter_empty"] => 1, dmap["hes6_protein"] => 1, dmap["her7_protein"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her7_promoter_with_76"]]) => [
                dmap["her7_promoter_with_76"] => -1, dmap["her7_promoter_empty"] => 1, dmap["hes6_protein"] => 1, dmap["her7_protein"] => 1
            ]
        ])
        start_time = time()
        simulate_discrete_runs(filename, n_examples_per_node, "fig.zinani.fully-coupled", params, bcs, config, 15000.0, 1000, discrete_ic, Dict{String,Float64}("temperature"=>temperature), Dict{String,String}("topo"=>topo_str))
        println("Done with Zinani fig with params:\n\ttemperature: ", temperature, "\n\ttype: ", state)
        println("Ran round in ", time() - start_time, " seconds")
    elseif state == "tangles-coupled"
        endpoint = her1_space + her1_len + her1_her7_space + her7_len + her7_space
        bcs = LinearBoundaryParameters(endpoint, false, false)
        config = DiscreteConfig([
            # Both genes are active, but only one of the mRNAs makes both proteins
            CoupledGene(time_rescaled_rates["mRNA_synthesis"], 1, her1_space + her1_len, her1_space,
                (discrete,_)->Float64(discrete[dmap["her1_promoter_empty"]])),
            CoupledGene(time_rescaled_rates["mRNA_synthesis"], 2, endpoint - (her7_space + her7_len), endpoint - her7_space,
                (discrete,_)->Float64(discrete[dmap["her7_promoter_empty"]])),
        ], 9, [
            ((discrete,_)->time_rescaled_rates["protein_synthesis"] * discrete[dmap["her1_mRNA"]]) => [dmap["her1_protein"] => 1],
            ((discrete,_)->time_rescaled_rates["protein_synthesis"] * discrete[dmap["her7_mRNA"]]) => [dmap["her7_protein"] => 1],
            ((discrete,_)->time_rescaled_rates["protein_degradation"] * discrete[dmap["her1_protein"]]) => [dmap["her1_protein"] => -1],
            ((discrete,_)->time_rescaled_rates["protein_degradation"] * discrete[dmap["her7_protein"]]) => [dmap["her7_protein"] => -1],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her1_promoter_empty"]] * discrete[dmap["her1_protein"]] * (discrete[dmap["her1_protein"]] - 1) / 2) => [
                dmap["her1_promoter_empty"] => -1, dmap["her1_protein"] => -2, dmap["her1_promoter_with_11"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her7_promoter_empty"]] * discrete[dmap["her1_protein"]] * (discrete[dmap["her1_protein"]] - 1) / 2) => [
                dmap["her7_promoter_empty"] => -1, dmap["her1_protein"] => -2, dmap["her7_promoter_with_11"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her1_promoter_empty"]] * discrete[dmap["hes6_protein"]] * discrete[dmap["her7_protein"]]) => [
                dmap["her1_promoter_empty"] => -1, dmap["hes6_protein"] => -1, dmap["her7_protein"] => -1, dmap["her1_promoter_with_76"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_association"] * discrete[dmap["her7_promoter_empty"]] * discrete[dmap["hes6_protein"]] * discrete[dmap["her7_protein"]]) => [
                dmap["her7_promoter_empty"] => -1, dmap["hes6_protein"] => -1, dmap["her7_protein"] => -1, dmap["her7_promoter_with_76"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her1_promoter_with_11"]]) => [
                dmap["her1_promoter_with_11"] => -1, dmap["her1_promoter_empty"] => 1, dmap["her1_protein"] => 2
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her7_promoter_with_11"]]) => [
                dmap["her7_promoter_with_11"] => -1, dmap["her7_promoter_empty"] => 1, dmap["her1_protein"] => 2
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her1_promoter_with_76"]]) => [
                dmap["her1_promoter_with_76"] => -1, dmap["her1_promoter_empty"] => 1, dmap["hes6_protein"] => 1, dmap["her7_protein"] => 1
            ],
            ((discrete,_)->time_rescaled_rates["dimer_dissociation"] * discrete[dmap["her7_promoter_with_76"]]) => [
                dmap["her7_promoter_with_76"] => -1, dmap["her7_promoter_empty"] => 1, dmap["hes6_protein"] => 1, dmap["her7_protein"] => 1
            ]
        ])
        start_time = time()
        simulate_discrete_runs(filename, n_examples_per_node, "fig.zinani.tangles-coupled", params, bcs, config, 15000.0, 1000, discrete_ic, Dict{String,Float64}("temperature"=>temperature), Dict{String,String}("topo"=>topo_str))
        println("Done with Zinani fig with params:\n\ttemperature: ", temperature, "\n\ttype: ", state)
        println("Ran round in ", time() - start_time, " seconds")
    end
end
println("Done with ALL simulations; node shutting down after ", time() - overall_start, " seconds!")
