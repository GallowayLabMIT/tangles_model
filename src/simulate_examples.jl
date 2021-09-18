using TanglesModel

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

base_rate = 1.0 / 120.0
n_examples = 5
for topo_multiplier in [0.1, 0.3, 0.5, 1.0],
    is_plasmid in [false, true],
    sc_dependent in [true],
    induction in exp10.(range(-2,0.5,length=5)),
    σ2 in vcat([0.0], exp10.(range(-1.5,-0.5,length=9)))

    start_time = time()

    bcs_2_gene = is_plasmid ? CircularBoundaryParameters(9616.0 * 0.34) : LinearBoundaryParameters(9616.0 * 0.34, false, false)
    genes_tandem_2 = [Gene(base_rate * induction, 1, 2329.0 * 0.34, 3372.0 * 0.34), Gene(base_rate, 2, 5067.0 * 0.34, 5810.0 * 0.34)]
    genes_convergent_2 = [Gene(base_rate * induction, 1, 2329.0 * 0.34, 3372.0 * 0.34), Gene(base_rate, 2, 4758.0 * 0.34, 4050.0 * 0.34)]
    genes_divergent_2 = [Gene(base_rate * induction, 1, 3614.0 * 0.34, 2555.0 * 0.34), Gene(base_rate, 2, 5067.0 * 0.34, 5810.0 * 0.34)]

    simulate_full_examples("output/julia_s2_example_sims.h5", n_examples, "2_gene.tandem",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent, σ2_coeff=σ2), bcs_2_gene,
            genes_tandem_2, 2, 15000.0)
    println("Done with 2_gene.tandem with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction, "\n\t σ2: ", σ2)
    simulate_full_examples("output/julia_s2_example_sims.h5", n_examples, "2_gene.convergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent, σ2_coeff=σ2), bcs_2_gene,
            genes_convergent_2, 2, 15000.0)
    println("Done with 2_gene.convergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction, "\n\t σ2: ", σ2)
    simulate_full_examples("output/julia_s2_example_sims.h5", n_examples, "2_gene.divergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent, σ2_coeff=σ2), bcs_2_gene,
            genes_divergent_2, 2, 15000.0)
    println("Done with 2_gene.convergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction, "\n\t σ2: ", σ2)

    bcs_3_gene = is_plasmid ? CircularBoundaryParameters(11238.0 * 0.34) : LinearBoundaryParameters(11238.0 * 0.34, false, false)
    genes_tandem_3 = [Gene(base_rate, 3, 2620.0 * 0.34, 3723.0 * 0.34), Gene(base_rate * induction, 1, 4225.0 * 0.34, 5003.0 * 0.34), Gene(base_rate, 2, 6682.0 * 0.34, 7443.0 * 0.34)]
    genes_convergent_3 = [Gene(base_rate, 3, 2620.0 * 0.34, 3723.0 * 0.34), Gene(base_rate * induction, 1, 4225.0 * 0.34, 5003.0 * 0.34), Gene(base_rate, 2, 6387.0 * 0.34, 5669.0 * 0.34)]
    genes_divergent_3 = [Gene(base_rate, 3, 3398.0 * 0.34, 2487.0 * 0.34), Gene(base_rate * induction, 1, 4917.0 * 0.34, 4151.0 * 0.34), Gene(base_rate, 2, 6682.0 * 0.34, 7443.0 * 0.34)]
    simulate_full_examples("output/julia_s2_example_sims.h5", n_examples, "3_gene.tandem",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent, σ2_coeff=σ2), bcs_3_gene,
            genes_tandem_3, 3, 15000.0)
    println("Done with 3_gene.tandem with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction, "\n\t σ2: ", σ2)
    simulate_full_examples("output/julia_s2_example_sims.h5", n_examples, "3_gene.convergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent, σ2_coeff=σ2), bcs_3_gene,
            genes_convergent_3, 3, 15000.0)
    println("Done with 3_gene.convergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction, "\n\t σ2: ", σ2)
    simulate_full_examples("output/julia_s2_example_sims.h5", n_examples, "3_gene.divergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent, σ2_coeff=σ2), bcs_3_gene,
            genes_divergent_3, 3, 15000.0)
    println("Done with 3_gene.divergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction, "\n\t σ2: ", σ2)
    
    println("Ran entire round in ", time() - start_time, " seconds")
end