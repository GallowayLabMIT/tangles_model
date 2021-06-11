using TanglesModel

function gen_sim_params(;
    topo_rate_factor::Float64=1.0,
    sc_dependent::Bool=DEFAULT_SIM_PARAMS.sc_dependent)
    return SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate * topo_rate_factor,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate,
        sc_dependent
    )
end

base_rate = 1.0 / 120.0
n_examples = 5
for topo_multiplier in exp10.(range(-1,1,length=5)),
    is_plasmid in [false, true],
    sc_dependent in [false, true],
    induction in exp10.(range(-2,0.5,length=10))

    start_time = time()

    bcs_2_gene = is_plasmid ? CircularBoundaryParameters(9616) : LinearBoundaryParameters(9616, false, false)
    genes_tandem_2 = [Gene(base_rate * induction, 1, 2329, 3372), Gene(base_rate, 2, 5067, 5810)]
    genes_convergent_2 = [Gene(base_rate * induction, 1, 2329, 3372), Gene(base_rate, 2, 4758, 4050)]
    genes_divergent_2 = [Gene(base_rate * induction, 1, 3614, 2555), Gene(base_rate, 2, 5067, 5810)]

    simulate_full_examples("output/julia_example_sims.h5", n_examples, "2_gene.tandem",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent), bcs_2_gene,
            genes_tandem_2, 2, 15000.0)
    println("Done with 2_gene.tandem with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction)
    simulate_full_examples("output/julia_example_sims.h5", n_examples, "2_gene.convergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent), bcs_2_gene,
            genes_convergent_2, 2, 15000.0)
    println("Done with 2_gene.convergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction)
    simulate_full_examples("output/julia_example_sims.h5", n_examples, "2_gene.divergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent), bcs_2_gene,
            genes_divergent_2, 2, 15000.0)
    println("Done with 2_gene.convergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction)

    bcs_3_gene = is_plasmid ? CircularBoundaryParameters(11238) : LinearBoundaryParameters(11238, false, false)
    genes_tandem_3 = [Gene(base_rate, 3, 2620, 3723), Gene(base_rate * induction, 1, 4225, 5003), Gene(base_rate, 2, 6682, 7443)]
    genes_convergent_3 = [Gene(base_rate, 3, 2620, 3723), Gene(base_rate * induction, 1, 4225, 5003), Gene(base_rate, 2, 6387, 5669)]
    genes_divergent_3 = [Gene(base_rate, 3, 3398, 2487), Gene(base_rate * induction, 1, 4917, 4151), Gene(base_rate, 2, 6682, 7443)]
    simulate_full_examples("output/julia_example_sims.h5", n_examples, "3_gene.tandem",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent), bcs_3_gene,
            genes_tandem_3, 3, 15000.0)
    println("Done with 3_gene.tandem with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction)
    simulate_full_examples("output/julia_example_sims.h5", n_examples, "3_gene.convergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent), bcs_3_gene,
            genes_convergent_3, 3, 15000.0)
    println("Done with 3_gene.convergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction)
    simulate_full_examples("output/julia_example_sims.h5", n_examples, "3_gene.divergent",
            gen_sim_params(topo_rate_factor=topo_multiplier,sc_dependent=sc_dependent), bcs_3_gene,
            genes_divergent_3, 3, 15000.0)
    println("Done with 3_gene.convergent with params:\n\ttopo: ", topo_multiplier, "\n\tis_plasmid: ", is_plasmid,
        "\n\tsc_dependent: ", sc_dependent, "\n\t induction: ", induction)
    
    println("Ran entire round in ", time() - start_time, " seconds")
end