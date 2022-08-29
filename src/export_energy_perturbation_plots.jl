using TanglesModel
using HDF5

function gen_sim_params(;
    topo_rate_factor::Float64=1.0,
    sc_dependent::Bool=DEFAULT_SIM_PARAMS.sc_dependent,
    σ2_coeff::Float64=0.0,
    torque_perturb::TanglesModel.TorqueFunctionPerturbation=NoTorqueFunctionPerturbation(),
    rnap_init_perturb::TanglesModel.RNAPInitPerturbation=NoRNAPInitPerturbation())
    return SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate * topo_rate_factor,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate,
        sc_dependent,
        σ2_coeff,
        OriginalTopoisomerase(),
        torque_perturb,
        rnap_init_perturb
    )
end

h5open("output/energy_perturbation_plots.h5", "cw") do h5
    i = 0
    for (tperturb, tperturb_name) in zip(
            [NoTorqueFunctionPerturbation(), PositiveSupercoilingBuffering(0.031)],
            ["original", "buffering"]
        ),
        (initperturb, initperturb_name) in zip(
            [NoRNAPInitPerturbation(), RNAPInitEnergyWell(-0.06, 0.125)],
            ["original", "energy_well"]
        )
        i += 1

        g = create_group(h5, "sweep" * lpad(i, 5, "0"))
        torque = zeros(100,1000)
        init_energy = zeros(100,1000)
        eval_points = Vector(range(-0.2,0.2,length=100))
        σ2_coeffs = vcat([0.0], exp10.(range(-3,-0.5,length=999)))


        for i in 1:length(eval_points),
            j in 1:length(σ2_coeffs)

            σ::Float64 = eval_points[i]
            α::Float64 = σ2_coeffs[j]
            sim_params = gen_sim_params(sc_dependent=true,σ2_coeff=α,
                torque_perturb=tperturb, rnap_init_perturb=initperturb)
            p = TanglesModel.TanglesParams(
                TanglesModel.InternalParameters(sim_params), TanglesModel.CircularBoundaryParameters(1000))
            torque[i,j] = TanglesModel.torque_response(σ, p)
            init_energy[i,j] = TanglesModel.polymerase_init_energy(σ, p, p.sim_params.rnap_init_perturbation)
        end
        g["sigma"] = eval_points
        g["s2_coeff"] = σ2_coeffs
        g["torque"] = torque
        g["initiation_energy"] = init_energy
        attributes(g)["torque_perturb"] = tperturb_name
        attributes(g)["rnap_perturb"] = initperturb_name
    end
end
