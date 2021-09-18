using TanglesModel
using HDF5

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

h5open("output/sigma_sweep.h5", "cw") do h5
    g = create_group(h5, "sweep")
    values = zeros(100,1000)
    eval_points = Vector(range(-0.2,0.2,length=100))
    σ2_coeffs = vcat([0.0], exp10.(range(-3,-0.5,length=999)))

    for i in 1:length(eval_points),
        j in 1:length(σ2_coeffs)

        σ::Float64 = eval_points[i]
        α::Float64 = σ2_coeffs[j]
        sim_params = gen_sim_params(sc_dependent=true,σ2_coeff=α)
        p = TanglesModel.TanglesParams(
            TanglesModel.InternalParameters(sim_params), TanglesModel.CircularBoundaryParameters(1000))
        values[i,j] = (TanglesModel.torque_response(σ, p) + p.sim_params.σ2_coeff * (σ / p.sim_params.σ_s)^2 * p.sim_params.τ_0) * 1.2 * 2.0 * π
    end
    g["sigma"] = eval_points
    g["s2_coeff"] = σ2_coeffs
    g["values"] = values
end