__precompile__()

module TanglesModel


using Base: UInt32
using MultiScaleArrays: length
using DifferentialEquations
using MultiScaleArrays
using HDF5
using Interpolations
import FiniteDiff

module StochasticToggle
include("StochasticToggleModel.jl")
export generate_jump_problem
end

export mRNA_Parameters, DEFAULT_mRNA_PARAMS
export RNAP_Parameters, DEFAULT_RNAP_PARAMS, DNA_Parameters, DEFAULT_DNA_PARAMS
export SimulationParameters, DEFAULT_SIM_PARAMS
export LinearBoundaryParameters, CircularBoundaryParameters, UncoupledGene, MultiUncoupledGene, CoupledGene, MultiCoupledGene, DiscreteConfig
export simulate_full_examples, simulate_summarized_runs, simulate_discrete_runs

function check_nonnegative(x, name)
    if x < 0
        error(name + " must be non-negative!")
    end
end

function check_positive(x, name)
    if x <= 0
        error(name + " must be positive!")
    end
end

struct mRNA_Parameters
    drag_coeff::Float64 # (pN nm^(drag_exponent - 1)) The base drag (e.g. α in α (nascent length)^exp)
    drag_exponent::Float64 # (unitless) The exponent in the drag equation.
    function mRNA_Parameters(drag_coeff,drag_exponent)
        check_positive(drag_coeff, "mRNA drag coefficient")
        new(drag_coeff, drag_exponent)
    end
end
DEFAULT_mRNA_PARAMS = mRNA_Parameters(1/20, 1)

struct RNAP_Parameters
    radius::Float64 # (nm) The width of a RNA polymerase.
    max_velocity::Float64 # (nm / s) The maximum speed of a RNA polymerase.
    stall_torque::Float64 # (pN nm) The applied torque at which the RNA polymerase slows down and stops.
    stall_width::Float64 # (pN nm) The width over which the polymerase stops, smoothing out the stall behavior.
    function RNAP_Parameters(radius, max_velocity, stall_torque, stall_width)
        check_nonnegative(radius, "RNAP radius")
        check_nonnegative(max_velocity, "RNAP max velocity")
        check_positive(stall_width, "Stall torque width")
        new(radius, max_velocity, stall_torque, stall_width)
    end
end
DEFAULT_RNAP_PARAMS = RNAP_Parameters(15, 20, 12, 3)

struct DNA_Parameters
    applied_force::Float64 # (pN) The assumed constant applied force to the DNA.
    twist_mobility::Float64 # (s pN nm) The DNA twist mobility.
    bend_persistance_length::Float64 # (pN) The stat-mech distance over which bend motion is correlated.
    twist_persistance_length::Float64 # (pN) The stat-mech distance over which twist motion is correlated.
    plectoneme_twist_persistance_length::Float64 # (pN) The stat-mech distance over which plectonemic twist motion is correlated.
    function DNA_Parameters(applied_force, twist_mobility, bend_plength, twist_plength, plectoneme_twist_plength)
        check_positive(applied_force, "Applied force")
        check_positive(twist_mobility, "Twist mobility")
        check_positive(bend_plength, "Bend persistance length")
        check_positive(twist_plength, "Twist persistance length")
        check_positive(plectoneme_twist_plength, "Plectoneme twist persistance length")
        new(applied_force, twist_mobility, bend_plength, twist_plength, plectoneme_twist_plength)
    end
end
DEFAULT_DNA_PARAMS = DNA_Parameters(1, 10, 50, 95, 24)

struct SimulationParameters
    mRNA_params::mRNA_Parameters
    RNAP_params::RNAP_Parameters
    DNA_params::DNA_Parameters
    temperature::Float64            # (K) The temperature over which to run the simulation.
    topoisomerase_rate::Float64     # (1 / sec) The base rate of topoisomerase activity.
    mRNA_degradation_rate::Float64  # (1 / sec) The base rate of mRNA degradation.
    sc_dependent::Bool              # If supercoiling-dependent initiation is used
    σ2_coeff::Float64               # The leading coefficent on the σ^2 term
end
DEFAULT_SIM_PARAMS = SimulationParameters(
    DEFAULT_mRNA_PARAMS, DEFAULT_RNAP_PARAMS, DEFAULT_DNA_PARAMS,
    298, 1 / 1200, 1 / 1200, false, 0.0
)

struct InternalParameters
    k_b::Float64    # (pN nm / K) Boltzmann constant.
    ω_0::Float64    # (1 / nm) Relaxed twist frequency.
    A::Float64      # (nm) DNA bend persistance length.
    C::Float64      # (nm) DNA twist persistance length.
    P::Float64      # (nm) DNA plectoneme twist persistence length.
    T::Float64      # (nm) Temperature.
    f::Float64      # (pN) Constant force applied.
    v_0::Float64    # (nm / s) Max speed of RNAP.
    τ_c::Float64    # (pN nm) Stall torque.
    τ_w::Float64    # (pN nm) Stall torque distribution width.
    r_rnap::Float64 # (nm) RNAP radius.
    χ::Float64      # (s pN nm) DNA twist mobility.
    η::Float64      # (pN nm^(α - 1)) mRNA drag coefficient.
    α::Float64      # (unitless) mRNA drag power-law exponent.
    τ_0::Float64    # (pN nm) Intermediate torque value.
    τ_s::Float64    # (pN nm) Critical stretched-phase torque.
    τ_p::Float64    # (pN nm) Critical plectonemic torque.
    σ_s::Float64    # (unitless) Critical stretched-phase supercoiling density.
    σ_p::Float64    # (unitless) Critical plectonemic-phase supercoiling density.
    topo_rate::Float64      # (1 / s) base topoisomerase activity rate
    mRNA_deg_rate::Float64  # (1 / s) base mRNA degradation rate
    sc_dependent::Bool      # If supercoiling dependent initiation is used
    σ2_coeff::Float64       # The leading coefficent on the σ^2 term
end

abstract type BoundaryParameters end

struct CircularBoundaryParameters <: BoundaryParameters
    length::Float64
end

struct LinearBoundaryParameters <: BoundaryParameters
    length::Float64
    left_is_free::Bool
    right_is_free::Bool
end

struct Promoter
    position::Float64
    base_rate::Float64
    supercoiling_dependent::Bool
end

abstract type Gene end

struct UncoupledGene <: Gene
    base_rate::Float64
    idx::UInt32
    start::Float64
    terminate::Float64
end

struct MultiUncoupledGene <: Gene
    base_rate::Float64
    idxes::Array{UInt32,1}
    start::Float64
    terminate::Float64
end

struct CoupledGene <: Gene
    base_rate::Float64
    idx::UInt32
    start::Float64
    terminate::Float64
    coupling_function
end

struct MultiCoupledGene <: Gene
    base_rate::Float64
    idxes::Array{UInt32,1}
    start::Float64
    terminate::Float64
    coupling_function
end

struct DiscreteConfig
    genes::Array{<:Gene}
    n_other_discrete::Int64
    discrete_reactions::Vector{Pair{Function, Vector{Pair{Int64, Int64}}}} # Pairs of the form (propensity_func(discrete, t) => Array{Pair{Int64, Int64}})
                                    # where the second array of pairs is the stoichiometry
    function DiscreteConfig(genes::Array{<:Gene})
        new(genes, 0, [])
    end
    function DiscreteConfig(genes::Array{<:Gene}, n_other_discrete::Int64, discrete_reactions::Vector{Pair{Function, Vector{Pair{Int64,Int64}}}})
        n_genes = length(genes)
        for (_, stoich_coeffs) in discrete_reactions
            for (species_id, stoich_coeff) in stoich_coeffs
                if species_id < 1 || species_id > (n_genes + n_other_discrete)
                    error("Invalid reaction passed in discrete_reactions!\n\tid:" + species_id + "\n\tstoich coeff:", stoich_coeff)
                end
            end
        end
        new(genes, n_other_discrete, discrete_reactions)
    end
end


function InternalParameters(sim_params::SimulationParameters)
    # Computes intermediate values, initializing a InternalParameters struct
    # from the better-documented SimulationParameters object.
    k_b = 1380649 / 100000000 # pN nm / K
    ω_0 = 1.85 # 1 / nm
    T = sim_params.temperature
    f = sim_params.DNA_params.applied_force
    A = sim_params.DNA_params.bend_persistance_length
    C = sim_params.DNA_params.twist_persistance_length
    P = sim_params.DNA_params.plectoneme_twist_persistance_length

    # Compute intermediates
    c = k_b * T * C * ω_0^2
    p = k_b * T * P * ω_0^2
    g = f - sqrt(k_b * T * f / A)
    cs = c * (1 - C / (4 * A) * sqrt(k_b * T / (A * f)))

    # Compute critical values
    τ_s = cs / ω_0
    τ_0 = sqrt(2 * p * g / (ω_0^2 * (1 - p / cs)))
    τ_p = p / ω_0
    σ_s = 1 / cs * sqrt(2 * p * g / (1 - p / cs))
    σ_p = 1 / p * sqrt(2 * p * g / (1 - p / cs))

    InternalParameters(
        1380649 / 100000000, # k_b
        1.85, # ω_0
        sim_params.DNA_params.bend_persistance_length, # A
        sim_params.DNA_params.twist_persistance_length, # C
        sim_params.DNA_params.plectoneme_twist_persistance_length, # P
        sim_params.temperature, # T
        sim_params.DNA_params.applied_force, # f
        sim_params.RNAP_params.max_velocity, # v_0
        sim_params.RNAP_params.stall_torque, # τ_c
        sim_params.RNAP_params.stall_width, # τ_w
        sim_params.RNAP_params.radius, # r_rnap
        sim_params.DNA_params.twist_mobility, # χ
        sim_params.mRNA_params.drag_coeff, # η
        sim_params.mRNA_params.drag_exponent, # α
        τ_0,
        τ_s,
        τ_p,
        σ_s,
        σ_p,
        sim_params.topoisomerase_rate,
        sim_params.mRNA_degradation_rate,
        sim_params.sc_dependent,
        sim_params.σ2_coeff
)
end

mutable struct TanglesArray <: DEDataArray{Float64,1}
    x::Array{Float64,1}
    discrete_components::Array{Int32,1}
    polymerase_direction::Array{Int64, 1}
    polymerase_stop::Array{Float64,1}
    polymerase_gene::Array{Array{UInt32,1}}
end


function get_sc_region(position::Float64, u::TanglesArray)::Int64
    loc::Int64 = 1
    num_polymerases::Int64 = convert(UInt32, (length(u) - 1) / 3)
    x = @view u.x[1:3:(end-1)]
    while loc <= num_polymerases && x[loc] < position
        loc += 1
    end
    return loc
end

function interp_twist(position::Float64, u::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, bcs::BoundaryParameters, ω0::Float64)
    return interp_twist(position, u.u, bcs, ω0)
end

function interp_twist(position::Float64, u::TanglesArray, bcs::LinearBoundaryParameters, ω0::Float64)
    # Returns (insert_idx, insert_twist)
    num_polymerases::Int64 = convert(UInt32, (length(u) - 1) / 3)

    if num_polymerases == 0
        return (1, 0.0)
    end

    x = @view u.x[1:3:(end-1)]
    ϕ = @view u.x[2:3:(end-1)]

    insert_location::Int64 = get_sc_region(position, u)

    σ_target = supercoiling_density(u, insert_location, bcs, ω0)
    # σ = Δϕ / (Δx (-ω0)), so Δϕ = -ω0 σ Δx

    if insert_location == 1
        # Special case, depends on BCs
        if bcs.left_is_free
            return (insert_location, ϕ[insert_location])
        end
        return (insert_location, σ_target * -ω0 * position)
    end

    if insert_location == num_polymerases && bcs.right_is_free
        # Special case, no twist needed on free right end
        return (insert_location, ϕ[insert_location])
    end
    # Calculate distance to previous polymerase and use that to set twist.
    return (insert_location,
        σ_target * -ω0 * (position - x[insert_location - 1]) + ϕ[insert_location - 1])
end

function interp_twist(position::Float64, u::TanglesArray, bcs::CircularBoundaryParameters, ω0::Float64)
    # Returns (insert_idx, insert_twist)
    num_polymerases::Int64 = convert(UInt32, (length(u) - 1) / 3)

    if num_polymerases == 0
        return (1, 0.0)
    end

    x = @view u.x[1:3:(end-1)]
    ϕ = @view u.x[2:3:(end-1)]

    insert_location::Int64 = get_sc_region(position, u)
    σ_target = supercoiling_density(u, insert_location, bcs, ω0)

    # σ = Δϕ / (Δx (-ω0)), so Δϕ = -ω0 σ Δx
    if insert_location == 1
        # Here, calculate this by the distance to the last polymerase.
        return (insert_location,
            ϕ[num_polymerases] + -ω0 * σ_target * (position + bcs.length - x[num_polymerases]))
    end
    return (insert_location,
            σ_target * -ω0 * (position - x[insert_location - 1]) + ϕ[insert_location - 1])
end

function supercoiling_density(
    u::TanglesArray,
    i::Int64,
    bcs::CircularBoundaryParameters,
    ω0::Float64)::Float64
    if length(u.x) == 1
        return 0
    end
    # Computes supercoiling density in N regions, accounting for the circular behavior
    # of the system by duplicating the first region density into the N + 1 slot used later.
    num_rnap = convert(Int32, (length(u.x) - 1) / 3)
    x = @view u[1:3:(end-1)]
    ϕ = @view u[2:3:(end-1)]
    if i == 1 || i == num_rnap + 1
        return (ϕ[1] - ϕ[end]) / (x[1] + bcs.length - x[end]) / -ω0
    end
    return (ϕ[i] - ϕ[i - 1]) / (x[i] - x[i - 1]) / -ω0
end

function supercoiling_density(
    u::TanglesArray,
    i::Int64,
    bcs::LinearBoundaryParameters,
    ω0::Float64)::Float64
    # Computes supercoiling density in N+1 regions, allowing for free or fixed left/right
    # boundary conditions

    if length(u.x) == 1
        return 0
    end

    num_rnap = convert(Int32, (length(u.x) - 1) / 3)
    x = @view u[1:3:(end-1)]
    ϕ = @view u[2:3:(end-1)]

    if i == 1
        if bcs.left_is_free
            return 0
        else
            return ϕ[1] / (x[1] * -ω0)
        end
    elseif i == num_rnap + 1
        if bcs.right_is_free
            return 0
        else
            return -ϕ[end] / (bcs.length - x[end]) / -ω0
        end
    end

    return (ϕ[i] - ϕ[i - 1]) / (x[i] - x[i - 1]) / -ω0
end

struct TanglesParams
    sim_params::InternalParameters
    bc_params::BoundaryParameters
end


function torque_response(σ::Float64, params::TanglesParams)::Float64
    # Computes the torque response at a specified supercoiling density
    abs_σ = abs(σ)
    # Three regime torque region, from:
    # 0   <= |σ| < σ_s => σ τ_s
    # σ_s <= |σ| < σ_p => sgn(σ) τ_0   -- phase transition regime
    # σ_p <= |σ|       => σ τ_p
    if abs_σ < params.sim_params.σ_s
        return params.sim_params.τ_s * σ
    elseif abs_σ < params.sim_params.σ_p
        return params.sim_params.τ_0 * sign(σ)
    end
    return params.sim_params.τ_p * σ
end

function polymerase_velocity(σ_b::Float64, σ_f::Float64, params::TanglesParams)
    # Compute the polymerase velocity given neighboring supercoiling densities.

    # τ_s = stall torque, τ_w = stall torque width
    stall_behind = (abs(torque_response(σ_b, params)) - params.sim_params.τ_s) / params.sim_params.τ_w
    stall_ahead =  (abs(torque_response(σ_f, params)) - params.sim_params.τ_s) / params.sim_params.τ_w

    # Restrict exponential argument to 10 (e.g. truncate the stall) to prevent overflows.
    return params.sim_params.v_0 / (
          (1 + exp(min(10.0, stall_behind)))
        * (1 + exp(min(10.0, stall_ahead)))
        )
end

function polymerase_termination_check(u::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, t, integrator)
    return polymerase_termination_check(u.u, t, integrator)
end

function polymerase_termination_check(u::TanglesArray, t, integrator)
    if length(u) > 1
        return min(((u.polymerase_stop - u[1:3:(end-1)]) .* u.polymerase_direction)...)
    end
    # Never trigger polymerase termination
    return 1
end

function terminate!(c::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, idx::UInt32)
    terminate!(c.u, idx)
end

function terminate!(c::TanglesArray, idx::UInt32)
    remove_idx = (1 + ((idx - 1) * 3)):(idx * 3)

    if remove_idx[end] > length(c)
        @warn "Data array is too small!"
    else
        for i in c.polymerase_gene[idx]
            c.discrete_components[i] += 1
        end
        deleteat!(c.x, remove_idx)
        deleteat!(c.polymerase_direction, idx)
        deleteat!(c.polymerase_stop, idx)
        deleteat!(c.polymerase_gene, idx)
    end
end

function _destroy_mRNA!(c::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, idx::UInt32)
    _destroy_mRNA!(c.u, idx)
end

function _destroy_mRNA!(c::TanglesArray, idx::UInt32)
    c.discrete_components[idx] -= 1
end

function degrade_mRNA!(integrator, idx::UInt32)
    for c in full_cache(integrator)
        _destroy_mRNA!(c, idx)
    end
end

function _update_discrete!(c::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, updates::Array{Pair{Int64, Int64}})
    _update_discrete!(c.u, updates)
end

function _update_discrete!(c::TanglesArray, updates::Array{Pair{Int64, Int64}})
    for (species_id, update_amount) in updates
        c.discrete_components[species_id] += update_amount
    end
end

function update_discrete!(integrator, updates::Array{Pair{Int64, Int64}})
    for c in full_cache(integrator)
        _update_discrete!(c, updates)
    end
end

function FiniteDiff.resize!(c::TanglesArray, i::Int64)
    resize!(c.x, i)
end


function extend_rnap!(c::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, position::Float64, insert_idx::UInt32, twist::Float64, genes::Array{UInt32,1}, terminate_end::Float64)
    extend_rnap!(c.u, position, insert_idx, twist, genes, terminate_end)
end

function extend_rnap!(c::TanglesArray, position::Float64, insert_idx::UInt32, twist::Float64, genes::Array{UInt32,1}, terminate_end::Float64)
    # Insert items in reverse order
    insert!(c.x,1 + 3 * (insert_idx - 1), 0.0)
    insert!(c.x,1 + 3 * (insert_idx - 1), twist)
    insert!(c.x,1 + 3 * (insert_idx - 1), position)
    insert!(c.polymerase_direction, insert_idx, terminate_end > position ? 1 : -1)
    insert!(c.polymerase_stop, insert_idx, terminate_end)
    insert!(c.polymerase_gene, insert_idx, genes)
end

function add_polymerase!(integrator, position::Float64, genes::Array{UInt32,1}, terminate_end::Float64)
    #print("Attempting to add polymerase. Current u:")
    #println(integrator.u)
    insert_idx::UInt32, twist::Float64 = interp_twist(position, integrator.u, integrator.p.bc_params, integrator.p.sim_params.ω_0)
    #println("Adding to full_cache vars")
    for c in full_cache(integrator)
        extend_rnap!(c, position, insert_idx, twist,  genes, terminate_end)
    end
    #println("Printing new cache...")
    #for c in full_cache(integrator)
    #    println(c)
    #end
    #println("Done adding!")
end

function mRNA_degradation_rate(u::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, p::TanglesParams, t, idx::UInt32)
    mRNA_degradation_rate(u.u, p, t, idx)
end
function mRNA_degradation_rate(u::TanglesArray, p::TanglesParams, t, idx::UInt32)
    return p.sim_params.mRNA_deg_rate * u.discrete_components[idx]
end

function polymerase_initiation_rate(u::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, p::TanglesParams, t, promoter::Promoter, coupling_func)::Float64
    return polymerase_initiation_rate(u.u, p, t, promoter, coupling_func)
end


function polymerase_initiation_rate(u::TanglesArray, p::TanglesParams, t, promoter::Promoter, coupling_func)::Float64
    # Initiation rate is zero if initiation site is occupied
    if length(u.x) == 1
        return promoter.base_rate * coupling_func(u.discrete_components, t)
    end

    if minimum(abs,u.x[1:3:end-1] .- promoter.position) < p.sim_params.r_rnap * 2
        return 0.0
    end

    σ::Float64 = supercoiling_density(u, get_sc_region(promoter.position, u), p.bc_params, p.sim_params.ω_0)
    energy::Float64 = (torque_response(σ, p) + p.sim_params.σ2_coeff * (σ / p.sim_params.σ_s)^2 * p.sim_params.τ_0) * 1.2 * 2.0 * π
    sc_rate_factor::Float64 = min(50.0,exp(-energy / (p.sim_params.k_b * p.sim_params.T)))

    return promoter.base_rate * coupling_func(u.discrete_components, t) * (promoter.supercoiling_dependent ? sc_rate_factor : 1.0)
end

function generate_jump(gene::UncoupledGene, sc_dependent::Bool)::VariableRateJump
    return VariableRateJump(
                (u,p,t) -> polymerase_initiation_rate(
                    u,p,t,
                    Promoter(gene.start,gene.base_rate,sc_dependent),(mRNA, t)->1.0),
                (int)   -> add_polymerase!(int, gene.start, [gene.idx], gene.terminate))
end

function generate_jump(gene::MultiUncoupledGene, sc_dependent::Bool)::VariableRateJump
    return VariableRateJump(
                (u,p,t) -> polymerase_initiation_rate(
                    u,p,t,
                    Promoter(gene.start,gene.base_rate,sc_dependent),(mRNA, t)->1.0),
                (int)   -> add_polymerase!(int, gene.start, gene.idxes, gene.terminate))
end

function generate_jump(gene::CoupledGene, sc_dependent::Bool)::VariableRateJump
    return VariableRateJump(
                (u,p,t) -> polymerase_initiation_rate(
                    u,p,t,
                    Promoter(gene.start,gene.base_rate,sc_dependent),gene.coupling_function),
                (int)   -> add_polymerase!(int, gene.start, [gene.idx], gene.terminate))
end

function generate_jump(gene::MultiCoupledGene, sc_dependent::Bool)::VariableRateJump
    return VariableRateJump(
                (u,p,t) -> polymerase_initiation_rate(
                    u,p,t,
                    Promoter(gene.start,gene.base_rate,sc_dependent),gene.coupling_function),
                (int)   -> add_polymerase!(int, gene.start, gene.idxes, gene.terminate))
end

function terminate_polymerase!(integrator)
    # Identify polymerase to be terminated
    #print("Attempting to terminate. Current u:")
    #println(integrator.u)
    p_idx::UInt32 = findmin(abs.(integrator.u.u[1:3:(end-1)] - integrator.u.u.polymerase_stop))[2]
    # Update gene lists
    remove_idx = (1 + ((p_idx - 1) * 3)):(p_idx * 3)
    #println("Deleting full_cache vars")
    # Update all internal caches
    for c in full_cache(integrator)
        terminate!(c, p_idx)
    end
    #println("Deleting non-user cache")
    deleteat_non_user_cache!(integrator,remove_idx)
    #println("Printing remaining cache...")
    #for c in full_cache(integrator)
    #    println(c)
    #end
    #println("Done terminating")
end

function internal_relax_supercoiling!(u::TanglesArray, n_rnap::Int64, start_idx::Int64, end_idx::Int64, bcs::CircularBoundaryParameters)
    x = @view u.x[1:3:end-1]
    ϕ = @view u.x[2:3:end-1]

    # If we cover all polymerases, relax all:
    if start_idx == 1 && end_idx == n_rnap
        ϕ[start_idx:end_idx] .= 0
        return
    end
    # Handle the edge cases properly by setting the x "boundaries"
    if start_idx == 1
        x_lower = x[end] - bcs.length
        ϕ_lower = ϕ[end]
        x_upper = x[end_idx + 1]
        ϕ_upper = ϕ[end_idx + 1]
    elseif end_idx == n_rnap
        x_lower = x[start_idx - 1]
        ϕ_lower = ϕ[start_idx - 1]
        x_upper = bcs.length + x[1]
        ϕ_upper = ϕ[1]
    else
        x_lower = x[start_idx - 1]
        ϕ_lower = ϕ[start_idx - 1]
        x_upper = x[end_idx + 1]
        ϕ_upper = ϕ[end_idx + 1]
    end
    for i in start_idx:end_idx
        α = (x_upper - x[i]) / (x_upper - x_lower)
        ϕ[i] = (ϕ_lower * α) + (ϕ_upper * (1 - α))
    end
end

function internal_relax_supercoiling!(u::TanglesArray, n_rnap::Int64, start_idx::Int64, end_idx::Int64, bcs::LinearBoundaryParameters)
    x = @view u.x[1:3:end-1]
    ϕ = @view u.x[2:3:end-1]

    # Handle special cases
    # Relax everything if either we cover all polymerases...
    if start_idx == 1 && end_idx == n_rnap
        ϕ[start_idx:end_idx] .= 0
        return
    end
    # or if we are at the edges and that edge is free
    if (start_idx == 1 && bcs.left_is_free) || (start_idx == n_rnap && bcs.right_is_free)
        ϕ[start_idx:end_idx] .= 0
        return
    end

    if start_idx == 1
        # We have at least one polymerase to the right to interpolate with
        for i in start_idx:end_idx
            ϕ[i] = ϕ[end_idx + 1] * x[i] / x[end_idx + 1]
        end
        return
    end

    if end_idx == n_rnap
        # Interpolate from the left
        for i in start_idx:end_idx
            ϕ[i] = ϕ[start_idx - 1] * (bcs.length - x[i]) / (bcs.length - x[start_idx - 1])
        end
        return
    end

    # Otherwise, interpolate between left and right values
    for i in start_idx:end_idx
        α = (x[end_idx + 1] - x[i]) / (x[end_idx + 1] - x[start_idx - 1])
        ϕ[i] = (ϕ[start_idx - 1] * α) + (ϕ[end_idx + 1] * (1 - α))
    end
end

function relax_supercoiling!(u::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, left::Float64, right::Float64, bcs::BoundaryParameters)
    relax_supercoiling!(u.u, left, right, bcs)
end

function relax_supercoiling!(u::TanglesArray, left::Float64, right::Float64, bcs::BoundaryParameters)
    num_polymerases::Int64 = convert(UInt32, (length(u) - 1) / 3)
    if num_polymerases == 0
        return
    end

    x = @view u.x[1:3:end-1]
    ϕ = @view u.x[2:3:end-1]

    start_idx = 1
    end_idx = 1
    for i in 1:num_polymerases
        rnap_pos::Float64 = u.x[1 + 3 * (start_idx - 1)]
        if rnap_pos < left
            start_idx += 1
        end
        if rnap_pos < right
            end_idx += 1
        end
    end
    # Correct for the off-by-one error
    end_idx -= 1
    internal_relax_supercoiling!(u, num_polymerases, start_idx,end_idx, bcs)
end


function relax_supercoiling!(integrator, left::Float64, right::Float64)
    for c in full_cache(integrator)
        relax_supercoiling!(c, left, right, integrator.p.bc_params)
    end
end


function tangles_derivatives!(du, u::TanglesArray, params::TanglesParams, t)
    # Arguments
    # ---------
    # du: preallocated, mutable derivative array.
    # u: A (3N + 1,) vector encoding (location, phi, mRNA_length)
    #    for each polymerase. The first variable is a dummy variable that does not change,
    #    but is needed to support the zero polymerase case.
    # p: named_tuple
    #   A tuple containing simulation parameters. In particular, it contains:
    #       sim_params: a SimulationParameters object
    #       bc_params: a BoundaryParameters object

    ns = 3 # number of states
    num_polymerases::Int64 = convert(Int64, (length(u) - 1) / ns)

    # Unpack useful constants
    ω0 = params.sim_params.ω_0
    α = params.sim_params.α
    η = params.sim_params.η
    χ = params.sim_params.χ

    for i = 1:num_polymerases
        x = u[1 + (ns * (i - 1))]
        ϕ = u[2 + (ns * (i - 1))]
        z = u[3 + (ns * (i - 1))]

        σ_b = supercoiling_density(u, i, params.bc_params, params.sim_params.ω_0)
        σ_f = supercoiling_density(u, i + 1, params.bc_params, params.sim_params.ω_0)
        # Compute polymerase direction and velocity
        v = u.polymerase_direction[i] * polymerase_velocity(σ_b, σ_f, params)
        τ = torque_response(σ_f, params) - torque_response(σ_b,params)

        drag = η * z^α
        dϕ = drag * v * ω0 / (χ + drag) .- τ / (χ + drag)

        du[1 + (ns * (i - 1))] = v
        du[2 + (ns * (i - 1))] = dϕ
        du[3 + (ns * (i - 1))] = abs(v)
    end
    du[end] = 0.0
end

function out_of_domain(u, p, t)
    if length(u.u) == 1
        return false
    end
    return any(diff(u.u[1:3:end-1]) .< 0)
end

function build_problem(
    sim_params::SimulationParameters,
    bcs::BoundaryParameters,
    dconfig::DiscreteConfig,
    t_end::Float64;
    tsteps::Int64=-1)
    build_problem(sim_params, bcs, dconfig, t_end, zeros(Int32,length(dconfig.genes)+ dconfig.n_other_discrete), tsteps=tsteps)
end

function build_problem(
    sim_params::SimulationParameters,
    bcs::BoundaryParameters,
    dconfig::DiscreteConfig,
    t_end::Float64,
    ics_discrete::Array{Int32,1};
    tsteps::Int64=-1)
    extra_save_positions = (tsteps == -1 ? (true,true) : (false,false))
    extra_tstops = (tsteps == -1 ? [] : range(0.0, t_end, length=tsteps))
    u0 = TanglesArray([0.0], ics_discrete, [], [], [])
    problem = ODEProblem(tangles_derivatives!, u0, [0.0, t_end], TanglesParams(
        InternalParameters(sim_params), bcs))
    termination_callback = ContinuousCallback(polymerase_termination_check, terminate_polymerase!,
                                save_positions=extra_save_positions)

    # Calculate intergenic regions
    sorted_genes::Array{Gene} = sort(dconfig.genes, by=(g::Gene)->min(g.start, g.terminate))
    # Duplicate gene if there is only one
    if length(sorted_genes) == 1
        push!(sorted_genes, sorted_genes[1])
    end
    topo_jumps = [
        ConstantRateJump(
            (u,p,t)->sim_params.topoisomerase_rate / (length(dconfig.genes) - 1.0),
            (int)->relax_supercoiling!(int,
                min(sorted_genes[i].start, sorted_genes[i].terminate),
                min(sorted_genes[i+1].start, sorted_genes[i+1].terminate)))
            for i in 1:(length(sorted_genes)-1)]
    jump_problem = JumpProblem(problem, Direct(),
        generate_jump.(dconfig.genes, sim_params.sc_dependent)...,
        topo_jumps...,
        [VariableRateJump(
            (u,p,t)->mRNA_degradation_rate(u,p,t,convert(UInt32,i)),
            (int)->degrade_mRNA!(int, convert(UInt32,i))) for i in 1:length(dconfig.genes)]...,
        [ConstantRateJump(
            (u,_,t)-> convert(Float64, propensity(u.u.discrete_components, t)),
            (int)->update_discrete!(int, updates)
            )
        for (propensity, updates) in dconfig.discrete_reactions]...
    ,save_positions=extra_save_positions)
    return () -> solve(jump_problem, Tsit5(), callback=termination_callback, maxiters=1e6, dtmax=10, isoutofdomain=out_of_domain,
                       saveat=extra_tstops, tstops=extra_tstops)
end

function write_bcs(group::HDF5.Group, bcs::LinearBoundaryParameters)
    attributes(group)["bcs.is_circular"] = 0
    attributes(group)["bcs.length"] = bcs.length
    attributes(group)["bcs.left_free"] = bcs.left_is_free
    attributes(group)["bcs.right_free"] = bcs.right_is_free
end
function write_bcs(group::HDF5.Group, bcs::CircularBoundaryParameters)
    attributes(group)["bcs.is_circular"] = 1
    attributes(group)["bcs.length"] = bcs.length
end

function write_h5_attributes(
    group::HDF5.Group,
    comment::String,
    dconfig::DiscreteConfig,
    sim_params::SimulationParameters,
    bcs::BoundaryParameters)
    attributes(group)["gene.start"] = [gene.start for gene in dconfig.genes]
    attributes(group)["gene.end"] = [gene.terminate for gene in dconfig.genes]
    attributes(group)["gene.base_rate"] = [gene.base_rate for gene in dconfig.genes]
    attributes(group)["rates.topo"] = sim_params.topoisomerase_rate
    attributes(group)["rates.mRNA_degradation"] = sim_params.mRNA_degradation_rate
    attributes(group)["rates.sc_dependent"] = sim_params.sc_dependent ? 1.0 : 0.0
    attributes(group)["coeff.sigma_squared"] = sim_params.σ2_coeff
    attributes(group)["coeff.mRNA_drag_exponent"] = sim_params.mRNA_params.drag_exponent
    attributes(group)["coeff.mRNA_drag_coeff"] = sim_params.mRNA_params.drag_coeff
    attributes(group)["rnap.max_velocity"] = sim_params.RNAP_params.max_velocity
    attributes(group)["rnap.stall_torque"] = sim_params.RNAP_params.stall_torque
    attributes(group)["rnap.stall_width"] = sim_params.RNAP_params.stall_width
    attributes(group)["git_status"] = try read(`git log -n1 --format=format:"%H"`, String) * (run(ignorestatus(`git diff-index --quiet HEAD --`)).exitcode == 0 ? "" : "-dirty") catch IOError; "unknown" end
    attributes(group)["comment"] = comment
    write_bcs(group, bcs)
end
function postprocess_to_h5(
    filename::String,
    solution,
    comment::String,
    dconfig::DiscreteConfig,
    sim_params::SimulationParameters,
    bcs::BoundaryParameters)
    n_genes = length(dconfig.genes)
    h5open(filename, "cw") do h5
        run_idx = 1
        while haskey(h5, "tangles_full_run." * lpad(run_idx, 6, "0"))
            run_idx += 1
        end
        g = create_group(h5, "tangles_full_run." * lpad(run_idx, 6, "0"))
        g["time"] = solution.t
        len = length(solution.u)
        width = maximum([length(u.u.x) for u in solution.u]) - 1

        rnap_loc::Matrix{Float64} = -ones(len, convert(Int,width / 3))
        ϕ::Matrix{Float64} = -ones(len, convert(Int,width / 3))
        z::Matrix{Float64} = -ones(len, convert(Int,width / 3))
        mRNA::Matrix{Int32} = zeros(len, n_genes)
        discrete_components::Matrix{Int32} = zeros(len, length(solution.u[1].u.discrete_components - n_genes))
        for (index, val) in enumerate(solution.u)
            x = val.u.x[1:end-1]
            n_polymerases = convert(Int, length(x) / 3)
            rnap_loc[index,1:n_polymerases] = x[1:3:end]
            ϕ[index,1:n_polymerases] = x[2:3:end]
            z[index,1:n_polymerases] = x[3:3:end]
            mRNA[index,:] = val.u.discrete_components[1:n_genes]
            discrete_components[index,:] = val.u.discrete_components[(n_genes+1):end]
        end
        g["rnap_location"] = rnap_loc
        g["phi"] = ϕ
        g["mRNA_length"] = z
        g["mRNA"] = mRNA
        g["discrete_components"] = discrete_components
        write_h5_attributes(g, comment, dconfig, sim_params, bcs)
    end
end

function simulate_full_examples(
    filename::String,
    n_simulations::Int64,
    comment::String,
    sim_params::SimulationParameters,
    bcs::BoundaryParameters,
    dconfig::DiscreteConfig,
    t_end::Float64)
    solver = build_problem(sim_params, bcs, dconfig, t_end)
    for _ in 1:n_simulations
        try
            postprocess_to_h5(filename, solver(), comment, dconfig, sim_params, bcs)
        catch err
            @warn "Solver failed!"
        end
    end
end

# Simulate mRNA concentrations over time, saving just these results
# instead of the full run results or the final mRNA summary.
function simulate_discrete_runs(
    filename::String,
    n_simulations::Int64,
    comment::String,
    sim_params::SimulationParameters,
    bcs::BoundaryParameters,
    dconfig::DiscreteConfig,
    t_end::Float64,
    t_steps::Int64,
    ics_discrete::Array{Int32,1},
    extra_metadata::Dict{String,Float64})
    solver = build_problem(sim_params, bcs, dconfig, t_end, ics_discrete, tsteps=t_steps)
    simulate_discrete_runs(
        solver, filename, n_simulations, comment,
        sim_params, bcs, dconfig, t_end,
        t_steps, extra_metadata)
end
function simulate_discrete_runs(
    filename::String,
    n_simulations::Int64,
    comment::String,
    sim_params::SimulationParameters,
    bcs::BoundaryParameters,
    dconfig::DiscreteConfig,
    t_end::Float64,
    t_steps::Int64,
    extra_metadata::Dict{String,Float64})
    solver = build_problem(sim_params, bcs, dconfig, t_end, tsteps=t_steps)
    simulate_discrete_runs(
        solver, filename, n_simulations, comment,
        sim_params, bcs, dconfig, t_end,
        t_steps, extra_metadata)
end

function simulate_discrete_runs(
    solver,
    filename::String,
    n_simulations::Int64,
    comment::String,
    sim_params::SimulationParameters,
    bcs::BoundaryParameters,
    dconfig::DiscreteConfig,
    t_end::Float64,
    t_steps::Int64,
    extra_metadata::Dict{String,Float64})

    n_genes = length(dconfig.genes)
    mRNA_results = zeros(Int32, n_simulations, n_genes, t_steps)
    discrete_results = zeros(Int32, n_simulations, dconfig.n_other_discrete, t_steps)
    for i = 1:n_simulations
        done = false
        n_repeats = 0
        sol = false
        while !done && n_repeats < 5
            sol = solver()
            n_repeats += 1
            done = (sol.retcode == :Success)
        end
        if n_repeats >= 5
            continue
        end
        timesteps = Interpolations.deduplicate_knots!(sol.t)
        discrete = hcat([sol.u[i].u.discrete_components for i in 1:length(sol)]...)
        interp_discrete = ConstantInterpolation((1:(n_genes + dconfig.n_other_discrete), timesteps), discrete)
        for gene_id = 1:n_genes
            mRNA_results[i,gene_id,:] = interp_discrete(gene_id,range(0.0, stop=t_end, length=t_steps))
        end
        for discrete_id in 1:dconfig.n_other_discrete
            discrete_results[i,discrete_id,:] = interp_discrete(discrete_id + n_genes, range(0.0, stop=t_end, length=t_steps))
        end
    end
    print(".")
    h5open(filename, "cw") do h5
        run_idx = 1
        while haskey(h5, "tangles_discrete_run." * lpad(run_idx, 6, "0"))
            run_idx += 1
        end
        g = create_group(h5, "tangles_discrete_run." * lpad(run_idx, 6, "0"))
        g["mRNA"] = mRNA_results
        g["discrete_components"] = discrete_results
        g["time"] = Vector(range(0.0, stop=t_end, length=t_steps))
        write_h5_attributes(g, comment, dconfig, sim_params, bcs)
        for (key, val) in extra_metadata
            attributes(g)[key] = val
        end
    end
end

function simulate_summarized_runs(
    filename::String,
    n_simulations::Int64,
    comment::String,
    sim_params::SimulationParameters,
    bcs::BoundaryParameters,
    dconfig::DiscreteConfig,
    t_end::Float64)
    solver = build_problem(sim_params, bcs, dconfig, t_end, tsteps=100)
    discrete_results = zeros(Int32, n_simulations, length(dconfig.genes) + dconfig.n_other_discrete)
    for i in 1:n_simulations
        try
            discrete_results[i,:] = solver().u[end].u.discrete_components
        catch err
            @warn "Solver failed!"
        end
    end
    h5open(filename, "cw") do h5
        run_idx = 1
        while haskey(h5, "tangles_summarized_run." * lpad(run_idx, 6, "0"))
            run_idx += 1
        end
        g = create_group(h5, "tangles_summarized_run." * lpad(run_idx, 6, "0"))
        g["final_mRNA"] = discrete_results[:,1:length(dconfig.genes)]
        g["final_discrete"] = discrete_results[:,(length(dconfig.genes)+1):end]
        write_h5_attributes(g, comment, dconfig, sim_params, bcs)
    end
end

# precompile hints
pc_linear_params = TanglesParams(
    InternalParameters(DEFAULT_SIM_PARAMS),
    LinearBoundaryParameters(5000, false, false))
pc_circular_params = TanglesParams(
    InternalParameters(DEFAULT_SIM_PARAMS),
    CircularBoundaryParameters(5000))
ODEProblem(tangles_derivatives!, TanglesArray([0.0], [], [], [], []), [0, 1000.0], pc_linear_params)
ODEProblem(tangles_derivatives!, TanglesArray([0.0], [], [], [], []), [0, 1000.0], pc_circular_params)
# Test uncoupled with multiple BCs
TanglesModel.build_problem(TanglesModel.DEFAULT_SIM_PARAMS, TanglesModel.LinearBoundaryParameters(5000, false, false), DiscreteConfig([TanglesModel.UncoupledGene(1/120.0, 1, 100, 1000), TanglesModel.UncoupledGene(1/120.0, 2, 3000, 2000)]), 1000.0)()
TanglesModel.build_problem(TanglesModel.DEFAULT_SIM_PARAMS, TanglesModel.CircularBoundaryParameters(5000), DiscreteConfig([TanglesModel.UncoupledGene(1/120.0, 1, 100, 1000), TanglesModel.UncoupledGene(1/120.0, 2, 3000, 2000)]), 1000.0)()
# Test coupled with linear BCs
TanglesModel.build_problem(TanglesModel.DEFAULT_SIM_PARAMS, TanglesModel.LinearBoundaryParameters(5000, false, false), DiscreteConfig([TanglesModel.CoupledGene(1/120.0, 1, 100, 1000, (mRNA,_)->100.0 / (100.0 + mRNA[1])), TanglesModel.CoupledGene(1/120.0, 2, 3000, 2000, (_,t)->1000.0/(1000.0 + t))]), 1000.0)()
# Test coupled, with extra reactions
TanglesModel.build_problem(TanglesModel.DEFAULT_SIM_PARAMS, TanglesModel.LinearBoundaryParameters(5000, false, false), DiscreteConfig(
    [
        TanglesModel.CoupledGene(1/120.0, 1, 100, 1000, (mRNA,_)->100.0 / (100.0 + mRNA[1])), TanglesModel.CoupledGene(1/120.0, 2, 3000, 2000, (_,t)->1000.0/(1000.0 + t))],
        2,
        [
            ((discrete,t)->10.0 / (10.0 + discrete[4])) => [3 => 1],
            ((discrete,t)->10.0 / (10.0 + discrete[3])) => [4 => 1],
            ((discrete,t)->discrete[3]) => [3 => -1],
            ((discrete,t)->discrete[4]) => [4 => -1],
        ]), 1000.0)()

end # module TanglesModel
