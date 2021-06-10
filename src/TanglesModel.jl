__precompile__()

module TanglesModel

using Base: UInt32
using MultiScaleArrays: length
using DifferentialEquations
using MultiScaleArrays
import FiniteDiff

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
    base_promoter_rate::Float64     # (1 / sec) The base rate of RNAP initiation.
    topoisomerase_rate::Float64     # (1 / sec) The base rate of topoisomerase activity.
    mRNA_degradation_rate::Float64  # (1 / sec) The base rate of mRNA degradation.
end
DEFAULT_SIM_PARAMS = SimulationParameters(
    DEFAULT_mRNA_PARAMS, DEFAULT_RNAP_PARAMS, DEFAULT_DNA_PARAMS,
    298, 1 / 120, 1 / 1200, 1 / 1200
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

struct Gene
    base_rate::Float64
    idx::UInt32
    start::Float64
    terminate::Float64
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
        σ_p
)
end

mutable struct TanglesArray <: DEDataArray{Float64,1}
    x::Array{Float64,1}
    mRNA::Array{Int32,1}
    polymerase_direction::Array{Int64, 1}
    polymerase_stop::Array{Float64,1}
    polymerase_gene::Array{UInt32,1}
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
    ω0::Float64)
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
    ω0::Float64)
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
        c.mRNA[c.polymerase_gene[idx]] += 1
        deleteat!(c.x, remove_idx)
        deleteat!(c.polymerase_direction, idx)
        deleteat!(c.polymerase_stop, idx)
        deleteat!(c.polymerase_gene, idx)
    end
end

function FiniteDiff.resize!(c::TanglesArray, i::Int64)
    resize!(c.x, i)
end


function extend_rnap!(c::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, position::Float64, insert_idx::UInt32, twist::Float64, gene::UInt32, terminate_end::Float64)
    extend_rnap!(c.u, position, insert_idx, twist, gene, terminate_end)
end

function extend_rnap!(c::TanglesArray, position::Float64, insert_idx::UInt32, twist::Float64, gene::UInt32, terminate_end::Float64)
    # Insert items in reverse order
    insert!(c.x,1 + 3 * (insert_idx - 1), 0.0)
    insert!(c.x,1 + 3 * (insert_idx - 1), twist)
    insert!(c.x,1 + 3 * (insert_idx - 1), position)
    insert!(c.polymerase_direction, insert_idx, terminate_end > position ? 1 : -1)
    insert!(c.polymerase_stop, insert_idx, terminate_end)
    insert!(c.polymerase_gene, insert_idx, gene)
end

function add_polymerase!(integrator, position::Float64, gene::UInt32, terminate_end::Float64)
    print("Attempting to add polymerase. Current u:")
    println(integrator.u)
    insert_idx::UInt32, twist::Float64 = interp_twist(position, integrator.u, integrator.p.bc_params, integrator.p.sim_params.ω_0)
    println("Adding to full_cache vars")
    for c in full_cache(integrator)
        extend_rnap!(c, position, insert_idx, twist,  gene, terminate_end)
    end
    println("Printing new cache...")
    for c in full_cache(integrator)
        println(c)
    end
    println("Done adding!")
end


function polymerase_initiation_rate(u::ExtendedJumpArray{Float64, 1, TanglesArray, Vector{Float64}}, p::TanglesParams, t, promoter::Promoter)::Float64
    return polymerase_initiation_rate(u.u, p, t, promoter)
end

function polymerase_initiation_rate(u::TanglesArray, p::TanglesParams, t, promoter::Promoter)::Float64
    # Initiation rate is zero if initiation site is occupied
    if length(u.x) == 1
        return promoter.base_rate
    end

    if minimum(abs,u.x[1:3:end-1] .- promoter.position) < p.sim_params.r_rnap * 2
        return 0.0
    end

    energy::Float64 = torque_response(supercoiling_density(
        u, get_sc_region(promoter.position, u), p.bc_params, p.sim_params.ω_0), p) * 1.2 * 2.0 * π
    sc_rate_factor::Float64 = exp(-energy / (p.sim_params.k_b * p.sim_params.T))
    
    return promoter.base_rate * (promoter.supercoiling_dependent ? sc_rate_factor : 1.0)
end

function generate_jump(gene::Gene, sc_dependent::Bool)::VariableRateJump
    return VariableRateJump(
                (u,p,t) -> polymerase_initiation_rate(u,p,t,Promoter(gene.start,gene.base_rate,sc_dependent)),
                (int)   -> add_polymerase!(int, gene.start, gene.idx, gene.terminate))
end


function terminate_polymerase!(integrator)
    # Identify polymerase to be terminated
    print("Attempting to terminate. Current u:")
    println(integrator.u)
    p_idx::UInt32 = findmin(abs.(integrator.u.u[1:3:(end-1)] - integrator.u.u.polymerase_stop))[2]
    # Update gene lists
    remove_idx = (1 + ((p_idx - 1) * 3)):(p_idx * 3)
    println("Deleting full_cache vars")
    # Update all internal caches
    for c in full_cache(integrator)
        terminate!(c, p_idx)
    end
    println("Deleting non-user cache")
    deleteat_non_user_cache!(integrator,remove_idx)
    println("Printing remaining cache...")
    for c in full_cache(integrator)
        println(c)
    end
    println("Done terminating")
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
    println(u)

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
        τ = torque_response(σ_f - σ_b,params)

        drag = η * z^α
        dϕ = drag * v * ω0 / (χ + drag) .- τ / (χ + drag)

        du[1 + (ns * (i - 1))] = v
        du[2 + (ns * (i - 1))] = dϕ
        du[3 + (ns * (i - 1))] = abs(v)
    end
    du[end] = 0.0
end

function out_of_domain(u, p, t)
    println("In OOD check")
    print(length(u.u))
    print(",")
    println(u.u)
    if length(u.u) == 1
        return false
    end
    return any(diff(u.u[1:3:end-1]) .< 0)
end

function build_problem(sim_params::SimulationParameters, bcs::BoundaryParameters, t_end::Float64)
    u0 = TanglesArray([0.0], [0], [], [], [])
    problem = ODEProblem(tangles_derivatives!, u0, [0.0, t_end], TanglesParams(
        InternalParameters(sim_params), bcs))
    termination_callback = ContinuousCallback(polymerase_termination_check, terminate_polymerase!,
                                save_positions=(true,true))
    jump_problem = JumpProblem(problem, Direct(), generate_jump(Gene(1/120, 1, 100, 500), false))
    return () -> solve(jump_problem, Tsit5(), callback=termination_callback, dtmax=10, isoutofdomain=out_of_domain)
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

tangles_derivatives!([0.0, 0.0, 0.0, 0.0], TanglesArray([100,0,0, 0], [0], [1], [500], [1]),pc_linear_params, 0)
tangles_derivatives!([0.0, 0.0, 0.0, 0.0], TanglesArray([100,0,0, 0], [0], [1], [500], [1]),pc_circular_params, 0)

end # module TanglesModel