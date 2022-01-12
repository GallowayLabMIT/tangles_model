using DiffEqJump
# Parameters are p = (n, K, α, β)
n(p) = p[1]
K(p) = p[2]
α(p) = p[3]
β(p) = p[4]
rate_A_creation(u,p,_) = α(p) * K(p) / (K(p) + u[2]^n(p))
rate_A_degradation(u,p,_) = β(p) * u[1]
rate_B_creation(u,p,_) = α(p) * K(p) / (K(p) + u[1]^n(p))
rate_B_degradation(u,p,_) = β(p) * u[2]

function A_creation_affect!(integrator)
    integrator.u[1] += 1
end

function B_creation_affect!(integrator)
    integrator.u[2] += 1
end

function A_degradation_affect!(integrator)
    integrator.u[1] -= 1
end

function B_degradation_affect!(integrator)
    integrator.u[2] -= 1
end

function generate_jump_problem(ic::Vector{Int64}, tspan::Tuple{Float64,Float64}, n::Float64, K::Float64, α::Float64, β::Float64)
    A_create_jump =  ConstantRateJump(rate_A_creation,    A_creation_affect!)
    B_create_jump =  ConstantRateJump(rate_B_creation,    B_creation_affect!)
    A_degrade_jump = ConstantRateJump(rate_A_degradation, A_degradation_affect!)
    B_degrade_jump = ConstantRateJump(rate_B_degradation, B_degradation_affect!)

    return JumpProblem(
        DiscreteProblem(ic, tspan, (n, K, α, β)),
        Direct(),
        A_create_jump, B_create_jump, A_degrade_jump, B_degrade_jump
    )
end