using DiffEqJump
using DifferentialEquations.EnsembleAnalysis
using TanglesModel.StochasticToggle
using Plots
using RandomNumbers

function adjust_ics(prob, i, _)
    ics = [i % 20, trunc(Int64, i / 20)]
    print(ics)
    generate_jump_problem(ics, (0.0, 100000.0), 5.0, 8.0, 1/160.0, 1/1200.0)
end
# Interesting regime: 5.0, 8.0, 1/20.0, 1/1200.0

prob = generate_jump_problem([5,5], (0.0,100000.0), 5.0, 8.0, 1/160.0, 1/1200.0)
ensemble_prob = EnsembleProblem(prob, prob_func=adjust_ics)
sol= solve(ensemble_prob, SSAStepper(), trajectories=400)

plot(componentwise_vectors_timepoint(sol, 22000)..., seriestype=:scatter, markercolor="gray", alpha=0.1, label="")

timepoint_3d = cat([
    cat(componentwise_vectors_timepoint(sol, i)..., dims=2) for i in range(0.0, stop=10000.0, length=50)]...,
    dims=3)

n_frames = 500
timepoints = range(0.0, stop=10000.0, length=n_frames)
solution_3d = permutedims(cat([
    cat(componentwise_vectors_timepoint(sol, i)..., dims=2) for i in range(0.0, stop=10000.0, length=n_frames)]...,
    dims=3),[3, 1, 2])
@gif for i ∈ 1:n_frames
    plot(componentwise_vectors_timepoint(sol, timepoints[i])...,
        xlims=[0,20],
        ylims=[0,20],
        seriestype=:scatter, markercolor="gray", alpha=0.1, label="")
    if i > 1
        plot!(
            solution_3d[max(1,i - 10):i,:,1],
            solution_3d[max(1,i-10):i,:,2],
            alpha=range(0.05, stop=0.0, length=(i - max(0, i - 11))),
            color="purple",
            label="")
    end
end
