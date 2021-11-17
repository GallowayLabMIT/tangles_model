using DiffEqJump
using DifferentialEquations.EnsembleAnalysis
using TanglesModel
using TanglesModel.StochasticToggle
using Plots
using Interpolations
using LinearAlgebra

#------------Simulation distributed over the input space---------------------------------------------
function adjust_ics(prob, i, _)
    ics = [i % 20, trunc(Int64, i / 20)]
    print(ics)
    generate_jump_problem(ics, (0.0, 100000.0), 5.0, 8.0, 1/160.0, 1/1200.0)
end

function animate_trajectories(ensemble_solution, n_frames, tspan, dot_colors, filename)
    timepoints = range(tspan[1], tspan[2], length=n_frames)
    solution_3d = permutedims(cat([
        cat(componentwise_vectors_timepoint(ensemble_solution, i)..., dims=2)
            for i in range(tspan[1], stop=tspan[2], length=n_frames)]...,
        dims=3),[3, 1, 2])
    summary_animation = @animate for i ∈ 1:n_frames
        plot(componentwise_vectors_timepoint(ensemble_solution, timepoints[i])...,
            xlims=[0,20],
            ylims=[0,20],
            seriestype=:scatter, markercolor=dot_colors, alpha=0.2, label="",
            title="Stochastic toggle switch trajectories")
        if i > 1
            plot!(
                solution_3d[max(1,i - 10):i,:,1],
                solution_3d[max(1,i-10):i,:,2],
                alpha=range(0.08, stop=0.0, length=(i - max(0, i - 11))),
                color="purple",
                label="")
        end
    end
    mp4(summary_animation,filename)
end

# Interesting regime: 5.0, 8.0, 1/20.0, 1/1200.0
prob = generate_jump_problem([5,5], (0.0,100000.0), 5.0, 8.0, 1/160.0, 1/1200.0)
ensemble_prob = EnsembleProblem(prob, prob_func=adjust_ics)
sol= solve(ensemble_prob, SSAStepper(), trajectories=400)
dot_colors = [(i % 20) <= trunc(Int64, i / 20) ? 1 : 2 for i in 1:400]

animate_trajectories(sol, 400, [0.0, 10000.0], dot_colors, "output/toggle/movies/stochastic_toggle_summary.mp4")

# Simulate starting in one of the basins
single_basin_prob = EnsembleProblem(generate_jump_problem([8,0], (0.0,100000.0), 5.0, 8.0, 1/160.0, 1/1200.0))
single_basin_sol = solve(single_basin_prob, SSAStepper(), trajectories=400)
animate_trajectories(single_basin_sol, 400, [0.0, 10000.0], ones(Int,400), "output/toggle/movies/stochastic_toggle_single_basin.mp4")
basin_curves = []
for n in 1.0:5.0
    basin_ensemble = EnsembleProblem(generate_jump_problem([8,0], (0.0,50000.0), n, 8.0, 1/160.0, 1/1200.0))
    basin_solution = solve(basin_ensemble, SSAStepper(), trajectories=5000)

    basin_percentage(t) = sum(componentwise_vectors_timepoint(basin_solution, t)[1] .> componentwise_vectors_timepoint(basin_solution, t)[2]) / 5000
    append!(basin_curves, vcat([basin_percentage(t) for t in range(0.0, stop=50000.0, length=200)]...))
end
basin_curves = convert(Matrix{Float64},reshape(basin_curves, (200,5)))

#====================================TANGLES=======================================================
==================================================================================================#
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
tangles_params = gen_sim_params(;sc_dependent=true, σ2_coeff=0.01)
#-----------------Simulate w/ TANGLES, with large intergenic distance (hopefully uncoupled)--------
large_bcs = LinearBoundaryParameters(100000, true, true)
large_gene_spacing = [
    CoupledGene(1 / 160.0, 1, 10000, 10500, (mRNA)->8.0 / (8.0 + mRNA[2]^5.0)),
    CoupledGene(1 / 160.0, 2, 90000, 90500, (mRNA)->8.0 / (8.0 + mRNA[1]^5.0))
]
low_coupling_problem = TanglesModel.build_problem(
    tangles_params, large_bcs, large_gene_spacing, 2, 10000.0, convert(Array{Int32,1},[5,5]))
tangles_soln = low_coupling_problem()
mRNA_vals = hcat([tangles_soln.u[i].u.mRNA for i in 1:length(tangles_soln)]...)
plot(transpose(mRNA_vals))

#-----------------Simulate w/ TANGLES, with small intergenic distance (introducing coupling)-------
small_bcs = LinearBoundaryParameters(7300.0 * 0.34, false, false)
n_trajectories = 100
n_samples = 200
max_t = 50000.0
convergent_results = zeros((5,n_samples,2, n_trajectories))
divergent_results = zeros((5,n_samples,2, n_trajectories))
for n = 1.0:5.0
    convergent_setup = TanglesModel.build_problem(
        tangles_params, small_bcs, [
            CoupledGene(1 / 160.0, 1, 2500.0 * 0.34, 3500.0 * 0.34, (mRNA)->8.0 / (8.0 + mRNA[2]^n)),
            CoupledGene(1 / 160.0, 2, 4800.0 * 0.34, 3800.0 * 0.34, (mRNA)->8.0 / (8.0 + mRNA[1]^n))
        ], 2, max_t, convert(Array{Int32, 1}, [8, 0])
    )
    divergent_setup = TanglesModel.build_problem(
        tangles_params, small_bcs, [
            CoupledGene(1 / 160.0, 1, 3500.0 * 0.34, 2500.0 * 0.34, (mRNA)->8.0 / (8.0 + mRNA[2]^n)),
            CoupledGene(1 / 160.0, 2, 3800.0 * 0.34, 4800.0 * 0.34, (mRNA)->8.0 / (8.0 + mRNA[1]^n))
        ], 2, max_t, convert(Array{Int32, 1}, [8, 0])
    )
    for i = 1:n_trajectories
        done = false
        sol = false
        while !done
            sol = convergent_setup()
            done = (sol.retcode == :Success)
        end
        timesteps = Interpolations.deduplicate_knots!(sol.t)
        mRNA = hcat([sol.u[i].u.mRNA for i in 1:length(sol)]...)
        interp_mRNA = ConstantInterpolation((1:2, timesteps), mRNA)
        convergent_results[trunc(Int,n),:,1,i] = interp_mRNA(1,range(0.0, stop=max_t, length=n_samples))
        convergent_results[trunc(Int,n),:,2,i] = interp_mRNA(2,range(0.0, stop=max_t, length=n_samples))
        print(".")
    end
    for i = 1:n_trajectories
        done = false
        sol = false
        while !done
            sol = divergent_setup()
            done = (sol.retcode == :Success)
        end
        timesteps = Interpolations.deduplicate_knots!(sol.t)
        mRNA = hcat([sol.u[i].u.mRNA for i in 1:length(sol)]...)
        interp_mRNA = ConstantInterpolation((1:2, timesteps), mRNA)
        divergent_results[trunc(Int,n),:,1,i] = interp_mRNA(1,range(0.0, stop=max_t, length=n_samples))
        divergent_results[trunc(Int,n),:,2,i] = interp_mRNA(2,range(0.0, stop=max_t, length=n_samples))
        print(".")
    end
end
convergent_summary = dropdims(sum(convergent_results[:,:,1,:] .> convergent_results[:,:,2,:],dims=3) / n_trajectories, dims=3)
divergent_summary = dropdims(sum(divergent_results[:,:,1,:] .> divergent_results[:,:,2,:],dims=3) / n_trajectories, dims=3)

#----------------Calculate the eigenvalues for a transition matrix)------------------
n_vals = 1.0:5.0
convergence_factor = zeros(length(n_vals))
for i =  1:length(n_vals)
    n = n_vals[i]
    for max_val = 10:10
        transition_matrix = zeros(((max_val + 1)^2, (max_val + 1)^2))
        ab_to_i(a,b) = 1 + a + (max_val + 1) * b
        # Fill in the transitions from each state in the stochastic transition matrix
        for a=0:max_val, b=0:max_val
            our_state = ab_to_i(a,b)
            # There are potentially four non-zero transitions.
            
            # Degradation terms
            if a > 0
                transition_matrix[our_state, ab_to_i(a-1,b)] = 1/1200.0 * a
            end
            if b > 0
                transition_matrix[our_state, ab_to_i(a,b-1)] = 1/1200.0 * b
            end
            # Creation terms
            if a < max_val
                transition_matrix[our_state, ab_to_i(a+1,b)] = 1/160.0 * 8.0 / (8.0 + b^n)
            end
            if b < max_val
                transition_matrix[our_state, ab_to_i(a,b+1)] = 1/160.0 * 8.0 / (8.0 + a^n)
            end
            # Central term
            transition_matrix[our_state, our_state] = 1 - sum(transition_matrix[our_state, :])
        end
        # Save the absolute value of the second highest eigenvalue 
        convergence_factor[i] = sort(abs.(eigvals(transition_matrix)))[end-1]
    end
end


#-------------------Plotting-------------------------
# Do all combined plotting
plot_tspan = range(0.0, stop=50000.0, length=200)
plot(0.45 .+ 0.55 .* hcat([c.^(1:50000) for c in convergence_factor]...),
    palette=palette(:viridis, 6), label=["n=1" "n=2" "n=3" "n=4" "n=5"], lw=3,
    xlims=(0,50000.0), ylims=(0.42, 1.0),
    xlabel="Time", ylabel="Fraction in basin", title="Eigenvalue analysis")
savefig("output/toggle/images/eigenvalue.pdf")
plot(plot_tspan,basin_curves,
    palette=palette(:viridis, 6), label=["n=1" "n=2" "n=3" "n=4" "n=5"], lw=3,
    xlims=(0,50000.0), ylims=(0.42, 1.0),
    xlabel="Time", ylabel="Fraction in basin", title="Stochastic")
savefig("output/toggle/images/stochastic.pdf")
plot(plot_tspan, transpose(convergent_summary),
    palette=palette(:viridis, 6), label=["n=1" "n=2" "n=3" "n=4" "n=5"], lw=3,
    xlims=(0,50000.0), ylims=(0.42, 1.0),
    xlabel="Time", ylabel="Fraction in basin", title="TANGLES: convergent orientation")
savefig("output/toggle/images/tangles_convergent.pdf")
plot(plot_tspan, transpose(divergent_summary),
    palette=palette(:viridis, 6), label=["n=1" "n=2" "n=3" "n=4" "n=5"], lw=3,
    xlims=(0,50000.0), ylims=(0.42, 1.0), title="TANGLES: divergent orientation")
savefig("output/toggle/images/tangles_divergent.pdf")
plot(0.45 .+ 0.55 .* hcat([c.^(1:50000) for c in convergence_factor]...),
    color=1, label="", lw=1,
    xlims=(0,50000.0), ylims=(0.42, 1.0))
plot!(0.45 .+ 0.55 .* (convergence_factor[1].^(1:50000)), color=1, label="Eigenvalue", lw=1)
plot!(plot_tspan,basin_curves[:,1], color=2, label="Stochastic", lw=3)
plot!(plot_tspan,basin_curves,
    color=2, label="", lw=3, xlims=(0,50000.0), ylims=(0.42, 1.0))
savefig("output/toggle/images/eigenvalue_stochastic_compare.pdf")

plot(plot_tspan,basin_curves,
    color=1, label="", lw=3, xlims=(0,50000.0), ylims=(0.42, 1.0))
plot!(plot_tspan,basin_curves[:,1], color=1, label="Stochastic", lw=3, xlims=(0,50000.0), ylims=(0.42, 1.0))
plot!(plot_tspan,transpose(convergent_summary)[:,1], color=2, label="Convergent", lw=3, xlims=(0,50000.0), ylims=(0.42, 1.0))
plot!(plot_tspan,transpose(divergent_summary)[:,1], color=3, label="Divergent", lw=3, xlims=(0,50000.0), ylims=(0.42, 1.0))
plot!(plot_tspan,transpose(convergent_summary),
    color=2, label="", lw=3, xlims=(0,50000.0), ylims=(0.42, 1.0))
plot!(plot_tspan,transpose(divergent_summary),
    color=3, label="", lw=3, xlims=(0,50000.0), ylims=(0.42, 1.0))
savefig("output/toggle/images/stochastic_convergent_divergent_compare.pdf")