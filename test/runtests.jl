import TanglesModel as TM
using Test

linear_bcs = TM.LinearBoundaryParameters(5000, false, false)
circular_bcs = TM.CircularBoundaryParameters(5000)

empty_u = TM.TanglesArray([0], [], [], [], [])
single_u = TM.TanglesArray([100, 1, 0, 0], [0], [1], [500], [1])
double_u = TM.TanglesArray([100, 1, 0, 300, 2, 200, 0], [0], [1], [500, 500], [1, 1])

@testset "Internals" begin
    @testset "supercoiling_density!" begin
        # NOTE: we use ω0 = 1.0 in these cases for simple, transparent test cases
        @test TM.supercoiling_density(empty_u, 1, linear_bcs, 1.0) ≈ 0
        @test TM.supercoiling_density(empty_u, 1, circular_bcs, 1.0) ≈ 0

        @test TM.supercoiling_density(single_u, 1, linear_bcs, 1.0) ≈ -1.0 / 100.0
        @test TM.supercoiling_density(single_u, 2, linear_bcs, 1.0) ≈ 1.0 / 4900.0
        # Circular BCS has no supercoiling density
        @test TM.supercoiling_density(single_u, 1, circular_bcs, 1.0) ≈ 0
        @test TM.supercoiling_density(single_u, 2, circular_bcs, 1.0) ≈ 0

        # Test the double polymerase supercoiling density
        @test TM.supercoiling_density(double_u, 1, linear_bcs, 1.0) ≈ -1.0 / 100.0
        @test TM.supercoiling_density(double_u, 2, linear_bcs, 1.0) ≈ -1.0 / 200.0
        @test TM.supercoiling_density(double_u, 3, linear_bcs, 1.0) ≈ 2.0 / 4700.0
        # Check free BCs
        @testset "free_bcs" begin
            @test TM.supercoiling_density(double_u, 1,
                TM.LinearBoundaryParameters(5000, true, false), 1.0) ≈ 0
            @test TM.supercoiling_density(double_u, 2,
                TM.LinearBoundaryParameters(5000, true, false), 1.0) ≈ -1.0 / 200.0
            @test TM.supercoiling_density(double_u, 3,
                TM.LinearBoundaryParameters(5000, true, false), 1.0) ≈ 2.0 / 4700.0

            @test TM.supercoiling_density(double_u, 1,
                TM.LinearBoundaryParameters(5000, false, true), 1.0) ≈ -1.0 / 100.0
            @test TM.supercoiling_density(double_u, 2,
                TM.LinearBoundaryParameters(5000, false, true), 1.0) ≈ -1.0 / 200.0
            @test TM.supercoiling_density(double_u, 3,
                TM.LinearBoundaryParameters(5000, false, true), 1.0) ≈ 0

            @test TM.supercoiling_density(double_u, 1,
                TM.LinearBoundaryParameters(5000, true, true), 1.0) ≈ 0
            @test TM.supercoiling_density(double_u, 2,
                TM.LinearBoundaryParameters(5000, true, true), 1.0) ≈ -1.0 / 200.0
            @test TM.supercoiling_density(double_u, 3,
                TM.LinearBoundaryParameters(5000, true, true), 1.0) ≈ 0
        end # testset Internals/supercoiling_density!/free_bcs
        @test TM.supercoiling_density(double_u, 1, circular_bcs, 1.0) == TM.supercoiling_density(
                    double_u, 3, circular_bcs, 1.0)
        @test TM.supercoiling_density(double_u, 1, circular_bcs, 1.0) ≈ 1.0 / 4800.0
        @test TM.supercoiling_density(double_u, 2, circular_bcs, 1.0) ≈ -1.0 / 200.0
    end # testset Internals/supercoiling_density!
    @testset "torque_response" begin
        internal_params = TM.InternalParameters(TM.DEFAULT_SIM_PARAMS)
        test_params = TM.TanglesParams(internal_params, linear_bcs)

        τ_0 = internal_params.τ_0
        τ_s = internal_params.τ_s
        τ_p = internal_params.τ_p
        σ_s = internal_params.σ_s
        σ_p = internal_params.σ_p

        @test TM.torque_response(0.0, test_params) ≈ 0.0
        @test TM.torque_response(σ_s / 2.0, test_params) ≈ τ_s * σ_s / 2.0
        @test TM.torque_response(-σ_s / 2.0, test_params) ≈ -τ_s * σ_s / 2.0
        @testset "constant_torque_regime" for x in 0.0:0.1:1.0
            σ = (x * σ_s + (1 - x) * σ_p)
            @test TM.torque_response(σ, test_params)  ≈ τ_0
            @test TM.torque_response(-σ, test_params) ≈ -τ_0
        end
        @test TM.torque_response(σ_p + .1, test_params) ≈ τ_p * (σ_p + 0.1)
        @test TM.torque_response(-σ_p - .1, test_params) ≈ -τ_p * (σ_p + 0.1)
    end # testset Internals/torque_response
    @testset "interp_twist" begin
        @testset "empty u's" begin
            @test all(TM.interp_twist(50.0, empty_u, linear_bcs, 1.0) .≈ (1, 0.0))
            @test all(TM.interp_twist(50.0, empty_u, circular_bcs, 1.0) .≈ (1, 0.0))
        end
        @testset "single linear u's" begin
        @test all(TM.interp_twist(50.0, single_u, linear_bcs, 1.0) .≈ (1, 0.5))
        @test all(TM.interp_twist(2550.0, single_u, linear_bcs, 1.0) .≈ (2, 0.5))
        end
        @testset "linear free ends" begin
            @test all(TM.interp_twist(50.0, single_u, TM.LinearBoundaryParameters(5000, true, false), 1.0) .≈ (1, 1.0))
            @test all(TM.interp_twist(50.0, single_u, TM.LinearBoundaryParameters(5000, false, true), 1.0) .≈ (1, 0.5))
            @test all(TM.interp_twist(2550.0, single_u, TM.LinearBoundaryParameters(5000, false, true), 1.0) .≈ (2, 1.0))
            @test all(TM.interp_twist(2550.0, single_u, TM.LinearBoundaryParameters(5000, true, false), 1.0) .≈ (2, 0.5))
            @test all(TM.interp_twist(50.0, single_u, TM.LinearBoundaryParameters(5000, true, true), 1.0) .≈ (1, 1.0))
            @test all(TM.interp_twist(2550.0, single_u, TM.LinearBoundaryParameters(5000, true, true), 1.0) .≈ (2, 1.0))
        end
        @testset "circular bcs" begin
            @test all(TM.interp_twist(50.0, double_u, circular_bcs, 1.0) .≈ (1, 1 + 50.0 / 4800.0))
            @test all(TM.interp_twist(200.0, double_u, circular_bcs, 1.0) .≈ (2, 1.5))
            @test all(TM.interp_twist(400.0, double_u, circular_bcs, 1.0) .≈ (3, 2 - 100.0 / 4800.0))
        end

    end # testset Internals/interp_twist
end #testset Internals/

@testset "Integration" begin
    problem      = TM.build_problem(TM.DEFAULT_SIM_PARAMS, TM.LinearBoundaryParameters(5000, false, false), [TM.Gene(1/120.0, 1, 100, 1000), TM.Gene(1/120.0, 2, 3000, 2000)], 2, 1000.0)
    circ_problem = TM.build_problem(TM.DEFAULT_SIM_PARAMS, TM.CircularBoundaryParameters(5000), [TM.Gene(1/120.0, 1, 100, 1000), TM.Gene(1/120.0, 2, 3000, 2000)], 2, 1000.0)
    problem()
    circ_problem()
end
