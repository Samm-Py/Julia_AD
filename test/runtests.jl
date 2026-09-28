using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "OTISensitivity.jl"))
using .OTISensitivity

@testset "Dual-number multiplication representation" begin
    a, b = [2., 3.], [7., 11.]
    @test oti_matrix(a)*oti_matrix(b) == oti_matrix([14., 43.])
    N = oti_matrix([0., 1.])
    @test N^2 == zeros(2, 2)
    @test_throws ArgumentError oti_matrix(Float64[])
end

@testset "Scalar exponential: value, first derivative and zero timestep" begin
    for dt in (0., 0.7), p in 0:1
        r = oti_response(reshape([-0.3], 1, 1), zeros(1, 1), [0.], [2.], [0.], dt; order=p)
        for j in 0:p
            exact = 2exp(-0.3dt)*(-0.3)^j
            @test r.derivatives[1, j+1] ≈ exact rtol=1e-11 atol=1e-13
            @test r.coefficients[1, j+1] ≈ exact rtol=1e-11 atol=1e-13
        end
    end
end

@testset "Thermal response and analytical first derivative" begin
    (; A, B, e, x, u) = thermal_example()
    # Independent nominal reference: diagonalize A and integrate each mode.
    λ, V = eigen(Symmetric(A))
    for forcing in (e, collect(range(-0.2, 0.3; length=10))), dt in (0., 0.5, 50., 100.)
        f = B*u + forcing
        integral = [abs(v) < 1e-14 ? dt : expm1(v*dt)/v for v in λ]
        nominal = V*(exp.(λ*dt).*(V'*x) + integral.*(V'*f))
        @test response(A, B, forcing, x, u, dt) ≈ nominal rtol=1e-11 atol=1e-10
        for p in (0, 1)
            r = oti_response(A, B, forcing, x, u, dt; order=p)
            @test r.derivatives[:, 1] ≈ nominal rtol=1e-11 atol=1e-10
            p == 1 && @test r.derivatives[:, 2] ≈ A*nominal + f rtol=1e-9 atol=1e-11
        end
    end
    @test_throws ArgumentError oti_response(A, B, e, x, u, 50.; order=-1)
end

@testset "Zero dynamics, nonzero forcing" begin
    A, B, e, x, u = zeros(2, 2), ones(2, 1), [2., 3.], [4., 5.], [7.]
    r = oti_response(A, B, e, x, u, 2.; order=1)
    @test r.derivatives[:, 1] ≈ x + 2(B*u+e)
    @test r.derivatives[:, 2] ≈ B*u+e
end

include(joinpath(@__DIR__, "..", "src", "SensitivityComparison.jl"))
using .SensitivityComparison: ctse_response, analytical_response

@testset "CTSE first derivatives and complex forcing integration" begin
    (; A, B, e, x, u) = thermal_example()
    for dt in (0., 0.5, 50., 100.), forcing in (e, collect(range(-0.2,0.3;length=10)))
        exact = analytical_response(A,B,forcing,x,u,dt; order=1)
        for h in (1e-8, 1e-20, 1e-30)
            r = ctse_response(A,B,forcing,x,u,dt; h)
            @test r.value ≈ exact[:,1] rtol=1e-11 atol=1e-10
            @test r.derivative ≈ exact[:,2] rtol=1e-10 atol=1e-11
        end
    end
    @test_throws ArgumentError ctse_response(A,B,e,x,u,50.;h=0.)
end
