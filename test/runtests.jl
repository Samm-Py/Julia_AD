using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "OTISensitivity.jl"))
using .OTISensitivity

@testset "OTI multiplication representation" begin
    a, b = [2., 3., 5.], [7., 11., 13.]
    @test oti_matrix(a)*oti_matrix(b) == oti_matrix([14., 43., 94.])
    N = oti_matrix([0., 1., 0.])
    @test N^3 == zeros(3, 3)
    @test_throws ArgumentError oti_matrix(Float64[])
end

@testset "Scalar exponential: factorial and zero timestep" begin
    for dt in (0., 0.7), p in 0:5
        r = oti_response(reshape([-0.3], 1, 1), zeros(1, 1), [0.], [2.], [0.], dt; order=p)
        for j in 0:p
            exact = 2exp(-0.3dt)*(-0.3)^j
            @test r.derivatives[1, j+1] ≈ exact rtol=1e-11 atol=1e-13
            @test r.coefficients[1, j+1] ≈ exact/factorial(j) rtol=1e-11 atol=1e-13
        end
    end
end

@testset "Thermal response and analytical derivatives" begin
    (; A, B, e, x, u) = thermal_example()
    # Independent nominal reference: diagonalize A and integrate each mode.
    λ, V = eigen(Symmetric(A))
    for forcing in (e, collect(range(-0.2, 0.3; length=10))), dt in (0., 0.5, 50., 100.)
        f = B*u + forcing
        integral = [abs(v) < 1e-14 ? dt : expm1(v*dt)/v for v in λ]
        nominal = V*(exp.(λ*dt).*(V'*x) + integral.*(V'*f))
        @test response(A, B, forcing, x, u, dt) ≈ nominal rtol=1e-11 atol=1e-10
        for p in (0, 1, 2, 4)
            r = oti_response(A, B, forcing, x, u, dt; order=p)
            @test r.derivatives[:, 1] ≈ nominal rtol=1e-11 atol=1e-10
            exact = A*nominal + f
            for j in 1:p
                @test r.derivatives[:, j+1] ≈ exact rtol=1e-9 atol=1e-11
                exact = A*exact
            end
        end
    end
    @test_throws ArgumentError oti_response(A, B, e, x, u, 50.; order=-1)
end

@testset "Zero dynamics, nonzero forcing" begin
    A, B, e, x, u = zeros(2, 2), ones(2, 1), [2., 3.], [4., 5.], [7.]
    r = oti_response(A, B, e, x, u, 2.; order=3)
    @test r.derivatives[:, 1] ≈ x + 2(B*u+e)
    @test r.derivatives[:, 2] ≈ B*u+e
    @test r.derivatives[:, 3:4] ≈ zeros(2, 2) atol=1e-12
end

include(joinpath(@__DIR__, "..", "src", "SensitivityComparison.jl"))
using .SensitivityComparison: ctse_response, finite_difference, fd_step, analytical_response

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

@testset "Finite differences: polynomial exactness and thermal agreement" begin
    for j in 1:4
        exact = factorial(6)/factorial(6-j)*2.0^(6-j)
        @test finite_difference(t -> [t^6],2.,j;h=0.1)[1] ≈ exact rtol=1e-9
    end
    (; A,B,e,x,u,dt) = thermal_example()
    exact = analytical_response(A,B,e,x,u,dt;order=4)
    f(t) = response(A,B,e,x,u,t)
    for j in 1:4
        @test finite_difference(f,dt,j;h=fd_step(dt,j)) ≈ exact[:,j+1] rtol=1e-4 atol=1e-10
    end
    @test_throws ArgumentError finite_difference(f,dt,0;h=1.)
    @test_throws ArgumentError finite_difference(f,dt,1;h=-1.)
end
