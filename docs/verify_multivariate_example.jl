# Verify Appendix A without changing the single-parameter thermal implementations.
using LinearAlgebra, Test
using ExponentialUtilities: exponential!

function main()
    # Total-degree-two monomials: 1, ε₁, ε₂, ε₁², ε₁ε₂, ε₂².
    basis = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    q = length(basis)
    function multiplication_matrix(direction)
        N = zeros(q,q)
        for (col,alpha) in enumerate(basis)
            target = (alpha[1]+direction[1],alpha[2]+direction[2])
            row = findfirst(==(target),basis)
            isnothing(row) || (N[row,col] = 1)
        end
        N
    end
    N1 = multiplication_matrix((1,0))
    N2 = multiplication_matrix((0,1))
    I6 = Matrix{Float64}(I,q,q)
    powers = [I6,N1,N2,N1^2,N1*N2,N2^2]
    represent(a) = sum(a[i]*powers[i] for i in 1:q)
    @testset "Multivariate OTI appendix" begin
        @test N1[2,1] == N1[4,2] == N1[5,3] == 1
        @test N2[3,1] == N2[5,2] == N2[6,3] == 1
        @test count(!iszero,N1) == count(!iszero,N2) == 3
        @test N1*N2 == N2*N1
        for k in 0:3
            @test N1^k*N2^(3-k) == zeros(q,q)
        end
        a,b = [2.,3.,5.,7.,11.,13.],[17.,19.,23.,29.,31.,37.]
        product = zeros(q)
        for (i,alpha) in enumerate(basis), (j,beta) in enumerate(basis)
            target = (alpha[1]+beta[1],alpha[2]+beta[2])
            row = findfirst(==(target),basis)
            isnothing(row) || (product[row] += a[i]*b[j])
        end
        @test represent(a)*represent(b) == represent(product)
        @test represent(a)[:,1] == a
        X,Y = I6+N1,N2
        W = exponential!(X*X+Y)
        expected = exp(1.0)*[1.,2.,1.,3.,2.,0.5]
        @test W[:,1] ≈ expected rtol=1e-12
        @test W ≈ represent(expected) rtol=1e-12
        @test W[2:3,1] ≈ exp(1.0)*[2.,1.] rtol=1e-12
        hessian = [2W[4,1] W[5,1]; W[5,1] 2W[6,1]]
        @test hessian ≈ exp(1.0)*[6. 2.;2. 1.] rtol=1e-12
        # At (0,0), the expansion of exp(xA+yB) through degree two is
        # I + xA + yB + x² A²/2 + xy(AB+BA)/2 + y² B²/2.
        A,B = [0. 1.;0. 0.],[0. 0.;1. 0.]
        @test A*B != B*A
        lifted = kron(N1,A)+kron(N2,B)
        E = exponential!(lifted)
        expected_blocks = [Matrix{Float64}(I,2,2),A,B,A*A/2,(A*B+B*A)/2,B*B/2]
        for i in 1:q
            @test E[2i-1:2i,1:2] ≈ expected_blocks[i] atol=1e-13 rtol=1e-12
        end
        @test E ≈ sum(kron(powers[i],expected_blocks[i]) for i in 1:q) atol=1e-13 rtol=1e-12
    end
end
main()
