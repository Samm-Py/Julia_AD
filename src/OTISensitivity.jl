module OTISensitivity
using LinearAlgebra
using ExponentialUtilities: exponential!
export oti_matrix, discretize_linear_dynamics, response, oti_response, thermal_example

"""Real multiplication representation of OTI Taylor coefficients [a₀, …, aₚ]."""
function oti_matrix(coefficients::AbstractVector{<:Real})
    isempty(coefficients) && throw(ArgumentError("at least one coefficient is required"))
    m = length(coefficients)
    C = zeros(eltype(coefficients), m, m)
    for j in 1:m, i in j:m
        C[i, j] = coefficients[i-j+1]
    end
    return C
end

function augmented_dynamics(A, B)
    n, nu = size(B)
    size(A) == (n, n) || throw(DimensionMismatch("A and B must have matching state dimensions"))
    return [A B Matrix{Float64}(I, n, n); zeros(nu+n, 2n+nu)]
end

function discretize_linear_dynamics(A, B, e, dt::Number)
    n, nu = size(B)
    length(e) == n || throw(DimensionMismatch("e must have one entry per state"))
    G = exponential!(augmented_dynamics(A, B) * dt)
    return G[1:n, 1:n], G[1:n, n+1:n+nu], G[1:n, n+nu+1:2n+nu] * e
end

"""One-step response with constant input and forcing."""
function response(A, B, e, x, u, dt)
    Ad, Bd, ed = discretize_linear_dynamics(A, B, e, dt)
    return Ad*x + Bd*u + ed
end

"""
    oti_response(A, B, e, x, u, dt; order=2)

Differentiate the one-step response with respect to dt, holding all other inputs
fixed. Column j+1 of `coefficients` is the jth Taylor coefficient; column j+1
of `derivatives` is the jth derivative (column 1 is the nominal response).
Uses a dense real block exponential and one nilpotent OTI direction.
"""
function oti_response(A, B, e, x, u, dt::Real; order::Integer=2)
    order >= 0 || throw(ArgumentError("order must be nonnegative"))
    n, nu = size(B)
    length(x) == length(e) == n || throw(DimensionMismatch("x and e must match A"))
    length(u) == nu || throw(DimensionMismatch("u must match B"))
    H = augmented_dynamics(A, B)
    seed = zeros(promote_type(Float64, typeof(dt)), order+1)
    seed[1] = dt
    order > 0 && (seed[2] = 1)
    # Coefficient-major ordering: nominal states, first-order states, etc.
    G = exponential!(kron(oti_matrix(seed), H))
    m = size(H, 1)
    z = vcat(x, u, e)
    coefficients = hcat([G[j*m+1:j*m+n, 1:m] * z for j in 0:order]...)
    derivatives = copy(coefficients)
    factor = one(eltype(derivatives))
    for j in 1:order
        factor *= j
        derivatives[:, j+1] .*= factor
    end
    return (; coefficients, derivatives)
end

"""Original ten-voxel conduction example, without plotting or global state."""
function thermal_example()
    n = 10
    T∞, l, cₚ, ρ, k = 300.0, 10e-3, 10000.0, 5000.0, 200.0
    rate = k / ρ / cₚ / l^2
    A = Matrix(Tridiagonal(fill(rate, n-1), fill(-2rate, n), fill(rate, n-1)))
    A[1, 1] = A[end, end] = -rate
    B = zeros(n, 1)
    B[1] = 1 / ρ / cₚ / l^3
    x = fill(T∞, n)
    x[end] = 3T∞
    return (; A, B, e=zeros(n), x, u=[500.0], dt=50.0)
end
end
