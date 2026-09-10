using Plots, LinearAlgebra, ExponentialUtilities

# Pulled from https://github.com/CMU-EMIT-Lab/ADDOPT.jl/blob/main/src/utilities/voxelized_conduction.jl
# Lines 1 - 13
function discretize_linear_dynamics(A, B, e, dt)
    nx, nu = size(B)
    H = vcat([A B I(nx)],
        zeros(nu + nx, nx + nu + nx)) # Dynamics matrix for combined system of x and u
    H .*= dt
    G = exponential!(H) # State transition matrix for combined system

    Ad = G[1:nx, 1:nx]
    Bd = G[1:nx, (nx+1):(nx+nu)]
    ed = G[1:nx, (nx+nu+1):(2nx+nu)] * e

    return Ad, Bd, ed
end

## Basic 1D Thermal Simulation Setup in State Space ##
# The units in this simulation are more for reference
# This problem isn't representative of anything in the real world. 

n = 10 # Number of voxels
T∞ = 300 # Starting temperature (K)
l = 10e-3 # Voxel size (m)

# Thermal Properties (Made up values for reference)
cₚ = 10000 # J/kgK
ρ = 5000 # kg/m³
k = 200 # W/mK
α = k / ρ / cₚ

P = 500 # Input power (Watts)

# Setting up matrices/vectors for an equation in the form of ẋ = Ax + Bu + e

xₖ₊₁ = zeros(n)

A = zeros(n, n)
xₖ = ones(n) * T∞
xₖ[end] = 3 * T∞ # High temperature at base initially

B = zeros(n, 1)
u = [P]

e = zeros(n) # No exogenous inputs

for i in 1:n 
    for j in 1:n 
        if i == j # Main Diagonal
            A[i, j] -= i == 1 || i == n ? α / l^2 : 2α / l^2 # No heat transfer out of boundaries
        end

        if j == i + 1 || j == i - 1 # Off diagonals
            A[i, j] += α / l^2
        end
    end
end

B[1] = 1 / ρ / cₚ / l^3

Δt = 50 # seconds

Ad, Bd, ed = discretize_linear_dynamics(A, B, e, Δt)

# Use of the discrete state space equation
xₖ₊₁ = Ad * xₖ + Bd * u + e

p1 = scatter(xₖ, xlabel="Voxel Number", ylabel="Temperature", label="Initial State")
scatter!(xₖ₊₁, label="Final State", title="1D Conduction Matrix Exponential Example", legend=:top)

display(p1)