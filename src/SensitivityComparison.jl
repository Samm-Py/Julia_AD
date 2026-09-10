module SensitivityComparison
using LinearAlgebra
include("OTISensitivity.jl")
using .OTISensitivity
export thermal_example, response, oti_response, ctse_response, finite_difference,
       analytical_response, write_method_csv, fd_step

"""Complex Taylor series expansion: nominal value and first derivative only."""
function ctse_response(A, B, e, x, u, dt::Real; h::Real=1e-20)
    isfinite(h) && h > 0 || throw(ArgumentError("h must be finite and positive"))
    z = response(A, B, e, x, u, complex(dt, h))
    return (; value=real.(z), derivative=imag.(z)/h)
end

"""Seven-point centered finite difference, using only real function values."""
function finite_difference(f, t::Real, order::Integer; h::Real)
    order in 1:4 || throw(ArgumentError("supported FD orders are 1 through 4"))
    isfinite(h) && h > 0 || throw(ArgumentError("h must be finite and positive"))
    nodes = -3:3
    # Exact rational moment equations: sum(w_k*k^q) = δ(q,order)*order!.
    moments = [(big(k)//big(1))^q for q in 0:6, k in nodes]
    rhs = zeros(Rational{BigInt}, 7)
    rhs[order+1] = factorial(order)
    weights = Float64.(moments \ rhs)
    return sum(w*f(t+k*h) for (w, k) in zip(weights, nodes))/h^order
end

# Balance O(h^accuracy) truncation against O(eps/h^order) roundoff.
function fd_step(t, order)
    order in 1:4 || throw(ArgumentError("supported FD orders are 1 through 4"))
    accuracy = order <= 2 ? 6 : 4
    return max(abs(t), 1.0)*eps(Float64)^(1/(order+accuracy))
end

"""Independent modal reference for the symmetric thermal A, then ODE recurrence."""
function analytical_response(A, B, e, x, u, dt; order=4)
    issymmetric(A) || throw(ArgumentError("modal reference requires symmetric A"))
    λ, V = eigen(Symmetric(A))
    integral = [abs(v) < 1e-14 ? dt : expm1(v*dt)/v for v in λ]
    nominal = V*(exp.(λ*dt).*(V'*x) + integral.*(V'*(B*u+e)))
    derivatives = zeros(length(x), order+1)
    derivatives[:, 1] = nominal
    for j in 1:order
        derivatives[:, j+1] = A*derivatives[:, j] + (j == 1 ? B*u+e : zero(e))
    end
    return derivatives
end

function write_method_csv(path, dt, derivatives)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "dt,voxel,", join(["d$j" for j in 0:size(derivatives,2)-1], ','))
        for i in axes(derivatives, 1)
            println(io, join(vcat(dt, i, derivatives[i, :]), ','))
        end
    end
end
end
