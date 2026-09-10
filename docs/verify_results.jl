using LinearAlgebra, Printf
include(joinpath(@__DIR__, "..", "src", "OTISensitivity.jl"))
using .OTISensitivity

function main()
    (; A, B, e, x, u, dt) = thermal_example()
    r = oti_response(A, B, e, x, u, dt; order=4)
    λ, V = eigen(Symmetric(A))
    integral = [abs(v) < 1e-14 ? dt : expm1(v*dt)/v for v in λ]
    y = V*(exp.(λ*dt).*(V'*x) + integral.*(V'*(B*u+e)))
    errors = Float64[]
    for j in 0:4
        push!(errors, norm(r.derivatives[:, j+1]-y, Inf))
        @assert isapprox(r.derivatives[:, j+1], y; rtol=1e-9, atol=1e-10)
        y = j == 0 ? A*y+B*u+e : A*y
    end
    open(joinpath(@__DIR__, "verification_results.tex"), "w") do io
        println(io, "At the example duration of 50 seconds, the maximum absolute errors")
        println(io, "over all ten voxels against the independent modal reference are:")
        println(io, raw"\begin{center}\begin{tabular}{cc}\toprule")
        println(io, raw"Derivative order $j$ & Maximum absolute error (K/s$^j$) \\ \midrule")
        for j in 0:4
            @printf(io, "%d & ", j)
            mantissa, exponent = split(@sprintf("%.3e", errors[j+1]), 'e')
            println(io, "\$", mantissa, raw"\times 10^{", parse(Int, exponent), raw"}$ ", repeat(string(Char(92)), 2))
        end
        println(io, raw"\bottomrule\end{tabular}\end{center}")
    end
    println("Maximum absolute errors by order 0–4: ", errors)
end
main()
