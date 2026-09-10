# Reproduce the table and step-size convergence figure in Section 1.3.1.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots, Printf

function main()
    out = joinpath(@__DIR__, "..", "results", "worked_example")
    mkpath(out)
    setprecision(BigFloat, 256) do
        exact_value = exp(BigFloat(1))
        exact_derivative = 2exact_value
        function evaluate(h)
            z = exp(complex(1.0, h)^2)
            value, derivative = real(z), imag(z)/h
            value_error = abs(BigFloat(value)-exact_value)
            derivative_error = abs(BigFloat(derivative)-exact_derivative)
            return (; h, value, derivative, value_error, derivative_error)
        end
        selected = evaluate.([1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-12, 1e-20])
        sweep = evaluate.(10.0 .^ range(-20, -1; length=240))
        # Check the leading truncation term independently in high precision.
        hb = BigFloat("1e-4")
        db = imag(exp(complex(BigFloat(1), hb)^2))/hb
        predicted_error = -(10exact_value/3)*hb^2
        @assert isapprox(db-exact_derivative, predicted_error; rtol=1e-7)
        @assert selected[end].derivative_error < 1e-14
        for (name, rows) in (("ctse_table.csv", selected), ("ctse_convergence.csv", sweep))
            open(joinpath(out, name), "w") do io
                println(io, "h,value,derivative,absolute_value_error,absolute_derivative_error")
                for r in rows
                    println(io, join((r.h,r.value,r.derivative,r.value_error,r.derivative_error), ','))
                end
            end
        end
        latex_sci(x) = begin
            mantissa, power = split(@sprintf("%.3e", Float64(x)), 'e')
            return "\$" * mantissa * raw"\times10^{" * string(parse(Int,power)) * "}\$"
        end
        open(joinpath(@__DIR__, "ctse_worked_results.tex"), "w") do io
            println(io, raw"\begin{tabular}{rrrr}\toprule")
            println(io, raw"Step $h$ & $\operatorname{Re}\widetilde f$ & $\operatorname{Im}\widetilde f/h$ & Absolute derivative error \\ \midrule")
            for r in selected
                @printf(io, "%s & %.12f & %.12f & %s ", latex_sci(r.h),r.value,r.derivative,latex_sci(r.derivative_error))
                println(io, repeat(string(Char(92)),2))
            end
            println(io,raw"\bottomrule\end{tabular}")
        end
        errors = Float64.([r.derivative_error for r in sweep])
        p = plot([r.h for r in sweep], errors;
            xscale=:log10,yscale=:log10,linewidth=2,
            label="CTSE (Float64)",xlabel="Imaginary step h",
            ylabel="Absolute first-derivative error",
            title="CTSE for f(t) = exp(t²) at t = 1",size=(950,540),
            legend=:topleft,left_margin=6Plots.mm,bottom_margin=4Plots.mm,
            ylims=(1e-17,1.0))
        theory_h = 10.0 .^ range(-9,-1;length=100)
        plot!(p,theory_h,Float64(10exact_value/3).*theory_h.^2;
            label="Leading truncation error: (10e/3)h²",linestyle=:dash,linewidth=2)
        scatter!(p,[r.h for r in selected],Float64.([r.derivative_error for r in selected]);
            label="Table entries",markersize=4)
        savefig(p,joinpath(out,"ctse_convergence.png"))
        println("Generated CTSE worked-example table, CSV data, and convergence plot.")
        println("Derivative reference (256-bit): ",exact_derivative)
        println("Absolute derivative error at h=1e-20: ",selected[end].derivative_error)
    end
end
main()
