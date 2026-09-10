# Reproduce the centered finite-difference error plot in Section 1.2.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots

function main()
    out = joinpath(@__DIR__, "..", "results", "worked_example")
    mkpath(out)
    setprecision(BigFloat, 256) do
        exact_value = exp(BigFloat(1))
        exact_derivative = 2exact_value
        steps = 10.0 .^ range(-16, -1; length=300)
        rows = map(steps) do h
            hb = BigFloat(h)
            d64 = (exp((1.0+h)^2)-exp((1.0-h)^2))/(2h)
            # Same mathematical stencil in high precision isolates truncation.
            db = (exp((BigFloat(1)+hb)^2)-exp((BigFloat(1)-hb)^2))/(2hb)
            total_error = abs(BigFloat(d64)-exact_derivative)
            truncation_error = abs(db-exact_derivative)
            rounding_error = abs(BigFloat(d64)-db)
            roundoff_scale = BigFloat(eps(Float64))*exact_value/hb
            leading_truncation = (10exact_value/3)*hb^2
            model = leading_truncation+roundoff_scale
            (; h,d64,total_error,truncation_error,rounding_error,roundoff_scale,model)
        end
        # The centered formula has leading error f'''(1)h^2/6 = (10e/3)h^2.
        hcheck = BigFloat("1e-4")
        dcheck = (exp((1+hcheck)^2)-exp((1-hcheck)^2))/(2hcheck)
        @assert isapprox(dcheck-exact_derivative,(10exact_value/3)*hcheck^2;rtol=1e-7)
        errors = Float64.([r.total_error for r in rows])
        best = argmin(errors)
        @assert errors[best] < errors[1] && errors[best] < errors[end]
        open(joinpath(out,"fd_convergence.csv"),"w") do io
            println(io,"h,derivative,total_absolute_error,truncation_error,rounding_error,roundoff_scale,heuristic_error_model")
            for r in rows
                println(io,join(Tuple(r),','))
            end
        end
        p=plot(steps,errors;xscale=:log10,yscale=:log10,
            label="Observed error (Float64)",linewidth=1.8,
            xlabel="Real step h",ylabel="Absolute first-derivative error",
            title="Centered finite differences: f(t) = exp(t²), t = 1",
            size=(1000,600),legend=:topright,
            left_margin=6Plots.mm,bottom_margin=4Plots.mm,ylims=(1e-13,100.0),
            xticks=10.0 .^ (-16:3:-1),yticks=10.0 .^ (-12:3:0))
        # Restrict the truncation curve to the visible scale; CSV retains all values.
        visible = [r for r in rows if r.truncation_error >= 1e-13]
        plot!(p,[r.h for r in visible],Float64.([r.truncation_error for r in visible]);
            label="Truncation error (~ h²)",linewidth=2,linestyle=:dash)
        plot!(p,steps,Float64.([r.roundoff_scale for r in rows]);
            label="Roundoff scale: eps*e/h (heuristic)",linewidth=2,linestyle=:dot)
        plot!(p,steps,Float64.([r.model for r in rows]);
            label="Sum of leading error scales",linewidth=1.7,linestyle=:dashdot,color=:black)
        scatter!(p,[steps[best]],[errors[best]];label="Smallest sampled error",markersize=5)
        savefig(p,joinpath(out,"fd_convergence.png"))
        println("Generated finite-difference convergence plot and CSV data.")
        println("Smallest sampled error: ",errors[best]," at h = ",steps[best])
    end
end
main()
