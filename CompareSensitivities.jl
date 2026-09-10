# Run headlessly by default; the unchanged original still executes its display call.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots, LinearAlgebra, Printf
include(joinpath(@__DIR__, "CTSEExample.jl"))
include(joinpath(@__DIR__, "CRExample.jl"))

module OriginalBaseline
include(joinpath(@__DIR__, "MatrixExponentialTest.jl"))
end

function compare_sensitivities()
    model = thermal_example()
    (; A, B, e, x, u, dt) = model
    output = joinpath(@__DIR__, "results")
    mkpath(output)
    cs = run_ctse(; model)
    cr = run_cr(; model, order=4).derivatives
    original = OriginalBaseline.xₖ₊₁
    @assert A == OriginalBaseline.A && B == OriginalBaseline.B
    @assert x == OriginalBaseline.xₖ && u == OriginalBaseline.u && e == OriginalBaseline.e
    @assert dt == OriginalBaseline.Δt
    ref = analytical_response(A, B, e, x, u, dt; order=4)
    f(t) = response(A, B, e, x, u, t)
    fd = hcat(f(dt), [finite_difference(f, dt, j; h=fd_step(dt,j)) for j in 1:4]...)
    write_method_csv(joinpath(output, "finite_difference.csv"), dt, fd)
    write_method_csv(joinpath(output, "analytical.csv"), dt, ref)
    write_method_csv(joinpath(output, "original.csv"), dt, reshape(original, :, 1))
    @assert isapprox(original, ref[:,1]; rtol=1e-11, atol=1e-10)
    @assert isapprox(cs.value, original; rtol=1e-11, atol=1e-10)
    @assert isapprox(cr[:,1], original; rtol=1e-11, atol=1e-10)
    @assert isapprox(cs.derivative, ref[:,2]; rtol=1e-10, atol=1e-11)
    @assert isapprox(cs.derivative, cr[:,2]; rtol=1e-10, atol=1e-11)
    for j in 1:4
        @assert isapprox(cr[:,j+1], ref[:,j+1]; rtol=1e-9, atol=1e-11)
        # FD is a numerical approximation, increasingly cancellation-limited at high orders.
        @assert isapprox(fd[:,j+1], ref[:,j+1]; rtol=1e-4, atol=1e-10)
        @assert isapprox(fd[:,j+1], cr[:,j+1]; rtol=1e-4, atol=1e-10)
    end
    @assert isapprox(fd[:,2], cs.derivative; rtol=1e-4, atol=1e-10)
    rows = []
    for (name, values, orders) in (("Original", reshape(original,:,1), 0:0),
        ("CTSE", hcat(cs.value,cs.derivative), 0:1), ("OTI-CR", cr, 0:4), ("FD", fd, 0:4))
        for j in orders
            err = norm(values[:,j+1]-ref[:,j+1], Inf)
            rel = err/max(norm(ref[:,j+1], Inf), eps())
            push!(rows, (name,j,err,rel))
        end
    end
    open(joinpath(output,"agreement.csv"), "w") do io
        println(io,"method,order,max_absolute_error,relative_infinity_error")
        foreach(row -> println(io,join(row,',')), rows)
    end

    p0 = plot(original, label="Original", linewidth=2, xlabel="Voxel", ylabel="Temperature (K)", title="Function value at Δt = 50 s")
    scatter!(p0, cs.value, label="CTSE", marker=:circle, markersize=6)
    scatter!(p0, cr[:,1], label="OTI–CR", marker=:xcross, markersize=5)
    scatter!(p0, fd[:,1], label="FD baseline", marker=:diamond, markersize=3)
    savefig(p0, joinpath(output,"function_values.png"))
    panels = []
    for j in 1:4
        p = plot(ref[:,j+1], label="Analytical", linewidth=2, xlabel="Voxel", ylabel="d$j T / dt$j (K/s^$j)", title="Derivative order $j")
        scatter!(p, cr[:,j+1], label="OTI–CR", marker=:xcross, markersize=5)
        scatter!(p, fd[:,j+1], label="FD", marker=:diamond, markersize=3)
        j == 1 && scatter!(p, cs.derivative, label="CTSE", marker=:circle, markersize=6, markercolor=:transparent, markerstrokecolor=:purple)
        push!(panels,p)
    end
    savefig(plot(panels..., layout=(2,2),size=(1200,850),left_margin=7Plots.mm,bottom_margin=5Plots.mm), joinpath(output,"derivatives.png"))

    # Step-size sweep exposes truncation at large h and cancellation at small h.
    steps = 10.0 .^ range(-6,1; length=65)
    sweep_panels = []
    open(joinpath(output,"step_sweep.csv"), "w") do io
        println(io,"method,order,h,relative_infinity_error")
        for j in 1:4
            scale = norm(ref[:,j+1], Inf)
            errors = [norm(finite_difference(f,dt,j;h)-ref[:,j+1],Inf)/scale for h in steps]
            for (h,err) in zip(steps,errors)
                println(io,join(("FD",j,h,err),','))
            end
            p = plot(steps, max.(errors,eps()), xscale=:log10, yscale=:log10, label="FD", xlabel="Step h (s)", ylabel="Relative infinity error", title="Derivative order $j")
            crerr = norm(cr[:,j+1]-ref[:,j+1],Inf)/scale
            hline!(p,[max(crerr,eps())],label="OTI–CR",linestyle=:dash)
            chosen = fd_step(dt,j)
            scatter!(p,[chosen],[max(norm(fd[:,j+1]-ref[:,j+1],Inf)/scale,eps())],label="Selected FD h",markersize=5)
            if j == 1
                cssteps = 10.0 .^ range(-30,1; length=65)
                cserrors = [norm(ctse_response(A,B,e,x,u,dt;h).derivative-ref[:,2],Inf)/scale for h in cssteps]
                for (h,err) in zip(cssteps,cserrors)
                    println(io,join(("CTSE",1,h,err),','))
                end
                plot!(p,cssteps,max.(cserrors,eps()),label="CTSE")
            end
            push!(sweep_panels,p)
        end
    end
    savefig(plot(sweep_panels...,layout=(2,2),size=(1200,850),left_margin=7Plots.mm,bottom_margin=5Plots.mm),joinpath(output,"step_convergence.png"))
    open(joinpath(@__DIR__,"docs","comparison_results.tex"),"w") do io
        println(io,raw"\begin{center}\begin{tabular}{lrrr}\toprule")
        println(io,raw"Method & Order & Max. absolute error & Relative infinity error \\ \midrule")
        for (name,j,err,rel) in rows
            @printf(io,"%s & %d & %.2e & %.2e ",name,j,err,rel)
            println(io,repeat(string(Char(92)),2))
        end
        println(io,raw"\bottomrule\end{tabular}\end{center}")
    end
    println("All comparison checks passed. Errors against independent analytical reference:")
    foreach(row -> @printf("%-8s order %d: absolute %.3e, relative %.3e\n",row...),rows)
    println("FD steps (seconds): ", [fd_step(dt,j) for j in 1:4])
    println("Saved CSV data and comparison plots to results/.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_sensitivities()
end
