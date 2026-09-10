if !isdefined(@__MODULE__, :SensitivityComparison)
    include(joinpath(@__DIR__, "src", "SensitivityComparison.jl"))
end
using .SensitivityComparison

function run_cr(; model=thermal_example(), order=4,
                output=joinpath(@__DIR__, "results", "cr.csv"))
    (; A, B, e, x, u, dt) = model
    result = oti_response(A, B, e, x, u, dt; order)
    write_method_csv(output, dt, result.derivatives)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cr()
    println("Saved OTI–CR values and derivatives through order four to results/cr.csv")
end
