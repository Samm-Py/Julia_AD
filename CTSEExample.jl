if !isdefined(@__MODULE__, :SensitivityComparison)
    include(joinpath(@__DIR__, "src", "SensitivityComparison.jl"))
end
using .SensitivityComparison

function run_ctse(; model=thermal_example(), h=1e-20,
                  output=joinpath(@__DIR__, "results", "ctse.csv"))
    (; A, B, e, x, u, dt) = model
    result = ctse_response(A, B, e, x, u, dt; h)
    write_method_csv(output, dt, hcat(result.value, result.derivative))
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_ctse()
    println("Saved CTSE values and first derivatives to results/ctse.csv")
end
