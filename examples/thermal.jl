# Thermal results in docs/ctse_dual_cr_examples.pdf. Run with: julia --project=. examples/thermal.jl

include(joinpath(@__DIR__, "..", "src", "SensitivityComparison.jl"))
using .SensitivityComparison
(; A, B, e, x, u, dt) = thermal_example()   # ten voxels, t = 50 s

# CTSE: evaluate the response at t + ih.
ctse = ctse_response(A, B, e, x, u, dt; h=1e-20)
# Dual CR: exponentiate kron(T, H) with T = [t 0; 1 t].
cr = oti_response(A, B, e, x, u, dt; order=1)
@show ctse.value ctse.derivative cr.derivatives[:, 2]
