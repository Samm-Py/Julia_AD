# Thermal timestep sensitivities

The quantities of interest are the ten final voxel temperatures. Differentiate
with respect to elapsed time `Δt`, holding the initial state, matrices, input,
and forcing fixed.

- `MatrixExponentialTest.jl`: restored original script, unchanged.
- `CTSEExample.jl`: complex Taylor series expansion (complex step); values and first derivatives.
- `CRExample.jl`: real Cauchy–Riemann representation of OTI numbers; values and derivatives through order four.
- `CompareSensitivities.jl`: runs the original and both methods, computes real finite differences and an independent analytical reference, checks agreement, and saves plots and CSV data.

Run from the repository root with the existing project environment:

```sh
julia --project=. CTSEExample.jl
julia --project=. CRExample.jl
julia --project=. CompareSensitivities.jl
julia --project=. test/runtests.jl
```

The comparison runs headlessly by default, including the original's display call.
The individual method scripts do not require plotting. Results go into `results/`:

- `function_values.png`: original, CTSE, OTI–CR, and unperturbed FD baseline.
- `derivatives.png`: orders 1–4, with CTSE only at order one.
- `step_convergence.png`: finite-difference step sweep, CTSE first-order step sweep, and OTI–CR error reference.
- `original.csv`, `ctse.csv`, `cr.csv`, `finite_difference.csv`, `analytical.csv`: per-voxel values and supported derivatives (`d0` is the value).
- `agreement.csv`: maximum absolute and relative infinity-norm errors against the analytical reference.
- `step_sweep.csv`: raw errors and step sizes.

The independent reference integrates the eigenmodes of the symmetric thermal
matrix, then differentiates using the ODE. Finite differences use seven-point
centered stencils with rationally computed weights and order-dependent steps.
They use only real function values. The step sweep exposes truncation and
cancellation errors; higher-order FD is less accurate than OTI–CR.

Shared numerical code is in `src/OTISensitivity.jl` and
`src/SensitivityComparison.jl`. The original's `+ e` update is retained exactly;
it is valid for its zero-forcing example. The new response routine uses `+ ed`
and is tested with nonzero forcing.

Read [the explanation](docs/oti_sensitivities.pdf) or edit
[its LaTeX source](docs/oti_sensitivities.tex). CTSE is introduced before OTI and
the real matrix representation. To regenerate the document and its result table:

```sh
julia --project=. CompareSensitivities.jl
julia --project=. docs/ctse_worked_example.jl
julia --project=. docs/fd_worked_example.jl
cd docs
pdflatex -interaction=nonstopmode -halt-on-error oti_sensitivities.tex
pdflatex -interaction=nonstopmode -halt-on-error oti_sensitivities.tex
```

`docs/verify_results.jl` remains a standalone OTI-only numerical check from the
first implementation; the current document uses the comparison-generated table.
The earlier `matrix_exponential_result.png` is a legacy output; current comparison
plots live in `results/`.

The CTSE worked example in Section 1.3.1 has its own reproducible table and
convergence figure. Run `julia --project=. docs/ctse_worked_example.jl`; CSV
data and the plot are saved under `results/worked_example/`. Derivative errors
compare Float64 complex evaluations against a 256-bit reference.

Section 1.2 uses the same scalar function for a centered finite-difference
truncation-versus-roundoff plot. Regenerate it with
`julia --project=. docs/fd_worked_example.jl`; the figure and CSV are in
`results/worked_example/`. The roundoff curve is an illustrative scale,
while observed and truncation errors use a 256-bit reference.

Appendix A describes multivariate OTI arithmetic with total-degree truncation,
commuting multiplication matrices, and extraction of mixed derivatives. Run
`julia --project=. docs/verify_multivariate_example.jl` to verify its examples.
The production thermal entry points remain single-parameter methods.
