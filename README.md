# Thermal timestep sensitivities

The quantities of interest are the ten final voxel temperatures of a
one-dimensional thermal model. Differentiate them with respect to the elapsed
time `Δt`, holding the initial state, matrices, input and forcing fixed. This
repository computes the first derivatives two ways: complex Taylor series
expansion (CTSE, the complex step) and dual numbers in real Cauchy–Riemann (CR)
matrix form.

Start with [worked examples of CTSE and dual numbers in CR form](ctse_dual_cr_examples.pdf)
([LaTeX source](ctse_dual_cr_examples.tex)). It walks through scalar
`sin(x)` and matrix `sin(X)` with each method, shows where CTSE breaks for
matrix functions, and reports both methods' results for the thermal model
against 256-bit references.

- `examples/`: each Julia calculation in the document as a runnable script,
  shown verbatim in the PDF: CTSE and dual CR for scalar and matrix `sin`, the
  hand-written dual type, and both methods on the thermal model.
- `src/OTISensitivity.jl`: the thermal model, its one-step response, and the CR
  representation (`oti_matrix`, `oti_response`; the document uses `order=1`,
  which is ordinary dual numbers).
- `src/SensitivityComparison.jl`: `ctse_response` and an independent modal
  reference, `analytical_response`.
- `docs/verify_ctse_dual_cr_examples.jl`: runs the scripts in `examples/`,
  checks every example in the document and regenerates its five comparison
  tables.
- `test/runtests.jl`: tests of the CTSE and first-order CR implementations.

Run from the repository root with the existing project environment:

```sh
julia --project=. examples/thermal.jl
julia --project=. test/runtests.jl
julia --project=. docs/verify_ctse_dual_cr_examples.jl
pdflatex -interaction=nonstopmode -halt-on-error ctse_dual_cr_examples.tex
pdflatex -interaction=nonstopmode -halt-on-error ctse_dual_cr_examples.tex
```
