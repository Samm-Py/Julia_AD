# Section 2 of ctse_dual_cr_examples.pdf: CTSE for the matrix sine of tH at t = 1.

using LinearAlgebra
H = [0.0 1.0; -1.0 0.0]
t0, h = 1.0, 1e-20
G = sin((t0 + im*h)*H)
value = real.(G)
derivative = imag.(G)/h
@show value derivative
