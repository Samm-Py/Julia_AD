# Section 3 of ctse_dual_cr_examples.pdf: sin(x) at x = 1 with a dual number in CR form.

using LinearAlgebra
x0 = 1.0
S = [x0 0.0; 1.0 x0]
Y = sin(S)
value = Y[1, 1]
derivative = Y[2, 1]
@show value derivative
