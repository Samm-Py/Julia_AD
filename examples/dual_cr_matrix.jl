# Section 4 of ctse_dual_cr_examples.pdf: the matrix sine of tH at t = 1 in dual CR form.

using LinearAlgebra
H = [0.0 1.0; -1.0 0.0]
T = [1.0 0.0; 1.0 1.0]
G = sin(kron(T, H))
value = G[1:2, 1:2]
derivative = G[3:4, 1:2]
@show value derivative
