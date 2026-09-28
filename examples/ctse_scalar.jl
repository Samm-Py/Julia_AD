# Section 1 of docs/ctse_dual_cr_examples.pdf: CTSE for f(x) = sin(x) at x = 1.

x0, h = 1.0, 1e-20
f(x) = sin(x)
y = f(x0 + im*h)
value = real(y)
derivative = imag(y)/h
@show value derivative
