# Section 3 of ctse_dual_cr_examples.pdf: sin(x) at x = 1 in dual arithmetic.

struct Dual <: Real
    a::Float64   # value
    b::Float64   # derivative
end
Base.sin(x::Dual) = Dual(sin(x.a), cos(x.a)*x.b)
y = sin(Dual(1.0, 1.0))
value, derivative = y.a, y.b
@show value derivative
