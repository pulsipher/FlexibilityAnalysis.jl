# Test the random variable macro with complex input
@randomvariable(m2, xs[i = 1:5], mean = 0)
means = ones(4)
@randomvariable(m2, x2s[i = 1:4], mean = means[i])
@randomvariable(m2, x3s[i = 1:3, j = 1:2], mean = 42)
means2 = ones(4, 4)
@randomvariable(m2, x4s[i = 1:4, j = 1:4], mean = means2[i, j])

# Check the recourse variable macro with vector and array
@recoursevariable(m2, ys[1:3])
@recoursevariable(m2, y2s[1:5, 1:6])

# Check the base JuMP macro
@variable(m2, zs[1:2])
@variable(m2, z2s[1:3, 1:5])
true
