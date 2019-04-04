# Try to make constraints with different variable combos
@constraint(m2, x + y - z <= 0)
@constraint(m2, thor[i = 1:5], xs[i] + sum(y2s[i, j] for j = 1:6) == 42)
@constraint(m2, z + 2 >= -23)
@constraint(m2, y + z == 23)
@constraint(m2, x + z <= 42)

# Ensure all the variables and constraints were added
size(JuMP.prepConstrMatrix(m2)) == (9, 84)
