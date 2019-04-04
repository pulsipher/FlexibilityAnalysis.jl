# Set the variables
means = [620; 388; 583; 313]
@randomvariable(m, T[i = 1:4], mean = means[i])
@recoursevariable(m, Qc)
@variable(m, x)

# Define the constraints
@constraint(m, 0.67Qc - 2T[2] - x >= -100)
@constraint(m, -250 - T[2] == x)
@constraint(m, 1388.5 + 0.5Qc - 0.75T[1] - T[2] - T[3] <= 0.0)
@constraint(m, 2044 + Qc - 1.5T[1] - 2T[2] - T[3] <= 0.0)
@constraint(m, 2830 + Qc - 1.5T[1] - 2T[2] - T[3] - 2T[4] <= 0.0)
@constraint(m, -Qc + 1.5T[1] + 2T[2] + T[3] + 3T[4] <= 3153)
true
