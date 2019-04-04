using JuMP, FlexibilityAnalysis, Gurobi

# Set the covariance matrix for the uncertain parameters
means = [4.0; 5.0]
bea = 1
covar = [2. bea;
         bea 3.]

# Setup the model
m = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))

# Define variables
@randomvariable(m, x, mean = means[1])
@randomvariable(m, y, mean = means[2])

# Define the constraints
@constraint(m, x + y - 14 <= 0.0)
@constraint(m, x - 2y - 2 <= 0.0)
@constraint(m, -x <= 0.0)
@constraint(m, -y <= 0.0)

# Define the uncertainty set
setuncertaintyset(m, :Ellipsoid, covar)

# Solve
status = solve(m, active_constr = true)
if status == :Optimal
    # Retrieve optimized data
    flexibility_index = getflexibilityindex(m)
    x_opt = getvalue(x)
    y_opt = getvalue(y)
    actives = getactiveconstraints(m)

    # Get the confidence level
    conf_lvl = getconfidencelevel(m)

    # Rank the inequality constraints
    rank_data = rankinequalities(m, active_constr = true)
end

# Estimate the value of SF
SF = findstochasticflexibility(m, use_flexibility_index = true)
