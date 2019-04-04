using JuMP
using FlexibilityAnalysis
using Gurobi
using LinearAlgebra

# Set the dimensions
n_lines = 2
n_dems = 3

# Setup the uncertainty set parameters
means = [0.; 60.; 10.]
covar = [80. 0 0; 0 80. 0; 0 0 120.]
box_dev = 2 * sqrt.(diag(covar))

# Specify the network details
line_cap = [35, 50]
gen_cap = 100

# Setup the model
m = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))

# Define variables
@randomvariable(m, d[i = 1:n_dems], mean = means[i])
@recoursevariable(m, a[1:n_lines])
@recoursevariable(m, g)

# Set the line capacity constraints
@constraint(m, [line = 1:n_lines], -line_cap[line] <= a[line])
@constraint(m, [line = 1:n_lines], a[line] <= line_cap[line])

# Set the generator capacity constraints
@constraint(m, 0.0 <= g)
@constraint(m, g <= gen_cap)

# Set the node balance constraints
@constraint(m, a[1] - d[1] == 0)
@constraint(m, -a[1] -a[2] + g - d[2] == 0)
@constraint(m, a[2] - d[3] == 0)

# Define the uncertainty set
setcovariance(m, covar)
setuncertaintyset(m, :Ellipsoid, only_positive = true)
# setuncertaintyset(m, :Hyperbox, [[box_dev]; [box_dev]], only_positive = true)
# setuncertaintyset(m, :PNorm, 1, only_positive = true)
# setuncertaintyset(m, :PNorm, 2, only_positive = true)
# setuncertaintyset(m, :PNorm, Inf, only_positive = true)

# Compute a center to replace the mean if desired
new_mean = findcenteredmean(m, center = :analytic, update_mean = true, only_positive = true)
updated_mean = getmean(m)

# Ensure the mean is valid
if ismeanfeasible(m)
    # Solve
    status = solve(m, U = 1000) # Change slack variable upper bound

    if status == :Optimal
        # Retrieve optimized data
        flexibility_index = getflexibilityindex(m)
        lines = getvalue(a)
        generator = getvalue(g)
        demands = getvalue(d)
        actives = getactiveconstraints(m)

        # Parse data
        data = getflexibilitydata(m)

        # Get the confidence level
        if data.uncertainty_set.name == :Ellipsoid
            conf_lvl = getconfidencelevel(m)
        end

        # Rank the inequality constraints
        rank_data = rankinequalities(m, max_ranks = 4, U = 1000)
    end

    # Estimate the SF index
    SF = findstochasticflexibility(m, use_vulnerability_model = true, use_flexibility_index = true, only_positive = true)
end
