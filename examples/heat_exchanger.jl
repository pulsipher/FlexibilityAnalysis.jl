using JuMP
using FlexibilityAnalysis
using Gurobi

# Setup the uncertainty set parameters
means = [620; 388; 583; 313]
a = 5
covar = [11.11 a a a; a 11.11 a a; a a 11.11 a; a a a 11.11]
box_dev = ones(4) * 10

# Setup the model
m = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))

# Define variables
@randomvariable(m, T[i = 1:4], mean = means[i])
@recoursevariable(m, Qc)
@variable(m, x)

# Define the constraints
@constraint(m, -100 - 0.67Qc + 2T[2] + x <= 0.0)
@constraint(m, -250 - T[2] == x)
@constraint(m, 1388.5 + 0.5Qc - 0.75T[1] - T[2] - T[3] <= 0.0)
@constraint(m, 2044 + Qc - 1.5T[1] - 2T[2] - T[3] <= 0.0)
@constraint(m, 2830 + Qc - 1.5T[1] - 2T[2] - T[3] - 2T[4] <= 0.0)
@constraint(m, -3153 - Qc + 1.5T[1] + 2T[2] + T[3] + 3T[4] <= 0.0)

# Define the uncertainty set
setcovariance(m, covar)
# setuncertaintyset(m, :Ellipsoid)
setuncertaintyset(m, :Hyperbox, [[box_dev]; [box_dev]])
# setuncertaintyset(m, :PNorm, 1)
# setuncertaintyset(m, :PNorm, 2)
# setuncertaintyset(m, :PNorm, Inf)

# Ensure the mean is valid
if ismeanfeasible(m)

    # Can estimate the SF without first solving the flexibility model
    SF = findstochasticflexibility(m, num_pts = 2000, use_vulnerability_model = true)

    # Solve
    status = solve(m, active_constr = true)

    if status == :Optimal
        # Retrieve optimized data
        flexibility_index = getflexibilityindex(m)
        temperatures = getvalue(T)
        cooling = getvalue(Qc)
        state = getvalue(x)
        actives = getactiveconstraints(m)

        # Parse data
        data = getflexibilitydata(m)

        # Get the confidence level
        if data.uncertainty_set.name == :Ellipsoid
            conf_lvl = getconfidencelevel(m)
        end

        # Rank the inequality constraints
        rank_data = rankinequalities(m, max_ranks = 3, active_constr = true)
    end
end
