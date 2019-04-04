using JuMP
using FlexibilityAnalysis
using Gurobi
using Pavito, Ipopt
using LinearAlgebra

# Set the dimensions
n_gens = 5
n_lines = 20
n_dems = 11

# Setup the uncertainty set parameters
β = 240.
covar = Matrix(I, n_dems, n_dems) * 1200.
covar[covar .== 0] = β
box_dev = ones(n_dems) * 2 * sqrt(covar[1])

# Specify the network details
line_cap = 100
gen_cap = [332; 140; 100; 100; 100]

# Setup the model
m = FlexibilityModel(solver = PavitoSolver(mip_solver = GurobiSolver(OutputFlag = 0),
                     cont_solver = IpoptSolver(print_level = 0), log_level = 0, mip_solver_drives = false))

# Define variables
@randomvariable(m, d[i = 1:n_dems], mean = 0) # Temperarily set the mean to 0
@recoursevariable(m, a[1:n_lines])
@recoursevariable(m, g[1:n_gens])

# Set the line capacity constraints
@constraint(m, [line = 1:n_lines], -line_cap <= a[line])
@constraint(m, [line = 1:n_lines], a[line] <= line_cap)

# Set the generator capacity constraints
@constraint(m, [gen = 1:n_gens], 0.0 <= g[gen])
@constraint(m, [gen = 1:n_gens], g[gen] <= gen_cap[gen])

# Set the node balance constraints
@constraint(m, g[1] - a[1] - a[6] == 0)
@constraint(m, a[1] + g[2] - sum(a[i] for i = [2; 4; 5]) - d[1] == 0)
@constraint(m, g[3] + a[2] - a[3] - d[2] == 0)
@constraint(m, sum(a[i] for i = [3; 4; 8]) - sum(a[i] for i = [7; 11]) - d[3] == 0)
@constraint(m, sum(a[i] for i = [5; 6; 7; 12]) - d[4] == 0)
@constraint(m, g[4] + sum(a[i] for i = [16; 18]) - sum(a[i] for i = [12; 19]) - d[5] == 0)
@constraint(m, a[9] - sum(a[i] for i = [8; 10]) == 0)
@constraint(m, g[5] - a[9] == 0)
@constraint(m, sum(a[i] for i = [10; 11]) - sum(a[i] for i = [13; 14]) - d[6] == 0)
@constraint(m, sum(a[i] for i = [13; 20]) - d[7] == 0)
@constraint(m, a[19] - a[20] - d[8] == 0)
@constraint(m, a[17] - a[18] - d[9] == 0)
@constraint(m, a[15] - sum(a[i] for i = [16; 17]) - d[10] == 0)
@constraint(m, a[14] - a[15] - d[11] == 0)

# Define the uncertainty set
setcovariance(m, covar)
# setuncertaintyset(m, :Ellipsoid, covar)
# setuncertaintyset(m, :Hyperbox, [[box_dev]; [box_dev]])
# setuncertaintyset(m, :PNorm, 1)
# setuncertaintyset(m, :PNorm, 2)
setuncertaintyset(m, :PNorm, Inf)

# Compute a center to replace the mean if desired
new_mean = findcenteredmean(m, center = :feasible, update_mean = true)
updated_mean = getmean(m)

# Solve
status = solve(m)

if status == :Optimal
    # Retrieve optimized data
    flexibility_index = getflexibilityindex(m)
    lines = getvalue(a)
    generators = getvalue(g)
    demands = getvalue(d)
    actives = getactiveconstraints(m)

    # Parse data
    data = getflexibilitydata(m)

    # Get the confidence level
    if data.uncertainty_set.name == :Ellipsoid
        conf_lvl = getconfidencelevel(m)
    end

    # Rank the inequality constraints
    rank_data = rankinequalities(m, max_ranks = 3)
end

# Estimate the value of SF
SF = findstochasticflexibility(m, use_flexibility_index = true)
