# Test the solvehook's options
@test_logs (:warn, "Problem requires use of active set constraint. Setting active_constr = true") solve(m) == :Optimal
@test solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 3.6) <= 1e-2
@test solve(m, active_constr = true, diag = true) == :Optimal && abs(getflexibilityindex(m) - 3.6) <= 1e-2
@test solve(m, active_constr = true, conic_δ = true) == :Optimal && abs(getflexibilityindex(m) - 3.6) <= 1e-2
@test solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 3.6) <= 1e-2
@test_logs (:warn, "The slack upper bound (U) value of 1000 is too low for this problem. Thus, U is being set to 10000.0 !") solve(m, active_constr = true, U = 1000) == :Optimal

# Test with other uncertainty sets
dev = ones(flex_data.numRVs) * 10
setuncertaintyset(m, :Hyperbox, [[dev]; [dev]])
@test solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 0.5) <= 1e-2
setuncertaintyset(m, :PNorm, 1)
@test solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 6.667) <= 1e-2
setuncertaintyset(m, :PNorm, 2)
@test solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 6.325) <= 1e-2
setuncertaintyset(m, :PNorm, Inf)
@test solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 5) <= 1e-2

# Check safeguards
dev = ones(flex_data.numRVs - 1) * 10
setuncertaintyset(m, :Hyperbox, [[dev]; [dev]])
@test_throws ErrorException solve(m, active_constr = true) == :Optimal
setuncertaintyset(m, :Ellipsoid, Matrix(I, 2, 2))
@test_throws ErrorException solve(m, active_constr = true) == :Optimal
setuncertaintyset(m, :Ellipsoid, covar)
@test solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 3.6) <= 1e-2

# Test the solution time
@test getsolutiontime(m) != nothing
