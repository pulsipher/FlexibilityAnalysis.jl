# Verify that using F works with both methods
SF = findstochasticflexibility(m, num_pts = 100, seed = 42)
@test abs(findstochasticflexibility(m, num_pts = 100, use_flexibility_index = true, seed = 42) - SF) <= 1e-5
@test abs(findstochasticflexibility(m, num_pts = 100, use_vulnerability_model = true, use_flexibility_index = true, seed = 42) - SF) <= 1e-5

# Verify using F works with hyperbos set
dev = ones(flex_data.numRVs) * 10
setuncertaintyset(m, :Hyperbox, [[dev]; [dev]])
solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 0.5) <= 1e-2
SF = findstochasticflexibility(m, num_pts = 100, seed = 42)
@test abs(findstochasticflexibility(m, num_pts = 100, use_flexibility_index = true, seed = 42) - SF) <= 1e-5
@test abs(findstochasticflexibility(m, num_pts = 100, use_vulnerability_model = true, use_flexibility_index = true, seed = 42) - SF) <= 1e-5

# Verify usign F works with PNorm set
setuncertaintyset(m, :PNorm, Inf)
solve(m, active_constr = true) == :Optimal && abs(getflexibilityindex(m) - 5) <= 1e-2
SF = findstochasticflexibility(m, num_pts = 100, seed = 42)
@test abs(findstochasticflexibility(m, num_pts = 100, use_flexibility_index = true, seed = 42) - SF) <= 1e-5
@test abs(findstochasticflexibility(m, num_pts = 100, use_vulnerability_model = true, use_flexibility_index = true, seed = 42) - SF) <= 1e-5
