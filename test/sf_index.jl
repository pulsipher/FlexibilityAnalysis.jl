# Compute the base SF index and verify same answer with culnerability model
SF = findstochasticflexibility(m, num_pts = 100, seed = 42)
@test abs(findstochasticflexibility(m, num_pts = 100, use_vulnerability_model = true, seed = 42) - SF) <= 1e-5
@test_logs (:warn, "Flexibility index hasn't yet been computed. Setting use_flexibility_index = false.") abs(findstochasticflexibility(m, num_pts = 100, use_flexibility_index = true, seed = 42) - SF) <= 1e-5

# Compute the SF with all positive values and verify both methods yield the same solution
SF2 = findstochasticflexibility(m, num_pts = 100, only_positive = true, seed = 42)
@test abs(findstochasticflexibility(m, num_pts = 100, use_vulnerability_model = true, only_positive = true, seed = 42) - SF2) <= 1e-5
