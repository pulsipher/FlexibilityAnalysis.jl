# Test for the correct answer
@test abs(getconfidencelevel(m) - 0.53722) <= 1e-3

# Check that it doesn't work for other sets and then reset the set
setuncertaintyset(m, :PNorm, 1)
@test_throws ErrorException abs(getconfidencelevel(m) - 0.53722) <= 1e-3
setuncertaintyset(m, :Ellipsoid)
