# Test that basic function works and that it doesn change the mean
@test all(findcenteredmean(m, only_positive = true) .>= 0)
@test all(getmean(m) .== means)

# Try different centers and changing the mean
@test isa(findcenteredmean(m, center = :analytic), Vector)
@test all(findcenteredmean(m, center = :analytic, only_positive = true, update_mean = true) .>= 0)
@test !all(getmean(m) .== means)

# Test invalid input
@test_throws ErrorException findcenteredmean(m, center = :Analytic)

# Reset the mean back to correct value
setmean(m, means)
@test all(getmean(m) .== means)
