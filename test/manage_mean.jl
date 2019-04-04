# Test the get/set mean functions
@test all(getmean(m) .== means)
setmean(m, ones(flex_data.numRVs))
@test all(getmean(m) .== ones(flex_data.numRVs))
@test_throws ErrorException setmean(m, ones(flex_data.numRVs - 1))
setmean(m, means)
@test getmean(T[1]) == 620
