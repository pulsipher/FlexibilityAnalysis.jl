# Check default
@test flex_data.uncertainty_set.name == :Ellipsoid && flex_data2.uncertainty_set.name == :Ellipsoid
@test FlexibilityAnalysis.getuncertaintyset(m).name == :Ellipsoid

# Check safeguards
@test_throws ErrorException setuncertaintyset(m, :Ellipsoid)
@test_throws ErrorException setuncertaintyset(m, :Ellipsoid, 42)
@test_throws ErrorException setuncertaintyset(m, :PNorm)
@test_throws ErrorException setuncertaintyset(m, :ellipsoid)
@test_throws ErrorException setuncertaintyset(m, :PNorm, 3)
@test_throws ErrorException setuncertaintyset(m, :Hyperbox, [[1; 1];[1; 1]])
@test_throws ErrorException setuncertaintyset(m, :Hyperbox, [[[1; 1]];[[1; 1; 1]]])

# Test each set type
a = 0
covar = [11.11 a a a; a 11.11 a a; a a 11.11 a; a a a 11.11]
setcovariance(m, covar)
@test flex_data.uncertainty_set.name == :Ellipsoid
setuncertaintyset(m2, :PNorm, 1)
@test flex_data2.uncertainty_set.name == :PNorm && flex_data2.uncertainty_set.p == 1
setuncertaintyset(m2, :PNorm, 2)
@test flex_data2.uncertainty_set.name == :PNorm && flex_data2.uncertainty_set.p == 2
setuncertaintyset(m2, :PNorm, Inf)
@test flex_data2.uncertainty_set.name == :PNorm && flex_data2.uncertainty_set.p == Inf
setuncertaintyset(m2, :Hyperbox, [[[1; 1]]; [[1; 1]]])
@test flex_data2.uncertainty_set.name == :Hyperbox && all(flex_data2.uncertainty_set.neg_dev .== ones(2))
setuncertaintyset(m2, :PNorm, 1, only_positive = true)
@test flex_data2.uncertainty_set.name == :PNorm && flex_data2.uncertainty_set.only_positive
setuncertaintyset(m2, :Hyperbox, [[[1; 1]]; [[1; 1]]], only_positive = true)
@test flex_data2.uncertainty_set.name == :Hyperbox && flex_data2.uncertainty_set.only_positive
setuncertaintyset(m2, :Ellipsoid, Matrix(I, flex_data2.numRVs, flex_data2.numRVs), only_positive = true)
@test flex_data2.uncertainty_set.name == :Ellipsoid && flex_data2.uncertainty_set.only_positive
