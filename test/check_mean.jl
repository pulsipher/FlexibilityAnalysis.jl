# Test the default and optional arguments
@test ismeanfeasible(m)
@test ismeanfeasible(m, toler = 1e-8)
@test ismeanfeasible(m, solver = GLPKSolverLP())

# Test an infeasible case
setmean(m, ones(flex_data.numRVs) * -50)
@test !ismeanfeasible(m)

# Reset the mean back to the original value
setmean(m, means)
@test ismeanfeasible(m)
