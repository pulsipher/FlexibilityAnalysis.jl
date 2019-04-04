# Test the default case
rank_data = rankinequalities(m, active_constr = true)
@test isa(rank_data, Vector{Dict})

# Test varying the max_ranks to work as expected
rank_data = rankinequalities(m, max_ranks = 1, active_constr = true)
@test isa(rank_data, Vector{Dict}) && length(rank_data) == 1
@test_throws ErrorException rankinequalities(m, max_ranks = 0, active_constr = true)
