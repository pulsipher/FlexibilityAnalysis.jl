# Test the getvalue function for each variable type
@test isa(getvalue(T), Vector)
@test isa(getvalue(Qc), Float64)
@test isa(getvalue(x), Float64)
