# Extract the data
flex_data = getflexibilitydata(m)
flex_data2 = getflexibilitydata(m2)
@test_throws ErrorException getflexibilitydata(Model())

# Check data type
@test isa(flex_data, FlexibilityAnalysis.FlexibilityData)
@test isa(flex_data2, FlexibilityAnalysis.FlexibilityData)

# Check F index extractor function
@test getflexibilityindex(m) == nothing

# Check prelimnary errors
@test_throws ErrorException getconfidencelevel(m)

# checking the model printing
print_to_string(f::Function, args...) = begin
    io = IOBuffer()
    f(io,args...)
    String(take!(io))
end
@test print_to_string(JuMP.show, m.ext[:FlexData].flexibility_constraints[1]) == "0.67*Qc + -2*T[2] + -x + 100 >= 0"
@test print_to_string(JuMP.show, m.ext[:FlexData].flexibility_constraints[4]) == "1*Qc + -1.5*T[1] + -2*T[2] + -1*T[3] + 2044 <= 0"
