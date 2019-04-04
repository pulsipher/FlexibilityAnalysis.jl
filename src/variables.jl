"""
    RandomVariable(m::Model, mean::Number, name::AbstractString)
Return a `RandomVariable` DataType for given flexibility model `m` given the `mean` and
the `name`. An anonymous JuMP variable is added directly to the flexibility model and its
is appended to `FlexibilityData.RVcols`.

**Arguments**
- `m::Model` The flexibility model.
- `mean::Number` The variable mean.
- `name::AbstractString` The variable name.
"""
function RandomVariable(m::Model, mean::Number, name::AbstractString)
    # Get the data and update it
    flex_data = getflexibilitydata(m)
    flex_data.numRVs += 1
    push!(flex_data.RVmeans, mean)
    push!(flex_data.RVnames, name)

    # Add an associated anonymous variable directly to the model and save its index
    @variable(m)
    index = length(m.colVal)
    push!(flex_data.RVcols, index)

    return RandomVariable(m, flex_data.numRVs)
end

"""
    getmean(variable::RandomVariable)
Return the mean corresponding to a particular `RandomVariable`. Currently this only accepts a single random
variable and vector variables are not accepted directly.

**Arguments**
- `variable::RandomVariable` The random variable name, must be a single variable.

```julia
julia> getmean(T[1])
620
```
"""
function getmean(variable::RandomVariable)
    flex_data = getflexibilitydata(variable.m)
    return flex_data.RVmeans[variable.idx]
end

"""
    RecourseVariable(m::Model, name::AbstractString)
Return a `RecourseVariable` DataType for given flexibility model `m` given the `name`.
An anonymous JuMP variable is added directly to the flexibility model and its is
appended to `FlexibilityData.recourse_cols`.

**Arguments**
- `m::Model` The flexibility model.
- `name::AbstractString` The variable name.
"""
function RecourseVariable(m::Model, name::AbstractString)
    # Get the data and update it
    flex_data = getflexibilitydata(m)
    flex_data.num_recourse_vars += 1
    push!(flex_data.recourse_names, name)

    # Add an associated anonymous variable directly to the model and save its index
    @variable(m)
    index = length(m.colVal)
    push!(flex_data.recourse_cols, index)

    return RecourseVariable(m, flex_data.num_recourse_vars)
end

"""
    JuMP.getvalue(v::FlexibilityVariable)
Return the value of the a flexibility variable this is an extension of `JuMP.getvalue`.

**Arguments**
- `v::FlexibilityVariable` The flexibility variable.

```julia
julia> getvalue(T)
4-element Array{Float64,1}:
 620.0
 388.0
 581.0
 319.0
```
"""
JuMP.getvalue(v::FlexibilityVariable) = getvalue(Variable(v.m, v.idx))

"""
    JuMP.linearindex(v::FlexibilityVariable)
Return the index of the a flexibility variable this is an extension of `JuMP.linearindex`.

**Arguments**
- `v::FlexibilityVariable` The flexibility variable, must be a single variable.

```julia
julia> linearindex(Qc)
1
```
"""
JuMP.linearindex(v::FlexibilityVariable) = v.idx
