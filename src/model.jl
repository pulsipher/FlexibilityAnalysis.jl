"""
    FlexibilityModel(; [solver = JuMP.UnsetSolver()])
Return a flexibility model object which extends a JuMP model object to contain [`FlexibilityData`](@ref) and implement
a custom solvehook. An appropriate solver should be specified in order solve the flexibility index problem. A
solver capable of handling MIQCPs is required for ellipsoidal and 2-norm uncertainty sets otherwise a MILP solver
can be used. This model is solved with `solve`, see [`solvehook`](@ref) for documention on the accepted keyword arguments.

**Arguments**
- `solver = JuMP.UnsetSolver()` The solver, should use an MIQCP, MINLP, or MILP solver as appropriate.

```julia
julia> m = FlexibilityModel(solver = GurobiSolver())
Feasibility problem with:
 * 0 linear constraints
 * 0 variables
Solver is Gurobi
```
"""
function FlexibilityModel(;solver = JuMP.UnsetSolver())
    m = Model(solver = solver)
    m.solvehook = solvehook     #solvehook is defined in solve.jl
    m.ext[:FlexData] = FlexibilityData(FlexibilityConstraint[], 0, Float64[], String[], Int[], 0, String[], Int[], EllipsoidalSet(), Matrix(undef, 0, 0), nothing, Int[], nothing)
    return m
end

"""
    getflexibilitydata(m::Model)
Return the [`FlexibilityData`](@ref) corresponding the flexibility model `m`. An error is thrown if `m` is a
regular JuMP model.

**Arguments**
- `m::Model` The flexibility model.

```julia
julia> getflexibilitydata(m)
FlexibilityAnalysis.FlexibilityData(JuMP.AbstractConstraint[], 0, Number[], AbstractString[], Int64[], 0, AbstractString[], Int64[], FlexibilityAnalysis.EllipsoidalSet(:Ellipsoid, false), Array{Number}(undef,0,0), nothing, Int64[], nothing)
```
"""
function getflexibilitydata(m::Model)
    if haskey(m.ext,:FlexData)
        return m.ext[:FlexData]
    else
        error("This functionality is only available for Flexibility models")
    end
end

"""
    getuncertaintyset(m::Model)
Return the current uncertainty set datatype in the flexibility modelas stored in `FlexibilityData.uncertaintyset`.

**Arguments**
- `m::Model` The flexibility model.
"""
getuncertaintyset(m::Model) =  getflexibilitydata(m).uncertainty_set

"""
    getcovariance(m::Model)
Return the current covariance matrix in the flexibility model as stored in `FlexibilityData.covariance`.

**Arguments**
- `m::Model` The flexibility model.

```julia
julia> getcovariance(m)
4×4 Array{Number,2}:
 11.11   0.0    0.0    0.0
  0.0   11.11   0.0    0.0
  0.0    0.0   11.11   0.0
  0.0    0.0    0.0   11.11
```
"""
getcovariance(m::Model) = getflexibilitydata(m).covariance

"""
    getmean(m::Model)
Return the current mean vector in the flexibility model as stored in `FlexibilityData.RVmeans`.

**Arguments**
- `m::Model` The flexibility model.

```julia
julia> getmean(m)
4-element Array{Number,1}:
 620
 388
 583
 313
```
"""
getmean(m::Model) = getflexibilitydata(m).RVmeans

"""
    getflexibilityindex(m::Model)
Return the current flexibility index in the flexibility model as stored in `FlexibilityData.flexibility_index`.

**Arguments**
- `m::Model` The flexibility model.

```julia
julia> getflexibilityindex(m)
3.5993764186390327
```
"""
getflexibilityindex(m::Model) = getflexibilitydata(m).flexibility_index

"""
    getsolutiontime(m::Model)
Return the solution time to compute the flexibility index as stored in the flexibility model as stored in `FlexibilityData.solution_time`.

**Arguments**
- `m::Model` The flexibility model.

```julia
julia> getsolutiontime(m)
0.00199127197265625
```
"""
getsolutiontime(m::Model) = getflexibilitydata(m).solution_time

"""
    getactiveconstraints(m::Model)
Return the current vector of active constraint indexes in the flexibility model as stored in `FlexibilityData.active_constraints`.

**Arguments**
- `m::Model` The flexibility model.

```julia
julia> getactiveconstraints(m)
2-element Array{Int64,1}:
 3
 6
```
"""
getactiveconstraints(m::Model) = getflexibilitydata(m).active_constraints

"""
    setcovariance(m::Model, covariance::Matrix)
Specify the covariance matrix `covariance` to be stored in the flexibility model `m`. This method verifies that
the matrix is symmetric positive semi-definite and writes it to `FlexibilityData.covariance`.

**Arguments**
- `m::Model` The flexibility model.
- `covariance::Matrix` The covariance matrix.

```julia
julia> setcovariance(m, [2 1; 1 2])
2×2 Array{Int64,2}:
 2  1
 1  2
```
"""
function setcovariance(m::Model, covariance::Matrix)
    # Run checks on the covairance matrix before using it
    n_row,ncol = size(covariance)
    n_row == ncol || error("Covariance matrix should be square")
    issymmetric(covariance) || error("Covariance matrix should be symmetric")
    all(eigen(covariance).values .>= 0) || error("Covariance matrix is not positive semi-definite")

    # Set the covariance matrix
    flex_data = getflexibilitydata(m)
    flex_data.covariance = covariance
end

"""
    setmean(m::Model, mean::Vector)
Specify the mean corresponding to `FlexibilityData.RVmeans` stored in the flexibility model `m`. This method verifies that
the length of the input `mean` matches the length of `FlexibilityData.RVmeans` before overwriting the current mean.

**Arguments**
- `m::Model` The flexibility model.
- `mean::Vector` The means of the random variables.

```julia
setmean(m, [2.3; 5])
2-element Array{Float64,1}:
 2.3
 5.0
```
"""
function setmean(m::Model, mean::Vector)
    flex_data = getflexibilitydata(m)
    if length(flex_data.RVmeans) == length(mean)
        flex_data.RVmeans = mean
    else
        error("The length of the specified mean doesn't match the current mean length.")
    end
end
