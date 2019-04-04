# Library

```@meta
CurrentModule = FlexibilityAnalysis
```

## Module

```@docs
FlexibilityAnalysis
```

## Functions/Methods

```@docs
FlexibilityModel
@randomvariable
@recoursevariable
getflexibilitydata
setcovariance
getcovariance
setmean
getmean(m::Model)
getmean(variable::RandomVariable)
setuncertaintyset
ismeanfeasible
findcenteredmean
getflexibilityindex
getconfidencelevel
getsolutiontime
getactiveconstraints
JuMP.getvalue(v::FlexibilityVariable)
rankinequalities
findstochasticflexibility
JuMP.linearindex(v::FlexibilityVariable)
```

## DataTypes

```@docs
FlexibilityVariable
RandomVariable
RecourseVariable
FlexibilityExpr
FlexibilityConstraint
AbstractUncertaintySet
EllipsoidalSet
HyperboxSet
PNormSet
FlexibilityData
```

## Internals

```@docs
getuncertaintyset
solvehook
MakeInputDict
AddSystemExpressions
RandomVariable(m::Model, mean::Number, name::AbstractString)
RecourseVariable(m::Model, name::AbstractString)
Base.show(io::IO, a::FlexibilityExpr)
JuMP.addconstraint(m::Model, constr::FlexibilityConstraint)
JuMP.show(io::IO,c::FlexibilityConstraint)
JuMP.constructconstraint!(flex_aff::FlexibilityExpr, sense::Symbol)
ComputeCenter
```

## Index

```@index
Pages = ["api.md"]
Module = ["FlexibilityAnalysis"]
Order = [:function, :type]
```
