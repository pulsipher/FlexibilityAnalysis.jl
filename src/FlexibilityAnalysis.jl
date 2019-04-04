"""
    FlexibilityAnalysis

A Julia package that implements a computational framework for quantifying and analyzing system flexibility.

The basic functionality emulates typical JuMP models to facilitate general computation of the flexibility
index problem. Thus, basic functionality uses [`FlexibilityModel`](@ref), [`@randomvariable`](@ref),
[`@recoursevariable`](@ref), `@variable`, `@constraint`, [`setuncertaintyset`](@ref), and `solve`.

**Methods/Macros**
- [`FlexibilityModel`](@ref)
- [`@randomvariable`](@ref)
- [`@recoursevariable`](@ref)
- [`getflexibilitydata`](@ref)
- [`setcovariance`](@ref)
- [`getcovariance`](@ref)
- [`getmean`](@ref)
- [`setmean`](@ref)
- [`setuncertaintyset`](@ref)
- [`ismeanfeasible`](@ref)
- [`findcenteredmean`](@ref)
- [`getflexibilityindex`](@ref)
- [`getsolutiontime`](@ref)
- [`getconfidencelevel`](@ref)
- [`getactiveconstraints`](@ref)
- [`rankinequalities`](@ref)
- [`findstochasticflexibility`](@ref)
"""
module FlexibilityAnalysis

    using Distributions
    using JuMP
    using LinearAlgebra
    using Random
    using Ipopt, Clp # Default solver

    # Define functions/macros to be readily accesible to the user
    export FlexibilityModel, @randomvariable, @recoursevariable, getflexibilitydata, setcovariance, getcovariance,
           getactiveconstraints, setuncertaintyset, getflexibilityindex, findcenteredmean, ismeanfeasible,
           getconfidencelevel, findstochasticflexibility, rankinequalities, getmean, setmean, getsolutiontime

    # Import all of the datatypes, methods, macros, and definitions
    include("datatypes.jl")
    include("model.jl")
    include("solve.jl")
    include("variables.jl")
    include("constraints.jl")
    include("macros.jl")
    include("uncertaintyset.jl")
    include("mean.jl")
    include("functions.jl")
    include("operators.jl")

end
