"""
    FlexibilityVariable <: JuMP.AbstractJuMPScalar
An abstract type to define new variable types.
"""
abstract type FlexibilityVariable <: JuMP.AbstractJuMPScalar end

"""
    RandomVariable <: FlexibilityVariable
A DataType for random variables.

**Fields**
- `m::Model` Flexibility model.
- `idx::Int` Index of variable in model.
"""
struct RandomVariable <: FlexibilityVariable
    m::Model
    idx::Int
end

"""
    RecourseVariable <: FlexibilityVariable
A DataType for recourse variables.

**Fields**
- `m::Model` Flexibility model.
- `idx::Int` Index of variable in model.
"""
struct RecourseVariable <: FlexibilityVariable
    m::Model
    idx::Int
end

"""
    FlexibilityExpr <: JuMP.GenericAffExpr
A `GenericAffExpr` that contains random and/or recourse variables.
"""
const FlexibilityExpr = JuMP.GenericAffExpr{AffExpr, Union{RandomVariable, RecourseVariable}}
FlexibilityExpr() = FlexibilityExpr(RandomVariable[], AffExpr[], AffExpr())


"""
    FlexibilityConstraint <: JuMP.AbstractConstraint
A constraint that contains random and/or recourse variables.

**Fields**
- `flex_expr::FlexibilityExpr` Constraint expression.
- `sense::Symbol` The cosntraint sense symbol `:(<=)` or `:(>=)` or `:(==)`.
"""
mutable struct FlexibilityConstraint <: JuMP.AbstractConstraint
    flex_expr::FlexibilityExpr
    sense::Symbol # :(<=) or :(>=) or :(==), right-hand side assumed to be zero
end

"""
    AbstractUncertaintySet
An abstract type to define new uncertainty set types.
"""
abstract type AbstractUncertaintySet end

"""
    EllipsoidalSet <: AbstractUncertaintySet
An ellipsoidal uncertainty set that will use the covariance matrix stored in [`FlexibilityData`](@ref).

**Fields**
- `name::Symbol` Name of the set which will be `:Ellipsoid`.
- `only_positive::Bool` An option to indicate if the set should be intersected with the set all positive numbers.
"""
mutable struct EllipsoidalSet <: AbstractUncertaintySet
    name::Symbol
    only_positive::Bool
end

"""
    HyperboxSet <: AbstractUncertaintySet
A hyperbox uncertainty set whose nomimal dimensions are determined by `neg_dev` and `pos_dev`.

**Fields**
- `name::Symbol` Name of the set which will be `:Hyperbox`.
- `neg_dev::Vector{Number}` A vector of the expected negative deviation of the random variables.
- `pos_dev::Vector{Number}` A vector of the expected positive deviation of the random variables.
- `only_positive::Bool` An option to indicate if the set should be intersected with the set all positive numbers.
"""
mutable struct HyperboxSet <: AbstractUncertaintySet
    name::Symbol
    neg_dev::Vector{Number}
    pos_dev::Vector{Number}
    only_positive::Bool
end

"""
    PNormSet <: AbstractUncertaintySet
A p-norm based uncertainty set based on a bounded p-norm.

**Fields**
- `name::Symbol` Name of the set which will be `:PNorm`.
- `p::Number` The value of p which can be 1, 2, or Inf.
- `only_positive::Bool` An option to indicate if the set should be intersected with the set all positive numbers.
"""
mutable struct PNormSet <: AbstractUncertaintySet
    name::Symbol
    p::Number
    only_positive::Bool
end

"""
    FlexibilityData

A DataType for storing the data necessary to manage the bookkeeping of the flexibility variables (`RandomVariable`
and `RecourseVariable`), the uncertainty set, and solution results.

**Fields**

- `flexibility_constraints::Vector{JuMP.AbstractConstraint}` Constraints that involve flexibility variables.
- `numRVs::Int` The number of `RandomVariable` that have been added to the model.
- `RVmeans::Vector{Number}` The means corresponding to each `RandomVariable`.
- `RVnames::Vector{AbstractString}` The symbolic name of each `RandomVariable`.
- `RVcols::Vector{Int}` The index of each `RandomVariable`.
- `num_recourse_vars::Int` The number of `RecourseVariable` that have been added to the model.
- `recourse_names::Vector{AbstractString}` The symbolic name of each `RecourseVariable`.
- `recourse_cols::Vector{Int}` The index of each `RecourseVariable`.
- `uncertainty_set::AbstractUncertaintySet` The uncertainty set DataType with all of the set specfic attributes.
- `covariance::Matrix{Number}` The covariance matrix.
- `flexibility_index::Union{Nothing, Number}` The flexibility index result obtained from solving the flexibility model.
- `active_constraints::Vector{Int}` The indexes of the active inequality constraints at the solution of the flexibility model.
- 'solution_time::Number' The solution time in seconds.
"""
mutable struct FlexibilityData
    flexibility_constraints::Vector{JuMP.AbstractConstraint}

    # Random variable data
    numRVs::Int
    RVmeans::Vector{Number}
    RVnames::Vector{AbstractString}
    RVcols::Vector{Int}

    # Recourse variable data
    num_recourse_vars::Int
    recourse_names::Vector{AbstractString}
    recourse_cols::Vector{Int}

    # Various formulation data/results
    uncertainty_set::AbstractUncertaintySet  #e.g. ellipsoidal, 1-norm, 2-norm, etc...
    covariance::Matrix{Number}
    flexibility_index::Union{Nothing, Number}
    active_constraints::Vector{Int}
    solution_time::Union{Nothing, Number}
end

# Set methods to extract covariance information if appropriate
EllipsoidalSet() = EllipsoidalSet(:Ellipsoid, false)
