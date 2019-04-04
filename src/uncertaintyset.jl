"""
    setuncertaintyset(m::Model, uncertainty_set::Symbol, [attribute = nothing; only_positive::Bool = false])
Specify the type of uncertainty set to be used and stored in `FlexibilityData.uncertainty_set` and provide the necessary
attribute. An ellipsoidal uncertainty set can be specified with the `:Ellipsoid` symbol and the corresponding covariance
matrix will need to input via `attribute` if it has not already been set with [`setcovariance`](@ref). A hyperbox uncertainty
set is specified with the `:Hyperbox` symbol and corresponding negative/positive deviation vectors need to be inputed via
`attribute` as a vector of vectors of the form `[[neg_dev]; [pos_dev]]`. Finally, a p-norm uncertainty set can be specified
with the `:PNorm` symbol and providing the corresponding `p` value via `attribute` where `p` can equal `1`, `2`, or `Inf`.

**Arguments**
- `m::Model` The flexibility model.
- `uncertainty_set::Symbol` The uncertainty set name.
- `attribute = nothing` The necessary atribute for the specified uncertainty set.

**Keyword Arguments**
- `only_positive::Bool = false` Indicate if the uncertainty set should be intersected with the set of all positive values.

```julia
julia> setuncertaintyset(m, :Ellipsoid)

julia> setuncertaintyset(m, :Ellipsoid, eye(4)* 11.11)

julia> setuncertaintyset(m, :Hyperbox, [[ones(4)]; [ones(4)]])
FlexibilityAnalysis.HyperboxSet(:Hyperbox, Number[1.0, 1.0, 1.0, 1.0], Number[1.0, 1.0, 1.0, 1.0], false)

julia> setuncertaintyset(m, :PNorm, 1)
FlexibilityAnalysis.PNormSet(:PNorm, 1, false)

julia> setuncertaintyset(m, :PNorm, 2, only_positive = true)
FlexibilityAnalysis.PNormSet(:PNorm, 2, true)

julia> setuncertaintyset(m, :PNorm, Inf)
FlexibilityAnalysis.PNormSet(:PNorm, Inf, false)
```
"""
function setuncertaintyset(m::Model, uncertainty_set::Symbol, attribute = nothing; only_positive::Bool = false)
    # Extract flexibility data
    flex_data = getflexibilitydata(m)
    attr_type = typeof(attribute)

    # Setup ellipdoial set if specified and run checks
    if uncertainty_set == :Ellipsoid
        if isa(attribute, Nothing) && length(flex_data.covariance) == 0
            error("Ellipsoidal set requires a covariance matrix, but one is not provided.")
        elseif !isa(attribute, Matrix) && !isa(attribute, Nothing)
            error("Expected ellipsoidal attribute to be covariance matrix of type Matrix, but got attribute of type $attr_type.")
        end
        flex_data.uncertainty_set = EllipsoidalSet()
        if !isa(attribute, Nothing)
            setcovariance(m, attribute)
        end
        flex_data.uncertainty_set.only_positive = only_positive; nothing

    # Setup hyperbox set if specified and run checks
    elseif uncertainty_set == :Hyperbox
        if !isa(attribute, Union{Vector{Vector{Float64}}, Vector{Vector{Float32}}, Vector{Vector{Int}}})
            error("Expected hyperbox attribute to be deviations of form [[neg_dev]; [pos_dev]], but got attribute of type $attr_type.")
        elseif length(attribute[1]) != length(attribute[2])
            error("The dimensions of pos_dev and neg_dev do not match.")
        end
        flex_data.uncertainty_set = HyperboxSet(uncertainty_set, attribute[1], attribute[2], only_positive)

    # Setup a P-norm set if specified and run checks
    elseif uncertainty_set == :PNorm
        if !isa(attribute, Number)
            error("Expected pnorm attribute to be power of type Number, but got attribute of type $attr_type.")
        elseif attribute != 1 && attribute != 2 && attribute != Inf
            error("$attribute is not a valid value of p. Currently, only 1, 2, and Inf are accepted.")
        end
        flex_data.uncertainty_set = PNormSet(uncertainty_set, attribute, only_positive)
    else
        error("$uncertainty_set is not a valid uncertainty set.")
    end
end
