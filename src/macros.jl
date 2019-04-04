using Base.Meta, Pkg
if Pkg.installed()["JuMP"] >= v"0.18.4"
    import JuMP.undef
end

"""
    @randomvariable(m, x, mean)
Defines a random variable using [`RandomVariable(m::Model, mean::Number, name::AbstractString)`](@ref)
and requires that a mean for the variable be provided. This can later be overwritten with [`setmean`](@ref).

**Arguments**
- `m::Model` The flexibility model.
- `x::Symbol` The variable name.
- `mean::Number` The variable mean

```julia
julia> @randomvariable(m2, w, mean = 42)

julia> @randomvariable(m2, ws[i = 1:4], mean = 42)
4-element Array{FlexibilityAnalysis.RandomVariable,1}:
 FlexibilityAnalysis.RandomVariable(Feasibility problem with:
 * 9 linear constraints
 * 89 variables

 julia> @randomvariable(m2, ws[i = 1:4], mean = [1; 2; 3; 4][i])
 4-element Array{FlexibilityAnalysis.RandomVariable,1}:
  FlexibilityAnalysis.RandomVariable(Feasibility problem with:
  * 9 linear constraints
  * 96 variables
```
"""
macro randomvariable(m, x, mean)
    m = esc(m)
    kwsymbol = VERSION < v"0.6.0-dev.1934" ? :kw : :(=) # changed by julia PR #19868
    if isexpr(mean,kwsymbol)
        @assert isexpr(mean,kwsymbol) && mean.args[1] == :mean
        mean = esc(mean.args[2])
    end

    if isa(x,Symbol)
        # easy case
        return quote
            $(esc(x)) = RandomVariable($m,$mean,$(string(x)))
            nothing
        end
    else
        if !isexpr(x,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end

        variable = gensym()
        refcall, idxvars, idxsets, idxpairs, condition = JuMP.buildrefsets(x, variable)
        varname = JuMP.getname(x)
        escvarname = esc(varname)

        varstr = :(string($(string(varname)),"["))
        for idxvar in idxvars
            push!(varstr.args,:(string($(esc(idxvar)))))
            push!(varstr.args,",")
        end
        deleteat!(varstr.args,length(varstr.args))
        push!(varstr.args,"]")

        code = :( $(refcall) = RandomVariable($m, $mean, $varstr ) )
        looped = JuMP.getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :RandomVariable)
        return quote
            $looped
            $escvarname = $variable
        end
    end
end

"""
    @recoursevariable(m, x)
Defines a recourse variable using [`RecourseVariable(m::Model, name::AbstractString)`](@ref).

**Arguments**
- `m::Model` The flexibility model.
- `x::Symbol` The variable name.

```julia
julia> @recoursevariable(m2, d)

julia> @recoursevariable(m2, ds[1:4])
4-element Array{FlexibilityAnalysis.RecourseVariable,1}:
 FlexibilityAnalysis.RecourseVariable(Feasibility problem with:
 * 9 linear constraints
 * 101 variables
```
"""
macro recoursevariable(m, x)
    m = esc(m)
    kwsymbol = VERSION < v"0.6.0-dev.1934" ? :kw : :(=) # changed by julia PR #19868

    if isa(x,Symbol)
        # easy case
        return quote
            $(esc(x)) = RecourseVariable($m,$(string(x)))
            nothing
        end
    else
        if !isexpr(x,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end

        variable = gensym()
        refcall, idxvars, idxsets, idxpairs, condition = JuMP.buildrefsets(x, variable)
        varname = JuMP.getname(x)
        escvarname = esc(varname)

        varstr = :(string($(string(varname)),"["))
        for idxvar in idxvars
            push!(varstr.args,:(string($(esc(idxvar)))))
            push!(varstr.args,",")
        end
        deleteat!(varstr.args,length(varstr.args))
        push!(varstr.args,"]")

        code = :( $(refcall) = RecourseVariable($m, $varstr ) )
        looped = JuMP.getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :RecourseVariable)
        return quote
            $looped
            $escvarname = $variable
        end
    end
end

"""
    JuMP.constructconstraint!(flex_aff::FlexibilityExpr, sense::Symbol)
Extends `JuMP.constructconstraint!` for `FlexibilityExpr` types.
"""
function JuMP.constructconstraint!(flex_aff::FlexibilityExpr, sense::Symbol)
    if sense == :(<=) || sense == :≤
        return FlexibilityConstraint(flex_aff, :(<=))
    elseif sense == :(>=) || sense == :≥
        return FlexibilityConstraint(flex_aff, :(>=))
    elseif sense == :(==)
        return FlexibilityConstraint(flex_aff, :(==))
    end
    error("Unrecognized constraint type $sense")
end
