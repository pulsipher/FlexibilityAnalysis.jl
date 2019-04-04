import Base: +, -, /, *, <=, >=

# Random Normal Variable
(+)(lhs::FlexibilityVariable,rhs::FlexibilityVariable) = FlexibilityExpr([lhs,rhs],[1.0,1.0],0.0)
(-)(lhs::FlexibilityVariable,rhs::FlexibilityVariable) = FlexibilityExpr([lhs,rhs],[1.0,-1.0],0.0)
(/)(lhs::FlexibilityVariable,rhs::Number) = FlexibilityExpr([lhs],[1/rhs],0.0)
(-)(rhs::FlexibilityVariable) = FlexibilityExpr([rhs],[-1.0],0.0)

for op in (:+, :-, :*)
    @eval begin
        ($op)(lhs::Number, rhs::FlexibilityVariable) = ($op)(lhs, FlexibilityExpr([rhs],[1.0],0.0))
        ($op)(lhs::Variable, rhs::FlexibilityVariable) = ($op)(convert(AffExpr,lhs),rhs)
        ($op)(lhs::Variable, rhs::FlexibilityExpr) = ($op)(convert(AffExpr,lhs),rhs)
    end
    if op == :-
        @eval ($op)(lhs::FlexibilityExpr, rhs::Variable) = (+)(lhs,-rhs)
        @eval ($op)(lhs::FlexibilityVariable, rhs::Union{Variable,Number}) = (+)(lhs,-rhs)
    else
        @eval ($op)(lhs::FlexibilityVariable, rhs::Union{Variable,Number}) = ($op)(rhs,lhs)
        @eval ($op)(lhs::FlexibilityExpr, rhs::Variable) = ($op)(rhs,lhs)
    end
end

# AffExpr
# AffExpr--FlexibilityVariable
(+)(lhs::AffExpr, rhs::FlexibilityVariable) = FlexibilityExpr([rhs],[convert(AffExpr,1.0)],lhs)
(-)(lhs::AffExpr, rhs::FlexibilityVariable) = FlexibilityExpr([rhs],[convert(AffExpr,-1.0)],lhs)
(*)(lhs::AffExpr, rhs::FlexibilityVariable) = FlexibilityExpr([rhs],[lhs],AffExpr())
(+)(lhs::FlexibilityVariable, rhs::AffExpr) = rhs+lhs
(-)(lhs::FlexibilityVariable, rhs::AffExpr) = (+)(lhs,-rhs)
(*)(lhs::FlexibilityVariable, rhs::AffExpr) = rhs*lhs

# AffExpr--FlexibilityExpr
(+)(lhs::AffExpr,rhs::FlexibilityExpr) = convert(FlexibilityExpr,lhs)+rhs
(+)(lhs::FlexibilityExpr,rhs::AffExpr) = rhs+lhs
(-)(lhs::AffExpr,rhs::FlexibilityExpr) = convert(FlexibilityExpr,lhs)-rhs
(-)(lhs::FlexibilityExpr,rhs::AffExpr) = lhs-convert(FlexibilityExpr,rhs)
Base.promote_rule(::Type{AffExpr},::Type{FlexibilityExpr}) = FlexibilityExpr
Base.promote_rule(::Type{Float64},::Type{FlexibilityExpr}) = FlexibilityExpr
Base.convert(::Type{FlexibilityExpr},a::AffExpr) = FlexibilityExpr(FlexibilityVariable[],AffExpr[],a)

# AffExpr--FlexibilityExpr
(*)(lhs::AffExpr,rhs::FlexibilityExpr) = FlexibilityExpr(rhs.vars, [lhs*c for c in rhs.coeffs], lhs*rhs.constant)
(*)(lhs::FlexibilityExpr,rhs::AffExpr) = rhs*lhs
