"""
    Base.show(io::IO, a::FlexibilityExpr)
Extend `Base.show` to print flexibility expressions.
"""
function Base.show(io::IO, a::FlexibilityExpr)
    # Setup empty case
    if length(a.vars) == 0
        return print(io, a.constant)
    end

    # Otherwise get data and parse it
    m = a.vars[1].m
    flex_data = getflexibilitydata(m)
    vars = a.vars
    names = []
    coeffs = []
    for i = 1:length(vars)
        push!(coeffs, a.coeffs[i])
        var = vars[i]
        if isa(var, RecourseVariable)
            push!(names, flex_data.recourse_names[var.idx])
        elseif isa(var, RandomVariable)
            push!(names, flex_data.RVnames[var.idx])
        end
    end

    # Print the expressions to the io
    strs = ["$(JuMP.aff_str(JuMP.REPLMode,coeffs[i], true))*$(names[i])" for i in 1:length(a.vars)]
    print(io,string(join(strs," + "), " + ", JuMP.aff_str(JuMP.REPLMode, a.constant, true)))
end

"""
    JuMP.addconstraint(m::Model, constr::FlexibilityConstraint)
Extend the `JuMP.addconstraint` function to handle `FlexibilityConstraint` types.
"""
function JuMP.addconstraint(m::Model, constr::FlexibilityConstraint)
    # Get the data
    flex_data = getflexibilitydata(m)
    push!(flex_data.flexibility_constraints, constr)

    # Add constraint to model here using the flexibility constraint information
    expr = constr.flex_expr
    sense = constr.sense

    #NOTE Look into why these coeffcients get conflated as expressions
    flex_coeffs = [expr.coeffs[i].constant for i = 1:length(expr.coeffs)]
    flex_vars = expr.vars

    # Parse the variable information
    flex_vars_jump = []
    for var in flex_vars
        if isa(var, RecourseVariable)
            jump_var = Variable(m, flex_data.recourse_cols[linearindex(var)])
            push!(flex_vars_jump, jump_var)
        elseif isa(var, RandomVariable)
            jump_var = Variable(m, flex_data.RVcols[linearindex(var)])
            push!(flex_vars_jump, jump_var)
        else
            error("Variable type $(typeof(var)) not recognized for constructing constraint")
        end
    end

    # Get non flexibility variable information
    state_coeffs = expr.constant.coeffs
    state_vars = expr.constant.vars
    constant = expr.constant.constant #should be the right hand side of the constraint

    # Add constraints to the model using the anonymous variables
    all_coeffs = vcat(flex_coeffs,state_coeffs)
    all_vars = vcat(flex_vars_jump,state_vars)
    if sense == :(<=)
        con_reference = @constraint(m, sum(all_coeffs[i]*all_vars[i] for i = 1:length(all_coeffs)) + constant <= 0)
    elseif sense == :(>=)
        con_reference = @constraint(m, sum(all_coeffs[i]*all_vars[i] for i = 1:length(all_coeffs)) + constant >= 0)
    elseif sense == :(==)
        con_reference = @constraint(m, sum(all_coeffs[i]*all_vars[i] for i = 1:length(all_coeffs)) + constant == 0)
    else
        error("Constraint sense $(sense) not recognized")
    end

    return ConstraintRef{JuMP.Model,FlexibilityConstraint}(m, length(flex_data.flexibility_constraints))
end

"""
    JuMP.show(io::IO,c::FlexibilityConstraint)
Extend the `JuMP.show` function to handle `FlexibilityConstraint` types.
"""
function JuMP.show(io::IO,c::FlexibilityConstraint)
    s = "$(string(c.flex_expr)) $(c.sense) 0"
    return print(io,s)
end
