"""
    ismeanfeasible(m::Model; [toler::Number = 1e-5, solver = Clp.ClpSolver()])
Returns a `Bool` indicating if the mean stored in `FlexibilityData.RVmeans` is feasible, meaning that it
lies inside the feasible region. This check is done using the so-called feasibility function.

**Arguments**
- `m::Model` The flexibility model.

**Keyword Arguments**
- `toler::Number = 1e-5` The numerical tolerance for checking the feasibility.
- `solver = Clp.ClpSolver()` The solver, any LP or NLP solver shoudl work.

```julia
julia> ismeanfeasible(m)
true

julia> ismeanfeasible(m, solver = GurobiSolver(OutputFlag = 0))
Academic license - for non-commercial use only
true
```
"""
function ismeanfeasible(m::Model; toler::Number = 1e-5, solver = Clp.ClpSolver())
    # Make the input dictionary
    input_dict = MakeInputDict(m)

    # Parse the necessary data dimensions
    n_z = input_dict["dims"][1]
    n_θ = input_dict["dims"][2]
    n_x = input_dict["dims"][3]
    n_f = input_dict["dims"][4]
    n_h = input_dict["dims"][5]

    # Prepare the JuMP model
    m_solve = Model(solver = solver)

    # Initialize variables
    @variable(m_solve, θ[1:n_θ])
    @variable(m_solve, u)
    if n_z != 0
        @variable(m_solve, z[1:n_z])
    end
    if n_x != 0
        @variable(m_solve, x[1:n_x])
    end

    # Set the objective
    @objective(m_solve, Min, u)

    # Fix θ and set values
    for i = 1:n_θ
        JuMP.fix(θ[i], input_dict["theta_nom"][i])
    end

    # Add the constraints
    exprs = AddSystemExpressions(m_solve, input_dict)
    @constraint(m_solve, f[j = 1:n_f], exprs[1][j] <= u)
    if n_h != 0
        @constraint(m_solve, h[i = 1:n_h], exprs[2][i] == 0)
    end

    # Solve and return
    status = solve(m_solve)
    return getvalue(u) <= toler
end

"""
    ComputeCenter(m::Model, center::Symbol, solver, toler::Number, only_positive::Bool)
Returns a center point that can be used to replace the mean if desired.

**Arguments**
- `m::Model` The flexibility model.
- `center::Symbol` Indicates the type of center, accepted arguments are `:feasible` and `:analytic`.
- `solver` The solver which must be an NLP solver for the analytic center.
- `toler::Number` The tolerance to check solution validity.
- `only_positive::Bool` Indicates if the center need by strictly positive.
"""
function ComputeCenter(m::Model, center::Symbol, solver, toler::Number, only_positive::Bool)
    # Check input
    if center != :analytic && center != :feasible
        error("Invalid argument center = $center, must be :analytic or :feasible")
    end

    # Make the input dictionary
    input_dict = MakeInputDict(m)

    # Parse the necessary data dimensions
    n_z = input_dict["dims"][1]
    n_θ = input_dict["dims"][2]
    n_x = input_dict["dims"][3]
    n_f = input_dict["dims"][4]
    n_h = input_dict["dims"][5]

    # Prepare the JuMP model
    m_solve = Model(solver = solver)

    # Initialize variables
    if only_positive
        @variable(m_solve, θ[1:n_θ] >= 0)
    else
        @variable(m_solve, θ[1:n_θ])
    end
    if center == :feasible
        @variable(m_solve, u)
    else
        @variable(m_solve, s[1:n_f] >= 0)
    end
    if n_z != 0
        @variable(m_solve, z[1:n_z])
    end
    if n_x != 0
        @variable(m_solve, x[1:n_x])
    end

    # Set objective
    if center == :feasible
        @objective(m_solve, Min, u)
    else
        @NLobjective(m_solve, Min, -sum(log(s[i]) for i = 1:n_f))
    end

    # Add constraints
    exprs = AddSystemExpressions(m_solve, input_dict)
    if center == :feasible
        @constraint(m_solve, f[j = 1:n_f], exprs[1][j] <= u)
    else
        @constraint(m_solve, f[j = 1:n_f], exprs[1][j] + s[j] == 0)
    end
    if n_h != 0
        @constraint(m_solve, h[i = 1:n_h], exprs[2][i] == 0)
    end

    # Solve and return
    status = solve(m_solve)
    if center == :feasible && getvalue(u) > toler
        upper_b = getvalue(u)
        @warn "Optimized mean is not feasible, can only achieve inequalities with upper bound u = $upper_b"
    end
    return getvalue(θ)
end

"""
    findcenteredmean(m::Model; [center::Symbol = :feasible, solver = Ipopt.IpoptSolver(print_level = 0), toler::Number = 1e-5, update_mean::Bool = false, only_positive::Bool = false])
Returns a center point based on the analytic or feasible center. The result can overwrite the mean stored in `FlexibilityData.RVmeans`
if desired. This is a wrapper function for [`ComputeCenter`](@ref).

**Arguments**
- `m::Model` The flexibility model.

**Keyword Arguments**
- `center::Symbol = :feasible` Indicates the type of center, accepted arguments are `:feasible` and `:analytic`.
- `solver = Ipopt.IpoptSolver(print_level = 0)` The solver which must be an NLP solver for the analytic center.
- `toler::Number = 1e-5` The tolerance to check solution validity.
- `update_mean::Bool = false` Indicates if the computed center should overwrite `FlexibilityData.RVmeans`.
- `only_positive::Bool = false` Indicates if the center need by strictly positive.

```julia
julia> findcenteredmean(m, only_positive = true)
4-element Array{Float64,1}:
 1684.74
   79.0718
  195.073
    0.0

julia> findcenteredmean(m, center = :analytic, update_mean = true)
4-element Array{Float64,1}:
  898.125
 -507.214
  594.544
  317.23
```
"""
function findcenteredmean(m::Model; center::Symbol = :feasible, solver = Ipopt.IpoptSolver(print_level = 0), toler::Number = 1e-5, update_mean::Bool = false, only_positive::Bool = false)
    center = ComputeCenter(m, center, solver, toler, only_positive) # included in functions.jl
    if update_mean
        flex_data = getflexibilitydata(m)
        flex_data.RVmeans = center
    end
    return center
end
