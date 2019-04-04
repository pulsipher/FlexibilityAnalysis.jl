"""
    rankinequalities(m::Model; max_ranks::Int = 5, suppress_warnings::Bool = true, U::Int = 10000, diag::Bool = false, active_constr::Bool = false, real_recourse_dim::Int = -1, conic_δ::Bool = false)
Returns ranking data in the form `Vector{Dict}` where each dictionary corresponds to a particular rank level and
contains the optimal flexibility_index, active constraint indexes, and flexibility model. The function will iteratively
solve copies of the flexibility model via `solve` where the prior active constraints are turned off in order rank the most limiting
constraints. The user can specify the maximum number of rank levels and the flexibility index problem will be repeatedly
solved until that maximum is acheived or the problem becomes unbounded, which occurs first.

**Arguments**
- `m::Model` The flexibility model.

**Keyword Arguments**
- `max_ranks::Int = 5` The maximum number of rank levels. 2
- `suppress_warnings::Bool = false` Indicates if solver warnings should be suppressed.
- `U::Union{Int, Float64} = 10000` The slack variable upper bound.
- `diag::Bool = false` Indicates whether or not to diagnonalize ellipsoidal uncertainty set (this is only active when an ellipsoidal set is used).
- `active_constr::Bool = false` Indicates if the optional active constraint should be used which enforces how many inequalities are active at the solution, this must be set to `true` for systems without control variables and/or contain state variables.
- `real_recourse_dim::Int = -1` The actual number of recourse variables in case state variables are included as recourse variables. This is mandatory if `active_constr = true` and no state variables are provided.
- `conic_δ::Bool = false` This should be set to `true` if a conic solver is used such as Pajarito.jl.

```julia
julia> rankinequalities(m, active_constr = true)
2-element Array{Dict,1}:
 Dict{String,Any}(Pair{String,Any}("flexibility_index", 3.59938),Pair{String,Any}("model", Feasibility problem with:
 * 6 linear constraints
 * 6 variables
Solver is Pavito),Pair{String,Any}("active_constraints", [3, 6]))
 Dict{String,Any}(Pair{String,Any}("flexibility_index", 9.58983),Pair{String,Any}("model", Feasibility problem with:
 * 6 linear constraints
 * 6 variables
Solver is Pavito),Pair{String,Any}("active_constraints", [1, 5]))
```
"""
function rankinequalities(m::Model; max_ranks::Int = 5, suppress_warnings::Bool = true, U::Int = 10000, diag::Bool = false, active_constr::Bool = false, real_recourse_dim::Int = -1, conic_δ::Bool = false)
    # Check the max_ranks
    if max_ranks < 1
        error("max_ranks set to $max_ranks, but must be an integer greater than 0.")
    end

    # Initialize data
    m_copy = deepcopy(m)
    ranked_data = Vector{Dict}(undef, max_ranks)
    ranked_constrs = []
    inactives = []
    counter = 0

    # Get the inequality information
    constr_bounds = JuMP.prepConstrBounds(m_copy)
    inequal_inds = findall(constr_bounds[1] .!= constr_bounds[2])

    # Iteratively solve the flexibility index problem and extract data
    for i = 1:max_ranks
        status = solve(m_copy, suppress_warnings = suppress_warnings, U = U, diag = diag, active_constr = active_constr, real_recourse_dim = real_recourse_dim, conic_δ = conic_δ, inactives = inactives)
        if status == :Optimal
            ranked_constrs = [findall(inequal_inds .== getactiveconstraints(m_copy)[j])[1] for j = 1:length(getactiveconstraints(m_copy))]
            inactives = unique([inactives; ranked_constrs])
            ranked_data[i] = Dict("flexibility_index" => getflexibilityindex(m_copy), "active_constraints" => getactiveconstraints(m_copy), "model" => m_copy)
            counter += 1
        else
            break
        end
        if length(inactives) == length(inequal_inds)
            break
        end
    end

    # Return the results
    return ranked_data[1:counter]
end

"""
    getconfidencelevel(m::Model)
Return the confidence level provided that the flexibility model `m` has been solved with an
ellipsoidal uncertainty set. This is equivalent to the quantile associated with the optimized
uncertainty set. Note that this assumes a multivariate Gaussian distribution with mean
`FlexibilityData.RVmeans` and covariance `FlexibilityData.covariance`.

**Arguments**
- `m::Model` The flexibility model.

```julia
julia> getconfidencelevel(m)
0.5370703369769008
```
"""
function getconfidencelevel(m::Model)
    flex_data = getflexibilitydata(m)
    if flex_data.uncertainty_set.name == :Ellipsoid && flex_data.flexibility_index != nothing
        return cdf(Chisq(flex_data.numRVs), flex_data.flexibility_index)
    elseif flex_data.flexibility_index == nothing
        error("Cannot compute the confidence level without first solving the flexibility model.")
    else
        error("Invalid uncertainty set, must solve the flexibility model with an ellipsoidal set.")
    end
end

"""
    findstochasticflexibility(m::Model; num_pts::Int = 10000, toler::Number = 1e-5, solver = Clp.ClpSolver(), only_positive::Bool = false, use_vulnerability_model::Bool = false, use_flexibility_index::Bool = false, seed::Int = -1)
Returns the estimated stochastic flexibility index that is evaluated via Monte Carlo sampling. At default
this estimation is carried out by evaluating the feasibility of each Monte Carlo sample. The samples are
generating from a multivariate Gaussian distribution with mean `FlexibilityData.RVmeans` and covariance
`FlexibilityData.covariance`. The vulnerability model also tests the feasibility of the samples, but does
so in one large optimization problem instead of evaluating each sample indiviually. The optimized flexibility
index can also be used to reduce the number of samples that need to be evaluated.

**Arguments**
- `m::Model` The flexibility model.

**Keyword Arguments**
- `num_pts::Int = 10000` Number of Monte Carlo samples.
- `toler::Number = 1e-5` The feasibility check tolerance.
- `solver = Clp.ClpSolver()` The solver, any LP or NLP solver should work.
- `only_positive::Bool = false` Indicates if only positive samples should be used.
- `use_vulnerability_model::Bool = false` Indicates if the vulnerability model should be used.
- `use_flexibility_index::Bool = false` Indicates if the optimal flexibility index should be used.
- `seed::Int = -1` Random seed for sample collection, any negative value will turn off the random seed.

```julia
julia> findstochasticflexibility(m)
0.9687

julia> findstochasticflexibility(m, use_vulnerability_model = true)
0.9705

julia> findstochasticflexibility(m, num_pts = 5000, use_flexibility_index = true)
0.973
```
"""
function findstochasticflexibility(m::Model; num_pts::Int = 10000, toler::Number = 1e-5, solver = Clp.ClpSolver(), only_positive::Bool = false, use_vulnerability_model::Bool = false, use_flexibility_index::Bool = false, seed::Int = -1)
    # Parse data
    flex_data = getflexibilitydata(m)

    if length(flex_data.covariance) == 0
        error("The covariance matrix needs to be set.")
    end

    # Load in the neceassary data
    input_dict = MakeInputDict(m)
    θ_nom = input_dict["theta_nom"]

    # Parse the necessary data dimensions
    n_z = input_dict["dims"][1]
    n_θ = input_dict["dims"][2]
    n_x = input_dict["dims"][3]
    n_f = input_dict["dims"][4]
    n_h = input_dict["dims"][5]

    # Setup a multivariate normal distribution and Monte Carlo samples
    if seed >= 0
        Random.seed!(seed)
    end
    d = MvNormal(Vector{Float64}(flex_data.RVmeans), Matrix{Float64}(getcovariance(m)))
    if only_positive
        samples = Array{Float64}(undef, flex_data.numRVs, num_pts)
        not_done = true
        counter = 1
        while not_done
            sample = rand(d, 1)
            if all(sample .>= 0)
                samples[:, counter] = sample
                counter += 1
            end
            if counter > num_pts
                not_done = false
            end
        end
    else
        samples = rand(d, num_pts)
    end

    # Determine which points are outside of the set and check for F
    F = getflexibilityindex(m)
    if use_flexibility_index && F == nothing
        @warn "Flexibility index hasn't yet been computed. Setting use_flexibility_index = false."
        use_flexibility_index = false
    elseif use_flexibility_index && flex_data.uncertainty_set.name == :Ellipsoid
        inv_covar = inv(getcovariance(m))
        outside_set = [(samples[:, k] - θ_nom)' * inv_covar * (samples[:, k] - θ_nom) > F for k = 1:num_pts]
    elseif use_flexibility_index && flex_data.uncertainty_set.name == :PNorm
        outside_set = [norm(samples[:, k] - θ_nom, flex_data.uncertainty_set.p) > F for k = 1:num_pts]
    elseif use_flexibility_index && flex_data.uncertainty_set.name == :Hyperbox
        outside_set = Vector{Bool}(undef, num_pts)
        for k = 1:num_pts
            all_inside = true
            for i = 1:flex_data.numRVs
                if samples[i, k] < θ_nom[i] - flex_data.uncertainty_set.neg_dev[i] * F || samples[i, k] > θ_nom[i] + flex_data.uncertainty_set.pos_dev[i] * F
                    all_inside = false
                end
            end
            outside_set[k] = !all_inside
        end
    end

    # Determine the points to test
    if use_flexibility_index
        test_pts = samples[:, outside_set]
        num_test_pts = size(test_pts)[2]
    else
        num_test_pts = num_pts
        test_pts = samples
    end

    # Stop if all the test points are feasible.
    if num_test_pts == 0
        return 1.0
    end

    # Setup the model to be solved
    m_solve = Model(solver = solver)

    if use_vulnerability_model
        # Initialize variables
        @variable(m_solve, θ[1:n_θ, 1:num_test_pts])
        @variable(m_solve, u[1:num_test_pts] >= 0)
        if n_z != 0
            @variable(m_solve, z[1:n_z, 1:num_test_pts])
        end
        if n_x != 0
            @variable(m_solve, x[1:n_x, 1:num_test_pts])
        end

        # Set the objective
        @objective(m_solve, Min, 1 / num_test_pts * sum(u[k] for k = 1:num_test_pts))

        # Add the constraints
        exprs = AddSystemExpressions(m_solve, input_dict, num_test_pts)
        @constraint(m_solve, f[j = 1:n_f, k = 1:num_test_pts], exprs[1][j, k] <= u[k])
        if n_h != 0
            @constraint(m_solve, h[i = 1:n_h, k = 1:num_test_pts], exprs[2][i, k] == 0)
        end
    else
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

        # Add the constraints
        exprs = AddSystemExpressions(m_solve, input_dict)
        @constraint(m_solve, f[j = 1:n_f], exprs[1][j] <= u)
        if n_h != 0
            @constraint(m_solve, h[i = 1:n_h], exprs[2][i] == 0)
        end
    end

    # Solve the model(s) and compute the SF
    if use_vulnerability_model
        for i = 1:n_θ
            for k = 1:num_test_pts
                JuMP.fix(θ[i, k], test_pts[i, k])
            end
        end
        status = solve(m_solve)
        if status == :Optimal
            return 1 - sum(getvalue(u) .>= toler) / num_pts
        else
            error("Problem not solved to optimality, cannot estimate SF.")
        end
    else
        infeasible_counter = 0
        feasible_results = zeros(Bool, num_test_pts)
        for k = 1:num_test_pts
            for i = 1:n_θ
                JuMP.fix(θ[i], test_pts[i, k])
            end
            status = solve(m_solve)
            if status == :Optimal
                if getvalue(u) <= toler
                    feasible_results[k] = true
                end
            else
                infeasible_counter += 1
            end
        end
        if infeasible_counter != 0
            @warn "Not all scenario subproblems not solved to optimality, estmation of SF might not be correct."
        end
        if use_flexibility_index
            return (sum(feasible_results) + num_pts - num_test_pts) / (num_pts - infeasible_counter)
        else
            return sum(feasible_results) / (num_pts - infeasible_counter)
        end
    end
end
