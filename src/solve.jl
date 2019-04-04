import MathProgBase

"""
    solvehook(m::Model; [suppress_warnings::Bool = false, U::Number = 10000, diag::Bool = false, active_constr::Bool = false, real_recourse_dim::Int = -1, conic_δ::Bool = false, inactives::Vector = []])
Returns the solution status to solving the flexibility model `m`. This solvehook what `JuMP.solve(::Model)` for flexibility models.
This solves the flexibility index problem using the variables and constraints specified in `m`.

**Arguments**
- `m::Model` The flexibility model.

**Keyword Arguments**
- `suppress_warnings::Bool = false` Indicates if solver warnings should be suppressed.
- `U::Number = 10000` The slack variable upper bound.
- `diag::Bool = false` Indicates if the ellipsoidal uncertainty set is diagonalized (this is only active when an ellipsoidal set is used).
- `active_constr::Bool = false` Indicates if the optional active constraint should be used which enforces how many inequalities are active at the solution, this must be set to `true` for systems without control variables and/or contain state variables.
- `real_recourse_dim::Int = -1` The actual number of recourse variables in case state variables are included as recourse variables. This is mandatory if `active_constr = true` and no state variables are provided.
- `conic_δ::Bool = false` This should be set to `true` if a conic solver is used such as Pajarito.jl.
- `inactives::Vector = []` The indexes of inequality constraints that should be turned off.
"""
function solvehook(m::Model; suppress_warnings::Bool = false, U::Number = 10000, diag::Bool = false, active_constr::Bool = false,
                   real_recourse_dim::Int = -1, conic_δ::Bool = false, inactives::Vector = [])

    # Pull in the constraint data from the JuMP model
    flex_data = getflexibilitydata(m)

    # Read in the data and extract it from the dictionary
    input_dict = MakeInputDict(m)
    fConsts = input_dict["fConsts"]
    fControls = input_dict["fControls"]
    fRandoms = input_dict["fRandoms"]
    fStates = input_dict["fStates"]
    hConsts = input_dict["hConsts"]
    hControls = input_dict["hControls"]
    hRandoms = input_dict["hRandoms"]
    hStates = input_dict["hStates"]
    θ_nom = input_dict["theta_nom"]
    n_z = input_dict["dims"][1]
    n_θ = input_dict["dims"][2]
    n_x = input_dict["dims"][3]
    n_f = input_dict["dims"][4]
    n_h = input_dict["dims"][5]
    if n_x != 0
        state_cols = input_dict["state_cols"]
    end
    inequal_inds = input_dict["inequal_inds"]

    # Check that U is sufficiently large
    if maximum(abs.(fConsts)) >= U
        U_new = 10 ^ ceil(log10(maximum(abs.(fConsts))))
        @warn "The slack upper bound (U) value of $U is too low for this problem. Thus, U is being set to $U_new !"
        U = U_new
    end

    # Inilialize the JuMP model to solve the MICP Problem
    m_solve = Model(solver = m.solver)

    # Initialize the variables
    @variable(m_solve, θ[1:n_θ])
    @variable(m_solve, δ >= 0)
    @variable(m_solve, s[1:n_f] >= 0)
    @variable(m_solve, y[1:n_f], Bin)
    if n_z != 0
        @variable(m_solve, z[1:n_z])
        @variable(m_solve, λ[1:n_f] >= 0)
        if n_h != 0
            @variable(m_solve, μ[1:n_h])
        end
    end
    if n_x != 0
        @variable(m_solve, x[1:n_x])
    end

    # Set objective function
    @objective(m_solve, Min, δ)

    # Turn off constraints that are set to be inactive
    if length(inactives) != 0
        for i = 1:length(inactives)
            JuMP.fix(y[inactives[i]], 0)
        end
    end

    # Add slack variable constraints
    @constraint(m_solve, sUyconst[j = 1:n_f], s[j] - U * (1 - y[j]) <= 0)

    # Add active set contraint
    if real_recourse_dim != -1
        active_constr = true
    end
    if active_constr || n_x != 0 || n_z == 0
        if !active_constr
            @warn "Problem requires use of active set constraint. Setting active_constr = true"
        end
        if real_recourse_dim == -1
            real_recourse_dim = n_z
            if n_h != 0 && n_x == 0
                @warn "real_recourse_dim not specified. Setting real_recourse_dim to $n_z"
            end
        end
        @constraint(m_solve, sum(y[j] for j = 1:n_f) == real_recourse_dim + 1)
    end

    # Add system constraints
    exprs = AddSystemExpressions(m_solve, input_dict)
    if length(inactives) != 0
        actives = deleteat!(Vector(1:n_f), sort(inactives))
        @constraint(m_solve, f[j = actives], exprs[1][j] + s[j] == 0)
    else
        @constraint(m_solve, f[j = 1:n_f], exprs[1][j] + s[j] == 0)
    end
    if n_h != 0
        @constraint(m_solve, h[i = 1:n_h], exprs[2][i] == 0)
    end

    # Add stationary conditions
    if n_z != 0
        @constraint(m_solve, sum(λ[j] for j = 1:n_f) == 1)
        if n_h != 0
            @constraint(m_solve, deriv[i = 1:n_z], sum(λ[j] * fControls[j, i] for j = 1:n_f) + sum(μ[j] * hControls[j, i] for j = 1:n_h) == 0)
        else
            @constraint(m_solve, deriv[i = 1:n_z], sum(λ[j] * fControls[j, i] for j = 1:n_f) == 0)
        end
    end

    # Add binary variable constraints
    if n_z != 0
        @constraint(m_solve, lyconst[j = 1:n_f], λ[j] - y[j] <= 0)
    end

    ## Setup the uncertainty set
    # Add hyperbox contraints if specified
    sqrt_δ = false
    if flex_data.uncertainty_set.name == :Hyperbox
        if length(flex_data.uncertainty_set.pos_dev) != n_θ
            error("The dimensions of random variables and the uncertainty set do not match.")
        end
        @constraint(m_solve, upper[i = 1:n_θ], θ[i] <= θ_nom[i] + flex_data.uncertainty_set.pos_dev[i] * δ)
        @constraint(m_solve, lower[i = 1:n_θ], θ[i] >= θ_nom[i] - flex_data.uncertainty_set.neg_dev[i] * δ)

    # Add p-norm constraints if specified
    elseif flex_data.uncertainty_set.name == :PNorm
        p = flex_data.uncertainty_set.p
        if p == 1
            @variable(m_solve, aux[1:n_θ] >= 0)
            @constraint(m_solve, abs_upper[i = 1:n_θ], aux[i] >= θ[i] - θ_nom[i])
            @constraint(m_solve, abs_lower[i = 1:n_θ], aux[i] >= -(θ[i] - θ_nom[i]))
            @constraint(m_solve, sum(aux[i] for i = 1:n_θ) <= δ)
        elseif p == 2
            @constraint(m_solve, sum((θ - θ_nom)[i]^2 for i = 1:n_θ) <= δ)
            sqrt_δ = true
        else
            @variable(m_solve, aux[1:n_θ] >= 0)
            @variable(m_solve, max_aux >= 0)
            @constraint(m_solve, abs_upper[i = 1:n_θ], aux[i] >= θ[i] - θ_nom[i])
            @constraint(m_solve, abs_lower[i = 1:n_θ], aux[i] >= -(θ[i] - θ_nom[i]))
            @constraint(m_solve, maxs[i = 1:n_θ], max_aux >= aux[i])
            @constraint(m_solve, max_aux <= δ)
        end

    # Add constraints for diagonalized region
    elseif diag && n_θ > 1
        covar = getcovariance(m)
        if n_θ != size(covar)[1]
            error("The dimensions of the covariance matrix and the random variables do not match.")
        end
        D, V = eigen(inv(covar))
        @variable(m_solve, w[1:n_θ])
        @constraint(m_solve, ws[i = 1:n_θ], w[i] == sqrt(D[i]) * sum(V[j, i] * (θ - θ_nom)[j] for j = 1:n_θ))
        if conic_δ
            @constraint(m_solve, sum(w[i] * w[i] for i = 1:n_θ) <= δ^2)
        else
            @constraint(m_solve, sum(w[i] * w[i] for i = 1:n_θ) <= δ)
        end

    # Add constraint for normal undiagonalized region
    elseif n_θ > 1
        covar = getcovariance(m)
        if n_θ != size(covar)[1]
            error("The dimensions of the covariance matrix and the random variables do not match.")
        end
        inv_covar = inv(covar)
        if conic_δ
            @constraint(m_solve, sum(sum((θ - θ_nom)[j] * inv_covar[j, i] * (θ - θ_nom)[i] for j = 1:n_θ) for i = 1:n_θ) <= δ^2)
        else
            @constraint(m_solve, sum(sum((θ - θ_nom)[j] * inv_covar[j, i] * (θ - θ_nom)[i] for j = 1:n_θ) for i = 1:n_θ) <= δ)
        end

    # Add constraint for speical case that there is only 1 random variable
    else
        covar = getcovariance(m)
        if n_θ != size(covar)[1]
            error("The dimensions of the covariance matrix and the random variables do not match.")
        end
        if isa(covar, Array)
            if conic_δ
                @constraint(m_solve, (θ[1] - θ_nom[1]) * (1 / covar[1]) * (θ[1] - θ_nom[1]) <= δ^2)
            else
                @constraint(m_solve, (θ[1] - θ_nom[1]) * (1 / covar[1]) * (θ[1] - θ_nom[1]) <= δ)
            end
        else
            if conic_δ
                @constraint(m_solve, (θ[1] - θ_nom) * (1 / covar) * (θ[1] - θ_nom) <= δ^2)
            else
                @constraint(m_solve, (θ[1] - θ_nom) * (1 / covar) * (θ[1] - θ_nom) <= δ)
            end
        end
    end

    # Enforce the intersection of the traditional set with the positive θ set
    if flex_data.uncertainty_set.only_positive
        @constraint(m_solve, pos_theta[i = 1:n_θ], θ[i] >= 0)
    end

    # Solve the model
    tic_time = @elapsed status = solve(m_solve, suppress_warnings = suppress_warnings, ignore_solve_hook = true)

    if status == :Optimal
        # Parse the optimized random variable values
        theta = getvalue(θ)
        for i = 1:n_θ
            setvalue(Variable(m, flex_data.RVcols[i]), theta[i])
        end

        # Parse the optimized recourse variable values
        if n_z != 0
            recourse = getvalue(z)
            for i = 1:n_z
                setvalue(Variable(m, flex_data.recourse_cols[i]), recourse[i])
            end
        end

        # Parse the optimized recourse variable values
        if n_x != 0
            state = getvalue(x)
            for i = 1:n_x
                setvalue(Variable(m, state_cols[i]), state[i])
            end
        end

        # Parse the active active constraints
        active_inds = findall(abs.(getvalue(y) - ones(n_f)) .<= 1e-4)
        flex_data.active_constraints = inequal_inds[active_inds]

        # Save the flexibility index
        if sqrt_δ
            flex_data.flexibility_index = sqrt(getobjectivevalue(m_solve))
        elseif flex_data.uncertainty_set.name == :Ellipsoid && conic_δ
            flex_data.flexibility_index = getobjectivevalue(m_solve)^2
        else
            flex_data.flexibility_index = getobjectivevalue(m_solve)
        end

        # Save the solution time
        try
            flex_data.solution_time = getsolvetime(m_solve)
        catch
            flex_data.solution_time = tic_time
        end
    end

    return status
end

"""
    MakeInputDict(m::Model)
Returns `input_dict::Dict` which contains the state space representation of the system equations stored in
the flexibility model. The resulting `input_dict` is used by most of the flexibility analysis functions.

**Arguments**
- `m::Model` The flexibility model.
"""
function MakeInputDict(m::Model)
    # Pull in the constraint data from the JuMP model
    flex_data = getflexibilitydata(m)
    constr_coeffs = JuMP.prepConstrMatrix(m)  #A coefficients
    constr_bounds = JuMP.prepConstrBounds(m)  #A lower and A upper

    # Determine which constraints are inequalities and equalities
    equal_inds = findall(constr_bounds[1] .== constr_bounds[2])
    inequal_inds = findall(constr_bounds[1] .!= constr_bounds[2])

    # Parse the problem dimensions
    n_f = length(inequal_inds)
    n_h = length(equal_inds)
    n_θ = flex_data.numRVs
    n_z = flex_data.num_recourse_vars
    n_x = size(constr_coeffs)[2] - n_θ - n_z

    # Check inequality directions and adjust as necessary
    reversed_inds = findall(constr_bounds[2] .== Inf)
    constr_coeffs[reversed_inds, :] *= -1
    constr_bounds[2][reversed_inds] = -constr_bounds[1][reversed_inds]
    constr_bounds[1][reversed_inds] .= -Inf

    # Parse the coefficent matrices needed for the flexibility index problem
    fConsts = -constr_bounds[2][inequal_inds]
    hConsts = -constr_bounds[2][equal_inds]
    fRandoms = constr_coeffs[inequal_inds, flex_data.RVcols]
    hRandoms = constr_coeffs[equal_inds, flex_data.RVcols]
    if n_z != 0
        fControls = constr_coeffs[inequal_inds, flex_data.recourse_cols]
        hControls = constr_coeffs[equal_inds, flex_data.recourse_cols]
    else
        fControls = 0
        hControls = 0
    end
    if n_x != 0
        num_vars = n_x + n_z + n_θ
        var_cols = collect(1:num_vars)
        state_cols = var_cols[BitArray(!any((y-> ==(y, var_cols[i])), flex_data.RVcols) && !any((y-> ==(y, var_cols[i])), flex_data.recourse_cols) for i = 1:num_vars)]
        fStates = constr_coeffs[inequal_inds, state_cols]
        hStates = constr_coeffs[equal_inds, state_cols]
    else
        fStates = 0
        hStates = 0
    end

    # Get the mean
    θ_nom = flex_data.RVmeans

    # Make the input dictionary
    input_dict = Dict("fConsts" => fConsts, "fControls" => fControls, "fRandoms" => fRandoms, "fStates" => fStates,
                      "hConsts" => hConsts, "hControls" => hControls, "hRandoms" => hRandoms, "hStates" => hStates,
                      "theta_nom" => θ_nom, "dims" => [n_z; n_θ; n_x; n_f; n_h], "inequal_inds" => inequal_inds)
    if n_x != 0
        input_dict["state_cols"] = state_cols
    end
    return input_dict
end

"""
    AddSystemExpressions(m::Model, input_dict::Dict, [num_scenarios::Int = 0])
Returns a vector of vectors where the first contains all the inequality expressions corresponding
to the inequalities defined in `input_dict` and the second contains all of the equality expressions
corresponding to the equalities defined in `input_dict`.

**Arguments**
- `m::Model` The flexibility model.
- `input_dict::Dict` Input dictionary as defined by [`MakeInputDict`](@ref).
- `num_scenarios::Int = 0` The number of scenerio subproblems, 0 turns off this feature.
"""
function AddSystemExpressions(m::Model, input_dict::Dict, num_scenarios::Int = 0)
    # Extract information
    fConsts = input_dict["fConsts"]
    fControls = input_dict["fControls"]
    fRandoms = input_dict["fRandoms"]
    fStates = input_dict["fStates"]
    hConsts = input_dict["hConsts"]
    hControls = input_dict["hControls"]
    hRandoms = input_dict["hRandoms"]
    hStates = input_dict["hStates"]
    n_z = input_dict["dims"][1]
    n_θ = input_dict["dims"][2]
    n_x = input_dict["dims"][3]
    n_f = input_dict["dims"][4]
    n_h = input_dict["dims"][5]

    # Make the expressions
    if num_scenarios == 0
        if n_z != 0
            if n_x != 0
                @expression(m, fexpr[j = 1:n_f], fConsts[j] + sum(fControls[j, i] * m[:z][i] for i = 1:n_z) + sum(fRandoms[j, i] * m[:θ][i] for i = 1:n_θ) +
                            sum(fStates[j, i] * m[:x][i] for i = 1:n_x))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h], hConsts[j] + sum(hControls[j, i] * m[:z][i] for i = 1:n_z) + sum(hRandoms[j, i] * m[:θ][i] for i = 1:n_θ) +
                                sum(hStates[j, i] * m[:x][i] for i = 1:n_x))
                end
            else
                @expression(m, fexpr[j = 1:n_f], fConsts[j] + sum(fControls[j, i] * m[:z][i] for i = 1:n_z) + sum(fRandoms[j, i] * m[:θ][i] for i = 1:n_θ))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h], hConsts[j] + sum(hControls[j, i] * m[:z][i] for i = 1:n_z) + sum(hRandoms[j, i] * m[:θ][i] for i = 1:n_θ))
                end
            end
        else
            if n_x != 0
                @expression(m, fexpr[j = 1:n_f], fConsts[j] + sum(fRandoms[j, i] * m[:θ][i] for i = 1:n_θ) + sum(fStates[j, i] * m[:x][i] for i = 1:n_x))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h], hConsts[j] + sum(hRandoms[j, i] * m[:θ][i] for i = 1:n_θ) + sum(hStates[j, i] * m[:x][i] for i = 1:n_x))
                end
            else
                @expression(m, fexpr[j = 1:n_f], fConsts[j] + sum(fRandoms[j, i] * m[:θ][i] for i = 1:n_θ))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h], hConsts[j] + sum(hRandoms[j, i] * m[:θ][i] for i = 1:n_θ))
                end
            end
        end
    else
        if n_z != 0
            if n_x != 0
                @expression(m, fexpr[j = 1:n_f, k = 1:num_scenarios], fConsts[j] + sum(fControls[j, i] * m[:z][i, k] for i = 1:n_z) + sum(fRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ) +
                            sum(fStates[j, i] * m[:x][i, k] for i = 1:n_x))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h, k = 1:num_scenarios], hConsts[j] + sum(hControls[j, i] * m[:z][i, k] for i = 1:n_z) + sum(hRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ) +
                                sum(hStates[j, i] * m[:x][i, k] for i = 1:n_x))
                end
            else
                @expression(m, fexpr[j = 1:n_f, k = 1:num_scenarios], fConsts[j] + sum(fControls[j, i] * m[:z][i, k] for i = 1:n_z) + sum(fRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h, k = 1:num_scenarios], hConsts[j] + sum(hControls[j, i] * m[:z][i, k] for i = 1:n_z) + sum(hRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ))
                end
            end
        else
            if n_x != 0
                @expression(m, fexpr[j = 1:n_f, k = 1:num_scenarios], fConsts[j] + sum(fRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ) + sum(fStates[j, i] * m[:x][i, k] for i = 1:n_x))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h, k = 1:num_scenarios], hConsts[j] + sum(hRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ) + sum(hStates[j, i] * m[:x][i, k] for i = 1:n_x))
                end
            else
                @expression(m, fexpr[j = 1:n_f, k = 1:num_scenarios], fConsts[j] + sum(fRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ))
                if n_h != 0
                    @expression(m, hexpr[j = 1:n_h, k = 1:num_scenarios], hConsts[j] + sum(hRandoms[j, i] * m[:θ][i, k] for i = 1:n_θ))
                end
            end
        end
    end

    # Return expressions
    if n_h != 0
        return [fexpr, hexpr]
    else
        return [fexpr, []]
    end
end
