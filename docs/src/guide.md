# User Guide
```@meta
CurrentModule = FlexibilityAnalysis
```

This page provides an overview of how to use FlexibilityAnalysis to analyze system flexibility.
Detailed explanations on the syntax of each method/function and datatype is
provided in [Library](@ref).

The package needs to be loaded along with [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) in the usual manner:

```julia
using FlexibilityAnalysis, JuMP
```

## Model Definition and Setup
### Flexibility Model Definition
The flexibility model is defined with the [`FlexibilityModel`](@ref) function
and the solver that will be used to solve the flexibility index problem should
be specified.

```julia
using Gurobi
m = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))
```
```julia
Feasibility problem with:
 * 0 linear constraints
 * 0 variables
Solver is Gurobi
```

Flexibility models are JuMP models that have been extended to include information
needed for flexibility analysis and incorporate [`solvehook`](@ref) which solves
the flexibility index problem.

!!! note
    A solver that is capable of solving MIQCPs if an ellipsoidal or 2-norm uncertainty
    set must be used. Otherwise, an MILP solver will work.

### Variable Definition
Now we can add variables to `m` using the [`@randomvariable`](@ref) macro for random
variables, the [`@recoursevariable`](@ref) macro for recourse/control variables,
the standard `@variable` JuMP macro to add state variables.

```julia
means = [620; 388; 583; 313]
@randomvariable(m, T[i = 1:4], mean = means[i])
@recoursevariable(m, Qc)
@variable(m, x)
```

In applications where it is not clear what variables are state/recourse variables,
all the nonrandom variables can simply be defined with the `@recoursevariable` macro.
The `@randomvariable` and `@recoursevariable` macros can define single variables
and/or arrays of variables if we append brackets to the variable as is done in JuMP.
For example

```julia
@recoursevariable(m, y[1:N, 1:M])
```
will create an `N` by `M` array of recourse variables.

However, these macros currently don't support specification of variable upper/lower
bounds in contrast to the traditional JuMP syntax. Similarly, bounds added via
the `@variable` macro will be ignored.

Please note that the `@randomvariable` macro requires that a mean be provided for
each random variable. This can be done in the three following ways:

```julia
@randomvariable(m, z], mean = 42) # create one variable with mean 42
@randomvariable(m, z[i = 1:N], mean = means[i]) # assign means via a predefined vector
@randomvariable(m, z[i = 1:N], mean = 42) # assign the same mean to each variable
```

The vector containing the means of all the random variables can later be changed
with the [`setmean`](@ref) function.

### Constraint Definition
Now we can add constraints to `m` via the `@constraint` macro as we normally would
with typical JuMP models.

```julia
@constraint(m, -100 - 0.67Qc + 2T[2] + x <= 0.0)
@constraint(m, -250 - T[2] == x)
@constraint(m, 0.5Qc - 0.75T[1] - T[2] - T[3] <= -1388.5)
@constraint(m, -Qc + 1.5T[1] + 2T[2] + T[3] >= 2044)
@constraint(m, Qc - 1.5T[1] - 2T[2] - T[3] - 2T[4] <= -2830)
@constraint(m, -Qc + 1.5T[1] + 2T[2] + T[3] + 3T[4] <= 3153)
```
```julia
JuMP.ConstraintRef{JuMP.Model,FlexibilityAnalysis.FlexibilityConstraint}(Feasibility problem with:
 * 6 linear constraints
 * 6 variables
Solver is Gurobi, 6)
```

!!! note
    Currently, only linear constraints can be used with FlexibilityAnalysis.

### Uncertainty Set Definition
Now we can specify which uncertainty set we would like to use. This is done via
the [`setuncertaintyset`](@ref) function. One of the 5 uncertainty set types
described in the [Uncertainty Set Characterization](@ref) section can specified.
For a flexibility model `m`, the set type is specified in the second argument with
one of the following symbols: `:Ellipsoid`, `:Hyperbox`, or `:PNorm`. The third
argument should contain whatever attribute is needed for that uncertainty set type.
Continuing the example above we would define an ellipsoidal set

```julia
covar = [11.11 0 0 0; 0 11.11 0 0; 0 0 11.11 0; 0 0 0 11.11]
setuncertaintyset(m, :Ellipsoid, covar)
```

where the required attribute is the covariance matrix. Note that the covariance
matrix must be symmetric positive semi-definite dimensionality that matches the
number of random variables.

The other sets could have instead been defined by:

```julia
box_dev = [10; 10; 10; 10]
setuncertaintyset(m, :Hyperbox, [[box_dev]; [box_dev]])
setuncertaintyset(m, :PNorm, 1)
setuncertaintyset(m, :PNorm, 2)
setuncertaintyset(m, :PNorm, Inf)
```

The hyperbox set requires that a vector of two vectors that correspond to the
negative and positive deviations be provided (these deviations are explained in
[Uncertainty Set Characterization](@ref)). The p-norm sets required that the value
of `p` be provided which can `1`, `2`, or `Inf`. A summary of the sets and their
required inputs is shown below.

| Uncertainty Set Type | Symbol       | Attribute                                |
|----------------------|--------------|------------------------------------------|
| Ellipsoidal          | `:Ellipsoid` | `covariance::Matrix`                     |
| Hyperbox             | `:Hyperbox`  | `[[neg_dev]; [pos_dev]]::Vector{Vector}` |
| P-Norm               | `:PNorm`     | `1`, `2`, or `Inf`                       |

The `setuncertaintyset` function also accepts the keyword argument `only_positive::Bool`
to indicate if the uncertainty set should be intersected with the set of all positive
real numbers ``\mathbb{R}_+^{n_{\theta}}``. For example, to define the set
``T_{ellip}(\delta) \cap \mathbb{R}_+^{n_{\theta}}`` we would call

```julia
setuncertaintyset(m, :Ellipsoid, covar, only_positive = true)
```

!!! note
    By default the uncertainty set is taken to be ellipsoidal. Thus, one
    need not call `setuncertaintyset` if an ellipsoidal set is desired, but the
    the covariance matrix must still be specified via the [`setcovariance`](@ref)
    function.

## Pre-Solution Methods
This section will outline methods/functions that are geared to be called before
`m` is solved, but can still be applied after it is solved.

### Mean Extraction/Manipulation
The mean corresponding to a particular random variable can be extracted via the
[`getmean(variable::RandomVariable)`](@ref) method. Note that this only works for
individual variables and that arrays of random variables are not valid input.

```julia
variable_mean = getmean(T[2])
```
```julia
388
```

The means of all the random variables in `m` can be extracted via the
[`getmean(m::Model)`](@ref) method. For the current example we have

```julia
current_mean = getmean(m)
```
```julia
4-element Array{Number,1}:
 620
 388
 583
 313
```

The means of all the random variables in `m` can be redefined using the
[`setmean`](@ref) method where the new means are passed as a vector in the second
argument.

```julia
setmean(m, [1; 1; 1; 1])
```
```julia
4-element Array{Int64,1}:
 1
 1
 1
 1
```

Note that the new means vector must match the length of the current means vector,
otherwise calling `setmean` will throw an error. We will reset means back to their
original values before continuing.

```julia
setmean(m, means)
```
```julia
4-element Array{Int64,1}:
 620
 388
 583
 313
```

As discussed in the [Uncertainty Set Characterization](@ref)) section, it is critical
that the means correspond to a feasible instance of the random variables. The
feasibility of the means can be tested with the [`ismeanfeasible`](@ref) function.
This tests the feasibility of the current mean for `m` using the feasibility
function and returns `true` if it is feasible or `false` otherwise. In our current
example we have

```julia
result = ismeanfeasible(m)
```
```julia
true
```

Thus, our current mean is feasible. By default the `ClpSolver` is used, but the
keyword `solver` can be used to assign another solver LP or NLP solver.

Finally, the [`findcenteredmean`](@ref) function can be used to compute the
analytic center or feasible center (which are described in
[Uncertainty Set Characterization](@ref)). This function returns the feasible
center by default, but can be changed via the `center::Symbol` keyword parameter
where `:feasible` denotes the feasible center and `:analytic` refers to the analytic
center. In our current example we have

```julia
centered_mean = findcenteredmean(m, center = :analytic)
```
```julia
4-element Array{Float64,1}:
  898.125
 -507.214
  594.544
  317.23
```

This center can be used to replace the mean by setting the `update_mean::Bool`
keyword parameter to `true`. Optionally, center can be constrained to be strictly
positive with the `only_positive:Bool` keyword parameter. The default solver is
the `IpoptSolver` since an NLP solver is required for the analytic center, but an
LP solver can be used to compute the feasible center.

### Covariance Extraction/Manipulation
The covariance matrix stored in `m` can be extracted via the [`getcovariance`](@ref)
method.

```julia
covar = getcovariance(m)
```
```julia
4×4 Array{Number,2}:
 11.11   0.0    0.0    0.0
  0.0   11.11   0.0    0.0
  0.0    0.0   11.11   0.0
  0.0    0.0    0.0   11.11
```

The covariance matrix can be set or changed using the [`setcovariance`](@ref)
function which requires the second argument to be `covariance::Matrix`.

```julia
setcovariance(m, covar)
```
```julia
4×4 Array{Number,2}:
 11.11   0.0    0.0    0.0
  0.0   11.11   0.0    0.0
  0.0    0.0   11.11   0.0
  0.0    0.0    0.0   11.11
```

Note that the specified covariance must be symmetric positive semi-definite, otherwise
`setcovariance` will throw an error. Also, it is important that the covariance matrix
appropriately match the number of random variables in `m`.

## Model Solution
Now that the flexibility model `m` is defined we can solve it (i.e., solve the
flexibility index problem). This is done simply by calling the JuMP `solve` function
associated with the model. Thus, for our current example we have

```julia
solve(m, active_constr = true)
```
```julia
:Optimal
```

A number of keyword arguments can be passed which are each described in the
documentation for [`solvehook`](@ref). Here the `active_constr::Bool` keyword
argument is used to turn on the active constraint which enforces the number of
active constraints at the solution (this is can be used for systems with linearly
independent inequalities). The `U::Number` keyword argument specifies the slack
upper bound which can be changed to improve the solution time for a particular
problem. Also, the `conic_δ::Bool` keyword argument can be used to when an MICP
solver is used such as [Pajarito.jl](https://github.com/JuliaOpt/Pajarito.jl).

## Post-Solution Methods

### Critical Point Extraction
Now that `m` is solved, the optimized values of the variables can be extracted
with the `getvalue` method as is normally done with JuMP models. With the current
example we have

```julia
temperatures = getvalue(T)
cooling = getvalue(Qc)
state = getvalue(x)
```

Note that this can be done with single variables and/or arrays of variables.

### Flexibility Index Information
The optimized value of the flexibility index stored in `m` can now be retrieved
by using the [`getflexibilityindex`](@ref) method.

```julia
flexibility_index = getflexibilityindex(m)
```
```julia
3.600355086286672
```

Similarly, the stored flexibility index can be used to obtain the confidence level
if an ellipsoidal uncertainty set was used. This can be calculated using the
[`getconfidencelevel`](@ref) function.

```julia
conf_lvl = getconfidencelevel(m)
```
```julia
0.5372159367269034
```

As discussed in the [Uncertainty Set Characterization](@ref) section, the
confidence level provides a lower bound on the stochastic flexibility index.

The indexes of the active constraints can be obtained via the
[`getactiveconstraints`](@ref) method.

```julia
actives = getactiveconstraints(m)
```
```julia
2-element Array{Int64,1}:
 3
 6
```

If desired, we can also directly extract all of the flexibility data associated
with the [`FlexibilityData`](@ref) type.

```julia
data = getflexibilitydata(m)
```
```julia
FlexibilityAnalysis.FlexibilityData(JuMP.AbstractConstraint[-0.67*Qc + 2*T[2] + x - 100 <= 0, -1*T[2] + -x - 250 == 0, 0.5*Qc + -0.75*T[1] + -1*T[2] + -1*T[3] + 1388.5 <= 0, -1*Qc + 1.5*T[1] + 2*T[2] + 1*T[3] + -2044 >= 0, 1*Qc + -1.5*T[1] + -2*T[2] + -1*T[3] + -2*T[4] + 2830 <= 0, -1*Qc + 1.5*T[1] + 2*T[2] + 1*T[3] + 3*T[4] + -3153 <= 0], 4, Number[620, 388, 583, 313], AbstractString["T[1]", "T[2]", "T[3]", "T[4]"], [1, 2, 3, 4], 1, AbstractString["Qc"], [5], FlexibilityAnalysis.EllipsoidalSet(:Ellipsoid, false), Number[11.11 0.0 0.0 0.0; 0.0 11.11 0.0 0.0; 0.0 0.0 11.11 0.0; 0.0 0.0 0.0 11.11], 3.600355086286672, [3, 6])
```

### Solution Time Extraction
The optimal solution time stored in `m` can be extracted using the
[`getsolutiontime`](@ref) method. This is extracted from the if it is supported,
otherwise it is determined using the `@elapsed` macro. With the current example we
have:

```julia
opt_time = getsolutiontime(m)
```
```julia
0.0034827596723466
```

## Analysis Methods

### Ranking Limiting Constraints
In the [Analysis Techniques](@ref) section we discussed how the flexibility index
problem can be used to rank inequality constraints that most limit system flexibility.
This can be done automatically using the [`rankinequalities`] function for a
flexibility model `m`. This will return a vector of type `Vector{Dict}` where each
dictionary contains the flexibility index, active constraint indexes, and optimized
flexibility model corresponding to a particular rank level. With the current example
we obtain

```julia
rank_data = rankinequalities(m, max_ranks = 3, active_constr = true)
```
```julia
2-element Array{Dict,1}:
 Dict{String,Any}(Pair{String,Any}("flexibility_index", 3.60036),Pair{String,Any}("model", Feasibility problem with:
 * 6 linear constraints
 * 6 variables
Solver is Gurobi),Pair{String,Any}("active_constraints", [3, 6]))
 Dict{String,Any}(Pair{String,Any}("flexibility_index", 9.93886e-8),Pair{String,Any}("model", Feasibility problem with:
 * 6 linear constraints
 * 6 variables
Solver is Gurobi),Pair{String,Any}("active_constraints", [1, 5]))
```

The keyword argument `max_ranks::Int = 5` specifies the maximum number of rank
levels, and the `m` will be iteratively solved without the previous active constraints
until the maximum number of ranks is accomplished or the problem becomes unbounded,
whichever occurs first. We also note that all of the same keyword arguments available
to the `solve` function are accessible here since `rankinequalities` is a wrapper
function for `solve`.

### Stochastic Flexibility Index
The stochastic flexibility index can be computed via Monte Carlo sampling using
the [`findstochasticflexibility`](@ref) function. The number of samples can be
specified with the `num_pts::Int = 10000` keyword argument.

```julia
SF = findstochasticflexibility(m, num_pts = 10000, use_vulnerability_model = true)
```
```julia
0.9634
```

By default each sample is evaluated individually, but all of them can be evaluated
simultaneously by setting the `use_vulnerability_model::Bool` to `true` as explained
in [Stochastic Flexibility Index Problem](@ref). Furthermore, the
`use_flexibility_index::Bool = false` keyword argument can be set to `true` to
use the optimized uncertainty set to reduce the number of MC samples that need
to be evaluated. The `only_positive::Bool = false` keyword argument enforces that
only positive MC samples are used.

!!! note
    By default the `ClpSolver` is used, but other solvers that support warm
    starts such as Gurobi and Cplex will improve performance if
    `use_vulnerability_model::Bool = false`. The solver can be changed with the
    `solver` keyword argument.
