var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#FlexibilityAnalysis.jl-1",
    "page": "Home",
    "title": "FlexibilityAnalysis.jl",
    "category": "section",
    "text": "A package for analyzing and quantifying the flexibility of complex systems.(Image: ellipsoid)"
},

{
    "location": "#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "FlexibilityAnalysis.jl (formerly FlexJuMP.jl) provides a computational framework to analyze and quantify system flexibility. It was originally developed as a JuMP extension to automate the setup and computation of the flexibility index problem with several different types of uncertainty sets. However, it is currently being expanded to compute a number of useful metrics and carry out helpful analyses. Currently, it capabilities include:Computing the flexibility index with ellipsoidal, hyperbox, or p-norm uncertainty sets\nChecking nominal point (mean) feasibility\nComputing the stochastic flexibility index via Monte Carlo sampling methods\nUsing the flexibility index to expedite the solution to the stochastic flexibility index\nComputing the analytic and feasible centers of the feasible region to better nominal points\nComputing the confidence level which provides a lower bound on the stochastic flexibility index\nRanking the most limiting system constraints via iterative solution of the flexibility index problemThese techniques are described in greater detail in Background.note: Note\nCurrently, FlexibilityAnalysis only accepts linear constraints and assumes the random variables to be multivariate Gaussian with a certain mean and covariance matrix specified by the user. Development is underway to allow for quadratic constraints and general nonlinear constraints.  "
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "FlexibilityAnalysis.jl is a registered Julia package and can be installed normally.using Pkg\nPkg.add(\"FlexibilityAnalysis\")The latest version of FlexibilityAnalysis.jl only supports Julia 1.0 and above. A version is still available for use with Julia 0.6 under the julia-0.6 branch in FlexJuMP.jl."
},

{
    "location": "#Quick-Start-1",
    "page": "Home",
    "title": "Quick Start",
    "category": "section",
    "text": "Below is a brief example of the high-level API, more explanation and examples are provided in User Guide and Examples.using FlexibilityAnalysis, JuMP\nusing Gurobi\n\n# Setup the uncertainty set parameters\nmeans = [620; 388; 583; 313]\ncovar = [11.11 0 0 0; 0 11.11 0 0; 0 0 11.11 0; 0 0 0 11.11]\n\n# Setup the model\nm = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))\n\n# Define variables\n@randomvariable(m, T[i = 1:4], mean = means[i])\n@recoursevariable(m, Qc)\n\n# Define the constraints\n@constraint(m, -0.67Qc + T[2] <= 350)\n@constraint(m, 0.5Qc - 0.75T[1] - T[2] - T[3] <= -1388.5)\n@constraint(m, -Qc + 1.5T[1] + 2T[2] + T[3] >= 2044)\n@constraint(m, Qc - 1.5T[1] - 2T[2] - T[3] - 2T[4] <= -2830)\n@constraint(m, -Qc + 1.5T[1] + 2T[2] + T[3] + 3T[4] <= 3153)\n\n# Define the uncertainty set\nsetuncertaintyset(m, :Ellipsoid, covar)\n\n# Solve\nstatus = solve(m)\n\n# Retrieve optimized data\nflexibility_index = getflexibilityindex(m)\nconf_lvl = getconfidencelevel(m)\ncritical_temperatures = getvalue(T)\ncritical_cooling = getvalue(Qc)\nactives_constrs = getactiveconstraints(m)\n\n# Print results\nprint(flexibility_index)3.6003393030116135"
},

{
    "location": "#Outline-1",
    "page": "Home",
    "title": "Outline",
    "category": "section",
    "text": "Pages = [\"home.md\", \"background.md\", \"guide.md\", \"examples.md\", \"api.md\"]"
},

{
    "location": "#Acknowledgements-1",
    "page": "Home",
    "title": "Acknowledgements",
    "category": "section",
    "text": "We acknowledge our support from the Department of Energy under grant DE-SC0014114."
},

{
    "location": "background/#",
    "page": "Background",
    "title": "Background",
    "category": "page",
    "text": "CurrentModule = FlexibilityAnalysis"
},

{
    "location": "background/#Background-1",
    "page": "Background",
    "title": "Background",
    "category": "section",
    "text": "This page provides a background on the theory and methods that form the computational framework that is implemented by FlexibilityAnalysis to quantify and analyze system flexibility. This is not a rigorous or complete review and more information is provided here"
},

{
    "location": "background/#System-Characterization-1",
    "page": "Background",
    "title": "System Characterization",
    "category": "section",
    "text": "We consider general systems that are subjected to random variation and attempt to counteract it via recourse (if there is recourse). We characterize these systems as having random variables boldsymboltheta in mathbbR^n_theta, recourse variables mathbfz in mathbbR^n_z, state variables mathbfx in mathbbR^n_x, equality constraints h_i(mathbfz mathbfx boldsymboltheta) = 0  i in I and inequality constraints f_j(mathbfz mathbfx boldsymboltheta) leq 0  j in J. The feasibility (meaning there exists recourse to satisfy the system constraints) of a particular instance of boldsymboltheta is evaluated via the feasibility function psi(boldsymboltheta) which isbeginequation\n	beginaligned\n		psi(boldsymboltheta) = min_mathbfz mathbfx u in mathbbR  u \n		textst  f_j(mathbfz mathbfx boldsymboltheta) leq u  j in J \n		 h_i(mathbfz mathbfx boldsymboltheta) = 0  i in I\n	endaligned\nendequationwhere boldsymboltheta is feasible if psi(boldsymboltheta) leq 0. This feasible function is what the ismeanfeasible and findstochasticflexibility functions use to evaluate feasibility of particular instances of boldsymboltheta. The set of all feasible instances of boldsymboltheta forms the feasible region Theta = boldsymbolthetapsi(boldsymboltheta) leq 0 which is illustrated below.(Image: feasible_region)"
},

{
    "location": "background/#Stochastic-Flexibility-Index-Problem-1",
    "page": "Background",
    "title": "Stochastic Flexibility Index Problem",
    "category": "section",
    "text": "The stochastic flexibility index SF is a flexibility metric for systems that are subjected to uncertain parameters (e.g., disturbances, physical parameters) that are modeled as random variables, and is defined as the probability of finding recourse to maintain feasible operation. This can be computed by integrating the probability density function pmathbbR^n_thetato mathbbR of the random parameters over the feasible region:beginequation\n	SF = int_boldsymbolthetainTheta p(boldsymboltheta) dboldsymboltheta = mathbbPleft( psi(boldsymboltheta) leq 0 right)\nendequationThis is illustrated below for a multivariate Gaussian distribution.(Image: sf_index)The SF index can be computed rigorously via Monte Carlo (MC) sampling. This is done by assessing the feasibility of each realization. Such an approach converges exponentially with the number of MC samples but typically requires a very large number of samples leading to scalability issues. One possible alternative is to estimate the stochastic vulnerability index SV = 1 - SF which can be determined by incorporating all of the MC samples in one large linear program:beginequation\nbeginaligned\nmin_mathbfz^k mathbfx^k y^k in mathbbR  frac1K sum_k = 1^K y^k \ntextst  f_j(mathbfz^k mathbfx^k boldsymboltheta^k) leq y^k  j in J   k = 1 K \n h_i(mathbfz^k mathbfx^k boldsymboltheta^k) = 0  i in I   k = 1K \n y^k geq 0  k = 1K\nendaligned\nendequationwhere K is the number of samples. The sampled value of SV is then obtained by:beginequation\nSV = frac1K sum_k = 1^K mathbb1_mathbbR_0(y^k*)\nendequationThe function findstochasticflexibility employs both of these techniques estimate the SF index via MC sampling."
},

{
    "location": "background/#Flexibility-Index-Problem-1",
    "page": "Background",
    "title": "Flexibility Index Problem",
    "category": "section",
    "text": "The so-called flexibility index problem provides a more scalable approach to quantifying flexibility. This approach seeks to identify the largest uncertainty set T(delta) (where delta in mathbbR_+ is a parameter that scales T) for which the system remains feasible. In other words, we seek to find the largest uncertainty set under which there exists recourse to recover feasibility. The flexibility index F is defined as:beginequation\n	beginaligned\n		 F =  max_deltainmathbbR_+ delta \n		textst  max_boldsymboltheta in T(delta) psi(boldsymboltheta) leq 0\n	endaligned\nendequationThis approach is illustrated with a traditional hyperbox set below.(Image: sf_index)This approach is deterministic in nature; consequently, it does not have a direct probabilistic interpretation (as the SF index does). By transforming the inner problem in terms of its first order KKT conditions we obtain the following mixed integer formulation which is valid for any compact T(delta):beginequation\n	beginaligned\n		F = min_delta mathbfz mathbfx boldsymboltheta lambda_j s_j y_j mu_j  delta \n		textst  f_j(mathbfz mathbfx boldsymboltheta)+ s_j = 0  j in J \n		 h_i(mathbfz mathbfx boldsymboltheta) = 0  i in I \n		 sum_j in J lambda_j = 1 \n		 sum_j in J lambda_j fracpartial f_j(mathbfz mathbfx boldsymboltheta)partial mathbfz + sum_i in I mu_i fracpartial h_i(mathbfz mathbfx boldsymboltheta)partial mathbfz= 0 \n		 s_j leq U(1 - y_j)  j in J \n		 lambda_j leq y_j  j in J \n		 boldsymboltheta in T(delta)  \n		 lambda_j s_j geq 0     y_j in 0 1  j in J\n	endaligned\nendequationThis formulation can be rather tedious to implement and thus it what originally motivated the development of FlexibilityAnalysis. This MIP is what the solvehook associated with the flexibility model m uses to compute the F index when the function solve is used."
},

{
    "location": "background/#Uncertainty-Set-Characterization-1",
    "page": "Background",
    "title": "Uncertainty Set Characterization",
    "category": "section",
    "text": "Any compact set can be used to compute the flexibility index. The table below shows the 5 types FlexibilityAnalysis currently employs. Here barboldsymboltheta is the nominal point (mean), V_boldsymboltheta is the covariance matrix, and Deltaboldsymboltheta^- Deltaboldsymboltheta^+ are maximum lower and upper deviations, respectively. In FlexibilityAnalysis these sets are specified via setuncertaintyset.Name Uncertainty Set\nEllipsoidal Norm T_ellip(delta) = boldsymboltheta  lVertboldsymboltheta - barboldsymbolthetarVert_V_boldsymboltheta^-1^2 leq delta \nHyperbox Set T_box(delta) = boldsymboltheta  barboldsymboltheta - deltaDeltaboldsymboltheta^- leq boldsymboltheta leq barboldsymboltheta + deltaDeltaboldsymboltheta^+\nell_infty Norm T_infty(delta) = boldsymboltheta  lVertboldsymboltheta - barboldsymbolthetarVert_infty leq delta \nell_1 Norm T_1(delta) = boldsymboltheta  lVertboldsymboltheta - barboldsymbolthetarVert_1 leq delta \nell_2 Norm T_2(delta) = boldsymboltheta  lVertboldsymboltheta - barboldsymbolthetarVert_2 leq delta We observe that T(delta) can be represented by any combination of the uncertainty sets in the table above or other sets (provided that the resulting combined set is compact). This provides the ability to incorporate physical knowledge on the uncertain parameters and/or to create a wider range of uncertainty set shapes. For instance, the positive (truncated) ellipsoidal set T_ellip+(delta) = T_ellip(delta) cap mathbbR_+^n_theta, where mathbbR_+^n_theta = boldsymboltheta  boldsymboltheta geq 0 can be used to represent parameters that are known to be non-negative (e.g., demands, prices). The in FlexibilityAnalysis the only_positive option enables the intersection of the sets mentioned above with the positive set mathbbR_+^n_theta.One advantage of using T_ellip(delta)  is that the associated flexibility index (denoted by F_ellip) can be used to obtain the confidence level:beginequation\n	alpha^* = fracgamma(fracn_theta2 fracF_ellip2)Gamma(fracn_theta2)\nendequationwhere gamma(cdot) and Gamma(cdot) are the incomplete and complete gamma functions. Interestingly, this confidence level provides a lower bound for the stochastic flexibility index SF (i.e., alpha^*leq SF). This thus provides an avenue to obtain a probabilistic interpretation of the deterministic flexibility index while avoiding MC sampling. This confidence level can be obtained with the getconfidencelevel function.For flexibility analysis it is critical that barboldsymboltheta be feasible, thus the function [ismeanfeasible] is used to verify the mean\'s feasibility. In some applications, specifying barboldsymboltheta might be difficult if not enough data is available or if a system is not routinely operated at any given point (e.g., a power grid). Thus, we introduce methods that can be employed to find well-centered nominal points.The analytic center barboldsymboltheta_ac can be defined by:beginequation\n	beginaligned\n		barboldsymboltheta_ac in  undersetboldsymboltheta mathbfz mathbfx smathrmargmax  sum_j in J log left(s_j right) \n		textst  f_j(mathbfz mathbfx boldsymboltheta) +s_jleq 0  j in J \n		 h_i(mathbfz mathbfx boldsymboltheta) = 0  i in I\n	endaligned\nendequationHere, s_jin mathbbR are slack variables. This problems pushes the constraints to the interior of the feasible region and gives the analytic center as the point that maximizes the geometric mean of the constraints (slack variables). The feasible center barboldsymboltheta_fc is the point that maximizes the worst-case slack variable and is given by:beginequation\n	beginaligned\n		barboldsymboltheta_fc in undersetboldsymboltheta mathbfz mathbfx smathrmargmax  s \n		textst  f_j(mathbfz mathbfx boldsymboltheta) + s leq 0  j in J \n		 h_i(mathbfz mathbfx boldsymboltheta) = 0  i in I\n	endaligned\nendequationIn FlexibilityAnalysis these centers can be computed using the findcenteredmean function."
},

{
    "location": "background/#Analysis-Techniques-1",
    "page": "Background",
    "title": "Analysis Techniques",
    "category": "section",
    "text": "The flexibility index is a valuable metric that can be used to compare system designs. The utility and interpretation of the index for a system will depend on the choice of T(delta). For instance, the index can be used to determine the confidence level alpha^* if the set T_ellip(delta) is used when boldsymboltheta  sim mathcalN(barboldsymboltheta V_boldsymboltheta). The index might not have a clear interpretation or utility for some choices of T(delta) but can still be useful to quantify improvements in flexibility from design modifications (retrofits). An example of comparing two systems with ellipsoidal and hyperbox sets is shown below.(Image: comparison)It is apparent that Design B has a larger feasible region Theta and thus has a larger SF index. Both sets reflect this behavior and could be used to make a comparison without having to compute the SF index. A key observation is that the choice of T(delta) significantly affects how it measures system flexibility. For instance, for Design B, the ellipsoidal and hyperbox sets identify different limiting constraints, and thus exhibit distinct behavior. In particular, index F_box is limited by constraint f_2 and would not change if constraints f_1, f_3, and/or f_4 were varied to increase the area of Theta, even though the overall system flexibility would be improved. On the other hand, index F_ellip is limited by constraint f_1. This conservative behavior parallels what is observed with the use of uncertainty sets in robust optimization, where such sets are used to optimize against the worst case scenario such that the solution is robust in the face of uncertainty, and the conservativeness of a solution depends on the chosen shape of the uncertainty set.The flexibility index problem can also be used to rank constraints that most limit system flexibility, since it implicitly identifies constraints that limit flexibility. Specifically, at the solution of the flexibility index problem, the binary variables y_j indicate which constraints are active and therefore limit flexibility. Thus, the flexibility index problem can be resolved by excluding this first set of limiting constraints, so that the next set of limiting constraints can be identified. This step can be repeated to rank system constraints to a desired extent (as long as the solution of the problem is bounded). This is done by noticing that each set of limiting constraints will have an associated flexibility index, which indicates how limiting a particular set of constraints is relative to the other sets. This methodology can also be used to identify and rank system components that limit flexibility in order to guide design improvements or quantify value of specific components. Limiting constraint information can also be used to quantify the impact of failure of system components (e.g., a production facility). This process is illustrated in the figure below.(Image: comparison)This shows how limiting (active) constraints are iteratively turned off to find subsequent limiting constraints. In FlexibilityAnalysis the function rankinequalities automatically carries out this iterative methodology. Also, for a particular solved flexibility model m, the active constraint indexes can be retrieved via the getactiveconstraints function.  "
},

{
    "location": "guide/#",
    "page": "User Guide",
    "title": "User Guide",
    "category": "page",
    "text": ""
},

{
    "location": "guide/#User-Guide-1",
    "page": "User Guide",
    "title": "User Guide",
    "category": "section",
    "text": "CurrentModule = FlexibilityAnalysisThis page provides an overview of how to use FlexibilityAnalysis to analyze system flexibility. Detailed explanations on the syntax of each method/function and datatype is provided in Library.The package needs to be loaded along with JuMP.jl in the usual manner:using FlexibilityAnalysis, JuMP"
},

{
    "location": "guide/#Model-Definition-and-Setup-1",
    "page": "User Guide",
    "title": "Model Definition and Setup",
    "category": "section",
    "text": ""
},

{
    "location": "guide/#Flexibility-Model-Definition-1",
    "page": "User Guide",
    "title": "Flexibility Model Definition",
    "category": "section",
    "text": "The flexibility model is defined with the FlexibilityModel function and the solver that will be used to solve the flexibility index problem should be specified.using Gurobi\nm = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))Feasibility problem with:\n * 0 linear constraints\n * 0 variables\nSolver is GurobiFlexibility models are JuMP models that have been extended to include information needed for flexibility analysis and incorporate solvehook which solves the flexibility index problem.note: Note\nA solver that is capable of solving MIQCPs if an ellipsoidal or 2-norm uncertainty set must be used. Otherwise, an MILP solver will work."
},

{
    "location": "guide/#Variable-Definition-1",
    "page": "User Guide",
    "title": "Variable Definition",
    "category": "section",
    "text": "Now we can add variables to m using the @randomvariable macro for random variables, the @recoursevariable macro for recourse/control variables, the standard @variable JuMP macro to add state variables.means = [620; 388; 583; 313]\n@randomvariable(m, T[i = 1:4], mean = means[i])\n@recoursevariable(m, Qc)\n@variable(m, x)In applications where it is not clear what variables are state/recourse variables, all the nonrandom variables can simply be defined with the @recoursevariable macro. The @randomvariable and @recoursevariable macros can define single variables and/or arrays of variables if we append brackets to the variable as is done in JuMP. For example@recoursevariable(m, y[1:N, 1:M])will create an N by M array of recourse variables.However, these macros currently don\'t support specification of variable upper/lower bounds in contrast to the traditional JuMP syntax. Similarly, bounds added via the @variable macro will be ignored.Please note that the @randomvariable macro requires that a mean be provided for each random variable. This can be done in the three following ways:@randomvariable(m, z], mean = 42) # create one variable with mean 42\n@randomvariable(m, z[i = 1:N], mean = means[i]) # assign means via a predefined vector\n@randomvariable(m, z[i = 1:N], mean = 42) # assign the same mean to each variableThe vector containing the means of all the random variables can later be changed with the setmean function."
},

{
    "location": "guide/#Constraint-Definition-1",
    "page": "User Guide",
    "title": "Constraint Definition",
    "category": "section",
    "text": "Now we can add constraints to m via the @constraint macro as we normally would with typical JuMP models.@constraint(m, -100 - 0.67Qc + 2T[2] + x <= 0.0)\n@constraint(m, -250 - T[2] == x)\n@constraint(m, 0.5Qc - 0.75T[1] - T[2] - T[3] <= -1388.5)\n@constraint(m, -Qc + 1.5T[1] + 2T[2] + T[3] >= 2044)\n@constraint(m, Qc - 1.5T[1] - 2T[2] - T[3] - 2T[4] <= -2830)\n@constraint(m, -Qc + 1.5T[1] + 2T[2] + T[3] + 3T[4] <= 3153)JuMP.ConstraintRef{JuMP.Model,FlexibilityAnalysis.FlexibilityConstraint}(Feasibility problem with:\n * 6 linear constraints\n * 6 variables\nSolver is Gurobi, 6)note: Note\nCurrently, only linear constraints can be used with FlexibilityAnalysis."
},

{
    "location": "guide/#Uncertainty-Set-Definition-1",
    "page": "User Guide",
    "title": "Uncertainty Set Definition",
    "category": "section",
    "text": "Now we can specify which uncertainty set we would like to use. This is done via the setuncertaintyset function. One of the 5 uncertainty set types described in the Uncertainty Set Characterization section can specified. For a flexibility model m, the set type is specified in the second argument with one of the following symbols: :Ellipsoid, :Hyperbox, or :PNorm. The third argument should contain whatever attribute is needed for that uncertainty set type. Continuing the example above we would define an ellipsoidal setcovar = [11.11 0 0 0; 0 11.11 0 0; 0 0 11.11 0; 0 0 0 11.11]\nsetuncertaintyset(m, :Ellipsoid, covar)where the required attribute is the covariance matrix. Note that the covariance matrix must be symmetric positive semi-definite dimensionality that matches the number of random variables.The other sets could have instead been defined by:box_dev = [10; 10; 10; 10]\nsetuncertaintyset(m, :Hyperbox, [[box_dev]; [box_dev]])\nsetuncertaintyset(m, :PNorm, 1)\nsetuncertaintyset(m, :PNorm, 2)\nsetuncertaintyset(m, :PNorm, Inf)The hyperbox set requires that a vector of two vectors that correspond to the negative and positive deviations be provided (these deviations are explained in Uncertainty Set Characterization). The p-norm sets required that the value of p be provided which can 1, 2, or Inf. A summary of the sets and their required inputs is shown below.Uncertainty Set Type Symbol Attribute\nEllipsoidal :Ellipsoid covariance::Matrix\nHyperbox :Hyperbox [[neg_dev]; [pos_dev]]::Vector{Vector}\nP-Norm :PNorm 1, 2, or InfThe setuncertaintyset function also accepts the keyword argument only_positive::Bool to indicate if the uncertainty set should be intersected with the set of all positive real numbers mathbbR_+^n_theta. For example, to define the set T_ellip(delta) cap mathbbR_+^n_theta we would callsetuncertaintyset(m, :Ellipsoid, covar, only_positive = true)note: Note\nBy default the uncertainty set is taken to be ellipsoidal. Thus, one need not call setuncertaintyset if an ellipsoidal set is desired, but the the covariance matrix must still be specified via the setcovariance function."
},

{
    "location": "guide/#Pre-Solution-Methods-1",
    "page": "User Guide",
    "title": "Pre-Solution Methods",
    "category": "section",
    "text": "This section will outline methods/functions that are geared to be called before m is solved, but can still be applied after it is solved."
},

{
    "location": "guide/#Mean-Extraction/Manipulation-1",
    "page": "User Guide",
    "title": "Mean Extraction/Manipulation",
    "category": "section",
    "text": "The mean corresponding to a particular random variable can be extracted via the getmean(variable::RandomVariable) method. Note that this only works for individual variables and that arrays of random variables are not valid input.variable_mean = getmean(T[2])388The means of all the random variables in m can be extracted via the getmean(m::Model) method. For the current example we havecurrent_mean = getmean(m)4-element Array{Number,1}:\n 620\n 388\n 583\n 313The means of all the random variables in m can be redefined using the setmean method where the new means are passed as a vector in the second argument.setmean(m, [1; 1; 1; 1])4-element Array{Int64,1}:\n 1\n 1\n 1\n 1Note that the new means vector must match the length of the current means vector, otherwise calling setmean will throw an error. We will reset means back to their original values before continuing.setmean(m, means)4-element Array{Int64,1}:\n 620\n 388\n 583\n 313As discussed in the Uncertainty Set Characterization) section, it is critical that the means correspond to a feasible instance of the random variables. The feasibility of the means can be tested with the ismeanfeasible function. This tests the feasibility of the current mean for m using the feasibility function and returns true if it is feasible or false otherwise. In our current example we haveresult = ismeanfeasible(m)trueThus, our current mean is feasible. By default the ClpSolver is used, but the keyword solver can be used to assign another solver LP or NLP solver.Finally, the findcenteredmean function can be used to compute the analytic center or feasible center (which are described in Uncertainty Set Characterization). This function returns the feasible center by default, but can be changed via the center::Symbol keyword parameter where :feasible denotes the feasible center and :analytic refers to the analytic center. In our current example we havecentered_mean = findcenteredmean(m, center = :analytic)4-element Array{Float64,1}:\n  898.125\n -507.214\n  594.544\n  317.23This center can be used to replace the mean by setting the update_mean::Bool keyword parameter to true. Optionally, center can be constrained to be strictly positive with the only_positive:Bool keyword parameter. The default solver is the IpoptSolver since an NLP solver is required for the analytic center, but an LP solver can be used to compute the feasible center."
},

{
    "location": "guide/#Covariance-Extraction/Manipulation-1",
    "page": "User Guide",
    "title": "Covariance Extraction/Manipulation",
    "category": "section",
    "text": "The covariance matrix stored in m can be extracted via the getcovariance method.covar = getcovariance(m)4×4 Array{Number,2}:\n 11.11   0.0    0.0    0.0\n  0.0   11.11   0.0    0.0\n  0.0    0.0   11.11   0.0\n  0.0    0.0    0.0   11.11The covariance matrix can be set or changed using the setcovariance function which requires the second argument to be covariance::Matrix.setcovariance(m, covar)4×4 Array{Number,2}:\n 11.11   0.0    0.0    0.0\n  0.0   11.11   0.0    0.0\n  0.0    0.0   11.11   0.0\n  0.0    0.0    0.0   11.11Note that the specified covariance must be symmetric positive semi-definite, otherwise setcovariance will throw an error. Also, it is important that the covariance matrix appropriately match the number of random variables in m."
},

{
    "location": "guide/#Model-Solution-1",
    "page": "User Guide",
    "title": "Model Solution",
    "category": "section",
    "text": "Now that the flexibility model m is defined we can solve it (i.e., solve the flexibility index problem). This is done simply by calling the JuMP solve function associated with the model. Thus, for our current example we havesolve(m, active_constr = true):OptimalA number of keyword arguments can be passed which are each described in the documentation for solvehook. Here the active_constr::Bool keyword argument is used to turn on the active constraint which enforces the number of active constraints at the solution (this is can be used for systems with linearly independent inequalities). The U::Number keyword argument specifies the slack upper bound which can be changed to improve the solution time for a particular problem. Also, the conic_δ::Bool keyword argument can be used to when an MICP solver is used such as Pajarito.jl."
},

{
    "location": "guide/#Post-Solution-Methods-1",
    "page": "User Guide",
    "title": "Post-Solution Methods",
    "category": "section",
    "text": ""
},

{
    "location": "guide/#Critical-Point-Extraction-1",
    "page": "User Guide",
    "title": "Critical Point Extraction",
    "category": "section",
    "text": "Now that m is solved, the optimized values of the variables can be extracted with the getvalue method as is normally done with JuMP models. With the current example we havetemperatures = getvalue(T)\ncooling = getvalue(Qc)\nstate = getvalue(x)Note that this can be done with single variables and/or arrays of variables."
},

{
    "location": "guide/#Flexibility-Index-Information-1",
    "page": "User Guide",
    "title": "Flexibility Index Information",
    "category": "section",
    "text": "The optimized value of the flexibility index stored in m can now be retrieved by using the getflexibilityindex method.flexibility_index = getflexibilityindex(m)3.600355086286672Similarly, the stored flexibility index can be used to obtain the confidence level if an ellipsoidal uncertainty set was used. This can be calculated using the getconfidencelevel function.conf_lvl = getconfidencelevel(m)0.5372159367269034As discussed in the Uncertainty Set Characterization section, the confidence level provides a lower bound on the stochastic flexibility index.The indexes of the active constraints can be obtained via the getactiveconstraints method.actives = getactiveconstraints(m)2-element Array{Int64,1}:\n 3\n 6If desired, we can also directly extract all of the flexibility data associated with the FlexibilityData type.data = getflexibilitydata(m)FlexibilityAnalysis.FlexibilityData(JuMP.AbstractConstraint[-0.67*Qc + 2*T[2] + x - 100 <= 0, -1*T[2] + -x - 250 == 0, 0.5*Qc + -0.75*T[1] + -1*T[2] + -1*T[3] + 1388.5 <= 0, -1*Qc + 1.5*T[1] + 2*T[2] + 1*T[3] + -2044 >= 0, 1*Qc + -1.5*T[1] + -2*T[2] + -1*T[3] + -2*T[4] + 2830 <= 0, -1*Qc + 1.5*T[1] + 2*T[2] + 1*T[3] + 3*T[4] + -3153 <= 0], 4, Number[620, 388, 583, 313], AbstractString[\"T[1]\", \"T[2]\", \"T[3]\", \"T[4]\"], [1, 2, 3, 4], 1, AbstractString[\"Qc\"], [5], FlexibilityAnalysis.EllipsoidalSet(:Ellipsoid, false), Number[11.11 0.0 0.0 0.0; 0.0 11.11 0.0 0.0; 0.0 0.0 11.11 0.0; 0.0 0.0 0.0 11.11], 3.600355086286672, [3, 6])"
},

{
    "location": "guide/#Solution-Time-Extraction-1",
    "page": "User Guide",
    "title": "Solution Time Extraction",
    "category": "section",
    "text": "The optimal solution time stored in m can be extracted using the getsolutiontime method. This is extracted from the if it is supported, otherwise it is determined using the @elapsed macro. With the current example we have:opt_time = getsolutiontime(m)0.0034827596723466"
},

{
    "location": "guide/#Analysis-Methods-1",
    "page": "User Guide",
    "title": "Analysis Methods",
    "category": "section",
    "text": ""
},

{
    "location": "guide/#Ranking-Limiting-Constraints-1",
    "page": "User Guide",
    "title": "Ranking Limiting Constraints",
    "category": "section",
    "text": "In the Analysis Techniques section we discussed how the flexibility index problem can be used to rank inequality constraints that most limit system flexibility. This can be done automatically using the [rankinequalities] function for a flexibility model m. This will return a vector of type Vector{Dict} where each dictionary contains the flexibility index, active constraint indexes, and optimized flexibility model corresponding to a particular rank level. With the current example we obtainrank_data = rankinequalities(m, max_ranks = 3, active_constr = true)2-element Array{Dict,1}:\n Dict{String,Any}(Pair{String,Any}(\"flexibility_index\", 3.60036),Pair{String,Any}(\"model\", Feasibility problem with:\n * 6 linear constraints\n * 6 variables\nSolver is Gurobi),Pair{String,Any}(\"active_constraints\", [3, 6]))\n Dict{String,Any}(Pair{String,Any}(\"flexibility_index\", 9.93886e-8),Pair{String,Any}(\"model\", Feasibility problem with:\n * 6 linear constraints\n * 6 variables\nSolver is Gurobi),Pair{String,Any}(\"active_constraints\", [1, 5]))The keyword argument max_ranks::Int = 5 specifies the maximum number of rank levels, and the m will be iteratively solved without the previous active constraints until the maximum number of ranks is accomplished or the problem becomes unbounded, whichever occurs first. We also note that all of the same keyword arguments available to the solve function are accessible here since rankinequalities is a wrapper function for solve."
},

{
    "location": "guide/#Stochastic-Flexibility-Index-1",
    "page": "User Guide",
    "title": "Stochastic Flexibility Index",
    "category": "section",
    "text": "The stochastic flexibility index can be computed via Monte Carlo sampling using the findstochasticflexibility function. The number of samples can be specified with the num_pts::Int = 10000 keyword argument.SF = findstochasticflexibility(m, num_pts = 10000, use_vulnerability_model = true)0.9634By default each sample is evaluated individually, but all of them can be evaluated simultaneously by setting the use_vulnerability_model::Bool to true as explained in Stochastic Flexibility Index Problem. Furthermore, the use_flexibility_index::Bool = false keyword argument can be set to true to use the optimized uncertainty set to reduce the number of MC samples that need to be evaluated. The only_positive::Bool = false keyword argument enforces that only positive MC samples are used.note: Note\nBy default the ClpSolver is used, but other solvers that support warm starts such as Gurobi and Cplex will improve performance if use_vulnerability_model::Bool = false. The solver can be changed with the solver keyword argument."
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": "CurrentModule = FlexibilityAnalysis"
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "This page provides two application examples that illustrate how FlexibilityAnalysis can be used and highlight some of its capabilities."
},

{
    "location": "examples/#Heat-Exchanger-Network-1",
    "page": "Examples",
    "title": "Heat Exchanger Network",
    "category": "section",
    "text": "This is a common benchmark example used in flexibility analysis literature. The heat exchanger network is characterized according to the diagram below.(Image: hx_network)The system constraints are:beginequation\n beginaligned\n  -350 K - 067 Q_c + T_3 leq 0 \n  13885 K + 05 Q_c - 075 T_1 - T_3 - T_5 leq 0 \n  2044 K + Q_c - 15 T_1 - 2 T_3 - T_5 leq 0 \n  2830 K + Q_c - 15 T_1 - 2 T_3 - T_5 - 2 T_8 leq 0 \n  -3153 K - Q_c + 15 T_1 + 2 T_3 + T_5 + 3 T_8 leq 0\n endaligned\nendequationwhere T_1, T_3, T_5, and T_8 denote the Gaussian parameters and Q_c denotes the recourse variable.  The mean and covariance matrix are given by:beginequation\n barboldsymboltheta =\n beginbmatrix\n   620  388  583  313\n   endbmatrix K\n        \n   V_boldsymboltheta =\n   beginbmatrix\n   1111  0  0  0 \n   0  1111  0  0 \n   0  0  1111  0 \n   0  0  0  1111\n endbmatrix K^2\nendequationThe parameter variance sigma_i^2 is taken to be 1111 K^2 which corresponds to bartheta_i pm 3 sigma_i, where 3 sigma_i is equated to 10 K in accordance to the pm 10 K variations historically reported.Now that we have formalized the system parameters and equations we can setup and solve a flexibility model. Let\'s first setup the model using a hyperbox uncertainty set.using FlexibilityAnalysis, JuMP\nusing Gurobi\n\n# Setup the uncertainty set parameters\nmeans = [620; 388; 583; 313]\ncovar = [11.11 0 0 0; 0 11.11 0 0; 0 0 11.11 0; 0 0 0 11.11]\nbox_dev = ones(4) * 10\n\n# Setup the model\nm = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))\n\n# Define variables\n@randomvariable(m, T[i = 1:4], mean = means[i])\n@recoursevariable(m, Qc)\n\n# Define the constraints\n@constraint(m, -350 - 0.67Qc + T[2] <= 0.0)\n@constraint(m, 1388.5 + 0.5Qc - 0.75T[1] - T[2] - T[3] <= 0.0)\n@constraint(m, 2044 + Qc - 1.5T[1] - 2T[2] - T[3] <= 0.0)\n@constraint(m, 2830 + Qc - 1.5T[1] - 2T[2] - T[3] - 2T[4] <= 0.0)\n@constraint(m, -3153 - Qc + 1.5T[1] + 2T[2] + T[3] + 3T[4] <= 0.0)\n\n# Define the uncertainty set\nsetuncertaintyset(m, :Hyperbox, [[box_dev]; [box_dev]]);Before solving let\'s verify that the means provide a feasible nominal point.isvalid = ismeanfeasible(m)trueSince the mean is feasible we can proceed to solve m using all the default settings and then extract the solution data.# Solve\n  status = solve(m)\n\n  if status == :Optimal\n      # Retrieve optimized data\n      flexibility_index = getflexibilityindex(m)\n      temperatures = getvalue(T)\n      cooling = getvalue(Qc)\n      actives = getactiveconstraints(m)\n\n      # Print the results\n      print(\"Flexibility Index:     \", round(flexibility_index, 4), \"\\n\")\n      print(\"Critical Temperatures: \", round.(temperatures, 2), \"\\n\")\n      print(\"Critical Cooling:      \", round(cooling, 2), \"\\n\")\n      print(\"Active Constraints:    \", actives)\n  endFlexibility Index:     0.5\nCritical Temperatures: [615.0, 383.0, 578.0, 318.0]\nCritical Cooling:      615.0\nActive Constraints:    [2, 5]This result indicates that the largest feasible hyperbox uncertainty set is given by 50% of its nominal scale, meaning only pm 5 K variations are feasible. Let\'s now specify the covariance matrix so we can estimate the stochastic flexibility index.setcovariance(m, covar)\nSF = findstochasticflexibility(m, use_vulnerability_model = true)0.9726Now let\'s change the uncertainty set to be ellipsoidal and resolve. Notice that we don\'t need to provide the covariance as an attribute because we already set it with setcovariance. We also will enable diagonalization of the ellipsoidal constraint by using the keyword argument diag:Bool.# Change the uncertainty set\nsetuncertaintyset(m, :Ellipsoid)\n\n# Solve\n  status = solve(m, diag = true)\n\n  if status == :Optimal\n      # Retrieve optimized data\n      flexibility_index = getflexibilityindex(m)\n      temperatures = getvalue(T)\n      cooling = getvalue(Qc)\n      actives = getactiveconstraints(m)\n      conf_lvl = getconfidencelevel(m)\n\n      # Print the results\n      print(\"Flexibility Index:     \", round(flexibility_index, 4), \"\\n\")\n      print(\"Confidence Level:      \", round(conf_lvl, 4), \"\\n\")\n      print(\"Critical Temperatures: \", round.(temperatures, 2), \"\\n\")\n      print(\"Critical Cooling:      \", round(cooling, 2), \"\\n\")\n      print(\"Active Constraints:    \", actives)\n  endFlexibility Index:     3.6004\nConfidence Level:      0.5372\nCritical Temperatures: [620.0, 388.0, 581.0, 319.0]\nCritical Cooling:      620.0\nActive Constraints:    [2, 5]We note that the confidence level provides a lower bound on the stochastic flexibility index as we expect."
},

{
    "location": "examples/#IEEE-14-Bus-Power-Network-1",
    "page": "Examples",
    "title": "IEEE-14 Bus Power Network",
    "category": "section",
    "text": "We now consider the IEEE 14-node power network. Which we model by performing balances at each node n in mathcalC and enforcing capacity constraints on the arcs a_k k in mathcalA and on the suppliers s_b b in mathcalS. The demands d_m m in mathcalD are assumed to be the uncertain parameters. The flexibility index problem those seek to identify the largest set of simultaneous demand withdrawals that the system can tolerate. The deterministic network model is given by:beginequation\n	sum_k in mathcalA_n^rec a_k - sum_k in mathcalA_n^snd a_k + sum_b in mathcalS_n s_b - sum_m in mathcalD_n d_m = 0    n in mathcalC\nendequationbeginequation\n	-a_k^C leq a_k leq a_k^C    k in mathcalA\nendequationbeginequation\n	0 leq s_b leq s_b^C    b in mathcalS\nendequationwhere mathcalA_n^rec denotes the set of receiving arcs at node n, mathcalA_n^snd denotes the set of sending arcs at n, mathcalS_n denotes the set of suppliers at n, mathcalD_n denotes the set of demands at n, a_k^C are the arc capacities, and s_b^C are the supplier capacities.This test case does not provide arc capacities so we enforce a capacity of 100 for all the arcs. A schematic of this system is provided below.(Image: ieee14)This system is subjected to a total of 10 uncertainty disturbances (the network demands). The demands are assumed to be boldsymboltheta sim mathcalN(barboldsymboltheta V_boldsymboltheta), where barboldsymboltheta = barboldsymboltheta_fc and V_boldsymboltheta = 1200 mathbbI. We select the set T_infty(delta). Thus, we setup flexibility model.using FlexibilityAnalysis, JuMP\nusing Gurobi, Pavito, Ipopt\n\n# Set the dimensions\nn_gens = 5\nn_lines = 20\nn_dems = 11\n\n# Setup the uncertainty set parameters\ncovar = eye(n_dems) * 1200.\n\n# Specify the network details\nline_cap = 100\ngen_cap = [332; 140; 100; 100; 100]\n\n# Setup the model\nm = FlexibilityModel(solver = PavitoSolver(mip_solver = GurobiSolver(OutputFlag = 0),\n                     cont_solver = IpoptSolver(print_level = 0), log_level = 0,\n                     mip_solver_drives = false))\n\n# Define variables\n@randomvariable(m, d[i = 1:n_dems], mean = 0) # Temperarily set the mean to 0\n@recoursevariable(m, a[1:n_lines])\n@recoursevariable(m, g[1:n_gens])\n\n# Set the line capacity constraints\n@constraint(m, [line = 1:n_lines], -line_cap <= a[line])\n@constraint(m, [line = 1:n_lines], a[line] <= line_cap)\n\n# Set the generator capacity constraints\n@constraint(m, [gen = 1:n_gens], 0.0 <= g[gen])\n@constraint(m, [gen = 1:n_gens], g[gen] <= gen_cap[gen])\n\n# Set the node balance constraints\n@constraint(m, g[1] - a[1] - a[6] == 0)\n@constraint(m, a[1] + g[2] - sum(a[i] for i = [2; 4; 5]) - d[1] == 0)\n@constraint(m, g[3] + a[2] - a[3] - d[2] == 0)\n@constraint(m, sum(a[i] for i = [3; 4; 8]) - sum(a[i] for i = [7; 11]) - d[3] == 0)\n@constraint(m, sum(a[i] for i = [5; 6; 7; 12]) - d[4] == 0)\n@constraint(m, g[4] + sum(a[i] for i = [16; 18]) - sum(a[i] for i = [12; 19]) - d[5] == 0)\n@constraint(m, a[9] - sum(a[i] for i = [8; 10]) == 0)\n@constraint(m, g[5] - a[9] == 0)\n@constraint(m, sum(a[i] for i = [10; 11]) - sum(a[i] for i = [13; 14]) - d[6] == 0)\n@constraint(m, sum(a[i] for i = [13; 20]) - d[7] == 0)\n@constraint(m, a[19] - a[20] - d[8] == 0)\n@constraint(m, a[17] - a[18] - d[9] == 0)\n@constraint(m, a[15] - sum(a[i] for i = [16; 17]) - d[10] == 0)\n@constraint(m, a[14] - a[15] - d[11] == 0)\n\n# Define the covariance and the uncertainty set\nsetcovariance(m, covar)\nsetuncertaintyset(m, :PNorm, Inf);The flexibility model m is now defined, so now we compute the feasible center and solve m.# Compute a center to replace the mean if desired\nnew_mean = findcenteredmean(m, center = :feasible, solver = GurobiSolver(OutputFlag = 0),\n                            update_mean = true)\nupdated_mean = getmean(m)\n\n# Solve\nstatus = solve(m)\n\nif status == :Optimal\n    # Retrieve optimized data\n    flexibility_index = getflexibilityindex(m)\n    lines = getvalue(a)\n    generators = getvalue(g)\n    demands = getvalue(d)\n    actives = getactiveconstraints(m)\n\n    # Print the results\n    print(\"Flexibility Index:   \", round(flexibility_index, 4), \"\\n\")\n    print(\"Active Constraints:  \", actives)\nendFlexibility Index:   26.3636\nActive Constraints:  [7, 14, 25, 33, 41, 42, 43, 44, 45]Now we will use the rankinequalities function to obtain a ranking of the most limiting components.# Rank the inequality constraints\nrank_data = rankinequalities(m, max_ranks = 3)3-element Array{Dict,1}:\n Dict{String,Any}(Pair{String,Any}(\"flexibility_index\", 26.3636),Pair{String,Any}(\"model\", Feasibility problem with:\n * 64 linear constraints\n * 36 variables\nSolver is Pavito),Pair{String,Any}(\"active_constraints\", [7, 14, 25, 33, 41, 42, 43, 44, 45]))\n Dict{String,Any}(Pair{String,Any}(\"flexibility_index\", 31.8181),Pair{String,Any}(\"model\", Feasibility problem with:\n * 64 linear constraints\n * 36 variables\nSolver is Pavito),Pair{String,Any}(\"active_constraints\", [21, 26, 47, 48, 49, 50]))\n Dict{String,Any}(Pair{String,Any}(\"flexibility_index\", 83.3332),Pair{String,Any}(\"model\", Feasibility problem with:\n * 64 linear constraints\n * 36 variables\nSolver is Pavito),Pair{String,Any}(\"active_constraints\", [1, 2, 3, 4, 10, 11, 12, 19, 36, 38, 46]))Thus, now we can identify which lines and generators most limit system flexibility. We also observe how the flexibility index increases with each rank level as is expected."
},

{
    "location": "api/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "api/#Library-1",
    "page": "Library",
    "title": "Library",
    "category": "section",
    "text": "CurrentModule = FlexibilityAnalysis"
},

{
    "location": "api/#FlexibilityAnalysis.FlexibilityAnalysis",
    "page": "Library",
    "title": "FlexibilityAnalysis.FlexibilityAnalysis",
    "category": "module",
    "text": "FlexibilityAnalysis\n\nA Julia package that implements a computational framework for quantifying and analyzing system flexibility.\n\nThe basic functionality emulates typical JuMP models to facilitate general computation of the flexibility index problem. Thus, basic functionality uses FlexibilityModel, @randomvariable, @recoursevariable, @variable, @constraint, setuncertaintyset, and solve.\n\nMethods/Macros\n\nFlexibilityModel\n@randomvariable\n@recoursevariable\ngetflexibilitydata\nsetcovariance\ngetcovariance\ngetmean\nsetmean\nsetuncertaintyset\nismeanfeasible\nfindcenteredmean\ngetflexibilityindex\ngetsolutiontime\ngetconfidencelevel\ngetactiveconstraints\nrankinequalities\nfindstochasticflexibility\n\n\n\n\n\n"
},

{
    "location": "api/#Module-1",
    "page": "Library",
    "title": "Module",
    "category": "section",
    "text": "FlexibilityAnalysis"
},

{
    "location": "api/#FlexibilityAnalysis.FlexibilityModel",
    "page": "Library",
    "title": "FlexibilityAnalysis.FlexibilityModel",
    "category": "function",
    "text": "FlexibilityModel(; [solver = JuMP.UnsetSolver()])\n\nReturn a flexibility model object which extends a JuMP model object to contain FlexibilityData and implement a custom solvehook. An appropriate solver should be specified in order solve the flexibility index problem. A solver capable of handling MIQCPs is required for ellipsoidal and 2-norm uncertainty sets otherwise a MILP solver can be used. This model is solved with solve, see solvehook for documention on the accepted keyword arguments.\n\nArguments\n\nsolver = JuMP.UnsetSolver() The solver, should use an MIQCP, MINLP, or MILP solver as appropriate.\n\njulia> m = FlexibilityModel(solver = GurobiSolver())\nFeasibility problem with:\n * 0 linear constraints\n * 0 variables\nSolver is Gurobi\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.@randomvariable",
    "page": "Library",
    "title": "FlexibilityAnalysis.@randomvariable",
    "category": "macro",
    "text": "@randomvariable(m, x, mean)\n\nDefines a random variable using RandomVariable(m::Model, mean::Number, name::AbstractString) and requires that a mean for the variable be provided. This can later be overwritten with setmean.\n\nArguments\n\nm::Model The flexibility model.\nx::Symbol The variable name.\nmean::Number The variable mean\n\njulia> @randomvariable(m2, w, mean = 42)\n\njulia> @randomvariable(m2, ws[i = 1:4], mean = 42)\n4-element Array{FlexibilityAnalysis.RandomVariable,1}:\n FlexibilityAnalysis.RandomVariable(Feasibility problem with:\n * 9 linear constraints\n * 89 variables\n\n julia> @randomvariable(m2, ws[i = 1:4], mean = [1; 2; 3; 4][i])\n 4-element Array{FlexibilityAnalysis.RandomVariable,1}:\n  FlexibilityAnalysis.RandomVariable(Feasibility problem with:\n  * 9 linear constraints\n  * 96 variables\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.@recoursevariable",
    "page": "Library",
    "title": "FlexibilityAnalysis.@recoursevariable",
    "category": "macro",
    "text": "@recoursevariable(m, x)\n\nDefines a recourse variable using RecourseVariable(m::Model, name::AbstractString).\n\nArguments\n\nm::Model The flexibility model.\nx::Symbol The variable name.\n\njulia> @recoursevariable(m2, d)\n\njulia> @recoursevariable(m2, ds[1:4])\n4-element Array{FlexibilityAnalysis.RecourseVariable,1}:\n FlexibilityAnalysis.RecourseVariable(Feasibility problem with:\n * 9 linear constraints\n * 101 variables\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getflexibilitydata",
    "page": "Library",
    "title": "FlexibilityAnalysis.getflexibilitydata",
    "category": "function",
    "text": "getflexibilitydata(m::Model)\n\nReturn the FlexibilityData corresponding the flexibility model m. An error is thrown if m is a regular JuMP model.\n\nArguments\n\nm::Model The flexibility model.\n\njulia> getflexibilitydata(m)\nFlexibilityAnalysis.FlexibilityData(JuMP.AbstractConstraint[], 0, Number[], AbstractString[], Int64[], 0, AbstractString[], Int64[], FlexibilityAnalysis.EllipsoidalSet(:Ellipsoid, false), Array{Number}(undef,0,0), nothing, Int64[], nothing)\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.setcovariance",
    "page": "Library",
    "title": "FlexibilityAnalysis.setcovariance",
    "category": "function",
    "text": "setcovariance(m::Model, covariance::Matrix)\n\nSpecify the covariance matrix covariance to be stored in the flexibility model m. This method verifies that the matrix is symmetric positive semi-definite and writes it to FlexibilityData.covariance.\n\nArguments\n\nm::Model The flexibility model.\ncovariance::Matrix The covariance matrix.\n\njulia> setcovariance(m, [2 1; 1 2])\n2×2 Array{Int64,2}:\n 2  1\n 1  2\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getcovariance",
    "page": "Library",
    "title": "FlexibilityAnalysis.getcovariance",
    "category": "function",
    "text": "getcovariance(m::Model)\n\nReturn the current covariance matrix in the flexibility model as stored in FlexibilityData.covariance.\n\nArguments\n\nm::Model The flexibility model.\n\njulia> getcovariance(m)\n4×4 Array{Number,2}:\n 11.11   0.0    0.0    0.0\n  0.0   11.11   0.0    0.0\n  0.0    0.0   11.11   0.0\n  0.0    0.0    0.0   11.11\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.setmean",
    "page": "Library",
    "title": "FlexibilityAnalysis.setmean",
    "category": "function",
    "text": "setmean(m::Model, mean::Vector)\n\nSpecify the mean corresponding to FlexibilityData.RVmeans stored in the flexibility model m. This method verifies that the length of the input mean matches the length of FlexibilityData.RVmeans before overwriting the current mean.\n\nArguments\n\nm::Model The flexibility model.\nmean::Vector The means of the random variables.\n\nsetmean(m, [2.3; 5])\n2-element Array{Float64,1}:\n 2.3\n 5.0\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getmean-Tuple{Model}",
    "page": "Library",
    "title": "FlexibilityAnalysis.getmean",
    "category": "method",
    "text": "getmean(m::Model)\n\nReturn the current mean vector in the flexibility model as stored in FlexibilityData.RVmeans.\n\nArguments\n\nm::Model The flexibility model.\n\njulia> getmean(m)\n4-element Array{Number,1}:\n 620\n 388\n 583\n 313\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getmean-Tuple{FlexibilityAnalysis.RandomVariable}",
    "page": "Library",
    "title": "FlexibilityAnalysis.getmean",
    "category": "method",
    "text": "getmean(variable::RandomVariable)\n\nReturn the mean corresponding to a particular RandomVariable. Currently this only accepts a single random variable and vector variables are not accepted directly.\n\nArguments\n\nvariable::RandomVariable The random variable name, must be a single variable.\n\njulia> getmean(T[1])\n620\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.setuncertaintyset",
    "page": "Library",
    "title": "FlexibilityAnalysis.setuncertaintyset",
    "category": "function",
    "text": "setuncertaintyset(m::Model, uncertainty_set::Symbol, [attribute = nothing; only_positive::Bool = false])\n\nSpecify the type of uncertainty set to be used and stored in FlexibilityData.uncertainty_set and provide the necessary attribute. An ellipsoidal uncertainty set can be specified with the :Ellipsoid symbol and the corresponding covariance matrix will need to input via attribute if it has not already been set with setcovariance. A hyperbox uncertainty set is specified with the :Hyperbox symbol and corresponding negative/positive deviation vectors need to be inputed via attribute as a vector of vectors of the form [[neg_dev]; [pos_dev]]. Finally, a p-norm uncertainty set can be specified with the :PNorm symbol and providing the corresponding p value via attribute where p can equal 1, 2, or Inf.\n\nArguments\n\nm::Model The flexibility model.\nuncertainty_set::Symbol The uncertainty set name.\nattribute = nothing The necessary atribute for the specified uncertainty set.\n\nKeyword Arguments\n\nonly_positive::Bool = false Indicate if the uncertainty set should be intersected with the set of all positive values.\n\njulia> setuncertaintyset(m, :Ellipsoid)\n\njulia> setuncertaintyset(m, :Ellipsoid, eye(4)* 11.11)\n\njulia> setuncertaintyset(m, :Hyperbox, [[ones(4)]; [ones(4)]])\nFlexibilityAnalysis.HyperboxSet(:Hyperbox, Number[1.0, 1.0, 1.0, 1.0], Number[1.0, 1.0, 1.0, 1.0], false)\n\njulia> setuncertaintyset(m, :PNorm, 1)\nFlexibilityAnalysis.PNormSet(:PNorm, 1, false)\n\njulia> setuncertaintyset(m, :PNorm, 2, only_positive = true)\nFlexibilityAnalysis.PNormSet(:PNorm, 2, true)\n\njulia> setuncertaintyset(m, :PNorm, Inf)\nFlexibilityAnalysis.PNormSet(:PNorm, Inf, false)\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.ismeanfeasible",
    "page": "Library",
    "title": "FlexibilityAnalysis.ismeanfeasible",
    "category": "function",
    "text": "ismeanfeasible(m::Model; [toler::Number = 1e-5, solver = Clp.ClpSolver()])\n\nReturns a Bool indicating if the mean stored in FlexibilityData.RVmeans is feasible, meaning that it lies inside the feasible region. This check is done using the so-called feasibility function.\n\nArguments\n\nm::Model The flexibility model.\n\nKeyword Arguments\n\ntoler::Number = 1e-5 The numerical tolerance for checking the feasibility.\nsolver = Clp.ClpSolver() The solver, any LP or NLP solver shoudl work.\n\njulia> ismeanfeasible(m)\ntrue\n\njulia> ismeanfeasible(m, solver = GurobiSolver(OutputFlag = 0))\nAcademic license - for non-commercial use only\ntrue\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.findcenteredmean",
    "page": "Library",
    "title": "FlexibilityAnalysis.findcenteredmean",
    "category": "function",
    "text": "findcenteredmean(m::Model; [center::Symbol = :feasible, solver = Ipopt.IpoptSolver(print_level = 0), toler::Number = 1e-5, update_mean::Bool = false, only_positive::Bool = false])\n\nReturns a center point based on the analytic or feasible center. The result can overwrite the mean stored in FlexibilityData.RVmeans if desired. This is a wrapper function for ComputeCenter.\n\nArguments\n\nm::Model The flexibility model.\n\nKeyword Arguments\n\ncenter::Symbol = :feasible Indicates the type of center, accepted arguments are :feasible and :analytic.\nsolver = Ipopt.IpoptSolver(print_level = 0) The solver which must be an NLP solver for the analytic center.\ntoler::Number = 1e-5 The tolerance to check solution validity.\nupdate_mean::Bool = false Indicates if the computed center should overwrite FlexibilityData.RVmeans.\nonly_positive::Bool = false Indicates if the center need by strictly positive.\n\njulia> findcenteredmean(m, only_positive = true)\n4-element Array{Float64,1}:\n 1684.74\n   79.0718\n  195.073\n    0.0\n\njulia> findcenteredmean(m, center = :analytic, update_mean = true)\n4-element Array{Float64,1}:\n  898.125\n -507.214\n  594.544\n  317.23\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getflexibilityindex",
    "page": "Library",
    "title": "FlexibilityAnalysis.getflexibilityindex",
    "category": "function",
    "text": "getflexibilityindex(m::Model)\n\nReturn the current flexibility index in the flexibility model as stored in FlexibilityData.flexibility_index.\n\nArguments\n\nm::Model The flexibility model.\n\njulia> getflexibilityindex(m)\n3.5993764186390327\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getconfidencelevel",
    "page": "Library",
    "title": "FlexibilityAnalysis.getconfidencelevel",
    "category": "function",
    "text": "getconfidencelevel(m::Model)\n\nReturn the confidence level provided that the flexibility model m has been solved with an ellipsoidal uncertainty set. This is equivalent to the quantile associated with the optimized uncertainty set. Note that this assumes a multivariate Gaussian distribution with mean FlexibilityData.RVmeans and covariance FlexibilityData.covariance.\n\nArguments\n\nm::Model The flexibility model.\n\njulia> getconfidencelevel(m)\n0.5370703369769008\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getsolutiontime",
    "page": "Library",
    "title": "FlexibilityAnalysis.getsolutiontime",
    "category": "function",
    "text": "getsolutiontime(m::Model)\n\nReturn the solution time to compute the flexibility index as stored in the flexibility model as stored in FlexibilityData.solution_time.\n\nArguments\n\nm::Model The flexibility model.\n\njulia> getsolutiontime(m)\n0.00199127197265625\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.getactiveconstraints",
    "page": "Library",
    "title": "FlexibilityAnalysis.getactiveconstraints",
    "category": "function",
    "text": "getactiveconstraints(m::Model)\n\nReturn the current vector of active constraint indexes in the flexibility model as stored in FlexibilityData.active_constraints.\n\nArguments\n\nm::Model The flexibility model.\n\njulia> getactiveconstraints(m)\n2-element Array{Int64,1}:\n 3\n 6\n\n\n\n\n\n"
},

{
    "location": "api/#JuMP.getvalue-Tuple{FlexibilityAnalysis.FlexibilityVariable}",
    "page": "Library",
    "title": "JuMP.getvalue",
    "category": "method",
    "text": "JuMP.getvalue(v::FlexibilityVariable)\n\nReturn the value of the a flexibility variable this is an extension of JuMP.getvalue.\n\nArguments\n\nv::FlexibilityVariable The flexibility variable.\n\njulia> getvalue(T)\n4-element Array{Float64,1}:\n 620.0\n 388.0\n 581.0\n 319.0\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.rankinequalities",
    "page": "Library",
    "title": "FlexibilityAnalysis.rankinequalities",
    "category": "function",
    "text": "rankinequalities(m::Model; max_ranks::Int = 5, suppress_warnings::Bool = true, U::Int = 10000, diag::Bool = false, active_constr::Bool = false, real_recourse_dim::Int = -1, conic_δ::Bool = false)\n\nReturns ranking data in the form Vector{Dict} where each dictionary corresponds to a particular rank level and contains the optimal flexibility_index, active constraint indexes, and flexibility model. The function will iteratively solve copies of the flexibility model via solve where the prior active constraints are turned off in order rank the most limiting constraints. The user can specify the maximum number of rank levels and the flexibility index problem will be repeatedly solved until that maximum is acheived or the problem becomes unbounded, which occurs first.\n\nArguments\n\nm::Model The flexibility model.\n\nKeyword Arguments\n\nmax_ranks::Int = 5 The maximum number of rank levels. 2\nsuppress_warnings::Bool = false Indicates if solver warnings should be suppressed.\nU::Union{Int, Float64} = 10000 The slack variable upper bound.\ndiag::Bool = false Indicates whether or not to diagnonalize ellipsoidal uncertainty set (this is only active when an ellipsoidal set is used).\nactive_constr::Bool = false Indicates if the optional active constraint should be used which enforces how many inequalities are active at the solution, this must be set to true for systems without control variables and/or contain state variables.\nreal_recourse_dim::Int = -1 The actual number of recourse variables in case state variables are included as recourse variables. This is mandatory if active_constr = true and no state variables are provided.\nconic_δ::Bool = false This should be set to true if a conic solver is used such as Pajarito.jl.\n\njulia> rankinequalities(m, active_constr = true)\n2-element Array{Dict,1}:\n Dict{String,Any}(Pair{String,Any}(\"flexibility_index\", 3.59938),Pair{String,Any}(\"model\", Feasibility problem with:\n * 6 linear constraints\n * 6 variables\nSolver is Pavito),Pair{String,Any}(\"active_constraints\", [3, 6]))\n Dict{String,Any}(Pair{String,Any}(\"flexibility_index\", 9.58983),Pair{String,Any}(\"model\", Feasibility problem with:\n * 6 linear constraints\n * 6 variables\nSolver is Pavito),Pair{String,Any}(\"active_constraints\", [1, 5]))\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.findstochasticflexibility",
    "page": "Library",
    "title": "FlexibilityAnalysis.findstochasticflexibility",
    "category": "function",
    "text": "findstochasticflexibility(m::Model; num_pts::Int = 10000, toler::Number = 1e-5, solver = Clp.ClpSolver(), only_positive::Bool = false, use_vulnerability_model::Bool = false, use_flexibility_index::Bool = false, seed::Int = -1)\n\nReturns the estimated stochastic flexibility index that is evaluated via Monte Carlo sampling. At default this estimation is carried out by evaluating the feasibility of each Monte Carlo sample. The samples are generating from a multivariate Gaussian distribution with mean FlexibilityData.RVmeans and covariance FlexibilityData.covariance. The vulnerability model also tests the feasibility of the samples, but does so in one large optimization problem instead of evaluating each sample indiviually. The optimized flexibility index can also be used to reduce the number of samples that need to be evaluated.\n\nArguments\n\nm::Model The flexibility model.\n\nKeyword Arguments\n\nnum_pts::Int = 10000 Number of Monte Carlo samples.\ntoler::Number = 1e-5 The feasibility check tolerance.\nsolver = Clp.ClpSolver() The solver, any LP or NLP solver should work.\nonly_positive::Bool = false Indicates if only positive samples should be used.\nuse_vulnerability_model::Bool = false Indicates if the vulnerability model should be used.\nuse_flexibility_index::Bool = false Indicates if the optimal flexibility index should be used.\nseed::Int = -1 Random seed for sample collection, any negative value will turn off the random seed.\n\njulia> findstochasticflexibility(m)\n0.9687\n\njulia> findstochasticflexibility(m, use_vulnerability_model = true)\n0.9705\n\njulia> findstochasticflexibility(m, num_pts = 5000, use_flexibility_index = true)\n0.973\n\n\n\n\n\n"
},

{
    "location": "api/#JuMP.linearindex-Tuple{FlexibilityAnalysis.FlexibilityVariable}",
    "page": "Library",
    "title": "JuMP.linearindex",
    "category": "method",
    "text": "JuMP.linearindex(v::FlexibilityVariable)\n\nReturn the index of the a flexibility variable this is an extension of JuMP.linearindex.\n\nArguments\n\nv::FlexibilityVariable The flexibility variable, must be a single variable.\n\njulia> linearindex(Qc)\n1\n\n\n\n\n\n"
},

{
    "location": "api/#Functions/Methods-1",
    "page": "Library",
    "title": "Functions/Methods",
    "category": "section",
    "text": "FlexibilityModel\n@randomvariable\n@recoursevariable\ngetflexibilitydata\nsetcovariance\ngetcovariance\nsetmean\ngetmean(m::Model)\ngetmean(variable::RandomVariable)\nsetuncertaintyset\nismeanfeasible\nfindcenteredmean\ngetflexibilityindex\ngetconfidencelevel\ngetsolutiontime\ngetactiveconstraints\nJuMP.getvalue(v::FlexibilityVariable)\nrankinequalities\nfindstochasticflexibility\nJuMP.linearindex(v::FlexibilityVariable)"
},

{
    "location": "api/#FlexibilityAnalysis.FlexibilityVariable",
    "page": "Library",
    "title": "FlexibilityAnalysis.FlexibilityVariable",
    "category": "type",
    "text": "FlexibilityVariable <: JuMP.AbstractJuMPScalar\n\nAn abstract type to define new variable types.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.RandomVariable",
    "page": "Library",
    "title": "FlexibilityAnalysis.RandomVariable",
    "category": "type",
    "text": "RandomVariable <: FlexibilityVariable\n\nA DataType for random variables.\n\nFields\n\nm::Model Flexibility model.\nidx::Int Index of variable in model.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.RecourseVariable",
    "page": "Library",
    "title": "FlexibilityAnalysis.RecourseVariable",
    "category": "type",
    "text": "RecourseVariable <: FlexibilityVariable\n\nA DataType for recourse variables.\n\nFields\n\nm::Model Flexibility model.\nidx::Int Index of variable in model.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.FlexibilityExpr",
    "page": "Library",
    "title": "FlexibilityAnalysis.FlexibilityExpr",
    "category": "type",
    "text": "FlexibilityExpr <: JuMP.GenericAffExpr\n\nA GenericAffExpr that contains random and/or recourse variables.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.FlexibilityConstraint",
    "page": "Library",
    "title": "FlexibilityAnalysis.FlexibilityConstraint",
    "category": "type",
    "text": "FlexibilityConstraint <: JuMP.AbstractConstraint\n\nA constraint that contains random and/or recourse variables.\n\nFields\n\nflex_expr::FlexibilityExpr Constraint expression.\nsense::Symbol The cosntraint sense symbol :(<=) or :(>=) or :(==).\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.AbstractUncertaintySet",
    "page": "Library",
    "title": "FlexibilityAnalysis.AbstractUncertaintySet",
    "category": "type",
    "text": "AbstractUncertaintySet\n\nAn abstract type to define new uncertainty set types.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.EllipsoidalSet",
    "page": "Library",
    "title": "FlexibilityAnalysis.EllipsoidalSet",
    "category": "type",
    "text": "EllipsoidalSet <: AbstractUncertaintySet\n\nAn ellipsoidal uncertainty set that will use the covariance matrix stored in FlexibilityData.\n\nFields\n\nname::Symbol Name of the set which will be :Ellipsoid.\nonly_positive::Bool An option to indicate if the set should be intersected with the set all positive numbers.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.HyperboxSet",
    "page": "Library",
    "title": "FlexibilityAnalysis.HyperboxSet",
    "category": "type",
    "text": "HyperboxSet <: AbstractUncertaintySet\n\nA hyperbox uncertainty set whose nomimal dimensions are determined by neg_dev and pos_dev.\n\nFields\n\nname::Symbol Name of the set which will be :Hyperbox.\nneg_dev::Vector{Number} A vector of the expected negative deviation of the random variables.\npos_dev::Vector{Number} A vector of the expected positive deviation of the random variables.\nonly_positive::Bool An option to indicate if the set should be intersected with the set all positive numbers.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.PNormSet",
    "page": "Library",
    "title": "FlexibilityAnalysis.PNormSet",
    "category": "type",
    "text": "PNormSet <: AbstractUncertaintySet\n\nA p-norm based uncertainty set based on a bounded p-norm.\n\nFields\n\nname::Symbol Name of the set which will be :PNorm.\np::Number The value of p which can be 1, 2, or Inf.\nonly_positive::Bool An option to indicate if the set should be intersected with the set all positive numbers.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.FlexibilityData",
    "page": "Library",
    "title": "FlexibilityAnalysis.FlexibilityData",
    "category": "type",
    "text": "FlexibilityData\n\nA DataType for storing the data necessary to manage the bookkeeping of the flexibility variables (RandomVariable and RecourseVariable), the uncertainty set, and solution results.\n\nFields\n\nflexibility_constraints::Vector{JuMP.AbstractConstraint} Constraints that involve flexibility variables.\nnumRVs::Int The number of RandomVariable that have been added to the model.\nRVmeans::Vector{Number} The means corresponding to each RandomVariable.\nRVnames::Vector{AbstractString} The symbolic name of each RandomVariable.\nRVcols::Vector{Int} The index of each RandomVariable.\nnum_recourse_vars::Int The number of RecourseVariable that have been added to the model.\nrecourse_names::Vector{AbstractString} The symbolic name of each RecourseVariable.\nrecourse_cols::Vector{Int} The index of each RecourseVariable.\nuncertainty_set::AbstractUncertaintySet The uncertainty set DataType with all of the set specfic attributes.\ncovariance::Matrix{Number} The covariance matrix.\nflexibility_index::Union{Nothing, Number} The flexibility index result obtained from solving the flexibility model.\nactive_constraints::Vector{Int} The indexes of the active inequality constraints at the solution of the flexibility model.\n\'solution_time::Number\' The solution time in seconds.\n\n\n\n\n\n"
},

{
    "location": "api/#DataTypes-1",
    "page": "Library",
    "title": "DataTypes",
    "category": "section",
    "text": "FlexibilityVariable\nRandomVariable\nRecourseVariable\nFlexibilityExpr\nFlexibilityConstraint\nAbstractUncertaintySet\nEllipsoidalSet\nHyperboxSet\nPNormSet\nFlexibilityData"
},

{
    "location": "api/#FlexibilityAnalysis.getuncertaintyset",
    "page": "Library",
    "title": "FlexibilityAnalysis.getuncertaintyset",
    "category": "function",
    "text": "getuncertaintyset(m::Model)\n\nReturn the current uncertainty set datatype in the flexibility modelas stored in FlexibilityData.uncertaintyset.\n\nArguments\n\nm::Model The flexibility model.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.solvehook",
    "page": "Library",
    "title": "FlexibilityAnalysis.solvehook",
    "category": "function",
    "text": "solvehook(m::Model; [suppress_warnings::Bool = false, U::Number = 10000, diag::Bool = false, active_constr::Bool = false, real_recourse_dim::Int = -1, conic_δ::Bool = false, inactives::Vector = []])\n\nReturns the solution status to solving the flexibility model m. This solvehook what JuMP.solve(::Model) for flexibility models. This solves the flexibility index problem using the variables and constraints specified in m.\n\nArguments\n\nm::Model The flexibility model.\n\nKeyword Arguments\n\nsuppress_warnings::Bool = false Indicates if solver warnings should be suppressed.\nU::Number = 10000 The slack variable upper bound.\ndiag::Bool = false Indicates if the ellipsoidal uncertainty set is diagonalized (this is only active when an ellipsoidal set is used).\nactive_constr::Bool = false Indicates if the optional active constraint should be used which enforces how many inequalities are active at the solution, this must be set to true for systems without control variables and/or contain state variables.\nreal_recourse_dim::Int = -1 The actual number of recourse variables in case state variables are included as recourse variables. This is mandatory if active_constr = true and no state variables are provided.\nconic_δ::Bool = false This should be set to true if a conic solver is used such as Pajarito.jl.\ninactives::Vector = [] The indexes of inequality constraints that should be turned off.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.MakeInputDict",
    "page": "Library",
    "title": "FlexibilityAnalysis.MakeInputDict",
    "category": "function",
    "text": "MakeInputDict(m::Model)\n\nReturns input_dict::Dict which contains the state space representation of the system equations stored in the flexibility model. The resulting input_dict is used by most of the flexibility analysis functions.\n\nArguments\n\nm::Model The flexibility model.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.AddSystemExpressions",
    "page": "Library",
    "title": "FlexibilityAnalysis.AddSystemExpressions",
    "category": "function",
    "text": "AddSystemExpressions(m::Model, input_dict::Dict, [num_scenarios::Int = 0])\n\nReturns a vector of vectors where the first contains all the inequality expressions corresponding to the inequalities defined in input_dict and the second contains all of the equality expressions corresponding to the equalities defined in input_dict.\n\nArguments\n\nm::Model The flexibility model.\ninput_dict::Dict Input dictionary as defined by MakeInputDict.\nnum_scenarios::Int = 0 The number of scenerio subproblems, 0 turns off this feature.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.RandomVariable-Tuple{Model,Number,AbstractString}",
    "page": "Library",
    "title": "FlexibilityAnalysis.RandomVariable",
    "category": "method",
    "text": "RandomVariable(m::Model, mean::Number, name::AbstractString)\n\nReturn a RandomVariable DataType for given flexibility model m given the mean and the name. An anonymous JuMP variable is added directly to the flexibility model and its is appended to FlexibilityData.RVcols.\n\nArguments\n\nm::Model The flexibility model.\nmean::Number The variable mean.\nname::AbstractString The variable name.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.RecourseVariable-Tuple{Model,AbstractString}",
    "page": "Library",
    "title": "FlexibilityAnalysis.RecourseVariable",
    "category": "method",
    "text": "RecourseVariable(m::Model, name::AbstractString)\n\nReturn a RecourseVariable DataType for given flexibility model m given the name. An anonymous JuMP variable is added directly to the flexibility model and its is appended to FlexibilityData.recourse_cols.\n\nArguments\n\nm::Model The flexibility model.\nname::AbstractString The variable name.\n\n\n\n\n\n"
},

{
    "location": "api/#Base.show-Tuple{IO,JuMP.GenericAffExpr{JuMP.GenericAffExpr{Float64,Variable},Union{RandomVariable, RecourseVariable}}}",
    "page": "Library",
    "title": "Base.show",
    "category": "method",
    "text": "Base.show(io::IO, a::FlexibilityExpr)\n\nExtend Base.show to print flexibility expressions.\n\n\n\n\n\n"
},

{
    "location": "api/#JuMP.addconstraint-Tuple{Model,FlexibilityAnalysis.FlexibilityConstraint}",
    "page": "Library",
    "title": "JuMP.addconstraint",
    "category": "method",
    "text": "JuMP.addconstraint(m::Model, constr::FlexibilityConstraint)\n\nExtend the JuMP.addconstraint function to handle FlexibilityConstraint types.\n\n\n\n\n\n"
},

{
    "location": "api/#Base.show-Tuple{IO,FlexibilityAnalysis.FlexibilityConstraint}",
    "page": "Library",
    "title": "Base.show",
    "category": "method",
    "text": "JuMP.show(io::IO,c::FlexibilityConstraint)\n\nExtend the JuMP.show function to handle FlexibilityConstraint types.\n\n\n\n\n\n"
},

{
    "location": "api/#JuMP.constructconstraint!-Tuple{JuMP.GenericAffExpr{JuMP.GenericAffExpr{Float64,Variable},Union{RandomVariable, RecourseVariable}},Symbol}",
    "page": "Library",
    "title": "JuMP.constructconstraint!",
    "category": "method",
    "text": "JuMP.constructconstraint!(flex_aff::FlexibilityExpr, sense::Symbol)\n\nExtends JuMP.constructconstraint! for FlexibilityExpr types.\n\n\n\n\n\n"
},

{
    "location": "api/#FlexibilityAnalysis.ComputeCenter",
    "page": "Library",
    "title": "FlexibilityAnalysis.ComputeCenter",
    "category": "function",
    "text": "ComputeCenter(m::Model, center::Symbol, solver, toler::Number, only_positive::Bool)\n\nReturns a center point that can be used to replace the mean if desired.\n\nArguments\n\nm::Model The flexibility model.\ncenter::Symbol Indicates the type of center, accepted arguments are :feasible and :analytic.\nsolver The solver which must be an NLP solver for the analytic center.\ntoler::Number The tolerance to check solution validity.\nonly_positive::Bool Indicates if the center need by strictly positive.\n\n\n\n\n\n"
},

{
    "location": "api/#Internals-1",
    "page": "Library",
    "title": "Internals",
    "category": "section",
    "text": "getuncertaintyset\nsolvehook\nMakeInputDict\nAddSystemExpressions\nRandomVariable(m::Model, mean::Number, name::AbstractString)\nRecourseVariable(m::Model, name::AbstractString)\nBase.show(io::IO, a::FlexibilityExpr)\nJuMP.addconstraint(m::Model, constr::FlexibilityConstraint)\nJuMP.show(io::IO,c::FlexibilityConstraint)\nJuMP.constructconstraint!(flex_aff::FlexibilityExpr, sense::Symbol)\nComputeCenter"
},

{
    "location": "api/#Index-1",
    "page": "Library",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"api.md\"]\nModule = [\"FlexibilityAnalysis\"]\nOrder = [:function, :type]"
},

]}
