# Test default model and model with other solver
m = FlexibilityModel(solver = PavitoSolver(mip_solver = GLPKSolverMIP(),
                     cont_solver = IpoptSolver(print_level = 0), log_level = 0, mip_solver_drives = false))
m2 = FlexibilityModel()
true
