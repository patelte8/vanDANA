from dolfin import PETScPreconditioner, PETScKrylovSolver, PETScOptions

# Form compiler parameters
FFC_parameters = {"representation": 'uflacs', "optimize": True, "cpp_optimize": True, "quadrature_degree": 5, "cpp_optimize_flags": "-O3"}

# Solver parameters
krylov_solvers=dict(
    monitor_convergence=False,
    report=False,
    error_on_nonconvergence=True,
    nonzero_initial_guess=True,
    maximum_iterations=300,
    absolute_tolerance=1e-8)

custom_newtons_solver = True
line_search_solver = False

# Solver dictionaries
tentative_velocity_solver=dict(
    solver_type='bicgstab',
    preconditioner_type='jacobi')

pressure_correction_solver=dict(
    solver_type='gmres',
    preconditioner_type='hypre_amg')

velocity_correction_solver=dict(
    solver_type='cg',
    preconditioner_type='jacobi')

energy_conservation_solver=dict(
    solver_type='bicgstab',
    preconditioner_type='jacobi')

solid_momentum_solver=dict(
    solver_type='bicgstab')                 # Use 'mumps' (direct solver), if solid is incompressible

if custom_newtons_solver == True:
    solid_momentum_solver.update(solver_type='bcgs')

# -----------------------------------------------------------------------------------------

piso_iterations = 2                         # no. of PISO iterations

# Define tentative_velocity_solver
precond = PETScPreconditioner(tentative_velocity_solver['preconditioner_type'])
u_solver = PETScKrylovSolver(tentative_velocity_solver['solver_type'], precond)
u_solver.parameters.update(krylov_solvers)

# Define pressure_correction_solver
precond = PETScPreconditioner(pressure_correction_solver['preconditioner_type'])
p_solver = PETScKrylovSolver(pressure_correction_solver['solver_type'], precond)
p_solver.parameters.update(krylov_solvers)
p_solver.set_reuse_preconditioner(True)

# Define velocity_correction_solver
precond = PETScPreconditioner(velocity_correction_solver['preconditioner_type'])
u_solver_c = PETScKrylovSolver(velocity_correction_solver['solver_type'], precond)
u_solver_c.parameters.update(krylov_solvers)
u_solver_c.set_reuse_preconditioner(True)

# Define energy_conservation_solver
precond = PETScPreconditioner(energy_conservation_solver['preconditioner_type'])
t_solver = PETScKrylovSolver(energy_conservation_solver['solver_type'], precond)
t_solver.parameters.update(krylov_solvers)


solid_displacement_parameters = {"newton_solver":{"linear_solver":solid_momentum_solver['solver_type'], "preconditioner":'hypre_amg', "report":True, \
                                                  "error_on_nonconvergence":True, "absolute_tolerance":1e-15, "relative_tolerance":1e-6, "maximum_iterations":20}}

# used only if its a custom newtons solver for compressible solid 
solid_displacement_custom_solver_parameters = {"absolute_tolerance":1e-15, "relative_tolerance":1e-6, "convergence_criterion":'residual', \
                                               "maximum_iterations":20, "report":True, "error_on_nonconvergence":True} #, "relaxation_parameter":1.0}

snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 10,
                                          "report": True,
                                          "line_search": "bt",
                                          "linear_solver": "bicgstab",
                                          "method": "newtonls",
                                          "absolute_tolerance": 1e-9,
                                          "relative_tolerance": 1e-7,
                                          "error_on_nonconvergence": True}}
