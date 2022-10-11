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

solid_displacement_parameters = {"newton_solver":{"linear_solver":'bicgstab',"preconditioner":'jacobi', "report":'True', "error_on_nonconvergence":'True',\
                                                  "absolute_tolerance":1e-12, "relative_tolerance":1e-6, "maximum_iterations":50}}

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



