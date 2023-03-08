from dolfin import *
from ufl import tensors, nabla_div
from .functions import *
from fenicstools import interpolate_nonmatching_mesh
from .solver_options import u_solver, p_solver, \
							u_solver_c
from .constitutive_eq import *
from .fem_stabilizations import *
import sys

sys.path.insert(0,  '..')
from user_inputs.user_parameters import *
from user_inputs.boundary_initial_conditions import constrained_domain
from utilities.read import *


PI = 3.14159265


class Fluid_problem:	
									
	def __init__(self, fluid_mesh, bool_stream):
		
		Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)
		if not problem_physics['viscous_dissipation']: Ec = 0.0

		mesh = fluid_mesh.mesh
		dim = mesh.geometry().dim()
		
		# Velocity components
		self.u_components = dim
		
		V  = FunctionSpace(mesh, 'P', fem_degree['velocity_degree'], constrained_domain = constrained_domain)		  # Fluid velocity   
		Q  = FunctionSpace(mesh, 'P', fem_degree['pressure_degree'], constrained_domain = constrained_domain)		  # Fluid pressure
		Z1 = VectorFunctionSpace(mesh, 'P', fem_degree['lagrange_degree'])                                            # Lagrange multiplier

		# --------------------------------

		Vp = VectorFunctionSpace(mesh, 'P', fem_degree['velocity_degree'])

		# Function assigner : scaler components to vector
		self.assigner_uv = FunctionAssigner(Vp, [V for ui in range(self.u_components)])

		# --------------------------------

		self.u1  = TrialFunction(V)
		self.v   = TestFunction(V)
		self.p   = TrialFunction(Q)
		self.q   = TestFunction(Q)
		self.Lm1 = [TrialFunction(Z1.sub(ui)) for ui in range(self.u_components)]

		variables = dict(); u_ = []; p_ = []
		for i in range(3):
			u_.insert(i, as_vector([Function(V) for ui in range(self.u_components)]))
			p_.insert(i, Function(Q))

		uv   = Function(Vp)		
		Lm_f = Function(Z1)
		Lm_f.vector().zero()

		variables.update(u_=u_, uv=uv, p_=p_, Lm_f=Lm_f)	
			
		vort, psi = Function(Q), Function(Q)
		vort.vector().zero(); psi.vector().zero()
		variables.update(vort=vort, psi=psi)

		# --------------------------------

		# Body force
		f = Constant((0,)*dim)
		if problem_physics['body_force'] == True: f = Constant(((1/(Fr*Fr))*f_dir(dim)))	
		
		# --------------------------------

		self.nx = tensors.unit_vector(0, dim)
		self.ny = tensors.unit_vector(1, dim)
		if dim == 3: self.nz = tensors.unit_vector(2, dim)

		self.A1 = None
		self.A2 = None
		self.A3 = None
		self.null_space = VectorSpaceBasis([])
		self.residual = [Function(V) for ui in range(self.u_components)]
		self.matrix = dict(Mij=None, Kij = None, Cij = None, \
						   Sij = [None for ui in range(self.u_components)], \
						   Bij = [None for ui in range(self.u_components)], \
						   Pij = [None for ui in range(self.u_components)], \
						   Yij = [None for ui in range(self.u_components)], \
						   b1_Ls1 = None, A1_SCW1 = None, b2 = None)
		
		self.f = f
		self.dim = dim
		self.variables = variables
		self.bool_stream = bool_stream
		self.h_f = CellDiameter(mesh)
		self.h_f_X = project(self.h_f, FunctionSpace(mesh, 'P', 1))
		
		self.F = [V, Q, Z1]
		self.dx = Measure("dx", domain=mesh)
		self.ds = Measure("ds", domain=mesh, subdomain_data=fluid_mesh.get_mesh_boundaries())	
		self.n = FacetNormal(mesh)
		if dim == 2: self.tang = as_vector([self.n[1], -self.n[0]])

		self.Re = Constant(Re)
		self.Fr = Constant(Fr)

		# Convection matrix
		self.u_ab = as_vector([Function(V) for ui in range(self.u_components)])    
		self.Cij = dot(dot(self.u_ab, nabla_grad(self.u1)), self.v)*self.dx

		# --------------------------------

	def pre_assemble(self, px, bcs, dt):

		d = self.matrix; Re = self.Re; dim = self.dim
		u1 = self.u1; v = self.v; p = self.p; q = self.q
		dx = self.dx; ds = self.ds; n = self.n; f = self.f
		
		d['Mij'] = self.A3 = assemble(dot(u1, v)*dx, tensor=d['Mij'])	                                     # Mass matrix 
		d['Kij'] = assemble(dot((0.5/Re)*nabla_grad(u1), nabla_grad(v))*dx, tensor=d['Kij'])      	   				# Viscous matrix
		
		for ui in range(self.u_components):
			d['Sij'][ui] = assemble(dot(p, v.dx(ui))*dx - dot(p*n[ui], v)*ds, tensor=d['Sij'][ui])      	# Pressure matrix
			d['Bij'][ui] = assemble(dot(f[ui], v)*dx, tensor=d['Bij'][ui])                                  # Body-force vector
			d['Pij'][ui] = assemble(dot(p.dx(ui), v)*dx, tensor=d['Pij'][ui])                               # Pressure-gradient matrix

		if problem_physics['solve_FSI'] == True:
			for ui in range(self.u_components):
				d['Yij'][ui] = assemble(dot(self.Lm1[ui], v)*dx, tensor=d['Yij'][ui])                       # Lagrange-multiplier matrix

		if time_control['adjustable_timestep'] == False:
			self.A1 = self.matrix['Kij'].copy()
			self.A1.axpy(1.0/float(dt), self.matrix['Mij'], True)
	    	
		self.A2 = assemble(dot(nabla_grad(p), nabla_grad(q))*dx)
		if bcs['pressure'] == []:
		    self.null_space = attach_nullspace(self.A2, px, self.F[1])    

		# Boundary conditions
		[bc.apply(self.A2) for bc in bcs['pressure']]	



	# Predict tentative velocity 	
	def residual_tentative_velocity(self, u_1, u_2, u_ab, p_, Lm_f, f, dt):
	
		for ui in range(self.u_components):			
			U = 0.5*(u_1[ui] + u_2[ui])	
			self.residual[ui] = (u_1[ui] - u_2[ui])/dt + dot(u_ab, nabla_grad(U)) + p_.dx(ui) - nabla_div((2/self.Re)*nabla_grad(U)) - f[ui] - Lm_f[ui]


	def assemble_tentative_velocity(self, u_, p_, Lm_f, dt):

		d = self.matrix; Re = self.Re; f = self.f
		u1 = self.u1; v = self.v; dx = self.dx
		u_ab = self.u_ab; residual = self.residual
		h_f = self.h_f; dx = self.dx; ds = self.ds

		if time_control['adjustable_timestep'] == False:
			A1 = self.A1.copy()
		else: 	
			A1 = Fluid_problem.optimized_lhs(self, dt)

		# Advecting velocity 
		for ui in range(self.u_components):
			u_ab[ui].vector().zero()
			u_ab[ui].vector().axpy(1.5, u_[1][ui].vector())
			u_ab[ui].vector().axpy(-0.5, u_[2][ui].vector())

		# Convective terms
		d['Cij'] = assemble(self.Cij, tensor=d['Cij'])
		A1.axpy(-0.5, d['Cij'], True)

		X1 = A1.copy(); X1.axpy(-2.0, self.matrix['Kij'], True); b1 = [None]*self.u_components
		for ui in range(self.u_components):
			b1[ui] = Fluid_problem.optimized_rhs(self, ui, X1, u_[1], p_[1])	
		
		A1.axpy(1.0, d['Cij'], True)

		# FSI lagrange multiplier
		if problem_physics['solve_FSI'] == True:
			for ui in range(self.u_components):
				b1[ui].axpy(1.0, d['Yij'][ui]*Lm_f.sub(ui).vector())
	
		# Residual vector
		Fluid_problem.residual_tentative_velocity(self, u_[1], u_[2], u_[2], p_[1], Lm_f, f, dt)
		
		# Stabilization terms	
		if stabilization_parameters['SUPG_NS'] == True:
			tau_supg = tau(alpha, u_[1], h_f, Re, dt); operator_supg = Pop(u_[1], v)
			for ui in range(self.u_components):
				d['b1_Ls1'] = assemble(tau_supg*dot(operator_supg, residual[ui])*dx, tensor=d['b1_Ls1'])
				b1[ui].axpy(-1.0, d['b1_Ls1'])

		if stabilization_parameters['crosswind_NS'] == True:
			R = as_vector([self.residual[ui] for ui in range(self.u_components)])
			d['A1_SCW1'] = assemble(inner(tau_cw(C_cw, u_[1], h_f, Re, R)*Pop_CW(u_[1], u1), nabla_grad(v))*dx, tensor=d['A1_SCW1'])
			A1.axpy(1.0, d['A1_SCW1'], True)

		return A1, b1	        

	def optimized_lhs(self, dt):      
	    
	    A = self.matrix['Kij'].copy()
	    A.axpy(1.0/float(dt), self.matrix['Mij'], True)

	    return A

	def optimized_rhs(self, ui, X1, u, p):

	    b = self.matrix['Bij'][ui].copy()
	    b.axpy(1.0, X1*u[ui].vector())
	    b.axpy(1.0, self.matrix['Sij'][ui]*p.vector())

	    return b
	
	def change_initial_guess(self, u):

		for i in range(self.dim):
			u[i].vector().zero()
		
	def solve_tentative_velocity(self, A, x, b, bcs):
	    
		for ui in range(self.u_components):
			[bc.apply(A, b[ui]) for bc in bcs[ui]]
			u_solver.solve(A, x[ui].vector(), b[ui])
			# solve(A, x[ui].vector(), b[ui], 'mumps')




	# Pressure correction
	def assemble_pressure_correction(self, u_, p_, Lm_f, dt):

		p = self.p; q = self.q; dx = self.dx 
		h_f = self.h_f; Re = self.Re; b2 = self.matrix['b2']
	
		L2 = dot(nabla_grad(p_[1]), nabla_grad(q))*dx - (1/dt)*divergence(u_[0], self.u_components)*q*dx
		
		if stabilization_parameters['PSPG_NS'] == True:
			R = as_vector([self.residual[ui] for ui in range(self.u_components)])
			L2 -= tau(alpha, u_[1], h_f, Re, dt)*dot(R, nabla_grad(q))*dx	
		
		b2 = assemble(L2, tensor=b2)
		return b2

	def solve_pressure_correction(self, x, b, bcs):
	    
		A = self.A2
		[bc.apply(b) for bc in bcs]
		if bcs == []:
		    self.null_space.orthogonalize(b)
		p_solver.solve(A, x.vector(), b)
		if bcs == []:
			normalize(x.vector())



	# Velocity correction	
	def assemble_velocity_correction(self, u_, p_, dt):
  
		b3 = [None]*self.u_components

		for ui in range(self.u_components):
			b3[ui] = self.matrix['Mij']*u_[0][ui].vector()
			b3[ui].axpy(-float(dt), self.matrix['Pij'][ui]*(p_[0].vector() - p_[1].vector())) 

		return b3

	def solve_velocity_correction(self, x, b, bcs):
		
		A = self.A3
		for ui in range(self.u_components):
			[bc.apply(A, b[ui]) for bc in bcs[ui]]
			u_solver_c.solve(A, x[ui].vector(), b[ui])


	
	# Post-processing functions
	def post_process_data(self, Mpi, u, p, t, tsp, text_file_handles):

		Re = self.Re; n = self.n; ds = self.ds

		# Compute drag, lift (note to self : written as per 3D sphere)
		traction = -1*dot(sigma(Re, u, p), n)
		drag = assemble(dot(traction, self.nx)*ds(4))/(0.5*(PI)/4)
		lift = assemble(dot(traction, self.ny)*ds(4))/(0.5*(PI)/4)

		Mpi.set_barrier()
		if Mpi.get_rank() == 0:
		    text_file_handles[0].write("{} {} {} {} {} {}" \
		        .format(t, "  ", drag, "  ", lift, "\n"))    

	def calc_vorticity_streamfunction(self, u, bcs):
		
		p = self.p; q = self.q; dx = self.dx
		vort = self.variables['vort']; psi = self.variables['psi']

		if self.bool_stream == True:
            
			# Compute vorticity by L2 projection
			a = p*q*dx
			L = (u.sub(0).dx(1) - u.sub(1).dx(0))*q*dx
			solve(a == L, vort, solver_parameters={'linear_solver': 'gmres',
                         'preconditioner': 'hypre_amg'})

			# Compute stream function : Laplacian(psi) = -vort
			a = inner(grad(p), grad(q))*dx
			L = vort*q*dx
			solve(a == L, psi, bcs, solver_parameters={'linear_solver': 'gmres',
                         'preconditioner': 'hypre_amg'})

		return vort, psi
