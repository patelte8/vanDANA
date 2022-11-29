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
		
		V  = VectorFunctionSpace(mesh, 'P', fem_degree['velocity_degree'], constrained_domain = constrained_domain)		  # Fluid velocity   
		Q  = FunctionSpace(mesh, 'P', fem_degree['pressure_degree'], constrained_domain = constrained_domain)		      	  # Fluid pressure
		Z1 = VectorFunctionSpace(mesh, 'P', fem_degree['lagrange_degree'])                                                  # Lagrange multiplier

		# --------------------------------

		self.u1 = TrialFunction(V)
		self.v = TestFunction(V)
		self.p = TrialFunction(Q)
		self.q = TestFunction(Q)
		self.Lm1 = TrialFunction(Z1) 

		variables = dict(); u_ = []; p_ = []
		for i in range(3):
			u_.insert(i, Function(V))
			p_.insert(i, Function(Q))

		Lm_f = Function(Z1)
		Lm_f.vector()[:] = 0.0
		variables.update(u_=u_, p_=p_, Lm_f=Lm_f)	
			
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
		self.matrix = dict(Mij=None, Kij = None, Sij = None, Bij = None, Pij = None, Yij = None, \
						A1_as1 = None, A1_SCW1 = None, A1_SLS1 = None, A1_Cij = None, BS = None, b1_Ls1 = None, b2 = None)
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
		self.u_ab = Function(V)    
		self.Cij = dot(dot(self.u_ab, nabla_grad(self.u1)), self.v)*self.dx

		# --------------------------------

	def pre_assemble(self, px, bcs, dt):

		d = self.matrix; Re = self.Re;
		u1 = self.u1; v = self.v; p = self.p; q = self.q
		dx = self.dx; ds = self.ds; n = self.n; f = self.f
		
		d['Mij'] = self.A3 = assemble(dot(u1, v)*dx, tensor=d['Mij'])                                       # Mass matrix 
		d['Kij'] = assemble(inner((1/Re)*epsilon(u1), epsilon(v))*dx - \
				   dot((1/Re)*nabla_grad(0.5*u1)*n, v)*ds, tensor=d['Kij'])            	   					# Viscous matrix
		d['Sij'] = assemble(inner(p*Identity(len(u1)), epsilon(v))*dx - dot(p*n, v)*ds, tensor=d['Sij'])    # Pressure matrix
		d['Bij'] = assemble(dot(f, v)*dx, tensor=d['Bij'])                                            		# Body-force vector
		d['Pij'] = assemble(dot(nabla_grad(p), v)*dx, tensor=d['Pij'])                                  	# Pressure-gradient matrix

		if problem_physics['solve_FSI'] == True:
			Lm1 = self.Lm1
			d['Yij'] = assemble(dot(Lm1, v)*dx, tensor=d['Yij'])                                           # Lagrange-multiplier matrix

		if time_control['adjustable_timestep'] == False:
			self.A1 = self.matrix['Kij'].copy()
			self.A1.axpy(1.0/float(dt), self.matrix['Mij'], True)
	    	
		self.A2 = assemble(dot(nabla_grad(p), nabla_grad(q))*dx)
		if bcs['pressure'] == []:
		    self.null_space = attach_nullspace(self.A2, px, self.F[1])    

		# Boundary conditions
		[bc.apply(self.A2) for bc in bcs['pressure']]
		[bc.apply(self.A3) for bc in bcs['velocity']]

		# Stabilization terms
		if stabilization_parameters['stab_LSIC_NS'] == True:
			SLS1 = tau_lsic(Re)*nabla_div(u1)*nabla_div(v)*dx
			d['A1_SLS1'] = assemble(SLS1, tensor=d['A1_SLS1'])

		if stabilization_parameters['stab_backflow_NS'] == True:
			d['BS'] = assemble(Constant(1e4)*dot((Identity(self.dim) - outer(n,n))*u1, v)*ds(2), tensor=d['BS'])	



	# Predict tentative velocity 	
	def residual_tentative_velocity(self, u_0, u_1, p_, Lm_f, dt):
	
		Re = self.Re; f = self.f

		U = 0.5*(u_0 + u_1)
		R = (u_0 - u_1)/dt + dot(u_1, nabla_grad(U)) - nabla_div(sigma(Re, U, p_)) - f - Lm_f
		
		return R	
	
	def assemble_tentative_velocity(self, u_, p_, Lm_f, dt):

		d = self.matrix; Re = self.Re; u_ab = self.u_ab
		u1 = self.u1; v = self.v; dx = self.dx
		h_f = self.h_f; dx = self.dx; ds = self.ds

		if time_control['adjustable_timestep'] == False:
			A1 = self.A1.copy()
		else: 	
			A1 = Fluid_problem.optimized_lhs(self, dt)

		# Advecting velocity 
		u_ab.vector().zero()
		u_ab.vector().axpy(1.5, u_[1].vector())
		u_ab.vector().axpy(-0.5, u_[2].vector())

		# Convective terms
		d['A1_Cij'] = assemble(self.Cij, tensor=d['A1_Cij'])
		A1.axpy(-0.5, d['A1_Cij'], True)

		b1 = Fluid_problem.optimized_rhs(self, A1, u_[1], p_[1])	
		A1.axpy(1.0, d['A1_Cij'], True)

		# FSI lagrange multiplier
		if problem_physics['solve_FSI'] == True:
			b1.axpy(1.0, d['Yij']*Lm_f.vector())
	
		# Stabilization terms 
		if stabilization_parameters['stab_SUPG_NS'] == True: 
			R = Fluid_problem.residual_tentative_velocity(self, u1, u_[1], p_[1], Lm_f, dt)
			S1 = tau(alpha, u_[1], h_f, Re, dt)*dot(R, Pop(u_[1], v))*dx
			as1 = lhs(S1); Ls1 = rhs(S1)	    
			d['A1_as1'] = assemble(as1, tensor=d['A1_as1'])
			d['b1_Ls1'] = assemble(Ls1, tensor=d['b1_Ls1'])
			A1.axpy(1.0, d['A1_as1'], True)
			b1.axpy(1.0, d['b1_Ls1'])

		if stabilization_parameters['stab_cross_NS'] == True:
			R = Fluid_problem.residual_tentative_velocity(self, u_[1], u_[2], p_[2], Lm_f, dt)
			SCW1 = tau_cw(C_cw, u_[1], h_f, Re, R)*inner(Pop_CW(u_[1], u1), nabla_grad(v))*dx
			d['A1_SCW1'] = assemble(SCW1, tensor=d['A1_SCW1'])
			A1.axpy(1.0, d['A1_SCW1'], True)

		if stabilization_parameters['stab_LSIC_NS'] == True:
		    A1.axpy(1.0, d['A1_SLS1'], True)

		if stabilization_parameters['stab_backflow_NS'] == True:
		    A1 += d['BS']

		return A1, b1	        

	def optimized_lhs(self, dt):      
	    
	    A1 = self.matrix['Kij'].copy()
	    A1.axpy(1.0/float(dt), self.matrix['Mij'], True)

	    return A1

	def optimized_rhs(self, AX, u, p):

	    X1 = AX.copy()
	    X1.axpy(-2.0, self.matrix['Kij'], True)

	    b1 = self.matrix['Bij'].copy()
	    b1.axpy(1.0, X1*u.vector())
	    b1.axpy(1.0, self.matrix['Sij']*p.vector())

	    return b1

	def solve_tentative_velocity(self, A, x, b, bcs):
	    
	    [bc.apply(A, b) for bc in bcs]
	    u_solver.solve(A, x.vector(), b)
	    # solve(A, x.vector(), b, 'mumps')



	# Pressure correction
	def assemble_pressure_correction(self, u_, p_, Lm_f, dt):

		p = self.p; q = self.q; dx = self.dx 
		h_f = self.h_f; Re = self.Re; b2 = self.matrix['b2']
	
		L2 = dot(nabla_grad(p_[1]), nabla_grad(q))*dx - (1/dt)*div(u_[0])*q*dx
		
		if stabilization_parameters['stab_PSPG_NS'] == True:
			R = Fluid_problem.residual_tentative_velocity(self, u_[0], u_[1], p_[1], Lm_f, dt)
			L2 -= tau(alpha, u_[0], h_f, Re, dt)*dot(R, nabla_grad(q))*dx	
		
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
  
		b3 = self.matrix['Mij']*u_[0].vector()   
		b3.axpy(-float(dt), self.matrix['Pij']*(p_[0].vector() - p_[1].vector())) 

		return b3

	def solve_velocity_correction(self, x, b, bcs):
		
		A = self.A3
		[bc.apply(b) for bc in bcs]
		u_solver_c.solve(A, x.vector(), b)


	
	# Post-processing functions
	def post_process_data(self, Mpi, u_, p_, t, tsp, text_file_handles):

		Re = self.Re; n = self.n; ds = self.ds

		# Compute drag, lift (note to self : written as per 3D sphere)
		traction = -1*dot(sigma(Re, u_[0], p_[0]), n)
		drag = assemble(dot(traction, self.nx)*ds(7))/(0.5*(PI)/4)
		lift = assemble(dot(traction, self.ny)*ds(7))/(0.5*(PI)/4)

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
