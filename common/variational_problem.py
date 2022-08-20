from dolfin import *
from ufl import tensors, nabla_div
from .functions import *
from .solver_options import u_solver, p_solver, \
							u_solver_c, t_solver
from .constitutive_eq import sigma
from .stabilizations import *
import sys

sys.path.insert(0,  '..')
from user_inputs.user_parameters import *
from user_inputs.problem_specific import perf, blood_perfusion

PI = 3.14159265

Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)
Pe = Re*Pr     
if not boolean_parameters.get('viscous_dissipation'): Ec = 0.0


class fluid_problem:	
									
	def __init__(self, S, boundaries, bool_stream):
		
		mesh = S[0].mesh()
		dim = mesh.geometry().dim()
		
		# --------------------------------

		self.u1 = TrialFunction(S[0])
		self.v = TestFunction(S[0])
		self.p = TrialFunction(S[1])
		self.q = TestFunction(S[1])

		variables = dict(); u_ = []; p_ = []
		for i in range(3):
			u_.insert(i, Function(S[0]))
			p_.insert(i, Function(S[1]))

		variables.update(u_=u_, p_=p_)	
		# LmTf_n = Function(S[1])
			
		vort, psi = Function(S[1]), Function(S[1])
		vort.vector()[:] = 0.0; psi.vector()[:] = 0.0
		variables.update(vort=vort, psi=psi)

		# --------------------------------

		# Body force
		f = Constant((0,)*dim)
		if boolean_parameters.get('body_force') == True: f = Constant(((1/(Fr*Fr))*f_dir(dim)))	
		
		# --------------------------------

		self.nx = tensors.unit_vector(0, dim)
		self.ny = tensors.unit_vector(1, dim)
		if dim == 3: self.nz = tensors.unit_vector(2, dim)

		self.A2 = None
		self.A3 = None
		self.null_space = VectorSpaceBasis([])
		self.matrix = dict(Mij=None, Kij = None, Sij = None, Bij = None, Pij = None, \
						A1_as1 = None, A1_SCW1 = None, A1_SLS1 = None, A1_Cij = None, \
						BS = None, b1_Ls1 = None, b1_Cij = None)
		self.f = f
		self.dim = dim
		self.variables = variables
		self.bool_stream = bool_stream
		self.h_f = CellDiameter(mesh)
		self.h_f_X = project(CellDiameter(mesh), S[1])
		
		self.dx = Measure("dx", domain=mesh)
		self.ds = Measure("ds", domain=mesh, subdomain_data=boundaries)	
		self.n = FacetNormal(mesh)
		if dim == 2: self.tang = as_vector([self.n[1], -self.n[0]])

		self.Re = Constant(Re)
		self.Fr = Constant(Fr)

		# --------------------------------

	def pre_assemble(self, px, FS, bcs):

		d = self.matrix; Re = self.Re
		u1 = self.u1; v = self.v; p = self.p; q = self.q
		dx = self.dx; ds = self.ds; n = self.n; f = self.f
		
		d['Mij'] = self.A3 = assemble(dot(u1, v)*dx, tensor=d['Mij'])                                       # Mass matrix 
		d['Kij'] = assemble(inner((1/Re)*epsilon(u1), epsilon(v))*dx - \
				   dot((1/Re)*nabla_grad(0.5*u1)*n, v)*ds, tensor=d['Kij'])            	   					# Viscous matrix
		d['Sij'] = assemble(inner(p*Identity(len(u1)), epsilon(v))*dx - dot(p*n, v)*ds, tensor=d['Sij'])    # Pressure matrix
		d['Bij'] = assemble(dot(f, v)*dx, tensor=d['Bij'])                                            		# Body-force vector
		d['Pij'] = assemble(dot(nabla_grad(p), v)*dx, tensor=d['Pij'])                                  	# Pressure-gradient matrix

		self.A2 = assemble(dot(nabla_grad(p), nabla_grad(q))*dx)
		if bcs['pressure'] == []:
		    self.null_space = attach_nullspace(self.A2, px, FS)

		# Boundary conditions
		[bc.apply(self.A2) for bc in bcs['pressure']]
		[bc.apply(self.A3) for bc in bcs['velocity']]

		# Stabilization terms
		if stabilization_parameters.get('stab_LSIC_NS') == True:
			SLS1 = tau_lsic(Re)*nabla_div(u1)*nabla_div(v)*dx
			d['A1_SLS1'] = assemble(SLS1, tensor=d['A1_SLS1'])

		if stabilization_parameters.get('stab_backflow_NS') == True:
			d['BS'] = assemble(Constant(1e4)*dot((Identity(self.dim) - outer(n,n))*u1, v)*ds(2), tensor=d['BS'])	



	# Predict tentative velocity 	
	def residual_tentative_velocity(self, u_0, u_1, p_, dt):
	
		Re = self.Re; f = self.f

		U = 0.5*(u_0 + u_1)
		R = (u_0 - u_1)/dt + dot(u_1, nabla_grad(U)) - nabla_div(sigma(Re, U, p_)) - f	
		
		return R	
	
	def assemble_tentative_velocity(self, u_, p_, dt):

		d = self.matrix; Re = self.Re
		u1 = self.u1; v = self.v; dx = self.dx
		h_f = self.h_f; dx = self.dx; ds = self.ds

		A1 = fluid_problem.optimized_lhs(self, dt)
		b1 = fluid_problem.optimized_rhs(self, A1, u_[1], p_[1])

		# Convective terms 
		u_ab = 1.5*u_[1] - 0.5*u_[2]
		Cij  = dot(dot(u_ab, nabla_grad(u_[1])), v)*dx                                       			
		Cij_imp = dot(dot(u_ab, nabla_grad(u1)), v)*dx
		d['A1_Cij'] = assemble(Cij_imp, tensor=d['A1_Cij'])
		d['b1_Cij'] = assemble(Cij, tensor=d['b1_Cij'])

		A1.axpy(0.5, d['A1_Cij'], True)
		b1.axpy(-0.5, d['b1_Cij'])

		# Stabilization terms 
		if stabilization_parameters.get('stab_SUPG_NS') == True: 
			R = fluid_problem.residual_tentative_velocity(self, u1, u_[1], p_[1], dt)
			S1 = tau(alpha, u_[1], h_f, Re, dt)*dot(R, Pop(u_[1], v))*dx
			as1 = lhs(S1); Ls1 = rhs(S1)	    
			d['A1_as1'] = assemble(as1, tensor=d['A1_as1'])
			d['b1_Ls1'] = assemble(Ls1, tensor=d['b1_Ls1'])
			A1.axpy(1.0, d['A1_as1'], True)
			b1.axpy(1.0, d['b1_Ls1'])

		if stabilization_parameters.get('stab_cross_NS') == True:
			R = fluid_problem.residual_tentative_velocity(self, u_[1], u_[2], p_[2], dt)
			SCW1 = tau_cw(C_cw, u_[1], h_f, Re, R)*inner(Pop_CW(u_[1], u1), nabla_grad(v))*dx
			d['A1_SCW1'] = assemble(SCW1, tensor=d['A1_SCW1'])
			A1.axpy(1.0, d['A1_SCW1'], True)

		if stabilization_parameters.get('stab_LSIC_NS') == True:
		    A1.axpy(1.0, d['A1_SLS1'], True)

		if stabilization_parameters.get('stab_backflow_NS') == True:
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



	# Pressure correction
	def assemble_pressure_correction(self, u_, p_, dt):

		p = self.p; q = self.q; dx = self.dx 
		h_f = self.h_f; Re = self.Re
	
		L2 = dot(nabla_grad(p_[1]), nabla_grad(q))*dx - (1/dt)*div(u_[0])*q*dx
		
		if stabilization_parameters.get('stab_PSPG_NS') == True:
			R = fluid_problem.residual_tentative_velocity(self, u_[0], u_[1], p_[1], dt)
			L2 -= tau(alpha, u_[0], h_f, Re, dt)*dot(R, nabla_grad(q))*dx	
		
		b2 = assemble(L2)
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
		        .format(t, "  ", drag, "  ", lift, "  "))    

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





class temperature_problem:	
									
	def __init__(self, S, boundaries):
		
		mesh = S[2].mesh()
		dim = mesh.geometry().dim()

		# --------------------------------

		self.Tp = TrialFunction(S[2])
		self.ttf = TestFunction(S[2])
		
		variables = dict(); T_ = []
		for i in range(4):
			T_.insert(i, Function(S[2]))

		variables.update(T_=T_)	

		# --------------------------------

		self.matrix = dict(Tij=None, Wij = None, BPij = None, BPj = None, \
					   A4_as4 = None, A4_SCW4 = None, A4_R1ij = None, \
					   b4_R2ij = None, b4_Ls4 = None, b4_Qij = None)

		self.dim = dim 
		self.variables = variables
		self.dx = Measure("dx", domain=mesh)
		self.ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
		self.area_sphere = assemble(1*self.ds(7))

		self.h_f  = CellDiameter(mesh)
		self.n 	  = FacetNormal(mesh)
		self.Re   = Constant(Re)
		self.Ec   = Constant(Ec)
		self.Pe   = Constant(Pe)
		self.perf = Constant(perf)

		# --------------------------------

	def pre_assemble(self):

		Tp = self.Tp; ttf = self.ttf
		d = self.matrix; Pe = self.Pe; dx = self.dx
				
		d['Tij']  = assemble(Tp*ttf*dx, tensor=d['Tij'])
		d['Wij']  = assemble((0.5/(Pe))*dot(nabla_grad(Tp), nabla_grad(ttf))*dx, tensor=d['Wij'])
		


	# Temperature equation 	
	def residual_energy_eq(self, T_0, T_1, u, dt):
	
		Re = self.Re; Ec = self.Ec; Pe = self.Pe

		R = (T_0 - T_1)/dt + dot(u, 0.5*nabla_grad(T_0 + T_1)) - (0.5/Pe)*nabla_div(nabla_grad(T_0) + nabla_grad(T_1)) - Qf(u, Ec, Re) #- LmTf_n	
		return R

	def assemble_temperature(self, T_, u, dt):

		d = self.matrix 
		Ec = self.Ec; Re = self.Re; Pe = self.Pe
		Tp = self.Tp; ttf = self.ttf; h_f = self.h_f
		dx = self.dx; ds = self.ds; pf = self.perf

		A4 = temperature_problem.optimized_lhs(self, dt)
		b4 = temperature_problem.optimized_rhs(self, A4, T_[1])

		# Advection terms
		R1ij = dot(u, nabla_grad(Tp))*ttf*dx				
		R2ij = dot(u, nabla_grad(T_[1]))*ttf*dx
		d['A4_R1ij'] = assemble(R1ij, tensor=d['A4_R1ij'])
		d['b4_R2ij'] = assemble(R2ij, tensor=d['b4_R2ij'])
		A4.axpy(0.5,  d['A4_R1ij'], True)
		b4.axpy(-0.5, d['b4_R2ij'])
		
		# Viscous dissipation source
		if boolean_parameters.get('viscous_dissipation') == True:
			Qij  = Qf(u, Ec, Re)*ttf*dx
			d['b4_Qij'] = assemble(Qij, tensor=d['b4_Qij'])
			b4.axpy(1.0, d['b4_Qij'])

		# Blood perfusion terms at boundary
		if blood_perfusion == True:
			BP  = pf*PFE(T_[1])*Tp*ttf*ds(1)
			BPb = pf*PFE(T_[1])*ttf*ds(1)
			d['BPij'] = assemble(BP, tensor=d['BPij'])
			d['BPj'] = assemble(BPb, tensor=d['BPj'])
			A4 += d['BPij']; b4 += d['BPj']
		
		# Stabilization terms
		if stabilization_parameters.get('stab_SUPG_HT') == True:
			R = temperature_problem.residual_energy_eq(self, Tp, T_[1], u, dt)
			S4 = tau(alpha, u, h_f, Pe, dt)*dot(R, Pop(u, ttf))*dx
			as4 = lhs(S4); Ls4 = rhs(S4)
			d['A4_as4'] = assemble(as4, tensor=d['A4_as4'])
			d['b4_Ls4'] = assemble(Ls4, tensor=d['b4_Ls4'])
			A4.axpy(1.0, d['A4_as4'], True)
			b4.axpy(1.0, d['b4_Ls4'])
				
		if stabilization_parameters.get('stab_cross_HT') == True:
			R = temperature_problem.residual_energy_eq(self, T_[1], T_[2], u, dt)
			SCW4 = inner(tau_cw(C_cw, T_[1], h_f, Pe, R)*Pop_CW(u, Tp), nabla_grad(ttf))*dx
			d['A4_SCW4'] = assemble(SCW4, tensor=d['A4_SCW4'])
			A4.axpy(1.0, d['A4_SCW4'], True)
		
		return A4, b4
		
	def optimized_lhs(self, dt):      
	    
	    F1 = self.matrix['Wij'].copy()
	    F1.axpy(1.0/float(dt), self.matrix['Tij'], True)

	    return F1

	def optimized_rhs(self, EZ, T):

	    E1 = EZ.copy()
	    E1.axpy(-2.0, self.matrix['Wij'], True)
	    b4 = E1*T.vector()
	    
	    return b4

	def solve_temperature(self, A, x, b, bcs):
	    
	    [bc.apply(A, b) for bc in bcs]
	    t_solver.solve(A, x.vector(), b)    		



	
	# Post-processing functions
	def post_process_data(self, Mpi, T_, t, text_file_handles):

		n = self.n; ds = self.ds
		
		# Compute average nusselt number
		nusselt = dot(nabla_grad(T_[0]), n)
		average_nusselt = assemble(nusselt*ds(7))/self.area_sphere

		Mpi.set_barrier()
		if Mpi.get_rank() == 0:
		    text_file_handles[0].write("{} {}".format(average_nusselt, "\n"))    



            
