from dolfin import *
from ufl import tensors, nabla_div
from .functions import *
from fenicstools import interpolate_nonmatching_mesh
from .solver_options import t_solver
from .constitutive_eq import *
from .fem_stabilizations import *
import sys

sys.path.insert(0,  '..')
from user_inputs.user_parameters import *
from user_inputs.boundary_initial_conditions import constrained_domain
from user_inputs.problem_specific import perf, blood_perfusion
from utilities.read import *


PI = 3.14159265



class Fluid_temperature_problem:	
									
	def __init__(self, fluid_mesh):
		
		Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)
		Pe = Re*Pr     
		if not problem_physics['viscous_dissipation']: Ec = 0.0

		mesh = fluid_mesh.mesh
		dim = mesh.geometry().dim()

		G  = FunctionSpace(mesh, 'P', fem_degree['temperature_degree'], constrained_domain = constrained_domain)           # Fluid temperature

		# --------------------------------

		self.Tp = TrialFunction(G)
		self.ttf = TestFunction(G)
		
		variables = dict(); T_ = []
		for i in range(3):
			T_.insert(i, Function(G))

		LmTf_ = Function(G)
		LmTf_.vector()[:] = 0.0	
		variables.update(T_=T_, LmTf_=LmTf_)	

		# --------------------------------

		self.A4 = None
		self.matrix = dict(Tij=None, Wij = None, BPij = None, BPj = None, \
					   A4_SCW4 = None, A4_R1ij = None, b4_R2ij = None, b4_Ls4 = None, b4_Qij = None)

		self.F = [G]
		self.dim = dim 
		self.variables = variables
		self.dx = Measure("dx", domain=mesh)
		self.ds = Measure("ds", domain=mesh, subdomain_data=fluid_mesh.get_mesh_boundaries())
		self.area_sphere = assemble(1*self.ds(4))

		self.h_f  = CellDiameter(mesh)
		self.n 	  = FacetNormal(mesh)
		self.Re   = Constant(Re)
		self.Ec   = Constant(Ec)
		self.Pe   = Constant(Pe)
		self.perf = Constant(perf)

		# --------------------------------

	def pre_assemble(self, dt):

		Tp = self.Tp; ttf = self.ttf
		d = self.matrix; Pe = self.Pe; dx = self.dx
				
		d['Tij']  = assemble(Tp*ttf*dx, tensor=d['Tij'])
		d['Wij']  = assemble((0.5/(Pe))*dot(nabla_grad(Tp), nabla_grad(ttf))*dx, tensor=d['Wij'])
		
		if time_control['adjustable_timestep'] == False:
			self.A4 = self.matrix['Wij'].copy()
			self.A4.axpy(1.0/float(dt), self.matrix['Tij'], True)

	# Temperature equation 	
	def residual_energy_eq(self, T_0, T_1, u, LmTf_, dt):
	
		Re = self.Re; Ec = self.Ec; Pe = self.Pe

		R = (T_0 - T_1)/dt + dot(u, 0.5*nabla_grad(T_0 + T_1)) - (0.5/Pe)*nabla_div(nabla_grad(T_0) + nabla_grad(T_1)) - Qf(u, Ec, Re) - LmTf_	
		return R

	def assemble_temperature(self, T_, u, LmTf_, dt):

		d = self.matrix 
		Ec = self.Ec; Re = self.Re; Pe = self.Pe
		Tp = self.Tp; ttf = self.ttf; h_f = self.h_f
		dx = self.dx; ds = self.ds; pf = self.perf

		if time_control['adjustable_timestep'] == False:
			A4 = self.A4.copy()
		else:	
			A4 = Fluid_temperature_problem.optimized_lhs(self, dt)

		b4 = Fluid_temperature_problem.optimized_rhs(self, A4, T_[1])

		# Advection terms
		R1ij = dot(u, nabla_grad(Tp))*ttf*dx				
		R2ij = dot(u, nabla_grad(T_[1]))*ttf*dx
		d['A4_R1ij'] = assemble(R1ij, tensor=d['A4_R1ij'])
		d['b4_R2ij'] = assemble(R2ij, tensor=d['b4_R2ij'])
		A4.axpy(0.5,  d['A4_R1ij'], True)
		b4.axpy(-0.5, d['b4_R2ij'])
		
		# Viscous dissipation source
		if problem_physics['viscous_dissipation'] == True:
			Qij  = Qf(u, Ec, Re)*ttf*dx
			d['b4_Qij'] = assemble(Qij, tensor=d['b4_Qij'])
			b4.axpy(1.0, d['b4_Qij'])

		# FSI temperature based lagrange multiplier
		if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
			b4.axpy(1.0, d['Tij']*LmTf_.vector())	

		# Blood perfusion terms at boundary
		if blood_perfusion == True:
			BP  = pf*PFE(T_[1])*Tp*ttf*ds(1)
			BPb = pf*PFE(T_[1])*ttf*ds(1)
			d['BPij'] = assemble(BP, tensor=d['BPij'])
			d['BPj'] = assemble(BPb, tensor=d['BPj'])
			A4 += d['BPij']; b4 += d['BPj']
		
		# Stabilization terms
		if stabilization_parameters['SUPG_HT'] == True:
			R = Fluid_temperature_problem.residual_energy_eq(self, T_[1], T_[2], u, LmTf_, dt)
			d['b4_Ls4'] = assemble(tau(alpha, u, h_f, Pe, dt)*dot(R, Pop(u, ttf))*dx, tensor=d['b4_Ls4'])
			b4.axpy(-1.0, d['b4_Ls4'])
				
		if stabilization_parameters['crosswind_HT'] == True:
			R = Fluid_temperature_problem.residual_energy_eq(self, T_[1], T_[2], u, LmTf_, dt)
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
	    # solve(A, x.vector(), b, 'mumps')    		



	
	# Post-processing functions
	def post_process_data(self, Mpi, T_, t, text_file_handles):

		n = self.n; ds = self.ds
		
		# Compute average nusselt number
		nusselt = dot(nabla_grad(T_[0]), n)
		average_nusselt = assemble(nusselt*ds(4))/self.area_sphere

		Mpi.set_barrier()
		if Mpi.get_rank() == 0:
		    text_file_handles[4].write("{} {} {} {}".format(t, "  ", average_nusselt, "\n"))    
