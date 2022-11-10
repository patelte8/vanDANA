from dolfin import *
from ufl import tensors, nabla_div
from .functions import *
import sys

sys.path.insert(0,  '..')
from user_inputs.user_parameters import *
from utilities.read import *


# Note: Lagrange mulitplier is the fictitious force in both domains 
class Lagrange_multiplier_problem:

	# Note to self : This problem is solved on the solid current configuration
	def __init__(self, solid_mesh):

		# --------------------------------

		mesh = solid_mesh.mesh
		dim = mesh.geometry().dim()

		Y  = VectorFunctionSpace(mesh, 'P', fem_degree['displacement_degree'])                                             # Solid displacement
		M  = FunctionSpace(mesh, 'P', fem_degree['pressure_degree'])                                                       # Solid pressure
		Z2 = VectorFunctionSpace(mesh, 'P', fem_degree['lagrange_degree'])                                                 # Lagrange multiplier

		self.Lm = TrialFunction(Z2)
		self.e  = TestFunction(Z2)

		variables = dict(); Lm_ = []; uf_ = []

		uf_, us_ = Function(Y), Function(Y)
		us_.vector()[:] = 0.0
		
		for i in range(2):
			Lm_.insert(i, Function(Z2))
			Lm_[i].vector()[:] = 0.0
		variables.update(Lm_=Lm_, uf_=uf_)	

		self.us_ = us_
		self.F = [Y, M, Z2]
		self.us_ = Function(Y)
		self.nx = tensors.unit_vector(0, dim)
		self.ny = tensors.unit_vector(1, dim)
		self.variables = variables
		self.dx = Measure("dx", domain=mesh)
		self.ds = Measure("ds", domain=mesh, subdomain_data=solid_mesh.get_mesh_boundaries())
		
		# --------------------------------

	def assemble_lagrange_multiplier(self, Lm_, us_, uf_, dt):

		e = self.e; Lm = self.Lm; dx = self.dx

		self.us_.vector()[:] = us_.vector().get_local()[:]

		a6 = assemble(dot(Lm, e)*dx)
		b6 = assemble((1/dt)*dot(self.us_ - uf_, e)*dx + dot(Lm_[1], e)*dx)

		return a6, b6

	def solve_lagrange_multiplier(self, A, x, b):

		solve(A, x.vector(), b, 'bicgstab', 'sor')


	# Compute drag and lift	(Note to self: written as per 2D cylinder)	
	def	post_process_data(self, Mpi, Lm, t, text_file_handles):

		dx = self.dx

		drag = 2*assemble(dot(Lm, self.nx)*dx)
		lift = 2*assemble(dot(Lm, self.ny)*dx)

		Mpi.set_barrier()
		if Mpi.get_rank() == 0:
		    text_file_handles[6].write("{} {} {} {} {} {}" \
		        .format(t, "  ", drag, "  ", lift, "\n"))	



class Solid_temperature_lagrange_multiplier_problem:

	# Note to self : This problem is solved on the solid current configuration
	def __init__(self, solid_mesh):

		# -------------------------------		

		rho, Spht, K, Ld, Sm = calc_non_dimensional_solid_properties(**physical_parameters, **characteristic_scales)
		Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)
		Pe = Re*Pr     
		if not problem_physics['viscous_dissipation']: Ec = 0.0

		mesh = solid_mesh.mesh
		dim = mesh.geometry().dim()

		M = FunctionSpace(mesh, 'P', fem_degree['lagrange_degree'])            		# Temperature based lagrange multiplier
		S = FunctionSpace(mesh, 'P', fem_degree['temperature_degree'])				# Solid Temperature
		
		# --------------------------------

		self.LmT = TrialFunction(M) 		      
		self.ls  = TestFunction(M)			 

		# --------------------------------

		variables = dict();  Ts_=[] ; LmTs_=[]

		for i in range(2):
			Ts_.insert(i, Function(S))
			LmTs_.insert(i, Function(M))
			LmTs_[i].vector()[:] = 0.0

		variables.update(Ts_=Ts_, LmTs_=LmTs_)	
		
		self.F = [M, S]
		self.Re  = Constant(Re)
		self.Ec  = Constant(Ec)
		self.Pe  = Constant(Pe)
		self.rho = Constant(rho)
		self.Spht = Constant(Spht)
		self.K = Constant(K)
		self.variables = variables
		self.n 	= FacetNormal(mesh)
		self.dx = Measure("dx", domain=mesh)
		self.ds = Measure("ds", domain=mesh, subdomain_data=solid_mesh.get_mesh_boundaries())
		
		# --------------------------------


	def assemble_solid_temperature_lagrange_multiplier(self, Ts_, uf_, dt):

		LmT = self.LmT; ls = self.ls; dx = self.dx

		a7 = assemble((-1)*LmT*ls*dx)
		b7 = ((self.rho*self.Spht)-1)*(1/dt)*(Ts_[0] - Ts_[1])*ls*dx + ((0.5*(self.K-1))/(self.Pe))*dot(nabla_grad(Ts_[0]) + nabla_grad(Ts_[1]), nabla_grad(ls))*dx 

		# Viscous dissipation source
		if problem_physics['viscous_dissipation'] == True:

			b7 += Qf(uf_, self.Ec, self.Re)*ls*dx
			
		return a7, assemble(b7)	

	def solve_solid_temperature_lagrange_multiplier(self, A, x, b):
	
		solve(A, x.vector(), b, 'mumps')

	def	post_process_data(self, Mpi, Ts, t, text_file_handles):

		ds = self.ds

		# Compute average nusselt number
		area_sphere = assemble(1*ds)
		nusselt = dot(nabla_grad(Ts), self.n)
		average_nusselt = assemble(nusselt*ds)/area_sphere

		Mpi.set_barrier()
		if Mpi.get_rank() == 0:
		    text_file_handles[8].write("{} {} {} {}".format(t, "  ", average_nusselt, "\n"))				