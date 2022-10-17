from dolfin import *
from ufl import tensors, nabla_div
from .functions import *
from fenicstools import interpolate_nonmatching_mesh
from .solver_options import solid_displacement_parameters, FFC_parameters
from .constitutive_eq import *
import sys

sys.path.insert(0,  '..')
from user_inputs.user_parameters import *
from user_inputs.boundary_initial_conditions import constrained_domain
from utilities.read import *


PI = 3.14159265



class Solid_problem:

	# Note to self : This problem is solved on the solid reference configuration
	def __init__(self, solid_mesh_R):

		rho, Spht, K, Ld, Sm = calc_non_dimensional_solid_properties(**physical_parameters, **characteristic_scales)
		Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)

		mesh = solid_mesh_R.mesh
		dim = mesh.geometry().dim()

		R1 = VectorElement('P', mesh.ufl_cell(), fem_degree.get('displacement_degree'))        
		T1 = FiniteElement('P', mesh.ufl_cell(), fem_degree.get('pressure_degree'))                          # Solid pressure 
		R  = FunctionSpace(mesh, R1)                                                                         # Solid displacement
		X  = FunctionSpace(mesh, MixedElement([R1, T1]))
		Z  = VectorFunctionSpace(mesh, 'P', fem_degree.get('lagrange_degree'))                               # Lagrange multiplier 

		# --------------------------------

		(h, j)  = TestFunctions(X)
		self.hc = TestFunction(R) 
		self.h  = h
		self.j  = j

		# --------------------------------

		variables = dict();  Dp_=[] 

		for i in range(3):
			Dp_.insert(i, Function(R))
		self.Lm_ = Function(Z)
		self.uf_ = Function(R)
		us_ = Function(R)
		mix = Function(X)
		ps_ = Function(FunctionSpace(mesh, T1))
		J_  = Function(FunctionSpace(mesh, T1))	

		variables.update(Dp_=Dp_, mix=mix, us_=us_, ps_=ps_, J_=J_)

		# Body force
		f = Constant((0,)*dim)
		if problem_physics.get('body_force') == True: f = Constant(((1/(Fr*Fr))*f_dir(dim)))	
		
		self.F = [R, X, Z]
		self.f = f
		self.n = FacetNormal(mesh)
		self.nx = tensors.unit_vector(0, dim)
		self.ny = tensors.unit_vector(1, dim)
		self.dim = dim
		self.Re = Constant(Re)
		self.Ld = Constant(Ld)
		self.Sm = Constant(Sm)
		self.rho = Constant(rho)
		self.variables = variables
		self.dx = Measure("dx", domain=mesh)
		self.ds = Measure("ds", domain=mesh, subdomain_data=solid_mesh_R.get_mesh_boundaries())
		
		# --------------------------------


	def assemble_solid_problem(self, compressible_solid, Dp_, mix, uf_, Lm_, dt):	

		rho = self.rho; Ld = self.Ld; Sm = self.Sm; j = self.j
		h = self.h; hc = self.hc; dx = self.dx; f = self.f

		self.uf_.vector()[:] = uf_.vector().get_local()[:]
		self.Lm_.vector()[:] = Lm_.vector().get_local()[:]

		# Define incompressible solid problem
		if compressible_solid == False:
			D_, ps_ = split(mix)
			a5 = rho*(1/(dt*dt))*dot(D_, h)*dx + inner(nabla_grad(h).T, stress_inc(Dp_[0] + D_, ps_, Sm))*dx + dot(J(F(Dp_[0] + D_))-1, j)*dx
			b5 = (rho-1)*(1/(dt*dt))*dot(Dp_[2], h)*dx + (rho-1)*dot(f, h)*dx      

			if problem_physics.get('solve_FSI') == True:
				b5 += (1/dt)*dot(self.uf_, h)*dx - dot(self.Lm_, h)*dx	

		# Define compressible solid problem	
		elif compressible_solid == True:
			a5 = rho*(1/(dt*dt))*dot(Dp_[1], hc)*dx + inner(nabla_grad(hc).T, stress_c(Dp_[0] + Dp_[1], Ld, Sm))*dx 
			b5 = (rho-1)*(1/(dt*dt))*dot(Dp_[2], hc)*dx + (rho-1)*dot(f, hc)*dx 

			if problem_physics.get('solve_FSI') == True:
				b5 += (1/dt)*dot(self.uf_, hc)*dx - dot(self.Lm_, hc)*dx
			
		a5 -= b5
		return a5	

	def solve_solid_displacement(self, compressible_solid, a5, Dp_, mix, ps_, p_, bcs):

		if compressible_solid == False:
		    solve(a5 == 0, mix, bcs, solver_parameters = solid_displacement_parameters,
		    							form_compiler_parameters = FFC_parameters)
		    (D, ps) = mix.split(deepcopy=True)

		    Dp_.vector()[:] = D.vector().get_local()[:]
		    ps_.vector()[:] = ps.vector().get_local()[:]         

		elif compressible_solid == True:
			solve(a5 == 0, Dp_, bcs, solver_parameters = solid_displacement_parameters,
										form_compiler_parameters = FFC_parameters)

			# Note to self: if it's a compressible solid, solid pressure is the same as fluid pressure
			ps_.assign(interpolate_nonmatching_mesh(p_, mix.sub(1).function_space().collapse()))

		# --------------------------------		
	
	# Initial guess for newtons solver
	def change_initial_guess(Dp_, mix):
		
		Dp_.vector().zero()
		assign(mix.sub(1), interpolate(Constant(0), mix.sub(1).function_space().collapse()))
		assign(mix.sub(0), interpolate(Expression(('0.0', '0.0', '0.0'), degree = 2), mix.sub(0).function_space().collapse()))	

	def compute_jacobian(self, J_, Dp):

		J_.assign(project(J(F(Dp)), J_.function_space()))


	# Compute drag and lift	(Note to self: written as per 2D cylinder)
	def	post_process_data(self, Mpi, u, p, t, text_file_handles):

		ds = self.ds; n = self.n

		traction = -1*dot(sigma(self.Re, u, p), n)
		drag = 2*assemble(dot(traction, self.nx)*ds)
		lift = 2*assemble(dot(traction, self.ny)*ds)

		Mpi.set_barrier()
		if Mpi.get_rank() == 0:
		    text_file_handles[5].write("{} {} {} {} {} {}".format(t, "  ", drag, "  ", lift, "\n")) 