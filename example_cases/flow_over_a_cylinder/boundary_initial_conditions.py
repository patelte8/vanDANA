from dolfin import DirichletBC, Constant, assign, Expression, interpolate, SubDomain, \
					MeshFunction, sqrt, DOLFIN_EPS, near
from .user_parameters import problem_physics
from .problem_specific import *				

constrained_domain = None

# Boundary conditions
def fluid_create_boundary_conditions(fluid_mesh, **V):

	boundaries = fluid_mesh.get_mesh_boundaries()

	# velocity
	bcu_left_x = DirichletBC(V['fluid'][0], Constant(1), boundaries, 1)
	bcu_cylinder_x = DirichletBC(V['fluid'][0], Constant(0), boundaries, 5)
	bcu_x = [bcu_left_x, bcu_cylinder_x]

	bcu_left_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_cylinder_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 5)
	bcu_y = [bcu_left_y, bcu_cylinder_y]

	bcu = [bcu_x, bcu_y]

	# pressure
	bcp = []

	# Streamfunction
	bcPSI = DirichletBC(V['fluid'][1], 0, boundaries, 5)

	bcs = dict(velocity = bcu, pressure = bcp, streamfunction = bcPSI)

	if problem_physics['solve_temperature'] == True:
		# temperature
		bcT_cylinder = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 5)
		bcT_left = DirichletBC(V['fluid_temp'][0], Constant(0), boundaries, 1)
		bcT = [bcT_cylinder, bcT_left]
		
		bcs.update(temperature = bcT)
			
	return bcs


def solid_create_boundary_conditions(solid_mesh, boundaries, dt, **V):

	# Note to self: Boundary conditions are for incremental displacement (delta D)

	# Solid
	if problem_physics['compressible_solid'] == False:
		bcx_cylinder = DirichletBC(V['solid'][1].sub(0), Constant((0, 0)), boundaries, 1)
	elif problem_physics['compressible_solid'] == True:
	    bcx_cylinder = DirichletBC(V['solid'][0], Constant((0, 0)), boundaries, 1)

	bcx = []  
	return bcx    


# Initial conditions
def fluid_create_initial_conditions(u_, p_, T_):

	# Velocity / pressure
	for i in range(3):
		u_[i][0].vector()[:] = 0.0
		u_[i][1].vector()[:] = 0.0
		p_[i].vector()[:] = 0.0

	# Temperature
	for i in range(3):
		T_[i].vector()[:] = 0.0
	

def solid_create_initial_conditions(Dp_, mix, dt):
	
	# Solid pressure (only defined for incompressible solid)
	assign(mix.sub(1), interpolate(Constant(0), mix.sub(1).function_space().collapse()))

	# Cumulative displacement
	Dp_[0].vector()[:] = 0.0 

	# Incremental displacement (delta D)
	Dp_[1].vector()[:] = 0.0 # V_init*dt
	Dp_[2].vector()[:] = 0.0 # V_init*dt
	assign(mix.sub(0), interpolate(Expression(('0.0', '0.0'), degree = 2), mix.sub(0).function_space().collapse()))

