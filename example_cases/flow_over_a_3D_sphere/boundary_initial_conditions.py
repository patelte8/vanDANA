from dolfin import DirichletBC, Constant, assign, Expression, interpolate, SubDomain, \
					MeshFunction, sqrt, DOLFIN_EPS, near
from .user_parameters import problem_physics
from .problem_specific import *				

constrained_domain = None

# Boundary conditions
def fluid_create_boundary_conditions(fluid_mesh, **V):

	boundaries = fluid_mesh.get_mesh_boundaries()

	# velocity
	bcu_inlet_x = DirichletBC(V['fluid'][0], Constant(1), boundaries, 1)
	bcu_sphere_x = DirichletBC(V['fluid'][0], Constant(0), boundaries, 4)
	bcu_wall_x = DirichletBC(V['fluid'][0], Constant(1), boundaries, 3)
	bcu_x = [bcu_inlet_x, bcu_sphere_x, bcu_wall_x]

	bcu_inlet_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_sphere_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 4)
	bcu_wall_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 3)
	bcu_y = [bcu_inlet_y, bcu_sphere_y, bcu_wall_y]

	bcu_inlet_z = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_sphere_z = DirichletBC(V['fluid'][0], Constant(0), boundaries, 4)
	bcu_wall_z = DirichletBC(V['fluid'][0], Constant(0), boundaries, 3)
	bcu_z = [bcu_inlet_z, bcu_sphere_z, bcu_wall_z]

	bcu = [bcu_x, bcu_y, bcu_z]

	# pressure
	bcp_outlet = DirichletBC(V['fluid'][1], Constant(0), boundaries, 2)
	bcp = [bcp_outlet]

	# Streamfunction
	bcPSI = DirichletBC(V['fluid'][1], 0, boundaries, 4)

	bcs = dict(velocity = bcu, pressure = bcp, streamfunction = bcPSI)

	if problem_physics['solve_temperature'] == True:
		# temperature
		bcT_sphere = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 4)
		bcT_inlet = DirichletBC(V['fluid_temp'][0], Constant(0), boundaries, 1)
		bcT = [bcT_sphere, bcT_inlet]
		
		bcs.update(temperature = bcT)
			
	return bcs


def solid_create_boundary_conditions(solid_mesh, boundaries, dt, **V):

	# Note to self: Boundary conditions are for incremental displacement (delta D)

	# Solid
	if problem_physics['compressible_solid'] == False:
		bcx_sphere = DirichletBC(V['solid'][1].sub(0), Constant((0, 0, 0)), boundaries, 1)
	elif problem_physics['compressible_solid'] == True:
	    bcx_sphere = DirichletBC(V['solid'][0], Constant((0, 0, 0)), boundaries, 1)

	bcx = []  
	return bcx    


# Initial conditions
def fluid_create_initial_conditions(u_, p_, T_):

	# Velocity / pressure
	for i in range(3):
		u_[i][0].vector()[:] = 0.0
		u_[i][1].vector()[:] = 0.0
		u_[i][2].vector()[:] = 0.0
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
	assign(mix.sub(0), interpolate(Expression(('0.0', '0.0', '0.0'), degree = 2), mix.sub(0).function_space().collapse()))

