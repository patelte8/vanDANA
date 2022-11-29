from dolfin import DirichletBC, Constant, assign, Expression, interpolate, SubDomain, \
					MeshFunction, sqrt, DOLFIN_EPS, near
from .user_parameters import problem_physics				

class PeriodicDomain(SubDomain):

    def inside(self, x, on_boundary):
        return bool(x[2] < DOLFIN_EPS and x[2] > -DOLFIN_EPS and on_boundary)

    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1]
        y[2] = x[2] - 1.0

constrained_domain = None

class RegionOfInterest(SubDomain):
    def inside(self,x,on_boundary):
        tol = 1e-6
        return sqrt(((x[0] - 2.0)*(x[0] - 2.0)) + ((x[1] - 2.0)*(x[1] - 2.0))) < 0.5 + tol


inflow_profile = Expression((('6.0*x[1]*(4.1 - x[1])/(4.1*4.1)', '0')), degree=2)

class Point_pressure(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 4.2) and near(x[1], 5.) and near(x[2], 2.)

# Boundary conditions
def fluid_create_boundary_conditions(fluid_mesh, **V):

	boundaries = fluid_mesh.get_mesh_boundaries()

	# velocity
	bcu_left = DirichletBC(V['fluid'][0], inflow_profile, boundaries, 1)
	bcu_bottom = DirichletBC(V['fluid'][0], Constant((0, 0)), boundaries, 2)
	bcu_top = DirichletBC(V['fluid'][0], Constant((0, 0)), boundaries, 4)
	bcu = [bcu_bottom, bcu_top, bcu_left]

	# pressure
	bcp_right = DirichletBC(V['fluid'][1], Constant(0), boundaries, 3)
	bcp = [bcp_right]

	# Streamfunction
	wall  = 'on_boundary'
	bcPSI = DirichletBC(V['fluid'][1], 0, wall)

	bcs = dict(velocity = bcu, pressure = bcp, streamfunction = bcPSI)

	if problem_physics.get('solve_temperature') == True:
		# temperature
		bcT_left = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 1)
		bcT_top = DirichletBC(V['fluid_temp'][0], Constant(0), boundaries, 4)
		bcT = [bcT_left, bcT_top]
		
		bcs.update(temperature = bcT)
			
	return bcs


def solid_create_boundary_conditions(solid_mesh_R, compressible_solid, dt, **V):

	# boundaries = solid_mesh_R.get_mesh_boundaries()
	cylinder = 0; Complement_cylinder = 1         
	mesh_part = MeshFunction("size_t", solid_mesh_R.mesh, 0, Complement_cylinder)     
	RegionOfInterest().mark(mesh_part, cylinder)
	subdomain_R = RegionOfInterest()

	# Note to self: Boundary conditions are for incremental displacement (delta D)

	# Solid
	if compressible_solid == False:
		bcx_cylinder = DirichletBC(V['solid'][1].sub(0), Constant((0, 0)), subdomain_R)
	elif compressible_solid == True:
	    bcx_cylinder = DirichletBC(V['solid'][0], Constant((0, 0)), subdomain_R)

	bcx = [bcx_cylinder]  
	return bcx    


# Initial conditions
def fluid_create_initial_conditions(u_, p_, T_):

	# init = Expression(('1.0', '0.0'), degree = 2)

	# Velocity / pressure
	for i in range(3):
		# assign(u_[i], interpolate(init, u_[i].function_space()))
		u_[i].vector()[:] = 0.0
		p_[i].vector()[:] = 0.0

	# temperature
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

