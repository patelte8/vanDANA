from dolfin import DirichletBC, Constant, assign, Expression, interpolate, SubDomain, \
					MeshFunction, sqrt, DOLFIN_EPS, near
from .user_parameters import problem_physics
from .problem_specific import *				

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

class Point_pressure(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 4.2) and near(x[1], 5.) and near(x[2], 2.)

# Boundary conditions
def fluid_create_boundary_conditions(fluid_mesh, **V):

	boundaries = fluid_mesh.get_mesh_boundaries()

	# velocity
	bcu_LIPV_x = DirichletBC(V['fluid'][0], inflow[1][0], boundaries, 3)
	bcu_LSPV_x = DirichletBC(V['fluid'][0], inflow[0][0], boundaries, 4)
	bcu_RIPV_x = DirichletBC(V['fluid'][0], inflow[3][0], boundaries, 5)
	bcu_RSPV_x = DirichletBC(V['fluid'][0], inflow[2][0], boundaries, 6)
	bcu_wall_x = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_CB_x   = DirichletBC(V['fluid'][0], Constant(0), boundaries, 7)	
	bcu_x = [bcu_LIPV_x, bcu_LSPV_x, bcu_RIPV_x, bcu_RSPV_x, bcu_wall_x, bcu_CB_x]

	bcu_LIPV_y = DirichletBC(V['fluid'][0], inflow[1][1], boundaries, 3)
	bcu_LSPV_y = DirichletBC(V['fluid'][0], inflow[0][1], boundaries, 4)
	bcu_RIPV_y = DirichletBC(V['fluid'][0], inflow[3][1], boundaries, 5)
	bcu_RSPV_y = DirichletBC(V['fluid'][0], inflow[2][1], boundaries, 6)
	bcu_wall_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_CB_y   = DirichletBC(V['fluid'][0], Constant(0), boundaries, 7)	
	bcu_y = [bcu_LIPV_y, bcu_LSPV_y, bcu_RIPV_y, bcu_RSPV_y, bcu_wall_y, bcu_CB_y]

	bcu_LIPV_z = DirichletBC(V['fluid'][0], inflow[1][2], boundaries, 3)
	bcu_LSPV_z = DirichletBC(V['fluid'][0], inflow[0][2], boundaries, 4)
	bcu_RIPV_z = DirichletBC(V['fluid'][0], inflow[3][2], boundaries, 5)
	bcu_RSPV_z = DirichletBC(V['fluid'][0], inflow[2][2], boundaries, 6)
	bcu_wall_z = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_CB_z   = DirichletBC(V['fluid'][0], Constant(0), boundaries, 7)	
	bcu_z = [bcu_LIPV_z, bcu_LSPV_z, bcu_RIPV_z, bcu_RSPV_z, bcu_wall_z, bcu_CB_z]

	bcu = [bcu_x, bcu_y, bcu_z]

	# pressure
	bcp_MV = DirichletBC(V['fluid'][1], Constant(0), boundaries, 2)
	bcp = [bcp_MV]

	# Streamfunction
	wall  = 'on_boundary'
	bcPSI = DirichletBC(V['fluid'][1], 0, wall)

	bcs = dict(velocity = bcu, pressure = bcp, streamfunction = bcPSI)

	if problem_physics['solve_temperature'] == True:
		# temperature
		bcT_LIPV = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 3)
		bcT_LSPV = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 4)
		bcT_RIPV = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 5)
		bcT_RSPV = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 6)
		bcT_CB   = DirichletBC(V['fluid_temp'][0], Constant(0), boundaries, 7)
		bcT = [bcT_LIPV, bcT_LSPV, bcT_RIPV, bcT_RSPV, bcT_CB]
		
		bcs.update(temperature = bcT)
			
	return bcs


def solid_create_boundary_conditions(solid_mesh, boundaries, dt, **V):

	# Note to self: Boundary conditions are for incremental displacement (delta D)

	# Solid
	if problem_physics['compressible_solid'] == False:
		bcx_cylinder = DirichletBC(V['solid'][1].sub(0), Constant((0, 0)), boundaries, 1) #, method="pointwise")
	elif problem_physics['compressible_solid'] == True:
	    bcx_cylinder = DirichletBC(V['solid'][0], Constant((0, 0)), boundaries, 1) #, method="pointwise")

	bcx = [bcx_cylinder]  
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
		T_[i].vector()[:] = 1.0
	

def solid_create_initial_conditions(Dp_, mix, dt):
	
	# Solid pressure (only defined for incompressible solid)
	assign(mix.sub(1), interpolate(Constant(0), mix.sub(1).function_space().collapse()))

	# Cumulative displacement
	Dp_[0].vector()[:] = 0.0 

	# Incremental displacement (delta D)
	Dp_[1].vector()[:] = 0.0 # V_init*dt
	Dp_[2].vector()[:] = 0.0 # V_init*dt
	assign(mix.sub(0), interpolate(Expression(('0.0', '0.0'), degree = 2), mix.sub(0).function_space().collapse()))

