from dolfin import DirichletBC, Constant

def create_boundary_conditions(boundaries, inflow, **V):

	# velocity
	bcu_LIPV = DirichletBC(V['fluid'][0], inflow[1], boundaries, 3)  	
	bcu_LSPV = DirichletBC(V['fluid'][0], inflow[0], boundaries, 4)  	
	bcu_RIPV = DirichletBC(V['fluid'][0], inflow[3], boundaries, 5)  	
	bcu_RSPV = DirichletBC(V['fluid'][0], inflow[2], boundaries, 6)  	
	bcu_wall = DirichletBC(V['fluid'][0], Constant((0, 0, 0)), boundaries, 1)
	bcu_CB   = DirichletBC(V['fluid'][0], Constant((0, 0, 0)), boundaries, 7)
	bcu = [bcu_LIPV, bcu_LSPV, bcu_RIPV, bcu_RSPV, bcu_wall, bcu_CB]

	# pressure
	bcp_MV   = DirichletBC(V['fluid'][1], Constant(0), boundaries, 2)
	bcp      = [bcp_MV]

	# temperature
	bcT_LIPV = DirichletBC(V['fluid'][2], Constant(1), boundaries, 3)
	bcT_LSPV = DirichletBC(V['fluid'][2], Constant(1), boundaries, 4)
	bcT_RIPV = DirichletBC(V['fluid'][2], Constant(1), boundaries, 5)
	bcT_RSPV = DirichletBC(V['fluid'][2], Constant(1), boundaries, 6)
	bcT_CB   = DirichletBC(V['fluid'][2], Constant(0), boundaries, 7)
	bcT = [bcT_LIPV, bcT_LSPV, bcT_RIPV, bcT_RSPV, bcT_CB]

	# Streamfunction
	wall  = 'on_boundary'
	bcPSI = DirichletBC(V['fluid'][1], 0, wall)

	return dict(velocity = bcu, pressure = bcp, temperature = bcT, streamfunction = bcPSI)


def create_initial_conditions(u_, p_, T_):

	# temperature
	for i in range(4):
		T_[i].vector()[:] = 1.0
	
	for i in range(3):
		u_[i].vector()[:] = 0.0
		p_[i].vector()[:] = 0.0
	