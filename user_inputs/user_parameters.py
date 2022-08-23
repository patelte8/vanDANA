from dolfin import sqrt, Constant
from ufl import tensors   

boolean_parameters = dict(
	
	restart = False,  				# Restart parameter
	viscous_dissipation = False,	# Heat release due to viscous gradients 
	body_force = False,      		# Gravitational force (uniform volumetric force)       	     
)		

# FEM stabilization
stabilization_parameters = dict(	

	# Navier-stokes
	stab_SUPG_NS = True,
	stab_PSPG_NS = True,
	stab_cross_NS = True,
	stab_LSIC_NS = True,
	stab_backflow_NS = True,

	# Energy-equation
	stab_SUPG_HT = True,
	stab_cross_HT = True
)

# Physical parameters    
physical_parameters = dict(
 
	g = 9.81,									# Gravity (m/s2)						  

	# Fluid 
	rho_f = 1060,								# Density (kg/m3)
	mu = 3.0*1e-6,								# Kinematic viscosity (m2/s)
	Spht_f = 3800,         						# Specific heat (J/kg.C)
	K_f = 0.52									# Thermal conductivity (W/m.C)
)

# Body force direction
def f_dir(dim):
	
	n = -1*tensors.unit_vector(1, dim)			# -ve y direction (by default)
	return n

# Characteristic scales
characteristic_scales = dict(
	
	Lsc = 0.02026,		            # m          
	Vsc = 0.093,         		    # m/s
	T0 = -1*52,						# lower_temp (C)
	Tm = 37							# higher_temp (c)
)

time_scale = characteristic_scales.get('Lsc')/characteristic_scales.get('Vsc')
characteristic_scales.update(Tsc = time_scale)

# Temporal control
time_control = dict(
				 C_no = 0.2, 					# Maximum possible Courant number
			   	 dt = 0.002,  					# Time-step
			   	 T = 1,							# Total runtime
			   	 adjustable_timestep = False 	# Calculate time-step using max Courant no. during runtime: used to accelerate temporal solution
			   )

# Nature of the problem
problem_physics = dict(
				  solve_temperature = True,		# Enter "True" is you want to solve for temperature
				  solve_FSI = False				# Enter "True" if you want to solve for fluid-structure interaction
				)

# Degree of variables
velocity_degree = 2
pressure_degree = 1
temperature_degree = 2

# Stabilization constants
alpha = Constant(0.85)                   	  	# SUPG/PSPG stabilization constant 
C_cw = Constant(0.7)                       		# Crosswind stabilization constant (quadratic elements: 0.35, for linear elements: 0.7)

# Non-dimensional numbers
def calc_non_dimensional_numbers(g, rho_f, mu, Spht_f, K_f, Vsc, Lsc, T0, Tm, Tsc):

	Re = (Vsc*Lsc)/mu            
	Pr = 27#mu/(K_f/(rho_f*Spht_f)) 
	Ec = (Vsc*Vsc)/(Spht_f*(Tm-T0))
	Fr = Vsc/sqrt(g*Lsc) 

	return Re, Pr, Ec, Fr

# If 2D problem?: Do u want to calculate stream function and vorticity! 
# Note to self: streamfunction is not defined for 3D. 
calc_stream_function = False

# File printing: control no. of time-steps
print_control = dict(
                  a = 40,  # for printing variables and restart files
                  b = 50,  # for post processing data/runtime_courant_no text file
                  c = 200  # for compute_time text file
                )


   
  


