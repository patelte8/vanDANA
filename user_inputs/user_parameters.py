from dolfin import sqrt, Constant
from ufl import tensors

restart = False									# Restart parameter

# Physics of the problem
# ---------------------------------------------------------------------
problem_physics = dict(
				  solve_temperature = True,		# enter "True" if you want to solve for temperature

				  solve_FSI = True,				# enter "True" if you want to solve for fluid-structure interaction
				  
				  compressible_solid = True,	# enter "True" if compressible: Also remember to specify compressibility (Ld) or possions ratio (nw)
				  								# enter "False" if incompressible

				  solid_material = 'neoHookean',# options: 1. neoHookean 2. linearelastic								

				  viscous_dissipation = False,	# Heat release due to viscous gradients 

				  body_force = False,      		# Gravitational force (uniform volumetric force)								 
				)
	
def f_dir(dim):									# Body force direction : -ve y direction (by default)
	
	n = -1*tensors.unit_vector(1, dim) 
	return n

interpolation_fx = 'phi4'						# Delta-function interpolation for FSI problems

# FEM stabilization and constants
# ---------------------------------------------------------------------
stabilization_parameters = dict(	

	# Navier-stokes
	SUPG_NS = False,							# explicit
	PSPG_NS = False,							# explicit		
	crosswind_NS = False,						# implicit

	# Energy-equation
	SUPG_HT = False,							# explicit
	crosswind_HT = False						# implicit
)

alpha = Constant(0.85)                   	  	# SUPG/PSPG stabilization constant 
C_cw = Constant(0.7)                       		# Crosswind stabilization constant (As per R Codina : quadratic elements: 0.35, for linear elements: 0.7)

# Physical parameters    
# ---------------------------------------------------------------------
physical_parameters = dict(
 
	g = 9.81,									# Gravity (m/s2)						  

	# Fluid 
	rho_f = 1,									# Density (kg/m3)
	nu = 1,										# Dynamic viscosity (kg/m.s)
	Spht_f = 1,         						# Specific heat (J/kg.C)
	K_f = 1,									# Thermal conductivity (W/m.C)

	# Solid
	rho_s = 10,									# Density (kg/m3)
	Sm = 0,										# Shear modulus (N/m2)
	Ld = 0,										# Compressibility (N/m2) ... only for neoHookean
	nw = 0.4,									# Poissons ratio ... only for linearelastic
	Spht_s = 0.11,								# Specific heat (J/kg.C)
	K_s = 1.2 									# Thermal conductivity (W/m.C)
)

def calc_non_dimensional_solid_properties(g, rho_f, nu, Spht_f, K_f, rho_s, Sm, Ld, nw, Spht_s, K_s, Lsc, Vsc, T0, Tm, Tsc):

	rho = rho_s/rho_f
	Spht = Spht_s/Spht_f
	K = K_s/K_f
	Ld = 2000 # Ld/(rho_f*Vsc*Vsc)
	Nw = nw
	Sm = 500 # Sm/(rho_f*Vsc*Vsc)
	
	return rho, Spht, K, Ld, Nw, Sm

# Characteristic scales
# ---------------------------------------------------------------------
characteristic_scales = dict(
	
	Lsc = 1,			            			# m          
	Vsc = 1,	         		    			# m/s
	T0 = -1*52,									# lower_temp (C)
	Tm = 37										# higher_temp (c)
)

# Temporal control
# ---------------------------------------------------------------------
time_control = dict(
				 C_no = 0.35, 					# Maximum possible Courant number
			   	 dt = 0.0025,  					# Time-step: constant throughout runtime if adjustable-timestep is "False"
			   	 T = 100,						# Total runtime
			   	 adjustable_timestep = True 	# Calculate time-step using max Courant no. during runtime: used to accelerate temporal solution
			   )

# FEM degree of variables
# ---------------------------------------------------------------------
fem_degree = dict(
				velocity_degree = 2,
				pressure_degree = 1,
				temperature_degree = 2, 
				displacement_degree = 2,
				lagrange_degree = 1
			   )

# Non-dimensional numbers
# ---------------------------------------------------------------------
def calc_non_dimensional_numbers(g, rho_f, nu, Spht_f, K_f, rho_s, Sm, Ld, nw, Spht_s, K_s, Lsc, Vsc, T0, Tm, Tsc):

	Re = 100#rho_f*(Vsc*Lsc)/nu            
	Pr = 1#(Spht_f*nu)/K_f 
	Ec = 0.4#(Vsc*Vsc)/(Spht_f*(Tm-T0))
	Fr = Vsc/sqrt(g*Lsc) 

	return Re, Pr, Ec, Fr

# Enter "True" if you want to post-process data
# ---------------------------------------------------------------------
post_process = True

# File printing / solid-remeshing control
# ---------------------------------------------------------------------
print_control = dict(
                  a = 40,   					# for printing variables and restart files
                  b = 20,  						# for post processing data
                  c = 20, 						# for simulation_wall_time text file
                  d = 5,   						# for remeshing solid current-configuration mesh		
                  e = 20    					# for runtime_tsp_courant_no_stats text file	
                )

# If 2D problem?: Do u want to calculate stream function and vorticity! # Note to self: streamfunction is not defined for 3D.
# --------------------------------------------------------------------- 
calc_stream_function = True
