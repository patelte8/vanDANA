from dolfin import *
from user_inputs import *
from common import *
from utilities import *
from mshr import *
import numpy as np
import array as arr
from fenicstools import *
from distutils.util import strtobool
import vtk_py3 as vtk_py3
import math, os, operator, copy, sys, io, json, vtk, matplotlib, cppimport, argparse, traceback
matplotlib.use('Agg')
from matplotlib import rc, pylab as plt


def vanDANA_solver(args):
	
	timer_total.start()

	restart = args.restart
	post_process = args.post_process
	print_control.update({"a": args.a, "b": args.b})
	fem_degree.update({"velocity_degree": args.velocity_degree, "displacement_degree": args.displacement_degree, "temperature_degree": args.temperature_degree})
	time_control.update({"T": args.T, "adjustable_timestep": args.adjustable_timestep})
	problem_physics.update({"solve_temperature": args.solve_temperature, "solve_FSI": args.solve_FSI, "solid_material": args.solid_material, "viscous_dissipation": args.viscous_dissipation})

	memory = MemoryUsage('Start')
	curr_dir = os.path.dirname(os.path.abspath(__file__)) + '/'
	remove_killvanDANA(curr_dir); remove_complete(curr_dir)	

	# MPI-initialize / terminal printing controls
	Mpi = MPI_Manage()

	blockPrint()
	if Mpi.get_rank() == 0: enablePrint() 

	# info(parameters, True)

	print(RED % "\nRestart : {}".format(str(restart)), flush = True)
	print(BLUE % "\nSolve scaler (temperature) / transport equation = {}".format(problem_physics['solve_temperature']), flush = True)
	print(BLUE % "Solve fluid-structure interactions = {}".format(problem_physics['solve_FSI']), flush = True)
	if problem_physics['solve_temperature'] == False: stabilization_parameters.update(SUPG_HT = False, crosswind_HT = False)
	print(BLUE % "\nFEM stabilizations = {}".format(str([k for k, v in stabilization_parameters.items() if v == True])), flush = True)

	time_scale = characteristic_scales['Lsc']/characteristic_scales['Vsc']
	characteristic_scales.update(Tsc = time_scale)

	# ---------------------------------------------------------------------------------

	# Calculate non-dimensional numbers
	Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)
	Pe = Re*Pr     
	if not problem_physics['viscous_dissipation'] : Ec = 0.0
	if problem_physics['solve_FSI'] or problem_physics['solve_temperature'] == True:
	    rho, Spht, K, Ld, Nw, Sm = calc_non_dimensional_solid_properties(**physical_parameters, **characteristic_scales)

	# ---------------------------------------------------------------------------------   

	# Read meshes
	mesh_path = path.join(curr_dir, "user_inputs/")
	fluid_mesh = get_mesh(Mpi.mpi_comm, mesh_path, "file_f.h5")
	hmax_f = Mpi.Max(fluid_mesh.mesh.hmax()); hmin_f = Mpi.Min(fluid_mesh.mesh.hmin())

	if problem_physics['solve_FSI'] == True:
	    mesh_path = path.join(curr_dir, "user_inputs/"); mesh_file = "file_s.h5"
	    if restart == True:
	    	mesh_path = path.join(curr_dir, "results/restart_variables/"); mesh_file = "solid_current_mesh.h5"

	    solid_mesh = get_mesh(Mpi.mpi_comm, mesh_path, mesh_file)
	    boundaries = solid_mesh.get_mesh_boundaries(); subdomains = solid_mesh.get_mesh_subdomains()
		
	    hmax_s = Mpi.Max(solid_mesh.mesh.hmax()); hmin_s = Mpi.Min(solid_mesh.mesh.hmin()) 

	# Problem dimension
	dim = fluid_mesh.mesh.geometry().dim()

	# ---------------------------------------------------------------------------------

	# Predict initial time-step
	if time_control['C_no'] > 0.8:	time_control.update(C_no = 0.8)
	if time_control['adjustable_timestep'] == True:
	    initial_time_step = round_decimals_down(0.2*hmin_f, 5)
	    time_control.update(dt = initial_time_step)     
	tsp = dt = time_control['dt']
	T = time_control['T']
	dt = Constant(dt)

	Mpi.set_barrier()
	print("\nFluid mesh specs | edge length: Max =",hmax_f, "; Min =",hmin_f, flush = True)
	print(GREEN % "\nReynolds number = {}".format(Re), flush = True)
	print(GREEN % "Froude number = {}".format(Fr), \
	      BLUE % "; considering body force = {}".format(problem_physics['body_force']), flush = True)
	if problem_physics['solve_temperature'] == True:
		print(GREEN % "Prandtl number = {}".format(Pr), flush = True)
		print(GREEN % "Eckert number = {}".format(Ec), \
		      BLUE % "; considering viscous dissipation = {}".format(problem_physics['viscous_dissipation']), flush = True)

	if problem_physics['solve_FSI'] == True:

	    print("\nSolid mesh specs | edge length: Max =",hmax_s, "; Min =",hmin_s, flush = True)
	    print(BLUE % "\nSolid behavior --- {} ; compressible = {}".format(problem_physics['solid_material'], problem_physics['compressible_solid']), flush = True)
	    print(GREEN % "Shear modulus = {}".format(Sm), flush = True)
	    if problem_physics['compressible_solid'] == True:
		    if problem_physics['solid_material'] == 'linearelastic':	print(GREEN % "Poissons ratio = {}".format(Nw), flush = True)
		    if problem_physics['solid_material'] == 'neoHookean':	print(GREEN % "Compressiblity = {}".format(Ld), flush = True)
	    print(RED % "\nThe following ratios are defined wrt fluid as the reference domain:", flush = True)
	    print(GREEN % "Density ratio = {}".format(rho), flush = True)
	    if problem_physics['solve_temperature'] == True:
		    print(GREEN % "Specific heat ratio = {}".format(Spht), flush = True)
		    print(GREEN % "Conductivity ratio = {}".format(K), flush = True)

	if restart == False: print(RED % "\nInitial time_step = {}".format(tsp), flush = True)
	print(RED % "Total time = {}".format(T), "\n", flush = True)

	# ---------------------------------------------------------------------------------
	                       
	# Create output folder
	result_folder = create_result_folder(curr_dir, restart, dim, calc_stream_function)

	# Initialize flow problem
	flow = Fluid_problem(fluid_mesh, result_folder.bool_stream); FS = dict(fluid = flow.F) 
	u_components = flow.u_components; u_ = flow.variables['u_'];  uv = flow.variables['uv']; assigner_uv = flow.assigner_uv 
	p_ = flow.variables['p_']; Lm_f = flow.variables['Lm_f']; vort = flow.variables['vort']; psi = flow.variables['psi']
	u_inner = flow.u_inner; p_inner = flow.p_inner

	# Initialize temperature problem
	flow_temp = Fluid_temperature_problem(fluid_mesh); FS.update(fluid_temp=flow_temp.F)
	T_ = flow_temp.variables['T_']; LmTf_ = flow_temp.variables['LmTf_']

	# Initialize solid problem
	if problem_physics['solve_FSI'] == True:
	    solid = Solid_problem(solid_mesh); FS.update(solid = solid.F)
	    Dp_ = solid.variables['Dp_']; mix = solid.variables['mix']; us_ = solid.variables['us_']; ps_ = solid.variables['ps_']; J_ = solid.variables['J_']
	    Mv = Function(VectorFunctionSpace(solid_mesh.mesh, 'P', fem_degree['displacement_degree']))

	# Initialize langrange multiplier problem
	if problem_physics['solve_FSI'] == True:
	    lagrange = Lagrange_multiplier_problem(solid_mesh); FS.update(lagrange = lagrange.F)
	    Lm_ = lagrange.variables['Lm_']; uf_ = lagrange.variables['uf_']

	# Initialize solid temperature based langrange multiplier problem
	    solid_temp = Solid_temperature_lagrange_multiplier_problem(solid_mesh); FS.update(solid_temp=solid_temp.F)
	    Ts_ = solid_temp.variables['Ts_']; LmTs_ = solid_temp.variables['LmTs_']

	variables = dict(flow = flow.variables)
	if problem_physics['solve_temperature'] == True: variables.update(flow_temp = flow_temp.variables)
	if problem_physics['solve_FSI'] == True:  variables.update(solid=solid.variables, lagrange=lagrange.variables)
	if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True: 
	    variables.update(solid_temp = solid_temp.variables)

	# ---------------------------------------------------------------------------------    

	# Initial conditions
	if restart == False:
	    fluid_create_initial_conditions(u_, p_, T_)

	    if problem_physics['solve_FSI'] == True:
	        solid_create_initial_conditions(Dp_, mix, dt)

	# Boundary conditions
	# cpp_code = compile_cpp_code(code)
	# RSPV_x = RIPV_x = LSPV_x = LIPV_x = CompiledExpression(cpp_code.Inflow_x(0, MeshFunction('size_t', fluid_mesh.mesh, 0)), degree = 2)
	# RSPV_y = RIPV_y = LSPV_y = LIPV_y = CompiledExpression(cpp_code.Inflow_y(0, MeshFunction('size_t', fluid_mesh.mesh, 0)), degree = 2)
	# RSPV_z = RIPV_z = LSPV_z = LIPV_z = CompiledExpression(cpp_code.Inflow_z(0, MeshFunction('size_t', fluid_mesh.mesh, 0)), degree = 2)
	inflow = [] #dict(x=[LSPV_x, LIPV_x, RSPV_x, RIPV_x], y=[LSPV_y, LIPV_y, RSPV_y, RIPV_y], z=[LSPV_z, LIPV_z, RSPV_z, RIPV_z])	
	bcs = fluid_create_boundary_conditions(fluid_mesh, inflow, **FS)

	if problem_physics['solve_FSI'] == True:
	    bcs.update(solid = solid_create_boundary_conditions(solid_mesh, boundaries, problem_physics['compressible_solid'], dt, **FS))

	# ---------------------------------------------------------------------------------    

	# Delta-interpolation (only required for FSI problems)
	if problem_physics['solve_FSI'] == True:
		fsi_interpolation = compile_cpp_code(fsi_interpolation_code)
		fsi_interpolation.create_bounding_box(solid_mesh.mesh)
		fsi_interpolation.calculate_fluid_mesh_size_h(fluid_mesh.mesh)
		fsi_interpolation.extract_dof_component_map_user(FS['fluid'][2], "F")
		fsi_interpolation.extract_dof_component_map_user(FS['lagrange'][0], "S")
		if problem_physics['solve_temperature'] == True:
		    fsi_interpolation.extract_dof_component_map_user(FS['fluid_temp'][0], "F")
		    fsi_interpolation.extract_dof_component_map_user(FS['solid_temp'][1], "S")

	# ---------------------------------------------------------------------------------        
	        
	# Pre-assemble matrices
	flow.pre_assemble(p_[0], bcs, dt)
	flow_temp.pre_assemble(dt)

	# Time
	t = 0

	piso_tol = 1e-3											# tolerance for PISO loop
	piso_iter = 2 											# no. of PISO interations
	recovering = False; no_consecutive_recovers = 0 		# recovery variables
	counters = create_counters(5)   						# enter number of counters required

	# Timer variables
	s1, s2, s3, s4, s5, s6, s7, si, sm, sr, s_dt = [0.0 for _ in range(11)]

	# --------------------------------------------------------------------------------- 

	# Output/write meshes
	pv1 = write_mesh(result_folder.folder, fluid_mesh.mesh, "fluid_mesh")
	pv1.write_mesh_boundaries(fluid_mesh.get_mesh_boundaries())

	if problem_physics['solve_FSI'] == True:
		timeseries = write_time_series(curr_dir, restart)

		filename = "solid_reference_mesh"
		if restart == True:	filename = "solid_restart_mesh"

		pv = write_mesh(result_folder.folder, solid_mesh.mesh, filename)
		pv.write_mesh_boundaries(boundaries)
		pv.write_mesh_subdomains(subdomains)

	# Output/write files
	files = ['u', 'p', 'T']
	text_files = ['flow_data', 'runtime_stats', 'restart', 'log_info', 'flow_temp_data']
	if problem_physics['solve_FSI'] == True:
	    files.extend(['Dp', 'us', 'ps', 'J', 'Lm'])
	    text_files.extend(['solid_data', 'lagrange_data', 'solid_mesh_quality'])
	    if problem_physics['solve_temperature'] == True:
	        files.extend(['Ts'])
	        text_files.extend(['solid_temp_data'])
	if result_folder.bool_stream == True:
	    files.extend(['vorticity', 'stream_function'])

	xdmf_file_handles, hdf5_file_handles = result_folder.create_files(files, Mpi.mpi_comm)
	text_file_handles = result_folder.create_text_files(text_files, Mpi.my_rank)
	result_folder.write_header_text_files(text_file_handles, Mpi.my_rank) 	

	# --------------------------------------------------------------------------------- 

	# Restart_variables
	restart_write_variables = dict(fluid = []); restart_read_variables = dict(fluid = [])
	for ui in range(u_components):	
		restart_write_variables['fluid'].extend([u_[0][ui]]); restart_read_variables['fluid'].extend([u_[1][ui]])
	for ui in range(u_components):	
		restart_write_variables['fluid'].extend([u_[1][ui]]); restart_read_variables['fluid'].extend([u_[2][ui]]) 
	restart_write_variables['fluid'].extend([p_[0], p_[1], vort, psi]); restart_write_variables.update(fluid_temp = [T_[0], T_[1]])
	restart_read_variables['fluid'].extend([p_[1], p_[2], vort, psi]); restart_read_variables.update(fluid_temp = [T_[1], T_[2]])

	if problem_physics['solve_FSI'] == True:
	    restart_write_variables.update(solid = [Dp_[0], Dp_[1], ps_], lagrange = [Lm_[0]])
	    restart_read_variables.update(solid = [Dp_[0], Dp_[2], ps_], lagrange = [Lm_[1]])       

	if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True: 
	    restart_write_variables.update(solid_temp = [Ts_[0]]); restart_write_variables['lagrange'].extend([LmTs_[0]]) 
	    restart_read_variables.update(solid_temp = [Ts_[1]]); restart_read_variables['lagrange'].extend([LmTs_[1]]) 

	if restart == True:
	    t, tsp = read_restart_files(result_folder.folder, Mpi.mpi_comm, text_file_handles[2], **restart_read_variables) 
	    dt = Constant(tsp)
	    print(RED % "Restart time_step = {}".format(tsp), flush = True)

	# ---------------------------------------------------------------------------------         

	# Calculate Total DOF's solved
	DOFS_variables = dict(velocity = [u_[0][ui] for ui in range(u_components)], pressure = [p_[0]])
	if problem_physics['solve_temperature'] == True: DOFS_variables.update(temperature = [T_[0]])
	if problem_physics['solve_FSI'] == True:
	    DOFS_variables.update(displacement = [Dp_[1]], lagrange_multiplier = [Lm_[0]])
	if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
	    DOFS_variables['lagrange_multiplier'].extend([LmTs_[0]])

	DOFS = Calc_total_DOF(Mpi, **DOFS_variables)
	print(GREEN % 'DOFs = {}'.format(DOFS), "\n", flush = True)

	# --------------------------------------------------------------------------------- 

	# Update temporal variables
	update = [u_, p_, T_]
	if problem_physics['solve_FSI'] == True:
	    update.extend([Dp_, Lm_])
	    if problem_physics['solve_temperature'] == True:
	        update.extend([Ts_, LmTs_])

	# Create progress bar
	LogLevel.ERROR; Mpi.set_barrier()

	initial_memory_use = Mpi.Sum(getMemoryUsage())
	print(RED % 'Total intitial memory usage for setting up & assembly of the problem = {} MB (RSS)'.format(initial_memory_use), "\n", flush = True)
	print(RED % 'Start Simulatons : t = {}'.format(t), "\n", flush = True)

	# ---------------------------------------------------------------------------------

	# Time loop
	try:

		while T > tsp and t < T:

			try: 	# loop for recoveries 

				timer_dt.start()
				inner_iter = 0
				
				if recovering == False:
				    update_counter(counters)

				    # Update current time
				    t += tsp					
				
				if problem_physics['solve_FSI'] == True:
					fsi_interpolation.create_bounding_box(solid_mesh.mesh)

				# ---------------------------------------------------------------------------------

				# Update boundary conditions : only if time-dependent
				# parabolic_profile.t = t; tim.t = t; num_cycle.cycle = int(t / t_period)     
				# for ui, value in inflow.items():     
				   #  inflow[ui][0].v = evaluate_boundary_val(param_LSPV); inflow[ui][1].v = evaluate_boundary_val(param_LIPV)
				   #  inflow[ui][2].v = evaluate_boundary_val(param_RSPV); inflow[ui][3].v = evaluate_boundary_val(param_RIPV)

				# ---------------------------------------------------------------------------------   

				if problem_physics['solve_FSI'] == True:
				    timer_si.start()
				    Lm_f.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, Lm_[1], FS['fluid'][2], interpolation_fx, "F"))
				    si += timer_si.stop()
				    
				timer_s1.start()
				# print(BLUE % "1: Predict tentative velocity step", flush = True)
				A1, b1 = flow.assemble_tentative_velocity(u_, p_, Lm_f, dt)
				try:
					flow.solve_tentative_velocity(A1, u_[0], b1, bcs['velocity'])
				except:
					print(BLUE % "changing initial guess for predicting tentative velocity", flush = True)
					flow.change_initial_guess(u_[0])
					flow.solve_tentative_velocity(A1, u_[0], b1, bcs['velocity'])
				s1 += timer_s1.stop()

				# PISO inner loop
				p_inner.assign(p_[1]); u_diff = 1e8
				while inner_iter < piso_iter:
					if u_diff > -1.0: # piso_tol:

						inner_iter += 1; u_diff = 0.0
						for ui in range(u_components):
							u_inner[ui].assign(u_[0][ui]) 
								
						timer_s2.start()
						# print(BLUE % "2: Pressure correction step", flush = True)
						b2 = flow.assemble_pressure_correction(u_[0], p_inner, dt)
						flow.solve_pressure_correction(p_[0], b2, bcs['pressure'])
						s2 += timer_s2.stop()

						timer_s3.start()
						# print(BLUE % "3: Velocity correction step", flush = True)
						b3 = flow.assemble_velocity_correction(u_[0], p_[0], p_inner, dt)
						flow.solve_velocity_correction(u_[0], b3, bcs['velocity'])
						s3 += timer_s3.stop()

						p_inner.assign(p_[0])
						for ui in range(u_components):
							u_inner[ui].vector().axpy(-1.0, u_[0][ui].vector())
							u_diff += u_inner[ui].vector().norm('l2')
						print("PISO loop {} : velocity error = {:.3e}".format(inner_iter, u_diff), flush = True)

				assigner_uv.assign(uv, [u_[0][ui] for ui in range(u_components)])

				# --------------------------------------------------------------------------------- 

				if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
				    timer_si.start()
				    LmTf_.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, LmTs_[1], FS['fluid_temp'][0], interpolation_fx, "F"))
				    si += timer_si.stop()

				timer_s4.start()
				# print(BLUE % "4: Energy conservation step", flush = True)
				if problem_physics['solve_temperature'] == True:
				    A4, b4 = flow_temp.assemble_temperature(T_, uv, LmTf_, dt)
				    flow_temp.solve_temperature(A4, T_[0], b4, bcs['temperature'])
				s4 += timer_s4.stop()	    

				# --------------------------------------------------------------------------------- 

				if problem_physics['solve_FSI'] == True:
				    timer_si.start()
				    uf_.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, uv, FS['lagrange'][0], interpolation_fx, "S"))
				    si += timer_si.stop()

				if problem_physics['solve_FSI'] == True:
					# Mapping to reference configuration
					timer_sm.start()
					Mv.assign(Dp_[0]); Mv.vector()[:]*=-1
					mapping(solid_mesh.mesh, Mv)
					sm += timer_sm.stop()

					timer_s5.start()
					# print(BLUE % "5: Solid momentum eq. step", flush = True)
					a5 = solid.assemble_solid_problem(problem_physics, Dp_, mix, uf_, Lm_[1], dt)
					try:
						solid.solve_solid_displacement(solid_mesh.mesh, problem_physics['compressible_solid'], a5, Dp_[1], mix, ps_, p_[0], bcs['solid'])
					except:
						print(BLUE % "changing initial guess for solid momentum eq.", flush = True)
						solid.change_initial_guess(Dp_[1], mix)	        		        	
						solid.solve_solid_displacement(solid_mesh.mesh, problem_physics['compressible_solid'], a5, Dp_[1], mix, ps_, p_[0], bcs['solid'])

					Dp_[0].vector().axpy(1.0, Dp_[1].vector())
					# solid.compute_jacobian(J_, Dp_[0])

					us_.vector().zero()
					us_.vector().axpy(1/float(dt), Dp_[1].vector())
					s5 += timer_s5.stop()

					# Mapping to current configuration
					timer_sm.start()
					Mv.vector()[:]*=-1
					mapping(solid_mesh.mesh, Mv)
					sm += timer_sm.stop()
				
				# --------------------------------------------------------------------------------- 

				timer_s6.start()
				# print(BLUE % "6: Lagrange multiplier (fictitious force) step", flush = True)
				if problem_physics['solve_FSI'] == True:
				    a6, b6 = lagrange.assemble_lagrange_multiplier(Lm_, us_, uf_, dt)
				    lagrange.solve_lagrange_multiplier(a6, Lm_[0], b6)
				s6 += timer_s6.stop()    

				# --------------------------------------------------------------------------------- 

				if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
				    timer_si.start()
				    Ts_[0].assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, T_[0], FS['solid_temp'][1], interpolation_fx, "S"))
				    si += timer_si.stop()

				timer_s7.start()
				# print(BLUE % "7: Solid temperature based lagrange multiplier step", flush = True)
				if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
				    a7, b7 = solid_temp.assemble_solid_temperature_lagrange_multiplier(Ts_, uf_, dt)
				    solid_temp.solve_solid_temperature_lagrange_multiplier(a7, LmTs_[0], b7)
				s7 += timer_s7.stop()

				recovering = False; no_consecutive_recovers = 0

			except Exception as e:
				
				recovering = True; no_consecutive_recovers += 1
				print(BLUE % 'error message : ', flush = True); traceback.print_exc(file=sys.stdout) #; print(e, flush = True)
				print(BLUE % "Recovering ... because vanDANA solver diverged --- at time : {} sec , corresponding timestep : {}".format(t, tsp), flush = True)

				tsp /= 2; t -= tsp
				dt = Constant(tsp)

				if no_consecutive_recovers > 5:
					print(BLUE % 'vanDANA solver - TERMINATED : t = {}, due to five consecutive recovers'.format(t), "\n", flush = True)
					break
				else:
					continue		

			else: 
		
				# Update progress on terminal
				print(' '*100 + "Progress : " + str((t/T)*100) + " %", flush = True)
				
				# Print output files
				if counters[0] >= print_control['a']:

					reset_counter(counters, 0); Mpi.set_barrier()
					print(BLUE % "File printing in progress --- Simulation run time : {} , Wall time elapsed : {} sec".format(t, timer_total.elapsed()[0]), flush = True) 

					vort, psi = flow.calc_vorticity_streamfunction(uv, bcs['streamfunction'])

					write_solution_files(restart, problem_physics, result_folder.bool_stream, t, xdmf_file_handles, hdf5_file_handles, **variables)

					# Write restart files    
					write_restart_files(result_folder.folder, Mpi, text_file_handles[2], t, tsp, **restart_write_variables)

					if problem_physics['solve_FSI'] == True:
						pv2 = write_mesh_H5(Mpi.mpi_comm, result_folder.folder, solid_mesh.mesh, "solid_current_mesh")
						pv2.write_mesh_H5_boundaries(boundaries)
						pv2.write_mesh_H5_subdomains(subdomains)
						pv2.hdf.close()
						
						timeseries.store(solid_mesh.mesh, t)

				# --------------------------------------------------------------------------------- 

				# Output post-processing data
				if post_process == True:
					if counters[1] >= print_control['b']:

						reset_counter(counters, 1)
						flow.post_process_data(Mpi, uv, p_[0], t, tsp, text_file_handles)
						if problem_physics['solve_temperature'] == True: 
							flow_temp.post_process_data(Mpi, T_, t, text_file_handles)
						if problem_physics['solve_FSI'] == True:
							solid.post_process_data(Mpi, us_, ps_, Dp_[0], t, text_file_handles)
							lagrange.post_process_data(Mpi, Lm_[0], uf_, t, dt, text_file_handles)
							if problem_physics['solve_temperature'] == True:	
								solid_temp.post_process_data(Mpi, Ts_[0], t, text_file_handles)

				# If required: calculate new time-step      
				if counters[4] >= print_control['e']:    

					reset_counter(counters, 4)
					tsp = calc_runtime_stats_timestep(Mpi, u_[0], u_components, u_diff, t, tsp, text_file_handles, fluid_mesh.mesh, hmin_f, flow.h_f_X, flow.Re, time_control)
					dt  = Constant(tsp)         

				# ---------------------------------------------------------------------------------     

				# Update previous solution
				update_variables(update, u_components, problem_physics)

				# Move mesh by delta D
				timer_sm.start()
				if problem_physics['solve_FSI'] == True:		
					ALE.move(solid_mesh.mesh, Dp_[1])
					solid_mesh.mesh.bounding_box_tree().build(solid_mesh.mesh)
				sm += timer_sm.stop()

				# Remeshing solid current-congifuration mesh
				timer_sr.start()
				if problem_physics['solve_FSI'] == True:
					if counters[3] >= print_control['d']:

						reset_counter(counters, 3)
						print(GREEN % "Remeshing solid current-congifuration mesh", flush = True)
						solid_mesh.mesh, ratio_min, ratio_max = mesh_smoothening(solid_mesh.mesh)   
						solid_mesh.mesh.bounding_box_tree().build(solid_mesh.mesh)
						if Mpi.get_rank() == 0:
							text_file_handles[7].write(f"{t:0,.10G}		{ratio_min:0,.10G}		{ratio_max:0,.10G}\n")
				sr += timer_sr.stop()

			finally: 

				# Timing tasks
				if counters[2] >= print_control['c']:
					
					reset_counter(counters, 2); Mpi.set_barrier() 
					if Mpi.get_rank() == 0:
						text_file_handles[3].truncate(0); text_file_handles[3].seek(0)
						text_file_handles[3].write("#Time		#Step_1			#Step_2			#Step_3			#Step_4			#Step_5			#Step_6			#Step_7			#Step_interpolation	#Step_move_mesh		#Step_remeshing\n")
						text_file_handles[3].write(f"{t:0,.10G}		{s1:0,.10G}		{s2:0,.10G}		{s3:0,.10G}		{s4:0,.10G}		{s5:0,.10G}		{s6:0,.10G}		{s7:0,.10G}		{si:0,.10G}		{sm:0,.10G}		{sr:0,.10G}\n\n")

				s_dt += timer_dt.stop()

				# Check for killvanDANA file
				if os.path.isfile(path.join(curr_dir, "results/killvanDANA")) == True:	
					print(RED % "--- killing vanDANA solver --- t = {}".format(t), "\n", flush = True)
					break

	# ---------------------------------------------------------------------------------     

	except Exception as e: print(BLUE % 'error message : ', flush = True); traceback.print_exc(file=sys.stdout) #; print(e, flush = True)		

	finally:

		if t >= T and Mpi.get_rank() == 0:
			print(BLUE % 'vanDANA solver - COMPLETED : t = {}'.format(t), "\n", flush = True)
			complete = io.TextIOWrapper(open(curr_dir + "results/complete", "wb", 0), write_through=True)
			complete.seek(0); complete.write("{}, T = {}".format("COMPLETED", T))
			complete.close()
		
		memory('Final memory use')
		print(RED % 'Total memory usage of solver = {} MB (RSS)'.format(str(memory.memory - initial_memory_use)), "\n", flush = True)    
		wall_time = timer_total.stop()

		Mpi.set_barrier() 
		if Mpi.get_rank() == 0: 
		    text_file_handles[3].write("{} {} {}".format("\n", "DOFs -->", json.dumps(DOFS)))
		    text_file_handles[3].write("{} {} {}".format("\n\n", "Total number of tasks : ", Mpi.size))
		    text_file_handles[3].write("{} {} {} {}".format("\n\n", "Total simulation wall time : ", wall_time, " sec"))
		    text_file_handles[3].write("{} {} {} {}".format("\n\n", "Total intitial memory usage for setting up & assembly of the problem : ", initial_memory_use, "MB (RSS)"))
		    text_file_handles[3].write("{} {} {} {} {}".format("\n\n", "Total memory usage of solver : ", str(memory.memory - initial_memory_use), "MB (RSS)", "\n\n\n"))
		    text_file_handles[3].write(timings(TimingClear.keep, [TimingType.wall]).str(True))

		print(RED % "Total simulation wall time : {} sec".format(wall_time), "\n", flush = True)

		for x in text_file_handles:
		    x.close()
		for y,z in hdf5_file_handles.items():
		    z.close(); del z	        

		if problem_physics['solve_FSI'] == True: 
			read_time_series(Mpi.mpi_comm, curr_dir)

		if restart == True:
			extract_hdf5_data_for_xdmf_visualization(Mpi.mpi_comm, curr_dir, result_folder.bool_stream, problem_physics, fem_degree)

	# --------------------------------------------------------------------------------- 
	
	
if __name__ == '__main__':
	
	# parsing arguments from command line
	parser = argparse.ArgumentParser(description = 'to append arguments from terminal')

	parser.add_argument('-restart', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=restart)
	parser.add_argument('-solve_temperature', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=problem_physics["solve_temperature"])
	parser.add_argument('-solve_FSI', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=problem_physics["solve_FSI"])
	parser.add_argument('-solid_material', type=str, metavar='', required=False, default=problem_physics["solid_material"])
	parser.add_argument('-viscous_dissipation', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=problem_physics["viscous_dissipation"])
	parser.add_argument('-velocity_degree', type=int, metavar='', required=False, default=fem_degree["velocity_degree"])
	parser.add_argument('-displacement_degree', type=int, metavar='', required=False, default=fem_degree["displacement_degree"])
	parser.add_argument('-temperature_degree', type=int, metavar='', required=False, default=fem_degree["temperature_degree"])
	parser.add_argument('-adjustable_timestep', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=time_control["adjustable_timestep"])
	parser.add_argument('-post_process', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=post_process)
	parser.add_argument('-T', type=float, metavar='', required=False, default=time_control["T"])
	parser.add_argument('-a', type=int, metavar='', required=False, default=print_control["a"])
	parser.add_argument('-b', type=int, metavar='', required=False, default=print_control["b"])

	# arguments are stored in "args"
	args = parser.parse_args()				

	# ----------------------------------------------------------------------------------

	vanDANA_solver(args)
