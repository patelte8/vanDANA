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
import math, os, operator, copy, sys, json, vtk, matplotlib, cppimport, argparse
matplotlib.use('Agg')
from matplotlib import rc, pylab as plt


def vanDANA_solver(args):
	
	restart = args.restart
	calc_stream_function = args.calc_stream_function
	print_control.update({"a": args.a})
	fem_degree.update({"velocity_degree": args.velocity_degree, "displacement_degree": args.displacement_degree, "temperature_degree": args.temperature_degree})
	time_control.update({"T": args.T, "adjustable_timestep": args.adjustable_timestep})
	problem_physics.update({"solve_temperature": args.solve_temperature, "solve_FSI": args.solve_FSI, "viscous_dissipation": args.viscous_dissipation})

	timer_total.start()
	memory = MemoryUsage('Start')
	curr_dir = os.path.dirname(os.path.abspath(__file__)) + '/'

	# MPI-initialize / terminal printing controls
	Mpi = MPI_Manage()

	blockPrint()
	if Mpi.get_rank() == 0: enablePrint() 

	# info(parameters, True)

	print(RED % "\nRestart : {}".format(str(restart)), flush = True)
	print(BLUE % "\nSolve scaler (temperature) / transport equation = {}".format(problem_physics.get('solve_temperature')), flush = True)
	print(BLUE % "Solve fluid-structure interactions = {}".format(problem_physics.get('solve_FSI')), flush = True)

	time_scale = characteristic_scales.get('Lsc')/characteristic_scales.get('Vsc')
	characteristic_scales.update(Tsc = time_scale)

	# ---------------------------------------------------------------------------------

	# Calculate non-dimensional numbers
	Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)
	Pe = Re*Pr     
	if not problem_physics.get('viscous_dissipation') : Ec = 0.0
	if problem_physics.get('solve_FSI') or problem_physics.get('solve_temperature') == True:
	    rho, Spht, K, Ld, Sm = calc_non_dimensional_solid_properties(**physical_parameters, **characteristic_scales)

	# ---------------------------------------------------------------------------------   

	# Read meshes
	mesh_path = path.join(curr_dir, "user_inputs/")
	fluid_mesh = get_mesh(mesh_path, "file_f.h5")
	hmax_f = Mpi.Max(fluid_mesh.mesh.hmax()); hmin_f = Mpi.Min(fluid_mesh.mesh.hmin())

	if problem_physics.get('solve_FSI') == True:
	    mesh_path = path.join(curr_dir, "user_inputs/"); mesh_file = "file_s.h5"

	    solid_mesh = get_mesh(mesh_path, mesh_file)
	    solid_mesh_R = get_mesh(mesh_path, mesh_file)
	    hmax_s = Mpi.Max(solid_mesh_R.mesh.hmax()); hmin_s = Mpi.Min(solid_mesh_R.mesh.hmin()) 

	# Problem dimension
	dim = fluid_mesh.mesh.geometry().dim()

	# ---------------------------------------------------------------------------------

	# Predict initial time-step
	if time_control.get('adjustable_timestep') == True:
	    initial_time_step = round_decimals_down(0.5*time_control.get('C_no')*hmin_f, 5)
	    time_control.update(dt = initial_time_step)     
	tsp = dt = time_control.get('dt')
	T = time_control.get('T')
	dt = Constant(dt)

	Mpi.set_barrier()
	print("\nFluid mesh specs | edge length: Max =",hmax_f, "; Min =",hmin_f, flush = True)
	print(GREEN % "\nReynolds number = {}".format(Re), flush = True)
	print(GREEN % "Prandtl number = {}".format(Pr), flush = True)
	print(GREEN % "Froude number = {}".format(Fr), \
	      BLUE % "; considering body force = {}".format(problem_physics.get('body_force')), flush = True)
	print(GREEN % "Eckert number = {}".format(Ec), \
	      BLUE % "; considering viscous dissipation = {}".format(problem_physics.get('viscous_dissipation')), flush = True)

	if problem_physics.get('solve_FSI') == True:

	    print("\nSolid mesh specs | edge length: Max =",hmax_s, "; Min =",hmin_s, flush = True)
	    print(GREEN % "\nStiffness = {}".format(Sm), flush = True)
	    print(GREEN % "Compressiblity = {}".format(Ld), \
	      BLUE % "; considering compressible solid = {}".format(problem_physics.get('compressible_solid')), "\n", flush = True)
	    print(RED % "The following ratios are defined wrt fluid as the reference domain:", flush = True)
	    print(GREEN % "Density ratio = {}".format(rho), flush = True)
	    print(GREEN % "Specific heat ratio = {}".format(Spht), flush = True)
	    print(GREEN % "Conductivity ratio = {}".format(K), flush = True)

	if restart == False: print(RED % "\nInitial time_step = {}".format(tsp), flush = True)
	print(RED % "Total time = {}".format(T), "\n", flush = True)

	# ---------------------------------------------------------------------------------
	                       
	# Create output folder
	result_folder = create_result_folder(curr_dir, restart, dim, calc_stream_function)

	# Initialize flow problem
	flow = Fluid_problem(fluid_mesh, result_folder.bool_stream); FS = dict(fluid = flow.F) 
	u_ = flow.variables['u_']; p_ = flow.variables['p_']; Lm_f = flow.variables['Lm_f']
	vort = flow.variables['vort']; psi = flow.variables['psi']

	# Initialize temperature problem
	flow_temp = Fluid_temperature_problem(fluid_mesh); FS.update(fluid_temp=flow_temp.F)
	T_ = flow_temp.variables['T_']; LmTf_ = flow_temp.variables['LmTf_']

	# Initialize solid problem
	if problem_physics.get('solve_FSI') == True:
	    solid = Solid_problem(solid_mesh_R); FS.update(solid = solid.F)
	    Dp_ = solid.variables['Dp_']; mix = solid.variables['mix']; us_ = solid.variables['us_']; ps_ = solid.variables['ps_']; J_ = solid.variables['J_']
	    Mv = Function(VectorFunctionSpace(solid_mesh.mesh, 'P', fem_degree.get('displacement_degree')))

	# Initialize langrange multiplier problem
	if problem_physics.get('solve_FSI') == True:
	    lagrange = Lagrange_multiplier_problem(solid_mesh); FS.update(lagrange = lagrange.F)
	    Lm_ = lagrange.variables['Lm_']; uf_ = lagrange.variables['uf_']

	# Initialize solid temperature based langrange multiplier problem
	    solid_temp = Solid_temperature_lagrange_multiplier_problem(solid_mesh); FS.update(solid_temp=solid_temp.F)
	    Ts_ = solid_temp.variables['Ts_']; LmTs_ = solid_temp.variables['LmTs_']

	variables = dict(flow = flow.variables)
	if problem_physics.get('solve_temperature') == True: variables.update(flow_temp = flow_temp.variables)
	if problem_physics.get('solve_FSI') == True:  variables.update(solid=solid.variables, lagrange=lagrange.variables)
	if problem_physics.get('solve_FSI') and problem_physics.get('solve_temperature') == True: 
	    variables.update(solid_temp = solid_temp.variables)

	# ---------------------------------------------------------------------------------    

	# Initial conditions
	if restart == False:
	    fluid_create_initial_conditions(u_, p_, T_)

	    if problem_physics.get('solve_FSI') == True:
	        solid_create_initial_conditions(Dp_, mix, dt)

	# Boundary conditions
	# cpp_code = compile_cpp_code(code).Inflow(0, MeshFunction('size_t', fluid_mesh.mesh, 0))
	# RSPV = RIPV = LSPV = LIPV = CompiledExpression(cpp_code, degree = 2)
	# inflow = [LSPV, LIPV, RSPV, RIPV]
	bcs = fluid_create_boundary_conditions(fluid_mesh, **FS)

	if problem_physics.get('solve_FSI') == True:
	    bcs.update(solid = solid_create_boundary_conditions(solid_mesh_R, problem_physics.get('compressible_solid'), dt, **FS))

	# ---------------------------------------------------------------------------------    

	# Delta-interpolation (only required for FSI problems)
	if problem_physics.get('solve_FSI') == True:
		fsi_interpolation = compile_cpp_code(fsi_interpolation_code)
		fsi_interpolation.extract_dof_component_map_user(FS['fluid'][2], "F")
		fsi_interpolation.extract_dof_component_map_user(FS['lagrange'][0], "S")
		if problem_physics.get('solve_temperature') == True:
		    fsi_interpolation.extract_dof_component_map_user(FS['fluid_temp'][0], "F")
		    fsi_interpolation.extract_dof_component_map_user(FS['solid_temp'][1], "S")

	# ---------------------------------------------------------------------------------        
	        
	# Pre-assemble matrices
	flow.pre_assemble(p_[0], bcs, dt)
	flow_temp.pre_assemble(dt)

	# Time
	t = 0
	# tim.t = t

	counters = create_counters(4)   # enter number of counters required

	# Timer variables
	s1, s2, s3, s4, s5, s6, s7, si, sm, sr, s_dt = [0.0 for _ in range(11)]

	# --------------------------------------------------------------------------------- 

	# Output/write meshes
	pv1 = write_mesh(result_folder.folder, fluid_mesh.mesh, "fluid_mesh")
	pv1.write_mesh_boundaries(fluid_mesh.get_mesh_boundaries())

	if problem_physics.get('solve_FSI') == True:
	    pvR = write_mesh(result_folder.folder, solid_mesh_R.mesh, "solid_reference_mesh")
	    pvR.write_mesh_boundaries(solid_mesh_R.get_mesh_boundaries())
	    pvR.write_mesh_subdomains(solid_mesh_R.get_mesh_subdomains())

	# Output/write files
	files = ['u', 'p', 'T']
	text_files = ['flow_data', 'runtime_stats', 'restart', 'log_info', 'flow_temp_data']
	if problem_physics.get('solve_FSI') == True:
	    files.extend(['Dp', 'us', 'ps', 'J', 'Lm'])
	    text_files.extend(['solid_data', 'lagrange_data', 'solid_mesh_quality'])
	    if problem_physics.get('solve_temperature') == True:
	        files.extend(['Ts'])
	        text_files.extend(['solid_temp_data'])
	if result_folder.bool_stream == True:
	    files.extend(['vorticity', 'stream_function'])

	xdmf_file_handles, hdf5_file_handles = result_folder.create_files(files, Mpi.mpi_comm)
	text_file_handles = result_folder.create_text_files(text_files, Mpi.my_rank)

	# --------------------------------------------------------------------------------- 

	# Restart_variables
	restart_write_variables = dict(fluid = [u_[0], u_[1], p_[0], p_[1], vort, psi], fluid_temp = [T_[0], T_[1]]) 
	restart_read_variables = dict(fluid = [u_[1], u_[2], p_[1], p_[2], vort, psi], fluid_temp = [T_[1], T_[2]])

	if problem_physics.get('solve_FSI') == True:
	    restart_write_variables.update(solid = [Dp_[0], Dp_[1], ps_], lagrange = [Lm_[0]])
	    restart_read_variables.update(solid = [Dp_[0], Dp_[2], ps_], lagrange = [Lm_[1]])        

	if problem_physics.get('solve_FSI') and problem_physics.get('solve_temperature') == True: 
	    restart_write_variables.update(solid_temp = [Ts_[0]]); restart_write_variables['lagrange'].extend([LmTs_[0]]) 
	    restart_read_variables.update(solid_temp = [Ts_[1]]); restart_write_variables['lagrange'].extend([LmTs_[1]]) 

	if restart == True:
	    t, tsp = read_restart_files(result_folder.folder, Mpi.mpi_comm, text_file_handles[2], **restart_read_variables) 
	    dt = Constant(tsp)
	    print(RED % "Restart time_step = {}".format(tsp), flush = True)

	# Move solid_mesh to restart position before starting simulation
	if restart and problem_physics.get('solve_FSI') == True:
	    print(RED % "Translate initial mesh to restart position\n", flush = True)
	    Mv.vector()[:] = np.zeros(len(Mv.vector().get_local()))
	    Mv.vector()[:] = Dp_[0].vector().get_local()[:]
	    ALE.move(solid_mesh.mesh, project(Mv, VectorFunctionSpace(solid_mesh.mesh, 'P', 1)))
	    solid_mesh.mesh.bounding_box_tree().build(solid_mesh.mesh)
	    lagrange.dx = dolfin.dx(solid_mesh.mesh); lagrange.ds = dolfin.ds(solid_mesh.mesh)
	    if problem_physics.get('solve_temperature') == True:
	        solid_temp.dx = dolfin.dx(solid_mesh.mesh); solid_temp.ds = dolfin.ds(solid_mesh.mesh)

	# ---------------------------------------------------------------------------------         

	# Calculate Total DOF's solved
	DOFS_variables = dict(velocity = [u_[0]], pressure = [p_[0]])
	if problem_physics.get('solve_temperature') == True: DOFS_variables.update(temperature = [T_[0]])
	if problem_physics.get('solve_FSI') == True:
	    DOFS_variables.update(displacement = [Dp_[1]], lagrange_multiplier = [Lm_[0]])
	if problem_physics.get('solve_FSI') and problem_physics.get('solve_temperature') == True:
	        DOFS_variables['lagrange_multiplier'].extend([LmTs_[0]])

	DOFS = Calc_total_DOF(Mpi, **DOFS_variables)
	print(GREEN % 'DOFs = {}'.format(DOFS), "\n", flush = True)

	# --------------------------------------------------------------------------------- 

	# Update temporal variables
	update = [u_, p_, T_]
	if problem_physics.get('solve_FSI') == True:
	    update.extend([Dp_, Lm_])
	    if problem_physics.get('solve_temperature') == True:
	        update.extend([Ts_, LmTs_])

	# Create progress bar
	LogLevel.ERROR; Mpi.set_barrier()

	initial_memory_use = Mpi.Sum(getMemoryUsage())
	print(RED % 'Total intitial memory usage for setting up & assembly of the problem = {} MB (RSS)'.format(initial_memory_use), "\n")
	print(RED % 'Start Simulatons', "\n", flush = True)

	# Time loop
	while T > tsp and t < T:
	    
	    timer_dt.start()
	    update_counter(counters)

	    # Update current time
	    t += tsp   

	    # Update boundary conditions
	    # tim.t = t; num_cycle.cycle = int(t / t_period)     
	    # inflow[0].v = evaluate_boundary_val(param_LSPV); inflow[1].v = evaluate_boundary_val(param_LIPV)
	    # inflow[2].v = evaluate_boundary_val(param_RSPV); inflow[3].v = evaluate_boundary_val(param_RIPV)

	    if problem_physics.get('solve_FSI') == True:
	        timer_si.start()
	        Lm_f.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, Lm_[1], FS['fluid'][2], interpolation_fx, "F"))
	        si += timer_si.stop()
	        
	    timer_s1.start()
	    # print(BLUE % "1: Predict tentative velocity step", flush = True)
	    A1, b1 = flow.assemble_tentative_velocity(u_, p_, Lm_f, dt)
	    flow.solve_tentative_velocity(A1, u_[0], b1, bcs['velocity'])
	    s1 += timer_s1.stop()

	    timer_s2.start()
	    # print(BLUE % "2: Pressure correction step", flush = True)
	    b2 = flow.assemble_pressure_correction(u_, p_, Lm_f, dt)
	    flow.solve_pressure_correction(p_[0], b2, bcs['pressure'])
	    s2 += timer_s2.stop()

	    timer_s3.start()
	    # print(BLUE % "3: Velocity correction step", flush = True)
	    b3 = flow.assemble_velocity_correction(u_, p_, dt)
	    flow.solve_velocity_correction(u_[0], b3, bcs['velocity'])
	    s3 += timer_s3.stop()

	    # --------------------------------------------------------------------------------- 

	    if problem_physics.get('solve_FSI') and problem_physics.get('solve_temperature') == True:
	        timer_si.start()
	        LmTf_.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, LmTs_[1], FS['fluid_temp'][0], interpolation_fx, "F"))
	        si += timer_si.stop()

	    timer_s4.start()
	    # print(BLUE % "4: Energy conservation step", flush = True)
	    if problem_physics.get('solve_temperature') == True:
	        A4, b4 = flow_temp.assemble_temperature(T_, u_[0], LmTf_, dt)
	        flow_temp.solve_temperature(A4, T_[0], b4, bcs['temperature'])
	    s4 += timer_s4.stop()	    

	    # --------------------------------------------------------------------------------- 

	    if problem_physics.get('solve_FSI') == True:
	        timer_si.start()
	        uf_.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, u_[0], FS['lagrange'][0], interpolation_fx, "S"))
	        si += timer_si.stop()

	    timer_s5.start()    
	    # print(BLUE % "5: Solid momentum eq. step", flush = True)    
	    if problem_physics.get('solve_FSI') == True:    
	        a5 = solid.assemble_solid_problem(problem_physics.get('compressible_solid'), Dp_, mix, uf_, Lm_[1], dt)
	        solid.solve_solid_displacement(problem_physics.get('compressible_solid'), a5, Dp_[1], mix, ps_, p_[0], bcs['solid'])

	        Dp_[0].vector().axpy(1.0, Dp_[1].vector())
	        solid.compute_jacobian(J_, Dp_[0])

	        us_.vector().zero()
	        us_.vector().axpy(1/float(dt), Dp_[1].vector())
	    s5 += timer_s5.stop()
	    
	    # --------------------------------------------------------------------------------- 

	    timer_s6.start()
	    # print(BLUE % "6: Lagrange multiplier (fictitious force) step", flush = True)
	    if problem_physics.get('solve_FSI') == True:
	        a6, b6 = lagrange.assemble_lagrange_multiplier(Lm_, us_, uf_, dt)
	        lagrange.solve_lagrange_multiplier(a6, Lm_[0], b6)
	    s6 += timer_s6.stop()    

	    # --------------------------------------------------------------------------------- 

	    if problem_physics.get('solve_FSI') and problem_physics.get('solve_temperature') == True:
	        timer_si.start()
	        Ts_[0].assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, T_[0], FS['solid_temp'][1], interpolation_fx, "S"))
	        si += timer_si.stop()

	    timer_s7.start()
	    # print(BLUE % "7: Solid temperature based lagrange multiplier step", flush = True)
	    if problem_physics.get('solve_FSI') and problem_physics.get('solve_temperature') == True:
	        a7, b7 = solid_temp.assemble_solid_temperature_lagrange_multiplier(Ts_, uf_, dt)
	        solid_temp.solve_solid_temperature_lagrange_multiplier(a7, LmTs_[0], b7)
	    s7 += timer_s7.stop()

	    # --------------------------------------------------------------------------------- 

	    # Print output files
	    if counters[0] > print_control['a']:

	        reset_counter(counters, 0); Mpi.set_barrier()
	        print(BLUE % "File printing in progress --- Simulation run time : {} , Wall time elapsed : {} sec".format(t, timer_total.elapsed()[0]), flush = True) 

	        vort, psi = flow.calc_vorticity_streamfunction(u_[0], bcs['streamfunction'])

	        write_solution_files(restart, problem_physics, result_folder.bool_stream, t, xdmf_file_handles, hdf5_file_handles, **variables)

	        if problem_physics.get('solve_FSI') == True:
	            pv2 = write_mesh(result_folder.folder, solid_mesh.mesh, "solid_current_mesh")
	            pv2.write_mesh_boundaries(solid_mesh.get_mesh_boundaries())
	            pv2.write_mesh_subdomains(solid_mesh.get_mesh_subdomains())

	        # Write restart files    
	        write_restart_files(result_folder.folder, Mpi, text_file_handles[2], t, tsp, **restart_write_variables)

	    # Update progress on terminal
	    print(' '*100 + "Progress : " + str((t/T)*100) + " %", flush = True)

	    # --------------------------------------------------------------------------------- 

	    # Print post-processing data / calculate new time-step if required 
	    if counters[1] > print_control['b']:

	        reset_counter(counters, 1)
	        # flow.post_process_data(Mpi, u_, p_, t, tsp, text_file_handles)
	        # flow_temp.post_process_data(Mpi, T_, t, text_file_handles)
	        solid.post_process_data(Mpi, us_, ps_, t, text_file_handles)
	        lagrange.post_process_data(Mpi, Lm_[0], t, text_file_handles)
	        solid_temp.post_process_data(Mpi, Ts_[0], t, text_file_handles)

	        tsp = calc_runtime_stats_timestep(Mpi, u_[0], t, tsp, text_file_handles, flow.h_f_X, flow.Re, time_control)
	        dt  = Constant(tsp)         

	    # ---------------------------------------------------------------------------------     

	    # Update previous solution
	    update_variables(update, problem_physics)

	    # Move mesh
	    timer_sm.start()
	    if problem_physics.get('solve_FSI') == True:
	        
	        Mv.vector()[:] = np.zeros(len(Mv.vector().get_local()))
	        Mv.vector()[:] = Dp_[1].vector().get_local()[:]
	        ALE.move(solid_mesh.mesh, project(Mv, VectorFunctionSpace(solid_mesh.mesh, 'P', 1)))
	        solid_mesh.mesh.bounding_box_tree().build(solid_mesh.mesh)
	        lagrange.dx = dolfin.dx(solid_mesh.mesh); lagrange.ds = dolfin.ds(solid_mesh.mesh)
	        if problem_physics.get('solve_temperature') == True:
	            solid_temp.dx = dolfin.dx(solid_mesh.mesh); solid_temp.ds = dolfin.ds(solid_mesh.mesh)
	    sm += timer_sm.stop()

	    # Remeshing solid current-congifuration mesh
	    timer_sr.start()
	    if problem_physics.get('solve_FSI') == True:
	        if counters[3] > print_control['d']:

	            reset_counter(counters, 3)
	            print(GREEN % "Remeshing solid current-congifuration mesh", flush = True)
	            solid_mesh.mesh, ratio_min, ratio_max = mesh_smoothening(solid_mesh.mesh)   
	            solid_mesh.mesh.bounding_box_tree().build(solid_mesh.mesh)
	            lagrange.dx = dolfin.dx(solid_mesh.mesh); lagrange.ds = dolfin.ds(solid_mesh.mesh)
	            if problem_physics.get('solve_temperature') == True:
	                solid_temp.dx = dolfin.dx(solid_mesh.mesh); solid_temp.ds = dolfin.ds(solid_mesh.mesh)
	            if Mpi.get_rank() == 0:
	                text_file_handles[7].write("{} {} {} {} {} {}".format(t, "  ", ratio_min, "  ", ratio_max, "\n"))
	    sr += timer_sr.stop()

	    # --------------------------------------------------------------------------------- 

	    # Timing tasks
	    if counters[2] > print_control['c']:
	    	
	        reset_counter(counters, 2); Mpi.set_barrier() 
	        if Mpi.get_rank() == 0:
	            text_file_handles[3].truncate(0); text_file_handles[3].seek(0)
	            text_file_handles[3].write(f"{t}    {s1}    {s2}    {s3}    {s4}    {s5}    {s6}    {s7}    {si}    {sm}    {sr}\n\n")
	    s_dt += timer_dt.stop()

	# ---------------------------------------------------------------------------------     

	memory('Final memory use ')
	print(RED % 'Total memory usage of solver = {} MB (RSS)'.format(str(memory.memory - initial_memory_use)), "\n", flush = True)    
	wall_time = timer_total.stop()

	Mpi.set_barrier() 
	if Mpi.get_rank() == 0: 
	    text_file_handles[3].write(timings(TimingClear.keep, [TimingType.wall]).str(True))
	    text_file_handles[3].write("{} {}".format("\n\n", json.dumps(DOFS)))
	    text_file_handles[3].write("{} {} {} {}".format("\n\n", "Total simulation wall time : ", wall_time, " sec"))

	print(BLUE % "Total simulation wall time : {} sec".format(wall_time), "\n", flush = True)

	for x in text_file_handles:
	    x.close()
	for y,z in hdf5_file_handles.items():
	    z.close(); del z

	# ---------------------------------------------------------------------------------     

	if restart == True:
		extract_hdf5_data_for_xdmf_visualization(Mpi.mpi_comm, curr_dir, result_folder.bool_stream, problem_physics, fem_degree)



if __name__ == '__main__':
	
	# parsing arguments from command line
	parser = argparse.ArgumentParser(description = 'to append arguments from terminal')

	parser.add_argument('-restart', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=restart)
	parser.add_argument('-solve_temperature', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=problem_physics["solve_temperature"])
	parser.add_argument('-solve_FSI', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=problem_physics["solve_FSI"])
	parser.add_argument('-viscous_dissipation', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=problem_physics["viscous_dissipation"])
	parser.add_argument('-velocity_degree', type=int, metavar='', required=False, default=fem_degree["velocity_degree"])
	parser.add_argument('-displacement_degree', type=int, metavar='', required=False, default=fem_degree["displacement_degree"])
	parser.add_argument('-temperature_degree', type=int, metavar='', required=False, default=fem_degree["temperature_degree"])
	parser.add_argument('-adjustable_timestep', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=time_control["adjustable_timestep"])
	parser.add_argument('-calc_stream_function', type=lambda x:bool(strtobool(x)), metavar='', required=False, default=calc_stream_function)
	parser.add_argument('-T', type=float, metavar='', required=False, default=time_control["T"])
	parser.add_argument('-a', type=int, metavar='', required=False, default=print_control["a"])

	# arguments are stored in "args"
	args = parser.parse_args()				

	# ----------------------------------------------------------------------------------

	vanDANA_solver(args)
