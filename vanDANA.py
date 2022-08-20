from dolfin import *
from user_inputs import *
from common import *
from utilities import *
from mshr import *
import numpy as np
import array as arr
import math, os, operator, sys, json
from ufl import Jacobian, nabla_div, tensors
from matplotlib import pylab as plt
from matplotlib import rc

timer_total.start()
memory = MemoryUsage('Start')
restart = boolean_parameters.get('restart')
curr_dir = os.path.dirname(os.path.abspath(__file__)) + '/'
inp_dir = path.join(curr_dir, "user_inputs/")

# MPI-initialize / terminal printing controls
Mpi = MPI_Manage()

blockPrint()
if Mpi.get_rank() == 0: enablePrint() 

# info(parameters, True)

# Read meshes
fluid_mesh = get_mesh(inp_dir, "file.h5")
mesh1 = fluid_mesh.mesh
boundaries1 = fluid_mesh.get_mesh_boundaries()
hmax_f = Mpi.Max(mesh1.hmax()); hmin_f = Mpi.Min(mesh1.hmin())

# Problem dimension
dim = mesh1.geometry().dim()

# Predict initial time-step
if time_control.get('variable_timestep') == True:
    time_control.update(dt = 0.5*time_control.get('C_no')*hmin_f)     
tsp = dt = time_control.get('dt')
T = time_control.get('T')
dt = Constant(dt)

# Create output folder
result_folder = create_result_folder(curr_dir, restart, dim, calc_stream_function)

Mpi.set_barrier()
print("\nMesh specs | edge length: Max =",hmax_f, "; Min =",hmin_f, flush = True)
print(GREEN % "\nReynolds number = {}".format(Re), flush = True)
print(GREEN % "Prandtl number = {}".format(Pr), flush = True)
print(GREEN % "Froude number = {}".format(Fr), \
      BLUE % "; considering body force = {}".format(boolean_parameters.get('body_force')), flush = True)
print(GREEN % "Eckert number = {}".format(Ec), \
      BLUE % "; considering viscous dissipation = {}".format(boolean_parameters.get('viscous_dissipation')), flush = True)
print(RED % "\ntime_step = {}".format(tsp), flush = True)
print(RED % "Total time = {}".format(T), "\n", flush = True)

# Define function spaces
V  = VectorFunctionSpace(mesh1, 'P', velocity_degree)		        # Fluid velocity on mesh1  
Q  = FunctionSpace(mesh1, 'P', pressure_degree)		                # Fluid pressure on mesh1
G  = FunctionSpace(mesh1, 'P', temperature_degree)                  # Fluid temperature on mesh1

FS = dict(fluid = [V, Q, G], solid = [])                   

# Initialize flow problem
flow = fluid_problem(FS['fluid'], boundaries1, result_folder.bool_stream)
u_ = flow.variables['u_']; p_ = flow.variables['p_']; vort = flow.variables['vort']; psi = flow.variables['psi']

# Initialize temperature problem
flow_temp = temperature_problem(FS['fluid'], boundaries1)
T_ = flow_temp.variables['T_'] 

# Initial conditions
create_initial_conditions(u_, p_, T_)

# Boundary conditions
cpp_code = compile_cpp_code(code).Inflow(0, MeshFunction('size_t', mesh1, 0))
RSPV = RIPV = LSPV = LIPV = CompiledExpression(cpp_code, degree = 2)
inflow = [LSPV, LIPV, RSPV, RIPV]
bcs = create_boundary_conditions(boundaries1, inflow, **FS)

# Pre-assemble matrices
flow.pre_assemble(p_[0], FS['fluid'][1], bcs)
flow_temp.pre_assemble()

# Time
t = 0
tim.t = t

counters = create_counters(3)   # enter number of counters required

# Write meshes
pv1 = write_mesh(result_folder.folder, mesh1, "fluid_mesh")
pv1.write_mesh_boundaries(boundaries1)

# Output file variables
files = ['u', 'p', 'T']
text_files = ['data', 'runtime_stats', 'restart', 'log_info']

if result_folder.bool_stream == True:
    files.extend(['vorticity', 'stream_function'])
xdmf_file_handles, hdf5_file_handles = result_folder.create_files(files, Mpi.mpi_comm)
text_file_handles = result_folder.create_text_files(text_files, Mpi.my_rank)

# Restart_variables
restart_write_variables = dict(fluid = [u_[0], u_[1], p_[0], p_[1] , T_[0], T_[1], T_[2], vort, psi], solid = [])
restart_read_variables = dict(fluid = [u_[1], u_[2], p_[1], p_[2], T_[1], T_[2], T_[3], vort, psi], solid = [])

if restart == True:
    t, tsp = read_restart_files(result_folder.folder, Mpi.mpi_comm, text_file_handles[2], **restart_read_variables) 
    dt = Constant(tsp)

# Calculate Total DOF's
DOFS_variables = dict(velocity = [u_[0]], pressure = [p_[0]], temperature = [T_[0]])
DOFS = Calc_total_DOF(Mpi, **DOFS_variables)
print(GREEN % 'DOFs = {}'.format(DOFS), "\n", flush = True)

initial_memory_use = Mpi.Sum(getMemoryUsage())
print(RED % 'Total intitial memory usage for setting up & assembly of the problem = {} MB (RSS)'.format(initial_memory_use), "\n")

# Create progress bar
LogLevel.ERROR; Mpi.set_barrier()
print(BLUE % 'Start Simulatons', "\n", flush = True)

while t < T:
    
    timer_dt.start()
    update_counter(counters)

    # Update current time
    t += tsp   

    # Update boundary conditions
    tim.t = t; num_cycle.cycle = int(t / t_period)     
    inflow[0].v = evaluate_boundary_val(param_LSPV); inflow[1].v = evaluate_boundary_val(param_LIPV)
    inflow[2].v = evaluate_boundary_val(param_RSPV); inflow[3].v = evaluate_boundary_val(param_RIPV)

    timer_s1.start()
    # print(BLUE % "1: Predict tentative velocity step", flush = True)
    A1, b1 = flow.assemble_tentative_velocity(u_, p_, dt)
    flow.solve_tentative_velocity(A1, u_[0], b1, bcs['velocity'])
    s1 += timer_s1.stop()

    timer_s2.start()
    # print(BLUE % "2: Pressure correction step", flush = True)
    b2 = flow.assemble_pressure_correction(u_, p_, dt)
    flow.solve_pressure_correction(p_[0], b2, bcs['pressure'])
    s2 += timer_s2.stop()

    timer_s3.start()
    # print(BLUE % "3: Velocity correction step", flush = True)
    b3 = flow.assemble_velocity_correction(u_, p_, dt)
    flow.solve_velocity_correction(u_[0], b3, bcs['velocity'])
    s3 += timer_s3.stop()

    timer_s4.start()
    # print(BLUE % "Step 7: Energy conservation step", flush = True)
    A4, b4 = flow_temp.assemble_temperature(T_, u_[0], dt)
    flow_temp.solve_temperature(A4, T_[0], b4, bcs['temperature'])
    s4 += timer_s4.stop()	    

    # Print output files
    if counters[0] > print_control['a']:

        reset_counter(counters, 0); Mpi.set_barrier()
        
        print(BLUE % "File printing in progress --- Simulation run time : {} , Wall time elapsed : {} sec".format(t, timer_total.elapsed()[0]), flush = True) 

        vort, psi = flow.calc_vorticity_streamfunction(u_[0], bcs['streamfunction'])

        write_solution_files(t, xdmf_file_handles, hdf5_file_handles, flow.variables, \
                             flow_temp.variables, restart, result_folder.bool_stream)

        # Write restart files    
        write_restart_files(result_folder.folder, Mpi, text_file_handles[2], \
                            t, tsp, **restart_write_variables)

    # Update progress on terminal
    print(' '*92 + "Progress : " + str((t/T)*100) + " %", flush = True)

    # Print post-processing data / calculate new time-step if required 
    if counters[1] > print_control['b']:

        reset_counter(counters, 1)
        
        flow.post_process_data(Mpi, u_, p_, t, tsp, text_file_handles)
        flow_temp.post_process_data(Mpi, T_, t, text_file_handles)

        tsp = calc_runtime_stats_timestep(Mpi, u_[0], t, tsp, text_file_handles, flow.h_f_X, flow.Re, time_control)
        dt  = Constant(tsp)         

    # Update previous solution
    update_variables(u_, p_, T_)

    # Timing tasks
    if counters[2] > print_control['c']:
    	
        reset_counter(counters, 2); Mpi.set_barrier() 
        if Mpi.get_rank() == 0:
            text_file_handles[3].truncate(0); text_file_handles[3].seek(0)
            text_file_handles[3].write(f"{t}    {s1}    {s2}    {s3}    {s4}\n\n")

    s_dt += timer_dt.stop()

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
for y in hdf5_file_handles:
    y.close()

if restart == True:
    merge_xdmf_files(curr_dir, "u", "file.h5" , velocity_degree, "velocity", "v")
    merge_xdmf_files(curr_dir, "p", "file.h5" , pressure_degree, "pressure", "s")
    merge_xdmf_files(curr_dir, "T", "file.h5" , temperature_degree, "temperature", "s")

    if result_folder.bool_stream == True:
        merge_xdmf_files(curr_dir, "vorticity", "file.h5" , pressure_degree, "vorticity", "s")
        merge_xdmf_files(curr_dir, "streamfunction", "file.h5" , pressure_degree, "streamfunction", "s")