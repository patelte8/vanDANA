from dolfin import Expression, UserExpression, CompiledExpression, \
					compile_cpp_code, MeshFunction, MPI
from .user_parameters import characteristic_scales, time_control
from .time_oscillating_bc import xyz, code
from scipy.interpolate import splrep, splev

import sys, os, cppimport
sys.path.insert(0,  '..')
from utilities.read import read_boundary_conditions, get_mesh

directory = os.path.dirname(os.path.abspath(__file__)) + '/'
fluid_mesh = get_mesh(MPI.comm_world, directory, "file_f.h5")

Lsc = characteristic_scales['Lsc']
Vsc = characteristic_scales['Vsc']
Tsc = Lsc/Vsc

# Cardiac cycle specifics
t_period = 60/63					# heartrate (sec)            
total_heart_cycles = 40   

# Blood perfusion rate
perf = 0
blood_perfusion = True

if blood_perfusion == True:
	perf = 0.85/60 						# ml/s/gm
	perf *= (0.14*Tsc)					# 14% of total coronary-artery perfusion

# Non-dimensionalizing time-period
t_period /= Tsc

time_control.update(T = total_heart_cycles*t_period)

# Expressions used during runtime
tim = Expression('t', t=0.0, degree=1)
num_cycle = Expression('cycle', cycle=0.0, degree=1)

# ---------------------------------------------------------------------------------  

# Read boundary conditions from csv file      
xdata, ydata = read_boundary_conditions(directory, 'RSPV.csv')
cs_RSPV = splrep(xdata,ydata,per=True)
xdata, ydata = read_boundary_conditions(directory, 'LSPV.csv')
cs_LSPV = splrep(xdata,ydata,per=True)
xdata, ydata = read_boundary_conditions(directory, 'RIPV.csv')
cs_RIPV = splrep(xdata,ydata,per=True)
xdata, ydata = read_boundary_conditions(directory, 'LIPV.csv')
cs_LIPV = splrep(xdata,ydata,per=True)
xdata, ydata = read_boundary_conditions(directory, 'temp.csv')
cs_temperatue = splrep(xdata,ydata,per=True)

param_LSPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 322.5, "Vsc": Vsc, "Tsc": Tsc, "func": cs_LSPV};
param_LIPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 209.9, "Vsc": Vsc, "Tsc": Tsc, "func": cs_LIPV};
param_RSPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 188.35, "Vsc": Vsc, "Tsc": Tsc, "func": cs_RSPV};
param_RIPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 437.6, "Vsc": Vsc, "Tsc": Tsc, "func": cs_RIPV};                     
param_temperature={"time": tim, "func": cs_temperatue, "Tsc": Tsc}

RSPV = RIPV = LSPV = LIPV = CompiledExpression(compile_cpp_code(xyz).Inflow(0, MeshFunction('size_t', fluid_mesh.mesh, 0)), degree = 2)
inflow = [LSPV, LIPV, RSPV, RIPV]

# ---------------------------------------------------------------------------------  

def evaluate_boundary_val(a):

    val = (splev((a['time'].t - a['nm'].cycle*a['period'])*a['Tsc'], a['func'])/a['Area'])/a['Vsc']
    return val

def time_varying_bc(tt):

	tim.t = tt; num_cycle.cycle = int(tt / t_period)
	inflow[0].v = evaluate_boundary_val(param_LSPV); inflow[1].v = evaluate_boundary_val(param_LIPV)
	inflow[2].v = evaluate_boundary_val(param_RSPV); inflow[3].v = evaluate_boundary_val(param_RIPV)

	pass