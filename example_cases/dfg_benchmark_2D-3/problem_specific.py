from dolfin import Expression, UserExpression, CompiledExpression
from .user_parameters import characteristic_scales, time_control
from scipy.interpolate import splrep, splev

import sys, os, cppimport
sys.path.insert(0,  '..')
from utilities.read import read_boundary_conditions

directory = os.path.dirname(os.path.abspath(__file__)) + '/'

Lsc = characteristic_scales['Lsc']
Vsc = characteristic_scales['Vsc']
Tsc = Lsc/Vsc

# Non-dimensionalizing quantities
b = 0.41/Lsc							# channel width
total_time = time_control['T']/Tsc 		# total run-time
time_control.update(T = total_time)

# Blood perfusion rate
perf = 0
blood_perfusion = False

if blood_perfusion == True:
	perf = 0.85/60 						# ml/s/gm
	perf *= (0.14*Tsc)					# 14% of total coronary-artery perfusion

# ---------------------------------------------------------------------------------  

# Expressions used during runtime
inflow_profile = Expression('6.0*sin(3.14159265*t/T)*x[1]*(b - x[1])/(b*b)', t=0, T=time_control['T'], b=b, degree=2)

def time_varying_bc(tt):

	inflow_profile.t = tt
	pass