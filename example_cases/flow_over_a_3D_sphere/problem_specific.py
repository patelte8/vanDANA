from dolfin import Expression, UserExpression, CompiledExpression
from .user_parameters import characteristic_scales, time_control
from scipy.interpolate import splrep, splev

import sys, os, cppimport
sys.path.insert(0,  '..')
from utilities.read import read_boundary_conditions

Lsc = characteristic_scales['Lsc']
Vsc = characteristic_scales['Vsc']
Tsc = Lsc/Vsc

# Blood perfusion rate
perf = 0
blood_perfusion = False

if blood_perfusion == True:
	perf = 0.85/60 						# ml/s/gm
	perf *= (0.14*Tsc)					# 14% of total coronary-artery perfusion

def time_varying_bc(tt):

	pass