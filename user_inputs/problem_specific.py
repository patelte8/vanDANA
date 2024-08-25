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

# ---------------------------------------------------------------------------------  

# Expressions used during runtime
parabolic_profile = Expression('6.0*x[1]*(4.1 - x[1])/(4.1*4.1)', degree=2)

def time_varying_bc(tt):

	pass

# ---------------------------------------------------------------------------------  

class Shear_modulus(UserExpression):

	def __init__(self, subdomains, Mat_0, Mat_1, **kwargs):
		super().__init__(**kwargs)
		self.subdomains = subdomains
		self.Mat_0 = Mat_0
		self.Mat_1 = Mat_1

	def eval_cell(self, values, x, cell):
		if self.subdomains[cell.index] == 2:
			values[0] = self.Mat_0

		else:
			values[0] = self.Mat_1

	def value_shape(self):
		return ()