from dolfin import Expression
from .user_parameters import characteristic_scales, time_control
from scipy.interpolate import splrep, splev

import sys, os
sys.path.insert(0,  '..')
from utilities.read import read_boundary_conditions

directory = os.path.dirname(os.path.abspath(__file__)) + '/'

Lsc = characteristic_scales['Lsc']
Vsc = characteristic_scales['Vsc']
Tsc = Lsc/Vsc

# # Cardiac cycle specifics
# t_period = 60/63					# heartrate (sec)            
# total_heart_cycles = 40   

# Blood perfusion rate
perf = 0
blood_perfusion = False

if blood_perfusion == True:
	perf = 0.85/60 						# ml/s/gm
	perf *= (0.14*Tsc)					# 14% of total coronary-artery perfusion

# # Non-dimensionalizing time-period
# t_period /= Tsc 

# time_control.update(T = 0.04)#total_heart_cycles*t_period)

# # Read boundary conditions from csv file      
# xdata, ydata = read_boundary_conditions(directory, 'RSPV.csv')
# cs_RSPV = splrep(xdata,ydata,per=True)
# xdata, ydata = read_boundary_conditions(directory, 'LSPV.csv')
# cs_LSPV = splrep(xdata,ydata,per=True)
# xdata, ydata = read_boundary_conditions(directory, 'RIPV.csv')
# cs_RIPV = splrep(xdata,ydata,per=True)
# xdata, ydata = read_boundary_conditions(directory, 'LIPV.csv')
# cs_LIPV = splrep(xdata,ydata,per=True)

# # Expressions used during runtime
# tim = Expression('t', t=0.0, degree=1)
# num_cycle = Expression('cycle', cycle=0.0, degree=1)

# param_LSPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 322.5, "Vsc": Vsc, "Tsc": Tsc, "func": cs_LSPV};
# param_LIPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 209.9, "Vsc": Vsc, "Tsc": Tsc, "func": cs_LIPV};
# param_RSPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 188.35, "Vsc": Vsc, "Tsc": Tsc, "func": cs_RSPV};
# param_RIPV={"time": tim, "nm": num_cycle, "period": t_period, "Area": 437.6, "Vsc": Vsc, "Tsc": Tsc, "func": cs_RIPV};                     

def evaluate_boundary_val(a):

    val = (splev((a['time'].t - a['nm'].cycle*a['period'])*a['Tsc'], a['func'])/a['Area'])/a['Vsc']
    return val

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
