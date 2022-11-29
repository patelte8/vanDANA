# __init__.py

from dolfin import parameters
from .constitutive_eq import *
from .solver_options import * 
from .functions import *
from .flow_variational_problem import *
from .solid_variational_problem import *
from .flow_temperature_variational_problem import *
from .lagrange_variational_problem import *
from .delta_interpolation import fsi_interpolation_code
from .fem_stabilizations import *
import sys, os, cppimport

# Optimization options for dolfin
parameters.update({ "linear_algebra_backend": "PETSc",
                    "form_compiler": FFC_parameters,
                    "std_out_all_processes": False})
