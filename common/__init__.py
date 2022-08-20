# __init__.py

from dolfin import parameters
from .constitutive_eq import *
from .solver_options import * 
from .functions import *
from .variational_problem import *
from .stabilizations import *

# Optimization options for dolfin
parameters.update({ "linear_algebra_backend":"PETSc",
                    "form_compiler":{"representation": 'uflacs',
                                     "optimize": True,
                                     "cpp_optimize": True,
                                     "quadrature_degree": 5,
                                     "cpp_optimize_flags": "-O3"},
                    "std_out_all_processes": False})