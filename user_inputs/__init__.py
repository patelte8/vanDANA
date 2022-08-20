# __init__.py

from .user_parameters import *
from .boundary_initial_conditions import *
from .problem_specific import *

PI = 3.14159265

# Calculate non-dimensional numbers
Re, Pr, Ec, Fr = calc_non_dimensional_numbers(**physical_parameters, **characteristic_scales)

Pe = Re*Pr     
if not boolean_parameters.get('viscous_dissipation') : Ec = 0.0

