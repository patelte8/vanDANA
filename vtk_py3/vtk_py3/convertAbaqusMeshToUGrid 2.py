########################################################################

import sys

from readAbaqusMesh import *
from writeUGrid     import *

########################################################################

name = sys.argv[1]

mesh = readAbaqusMesh(name + ".inp", "hex")
writeUGrid(mesh, name + ".vtk")
