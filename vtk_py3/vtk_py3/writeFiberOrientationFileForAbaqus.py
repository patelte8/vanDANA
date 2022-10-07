########################################################################

import sys

from readUGrid import *

########################################################################

def writeFiberOrientationFileForAbaqus(mesh,
                                       fiber_orientation_file_name,
                                       eF_field_name="eF",
                                       eS_field_name="eS",
                                       sep=", ",
                                       verbose=True):
    if (verbose): print "*** writeFiberOrientationFileForAbaqus ***"

    orientation_file = open(fiber_orientation_file_name, "w")
    orientation_file.write(", 1., 0., 0., 0., 1., 0." + "\n")

    nb_cells = mesh.GetNumberOfCells()

    eF_array = mesh.GetCellData().GetArray(eF_field_name)
    eS_array = mesh.GetCellData().GetArray(eS_field_name)

    for num_cell in range(nb_cells):
        eF = eF_array.GetTuple(num_cell)
        eS = eS_array.GetTuple(num_cell)

        line = str(num_cell+1)
        for k in range(3): line += sep + str(eF[k])
        for k in range(3): line += sep + str(eS[k])
        line += "\n"
        orientation_file.write(line)

    orientation_file.close()

if (__name__ == "__main__"):
    assert (len(sys.argv) in [2,3]), "Number of arguments must be 1 or 2. Aborting."
    if (len(sys.argv) == 2):
        mesh_file_name              = sys.argv[1] + "-Mesh.vtk"
        fiber_orientation_file_name = sys.argv[1] + "-Orientation.inp"
    elif (len(sys.argv) == 3):
        mesh_file_name              = sys.argv[1]
        fiber_orientation_file_name = sys.argv[2]
    mesh = readUGrid(mesh_file_name)
    writeFiberOrientationFileForAbaqus(mesh, fiber_orientation_file_name)
