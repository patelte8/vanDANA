########################################################################

import numpy
import vtk
import vtk_py3

#from mat_vec_tools    import *
from vtk_py3.createFloatArray import *

########################################################################

def addFieldPrincipalDirections(ugrid_mesh,
                                field_name,
                                field_type="cell",
                                field_storage="vec",
                                verbose=True):

    if (verbose): print('*** addFieldPrincipalDirections ***')

    assert (field_type in ["point", "cell"]), "\"field_type\" must be \"point\" or \"cell\". Aborting."
    assert (field_storage in ["vec", "Cmat", "Fmat"]), "\"field_storage\" must be \"vec\", \"Cmat\" or \"Fmat\". Aborting."

    if   (field_type == "cell" ): nb_cells = ugrid_mesh.GetNumberOfCells()
    elif (field_type == "point"): nb_cells = ugrid_mesh.GetNumberOfPoints()

    if   (field_type == "cell" ): field = ugrid_mesh.GetCellData().GetArray(field_name)
    elif (field_type == "point"): field = ugrid_mesh.GetPointData().GetArray(field_name)

    field_Lmin = createFloatArray(field_name+'_Lmin', 1, nb_cells)
    field_Lmid = createFloatArray(field_name+'_Lmid', 1, nb_cells)
    field_Lmax = createFloatArray(field_name+'_Lmax', 1, nb_cells)

    field_Vmin = createFloatArray(field_name+'_Vmin', 3, nb_cells)
    field_Vmid = createFloatArray(field_name+'_Vmid', 3, nb_cells)
    field_Vmax = createFloatArray(field_name+'_Vmax', 3, nb_cells)

    for num_cell in range(nb_cells):
        if   (field_storage == "vec" ): matrix = vec_col_to_mat_sym(field.GetTuple(num_cell))
        elif (field_storage == "Cmat"): matrix = numpy.reshape(field.GetTuple(num_cell), (3,3), order='C')
        elif (field_storage == "Fmat"): matrix = numpy.reshape(field.GetTuple(num_cell), (3,3), order='F')

        if (numpy.linalg.norm(matrix) > 1e-6):
            #if (verbose): print 'num_cell =', num_cell

            val, vec = numpy.linalg.eig(matrix)
            #if (verbose): print 'val =', val
            #if (verbose): print 'vec =', vec
            idx = val.argsort()
            val = val[idx]
            vec = vec[:,idx]
            #if (verbose): print 'val =', val
            #if (verbose): print 'vec =', vec

            matrix_Lmin = [val[0]]
            matrix_Lmid = [val[1]]
            matrix_Lmax = [val[2]]

            matrix_Vmin = vec[:,0]
            matrix_Vmid = vec[:,1]
            matrix_Vmax = vec[:,2]
        else:
            matrix_Lmin = [0.]
            matrix_Lmid = [0.]
            matrix_Lmax = [0.]
            matrix_Vmin = [0.]*3
            matrix_Vmid = [0.]*3
            matrix_Vmax = [0.]*3

        field_Lmin.InsertTuple(num_cell, matrix_Lmin)
        field_Lmid.InsertTuple(num_cell, matrix_Lmid)
        field_Lmax.InsertTuple(num_cell, matrix_Lmax)
        field_Vmin.InsertTuple(num_cell, matrix_Vmin)
        field_Vmid.InsertTuple(num_cell, matrix_Vmid)
        field_Vmax.InsertTuple(num_cell, matrix_Vmax)

    if (field_type == "cell"):
        ugrid_mesh.GetCellData().AddArray(field_Lmin)
        ugrid_mesh.GetCellData().AddArray(field_Lmid)
        ugrid_mesh.GetCellData().AddArray(field_Lmax)
        ugrid_mesh.GetCellData().AddArray(field_Vmin)
        ugrid_mesh.GetCellData().AddArray(field_Vmid)
        ugrid_mesh.GetCellData().AddArray(field_Vmax)
    elif (field_type == "point"):
        ugrid_mesh.GetPointData().AddArray(field_Lmin)
        ugrid_mesh.GetPointData().AddArray(field_Lmid)
        ugrid_mesh.GetPointData().AddArray(field_Lmax)
        ugrid_mesh.GetPointData().AddArray(field_Vmin)
        ugrid_mesh.GetPointData().AddArray(field_Vmid)
        ugrid_mesh.GetPointData().AddArray(field_Vmax)

if (__name__ == "__main__"):
    assert (len(sys.argv) in [3,4,5]), "Number of arguments must be 2, 3 or 4. Aborting."
    ugrid_mesh = readUGrid(sys.argv[1])
    if (len(sys.argv) == 2):
        addFieldPrincipalDirections(ugrid_mesh, sys.argv[2])
    elif (len(sys.argv) == 3):
        addFieldPrincipalDirections(ugrid_mesh, sys.argv[2], sys.argv[3])
    elif (len(sys.argv) == 4):
        addFieldPrincipalDirections(ugrid_mesh, sys.argv[2], sys.argv[3], sys.argv[4])
    writeUGrid(ugrid_mesh, sys.argv[1])
