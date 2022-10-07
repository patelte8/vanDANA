########################################################################

import vtk

########################################################################

def readAbaqusMesh(mesh_filename,
                   elem_types="all",
                   verbose=True):

    if (verbose): print "*** readAbaqusMesh ***"

    points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()

    mesh_file = open(mesh_filename, "r")

    context = ""
    for line in mesh_file:
        if (line[-1:] == "\n"): line = line[:-1]
        #if (verbose): print "line =", line

        if line.startswith("**"): continue

        if (context == "reading nodes"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                points.InsertNextPoint([float(coord) for coord in splitted_line[1:4]])

        if (context == "reading elems"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                assert (len(splitted_line) == 1+cell_nb_points), "Wrong number of elements in line. Aborting."
                for num_point in range(cell_nb_points): cell.GetPointIds().SetId(num_point, int(splitted_line[1+num_point])-1)
                cell_array.InsertNextCell(cell)

        if line.startswith("*NODE"):
            context = "reading nodes"
        if line.startswith("*ELEMENT"):
            if ("TYPE=F3D4" in line) and (("quad" in elem_types) or ("all" in elem_types)):
                context = "reading elems"
                cell_vtk_type = vtk.VTK_QUAD
                cell_nb_points = 4
                cell = vtk.vtkQuad()
            elif ("TYPE=C3D4" in line) and (("tet" in elem_types) or ("all" in elem_types)):
                context = "reading elems"
                cell_vtk_type = vtk.VTK_TETRA
                cell_nb_points = 4
                cell = vtk.vtkTetra()
            elif ("TYPE=C3D8" in line) and (("hex" in elem_types) or ("all" in elem_types)):
                context = "reading elems"
                cell_vtk_type = vtk.VTK_HEXAHEDRON
                cell_nb_points = 8
                cell = vtk.vtkHexahedron()
            else:
                print "Warning: element type not taken into account."

    mesh_file.close()

    if (verbose): print "Creating UGrid..."

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.SetCells(cell_vtk_type, cell_array)

    if (verbose): print "nb_cells = " + str(ugrid.GetNumberOfCells())

    return ugrid
