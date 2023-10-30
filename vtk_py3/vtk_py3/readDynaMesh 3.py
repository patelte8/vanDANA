########################################################################

import vtk

########################################################################

def readDynaMesh(lsdyna_mesh_file_name, cell_type='hexahedron', verbose=True, filetype='post'):

    if (verbose): print '*** readDynaMesh ***'

    points = vtk.vtkPoints()

    if (cell_type == 'vertex'):
        cell_vtk_type = vtk.VTK_VERTEX
        cell = vtk.vtkVertex()
        cell_array = vtk.vtkCellArray()
    elif (cell_type == 'hexahedron'):
        cell_vtk_type = vtk.VTK_HEXAHEDRON
        cell = vtk.vtkHexahedron()
        cell_array = vtk.vtkCellArray()
	matid_array = vtk.vtkIntArray()
	matid_array.SetName("matid")
	fiber_array = vtk.vtkFloatArray()
	fiber_array.SetNumberOfComponents(3)
	fiber_array.SetName("fvec")
	sheet_array = vtk.vtkFloatArray()
	sheet_array.SetNumberOfComponents(3)
	sheet_array.SetName("svec")
    else:
        print 'Wrong cell type. Aborting.'
        exit()

    if (verbose): print 'Reading Dyna mesh file...'

    mesh_file = open(lsdyna_mesh_file_name, 'r')

    context = ''
    isfiber = False
    issheet = False
    for line in mesh_file:
        if (line[-1:] == '\n'): line = line[:-1]
        #if (verbose): print 'line =', line

        if line.startswith('$'): continue

        if (context == 'reading nodes'):
            if line.startswith('*'):
                context = ''
            else:
                splitted_line = line.split(',')
                points.InsertNextPoint([float(coord) for coord in splitted_line[1:4]])
                if (cell_type == 'vertex'):
                    cell.GetPointIds().SetId(0, points.GetNumberOfPoints()-1)
                    cell_array.InsertNextCell(cell)

        if (context == 'reading elems'):
            if line.startswith('*'):
                context = ''
            else:
                splitted_line = line.split(',')

		if(filetype == 'post'):
                	if (len(splitted_line) == 3 and isfiber == True):
			    fiber_array.InsertNextTuple3(float(splitted_line[0]), float(splitted_line[1]), float(splitted_line[2]))
			    isfiber = False; issheet = True;
			elif (len(splitted_line) == 3 and issheet == True):
			    sheet_array.InsertNextTuple3(float(splitted_line[0]), float(splitted_line[1]), float(splitted_line[2]))
			    isfiber = False; issheet = False;

                	if (len(splitted_line) != 3):
			    matid_array.InsertNextValue(int(splitted_line[1]))
                	    for num_node in range(8):
                	        cell.GetPointIds().SetId(num_node, int(splitted_line[2+num_node])-1)
			    isfiber = True;
                	    cell_array.InsertNextCell(cell)

		elif(filetype == 'pre'):

			matid_array.InsertNextValue(int(splitted_line[1]))
                	for num_node in range(8):
                	    cell.GetPointIds().SetId(num_node, int(splitted_line[2+num_node])-1)
			isfiber = True;
                	cell_array.InsertNextCell(cell)


        if line.startswith('*NODE'): context = 'reading nodes'
        if line.startswith('*ELEMENT_SOLID'): context = 'reading elems'

    mesh_file.close()

    if (verbose): print 'Creating VTK mesh...'

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.SetCells(cell_vtk_type, cell_array)
    ugrid.GetCellData().AddArray(matid_array)
    if(filetype == 'post'):
    	ugrid.GetCellData().AddArray(fiber_array)
    	ugrid.GetCellData().AddArray(sheet_array)

    return ugrid
