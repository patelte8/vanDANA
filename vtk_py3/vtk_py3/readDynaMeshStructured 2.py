########################################################################

import vtk
import math

########################################################################

def readDynaMeshStructured(lsdyna_mesh_file_name, cell_type='hexahedron', verbose=True):

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
    else:
        print 'Wrong cell type. Aborting.'
        exit()

    if (verbose): print 'Reading Dyna mesh file...'

    cnt = 1;

    mesh_file = open(lsdyna_mesh_file_name, 'r')

    context = ''
    for line in mesh_file:
        if (line[-1:] == '\n'): line = line[:-1]
        #if (verbose): print 'line =', line

        if line.startswith('$'): continue

        if (context == 'reading nodes'):
            if line.startswith('*'):
                context = ''
            else:
                #splitted_line = line.split(',')
		splitted_line = []
		splitted_line.append(line[0:8].replace(" ",""))
		splitted_line.append(line[8:24].replace(" ",""))
		splitted_line.append(line[25:25+16].replace(" ",""))
		splitted_line.append(line[25+17:50+17].replace(" ",""))
		print splitted_line

                points.InsertNextPoint([float(coord) for coord in splitted_line[1:4]])
                if (cell_type == 'vertex'):
                    cell.GetPointIds().SetId(0, points.GetNumberOfPoints()-1)
                    cell_array.InsertNextCell(cell)

        if (context == 'reading elems'):
            if line.startswith('*'):
                context = ''
            else:
		if(cnt%2 == 0):

			splitted_line = []
			splitted_line.append(line[0:8].replace(" ",""))
			splitted_line.append(line[9:16].replace(" ",""))
			splitted_line.append(line[17:24].replace(" ",""))
			splitted_line.append(line[25:32].replace(" ",""))
			splitted_line.append(line[33:40].replace(" ",""))
			splitted_line.append(line[41:48].replace(" ",""))
			splitted_line.append(line[49:56].replace(" ",""))
			splitted_line.append(line[57:64].replace(" ",""))
                	if (len(splitted_line) == 3): continue
                	if (cell_type == 'hexahedron'):
                	    for num_node in range(8):
                	        cell.GetPointIds().SetId(num_node, int(splitted_line[num_node])-1)
                	cell_array.InsertNextCell(cell)

		cnt = cnt + 1;

        if line.startswith('*NODE'): context = 'reading nodes'
        if line.startswith('*ELEMENT_SOLID'): context = 'reading elems'

    mesh_file.close()

    if (verbose): print 'Creating VTK mesh...'

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.SetCells(cell_vtk_type, cell_array)

    return ugrid
