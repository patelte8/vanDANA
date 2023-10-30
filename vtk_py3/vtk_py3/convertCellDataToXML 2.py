########################################################################

import vtk

########################################################################

def convertCellDataToXML(celldata, xmlfilename):

	xmlfile = open(xmlfilename, "w")
	print >>xmlfile, "<?xml version='1.0'?>"
	print >>xmlfile, "<dolfin xmlns:dolfin='http://fenicsproject.org'>"
	print >>xmlfile, "<function_data size='"+str(celldata.GetNumberOfTuples())+"'>"

	cnt = 0
	for cellid in range(0,celldata.GetNumberOfTuples()):
		data = celldata.GetTuple(cellid)[0]
		print >>xmlfile, "<dof index='" + str(cnt) + "' value='" + str(data) + "' cell_index='" + str(cellid) + "' cell_dof_index='0' />"
		cnt += 1

	print >>xmlfile, "</function_data>"
	print >>xmlfile, "</dolfin>"

