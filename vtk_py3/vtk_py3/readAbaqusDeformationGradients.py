########################################################################

import numpy
import vtk

########################################################################

def readAbaqusDeformationGradients(data_file_name,
                                   verbose=True):

    if (verbose): print '*** readAbaqusDeformationGradients ***'

    F_array = createFloatArray("F", 9)

    data_file = open(data_file_name, 'r')
    context = ""
    num_cell = 0
    for line in data_file:
        if (context == "reading deformation gradients"):
            #print line
            if ("MAXIMUM" in line):
                context = ""
                continue
            if ("OR" in line):
                splitted_line = line.split()
                assert (int(splitted_line[0]) == num_cell+1), "Wrong element number. Aborting."
                F_list = [float(splitted_line[ 3]), float(splitted_line[ 6]), float(splitted_line[7]),
                          float(splitted_line[ 9]), float(splitted_line[ 4]), float(splitted_line[8]),
                          float(splitted_line[10]), float(splitted_line[11]), float(splitted_line[5])]
                F_array.InsertNextTuple(F_list)
                num_cell += 1

        if (line == "    ELEMENT  PT FOOT-       DG11        DG22        DG33        DG12        DG13        DG23        DG21        DG31        DG32    \n"):
            context = "reading deformation gradients"

    data_file.close()

    if (verbose): print "nb_tuples = " + str(F_array.GetNumberOfTuples())

    return F_array





