########################################################################

import numpy
import vtk

########################################################################

def readAbaqusStress(data_file_name,
                     verbose=True):

    if (verbose): print '*** readAbaqusStress ***'

    s_array = createFloatArray("", 6)

    data_file = open(data_file_name, 'r')
    context = ""
    num_cell = 0
    for line in data_file:
        if (context == "reading stresses"):
            #print line
            if ("MAXIMUM" in line):
                context = ""
                continue
            if ("OR" in line):
                splitted_line = line.split()
                assert (int(splitted_line[0]) == num_cell+1), "Wrong element number. Aborting."
                s_list = [float(splitted_line[k]) for k in range(3,9)]
                s_array.InsertNextTuple(s_list)
                num_cell += 1

        if (line == "    ELEMENT  PT FOOT-       S11         S22         S33         S12         S13         S23     \n"):
            context = "reading stresses"

    data_file.close()

    if (verbose): print "nb_tuples = " + str(s_array.GetNumberOfTuples())

    return s_array





