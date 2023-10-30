#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def myPrint(
    verbose,
    string):

    if not hasattr(myPrint, "initialized"):
        myPrint.initialized = True
        myPrint.verbose_ini = verbose
    if (verbose > 0): print (myPrint.verbose_ini - verbose) * 4 * " " + string
