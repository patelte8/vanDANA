########################################################################

import numpy

########################################################################

def vec_col_to_mat_sym(vec):
    if (numpy.size(vec,0) == 3):
        return numpy.array([[vec[0], vec[2]],
                            [vec[2], vec[1]]])
    elif (numpy.size(vec,0) == 6):
        return numpy.array([[vec[0], vec[3], vec[4]],
                            [vec[3], vec[1], vec[5]],
                            [vec[4], vec[5], vec[2]]])
#        return numpy.array([[vec[0], vec[3], vec[5]],
#                            [vec[3], vec[1], vec[4]],
#                            [vec[5], vec[4], vec[2]]])

    else:
        print 'Wrong vector dimension in vec_col_to_mat_sym. Aborting.'
        exit()

def mat_sym_to_vec_col(mat):
    if (numpy.size(mat,0) == 2) and (numpy.size(mat,1) == 2):
        return numpy.array([mat[0,0], mat[1,1], mat[0,1]])
    elif (numpy.size(mat,0) == 3) and (numpy.size(mat,1) == 3):
        return numpy.array([mat[0,0], mat[1,1], mat[2,2], mat[0,1], mat[0,2], mat[1,2]])
        #return numpy.array([mat[0,0], mat[1,1], mat[2,2], mat[0,1], mat[1,2], mat[0,2]])
    else:
        print 'Wrong matrix dimension in mat_sym_to_vec_col. Aborting.'
        exit()
