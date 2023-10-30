import pyquaternion as pyq
import numpy as np 
import pdb

def bislerp(Qa, Qb, t):
	
    Qa_M = mat3([Qa[0,0], Qa[0,1], Qa[0,2], Qa[1,0], Qa[1,1], Qa[1,2], Qa[2,0], Qa[2,1], Qa[2,2]])
    Qb_M = mat3([Qb[0,0], Qb[0,1], Qb[0,2], Qb[1,0], Qb[1,1], Qb[1,2], Qb[2,0], Qb[2,1], Qb[2,2]])
    qa = quat(Qa_M)
    qb = quat(Qb_M)
    
    val = np.zeros(8)
    quat_i = quat(0,1,0,0)
    quat_j = quat(0,0,1,0)
    quat_k = quat(0,0,0,1)
    quat_array = [qa, -qa, qa*quat_i, -qa*quat_i, qa*quat_j, -qa*quat_j, qa*quat_k, -qa*quat_k] 
    cnt = 0
    for qt in quat_array:
        val[cnt] = qt.dot(qb)
        cnt = cnt + 1
	
    qt = quat_array[val.argmax(axis=0)]	
	
    if(t < 0):
        t = 0.0
	
    qm = slerp(t, qt, qb)
    qm = qm.normalize()
    Qm_M = qm.toMat3()
    Qm = [[Qm_M[0,0], Qm_M[0,1], Qm_M[0,2]], [Qm_M[1,0], Qm_M[1,1], Qm_M[1,2]], [Qm_M[2,0], Qm_M[2,1], Qm_M[2,2]]]

    return Qm


def my_bislerp(Qa, Qb, t):
	
    #Qa_M = mat3([Qa[0,0], Qa[0,1], Qa[0,2], Qa[1,0], Qa[1,1], Qa[1,2], Qa[2,0], Qa[2,1], Qa[2,2]])
    #Qb_M = mat3([Qb[0,0], Qb[0,1], Qb[0,2], Qb[1,0], Qb[1,1], Qb[1,2], Qb[2,0], Qb[2,1], Qb[2,2]])
    #qa = quat(Qa_M)
    #qb = quat(Qb_M)

    qa = pyq.Quaternion(matrix=Qa)
    qb = pyq.Quaternion(matrix=Qb)


    val = np.zeros(8)
    #quat_i = quat(0,1,0,0)
    #quat_j = quat(0,0,1,0)
    #quat_k = quat(0,0,0,1)

    quat_i = pyq.Quaternion(0,1,0,0)
    quat_j = pyq.Quaternion(0,0,1,0)
    quat_k = pyq.Quaternion(0,0,0,1)

    quat_array = [qa, -qa, qa*quat_i, -qa*quat_i, qa*quat_j, -qa*quat_j, qa*quat_k, -qa*quat_k] 

    
    #print('My Quarternion NUMPY array {arr}'.format(arr=np.array(qa)))
    #print('My Quarternion NUMPY array {arr}'.format(arr=np.array([qa[0], qa[1], qa[2], qa[3] ]) ) )


    cnt = 0
    qb_arr = np.array([qb[0], qb[1], qb[2], qb[3]])
    for qt in quat_array:
        qt_arr = np.array([qt[0], qt[1], qt[2], qt[3]])
        val[cnt] = np.dot(qt_arr, qb_arr)
        cnt = cnt + 1

    print('calculated val array is :{arr}'.format(arr=val))
    qt = quat_array[val.argmax(axis=0)]	
    print('calculated val array is :{arr}'.format(arr=qt))



    if(t < 0):
        t = 0.0
    
    #qm = slerp(t, qt, qb)
    qm = pyq.Quaternion.slerp(qt, qb, t)
    qm = qm.normalised
    Qm_M = qm.rotation_matrix

    print('calculated val matrix is :{mat}'.format(mat=Qm_M))

    #Qm = [[Qm_M[0,0], Qm_M[0,1], Qm_M[0,2]], [Qm_M[1,0], Qm_M[1,1], Qm_M[1,2]], [Qm_M[2,0], Qm_M[2,1], Qm_M[2,2]]]

    return Qm_M


def run_demo():
    #https://github.com/KieranWynn/pyquaternion

    # Create a quaternion representing a rotation of +90 degrees about positive y axis.
    my_quaternion = pyq.Quaternion(axis=[0, 1, 0], degrees=90)

    my_vector = [0, 0, 4]
    my_rotated_vector = my_quaternion.rotate(my_vector)

    print('\nBasic Rotation')
    print('--------------')
    print('My Vector: {}'.format(my_vector))
    print('Performing rotation of {angle} deg about {axis}'.format(angle=my_quaternion.degrees, axis=my_quaternion.axis))
    print('My Rotated Vector: {}'.format(my_rotated_vector))

    # Create another quaternion representing no rotation at all
    null_quaternion = pyq.Quaternion(axis=[0, 1, 0], angle=0)

    print('\nInterpolated Rotation')
    print('---------------------')

    # The following will create a sequence of 9 intermediate quaternion rotation objects
    for q in pyq.Quaternion.intermediates(null_quaternion, my_quaternion, 9, include_endpoints=True):
        my_interpolated_point = q.rotate(my_vector)
        print('My Interpolated Point: {point}\t(after rotation of {angle} deg about {axis})'.format(
            point=my_interpolated_point, angle=round(q.degrees, 4), axis=q.axis
        ))
    
    print('Done!')

def run_tests():
    rotation = np.eye(3)
    q = pyq.Quaternion(matrix=rotation) # Using 3x3 rotation matrix
    print('My first Quaternion from matrix is: {quat}\t with angle {angle} and axis {axis}'.format(
            quat=q, angle=q.degrees, axis=q.axis))

    q1 = pyq.Quaternion.random()
    print('My first Quaternion from matrix is: {quat}\t with angle {angle} and axis {axis}'.format(
            quat=q1, angle=q1.degrees, axis=q1.axis))

    q2 = pyq.Quaternion.random()
    print('My second random Quaternion is: {quat}\t with angle {angle} and axis {axis}'.format(
            quat=q2, angle=q2.degrees, axis=q2.axis))

    q1_q2 = q1*q2
    print('My multiplied Quaternion is: {quat}\t with angle {angle} and axis {axis}'.format(
            quat=q1_q2, angle=q1_q2.degrees, axis=q1_q2.axis))

    q_interp = pyq.Quaternion.slerp(q1, q2, amount = 0.5)
    print('My SLERP interpiolated Quaternion is: {quat}\t with angle {angle} and axis {axis}'.format(
            quat=q_interp, angle=q_interp.degrees, axis=q_interp.axis)) 

    q_interp_normalzied = q_interp.normalised 
    print('My normalized interpolated Quaternion is: {quat}\t with angle {angle} and axis {axis}'.format(
            quat=q_interp_normalzied, angle=q_interp_normalzied.degrees, axis=q_interp_normalzied.axis)) 

    q_interp_matrix = q_interp.rotation_matrix 
    print('My interpolated Quaternion in matrix form is: \n {quat}'.format(
            quat=q_interp_matrix))


if __name__ == "__main__": 

    #run_demo()
    #run_tests()

    q1 = pyq.Quaternion.random()
    q1mat = q1.rotation_matrix 
    print('My roated Quaternion in matrix form is: \n {quat}'.format(
            quat=q1mat))

    q2 = pyq.Quaternion.random()
    q2mat = q2.rotation_matrix 
    print('My rotated Quaternion in matrix form is: \n {quat}'.format(
            quat=q2mat))

    q_interp_mat = my_bislerp(q1mat, q2mat, 0.3)
    print('My interpolated Quaternion in matrix form is: \n {quat}'.format(
            quat=q_interp_mat))

    pdb.set_trace()
