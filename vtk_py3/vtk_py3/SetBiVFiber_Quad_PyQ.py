########################################################################
import argparse
from mpi4py import MPI as pyMPI
import numpy as np
#from thLib import quat 
import sys
import os
from dolfin import *
from sympy import Symbol, nsolve 
import math
from scipy import linalg
from scipy import optimize

import vtk_py3 as vtk_py

import pyquaternion as pyq
import pdb

import warnings
from ffc.quadrature.deprecation \
import QuadratureRepresentationDeprecationWarning


########################################################################


#def SetBiVFiber_Quad(mesh, boundaries, mtvfiberfilename, mtvsheetfilename, al_endo, al_epi, b_endo, b_epi, isrotatept, isreturn, outfilename, param, outdirectory="./", epiid=2, rvid=1, lvid=3, degree=2):
def SetBiVFiber_Quad_PyQ(param):
        '''
        What: This function uses pyquaternion lib instead of cgkit for quarternion calculations
        What more: pip install (http://kieranwynn.github.io/pyquaternion/) 

        Where: my_bislerp() instead of bislerp() in SetBiVFiber_Quad

        How: vtk_py.SetBiVFiber_J( ... )
        
        Plaussible issues:  matrix to quaternion is not unique, topology is trickier
                            http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/tim.htm 
                            http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/ 

        What else: comment out all other problematic lines 
                   for eg., Sdof_coordinates 
        
        '''

        #pdb.set_trace()
        
        #assert "mesh" in param,  "No mesh is found"
        #assert "facetboundaries" in param,  "No facetboundaries are found"
        warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)

        mesh = param["mesh"]
        boundaries = param["facetboundaries"]
        matid = param["matid"] if ("matid" in param) else None
        
        epiid = param["epiid"] if ("epiid" in param) else 2
        lvid = param["lvid"] if ("lvid" in param) else 3
        rvid = param["rvid"] if ("rvid" in param) else 1
        outdirectory = param["outdirectory"] if ("outdirectory" in param) else "./"
        isrotatept = param["isrotatept"] if ("isrotatept" in param) else False
        isreturn = param["isreturn"] if ("isreturn" in param) else True
        isLV = param["isLV"] if ("isLV" in param) else False
        isscaling = param["isscaling"] if ("isscaling" in param) else False
        outfilename = param["outfilename"] if ("outfilename" in param) else "fiber"
        degree = param["degree"] if ("degree" in param) else 2
        mtvfiberfilename = param["mtvfiberfilename"] if ("mtvfiberfilename" in param) else None
        mtvsheetfilename = param["mtvsheetfilename"] if ("mtvsheetfilename" in param) else None
        
        lvmatid = param["LV_matid"] if ("LV_matid" in param) else 0 
        rvmatid = param["RV_matid"] if ("RV_matid" in param) else 0 
        septummatid = param["Septum_matid"] if ("Septum_matid" in param) else 0 
        
        lv_fiber_angle = np.array(param["LV_fiber_angle"])*math.pi/180 if ("LV_fiber_angle" in param) else np.array([60, -60])*math.pi/180
        rv_fiber_angle = np.array(param["RV_fiber_angle"])*math.pi/180 if ("RV_fiber_angle" in param) else np.array([60, -60])*math.pi/180
        septum_fiber_angle = np.array(param["Septum_fiber_angle"])*math.pi/180 if ("Septum_fiber_angle" in param) else np.array([60, -60])*math.pi/180
        
        lv_sheet_angle = np.array(param["LV_sheet_angle"])*math.pi/180 if ("LV_sheet_angle" in param) else np.array([0.1, -0.1])*math.pi/180
        rv_sheet_angle = np.array(param["RV_sheet_angle"])*math.pi/180 if ("RV_sheet_angle" in param) else np.array([0.1, -0.1])*math.pi/180
        septum_sheet_angle = np.array(param["Septum_sheet_angle"])*math.pi/180 if ("Septum_sheet_angle" in param) else np.array([0.1, -0.1])*math.pi/180
        
        
        # Set Fiber angle
        #al_endo = al_endo*math.pi/180; b_endo = b_endo*math.pi/180;
        #al_epi = al_epi*math.pi/180;  b_epi = b_epi*math.pi/180;
        
        comm = MPI.comm_world #mesh.mpi_comm().tompi4py()
        
        minz = comm.allreduce(np.amin(mesh.coordinates()[:,2]), op=pyMPI.MIN)
        maxz = comm.allreduce(np.amax(mesh.coordinates()[:,2]), op=pyMPI.MAX)
        
        
        # Boundary condition for Laplace equation from the bottom to the top 
        
        def right_boundary(x):
                return x[2] < minz+1.0
        
        # Function that solves the Laplace equation and and calculates the gradient of the solution
        def lap(Me, boud, par,op, filename):
        
                set_log_level(40)
                parameters["form_compiler"]["quadrature_degree"] = degree
                parameters["form_compiler"]["representation"] = 'quadrature'
        
                epi=epiid; rv=rvid; lv=lvid;
                # Create mesh and define function space
                mesh = Me
                V = FunctionSpace(mesh, FiniteElement("Lagrange", mesh.ufl_cell(), 2))
                SQuad = FunctionSpace(mesh, FiniteElement("Quadrature", mesh.ufl_cell(), degree=degree, quad_scheme="default"))
        
                y = mesh.coordinates()
                sa = y.shape
                ttt = y[:, 2].min()
                cord =  np.argwhere(y == ttt); 
                point = y[cord[0,0],:]
                #print  sa, cord, point
                u0 = Constant(0.0)
                if op == 1:
                        bc =[ DirichletBC(V, u0, right_boundary),  DirichletBC(V, 1, boud, 4)] 
                # Define boundary conditions
                else:
                        bc =  [DirichletBC(V, par[0], boud, lv),  DirichletBC(V, par[1], boud, epi ), DirichletBC(V, par[2], boud, rv)]
                
                # Define variational problem
                u = TrialFunction(V)
                v = TestFunction(V)
                f = Constant(0)
                a = inner(nabla_grad(u), nabla_grad(v))*dx
                L = f*v*dx
        
                # Compute solution
                u = Function(V)
                solve(a == L, u, bc)#, solver_parameters={"linear_solver": "petsc"})
                u_a = u.vector().get_local()
                u_a_quad = project(u, SQuad).vector().get_local()#, solver_type="petsc").vector().array()
        
        
                # Compute gradient
                V_g = FunctionSpace(mesh, VectorElement("Lagrange", mesh.ufl_cell(), 2))
                VQuad = FunctionSpace(mesh, VectorElement("Quadrature", mesh.ufl_cell(), degree=degree, quad_scheme="default"))
                v = TestFunction(V_g)
                w = TrialFunction(V_g)
                
                a = inner(w, v)*dx
                L = inner(grad(u), v)*dx
                grad_u = Function(V_g)
                solve(a == L, grad_u)#, solver_parameters={"linear_solver": "petsc"})
                grad_u_quad = project(grad_u, VQuad)#, solver_type="petsc")
                grad_u.rename('grad_u', 'continuous gradient field')
                grad_ua_quad =  grad_u_quad.vector().get_local()#.array()
                
        
                #plot(u, interactive=True)
        
                # Dump solution to file in VTK format
                #newfilename = filename + '.pvd'
                #file1= File(newfilename)
                #file1<< u
        
                #newfilename = filename + 'grad.pvd'
                #file2 = File(newfilename)
                #file2 << grad_u
        
        #        Hold plot
                #interactive()
                return u_a_quad, grad_ua_quad
        
        # Function to calculate the lenght of the vector  
        
        def len_vec(vec):
                l = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                return l
        
        
        
        ##################################################################################################
        # Function to calculate the orthogonal axis using the gradients calculated in the previous step
        ###################################################################################################
        def L1(e0, v0):
        
                L1x = (v0[0] - (e0[0]*v0[0] + e0[1]*v0[1] + e0[2]*v0[2])*e0[0])**2 
                L1y = (v0[1] - (e0[0]*v0[0] + e0[1]*v0[1] + e0[2]*v0[2])*e0[1])**2 
                L1z = (v0[2] - (e0[0]*v0[0] + e0[1]*v0[1] + e0[2]*v0[2])*e0[2])**2 
        
                return (L1x + L1y + L1z)**0.5;
        
        def P(e0, v0):
        
                P1x = v0[0] - (e0[0]*v0[0] + e0[1]*v0[1] + e0[2]*v0[2])*e0[0]
                P1y = v0[1] - (e0[0]*v0[0] + e0[1]*v0[1] + e0[2]*v0[2])*e0[1]
                P1z = v0[2] - (e0[0]*v0[0] + e0[1]*v0[1] + e0[2]*v0[2])*e0[2]
        
                return [P1x, P1y, P1z];
        
        def function_e0(e0, v0, e1):
        
                f = [L1(e0,v0)*e0[0] - (e1[1]*P(e0,v0)[2] - e1[2]*P(e0,v0)[1]),
                     L1(e0,v0)*e0[1] - (e1[2]*P(e0,v0)[0] - e1[0]*P(e0,v0)[2]),
                     L1(e0,v0)*e0[2] - (e1[0]*P(e0,v0)[1] - e1[1]*P(e0,v0)[0])]
        
                return f;
        
        
        def axisf(vec1, vec2):
                
                len_vec2 = len_vec(vec2); len_vec1 = len_vec(vec1);
                e1 = vec1/len_vec1; ini_e2 = vec2/len_vec2;
                ini_e0 = np.cross(e1, ini_e2)
                
        
                e0 = np.zeros(3); e2 = np.zeros(3);     
        
                # Solve using numerical 
                def function_wrap(e0):
                        return function_e0(e0, vec2, e1)
        
                sol = optimize.root(function_wrap, [ini_e0[0], ini_e0[1], ini_e0[2]], method='hybr')
                e0[0] = sol.x[0];  e0[1] = sol.x[1]; e0[2] = sol.x[2]; 
                e2 = P(e0, vec2);
                len_e2 = len_vec(e2)    
                e2[0] = e2[0]/len_e2; e2[1] = e2[1]/len_e2;  e2[2] = e2[2]/len_e2; 
        
                Q = np.array([[e0[0], e1[0], e2[0] ],[ e0[1], e1[1], e2[1]],[ e0[2], e1[2], e2[2]]])
        
        
                return Q
        
        
        #################################################################################################
        #Function to ensure the that the axis are orthogonal 
        #################################################################################################
        def orto(Q, cnt):
                e0 = np.zeros(3); e1 = np.zeros(3); e2 = np.zeros(3);
                e0[0] = Q[0][0]; e0[1] = Q[1][0]; e0[2] = Q[2][0];
                e1[0] = Q[0][1]; e1[1] = Q[1][1]; e1[2] = Q[2][1];
                e2[0] = Q[0][2]; e2[1] = Q[1][2]; e2[2] = Q[2][2];
                e0x = Symbol('e0x'); e0y = Symbol('e0y') ; e0z = Symbol('e0z'); e2x = Symbol('e2x'); e2y = Symbol('e2y') ; e2z = Symbol('e2z');
        
        
                try:
        
                        aa = nsolve(  [ e0x - e1[1]*e2z + e1[2]*e2y, e0y - e1[2]*e2x + e1[0]*e2z, e0z - e1[0]*e2y + e1[1]*e2x, e1[0] - e2y*e0z + e2z*e0y, e1[1] - e2z*e0x + e2x*e0z, e1[2] - e2x*e0y + e2y*e0x  ], [ e0x, e0y, e0z, e2x, e2y, e2z ],[e0[0], e0[1], e0[2], e2[0], e2[1],e2[2] ])
        
        
                        
                        e0[0] = aa[0]; e0[1] = aa[1]; e0[2] = aa[2];
                        e2[0] = aa[3]; e2[1] = aa[4]; e2[2] = aa[5];
                        Qa = np.array([[e0[0], e1[0], e2[0] ],[ e0[1], e1[1], e2[1]],[ e0[2], e1[2], e2[2]]])
        
                except ZeroDivisionError as detail:
                        print('Handling run-time error:', detail)
                        Qa = Q
                        f = quat.rotmat2quat(Qa)
                        outfile = open("quat.txt", "w")
                        print  >>outfile, cnt, " ", f 
                return Qa
        
        
        def orto3(Q, cnt):
        
                Q2 = np.dot(np.transpose(Q),Q)
                Q2inv = np.linalg.inv(Q2)
                sqrtQ2 = linalg.sqrtm(Q2inv)
                Qa = np.dot(Q,sqrtQ2)
        
                return Qa
        
        #################################################################################################
        #Function to ensure the that the axis are orthogonal but this one is not working prorperly
        #################################################################################################
        
        def orto2(Q):
                e0 = np.zeros(3); e1 = np.zeros(3); e2 = np.zeros(3);
                e0[0] = Q[0][0]; e0[1] = Q[1][0]; e0[2] = Q[2][0];
                e1[0] = Q[0][1]; e1[1] = Q[1][1]; e1[2] = Q[2][1];
                e2[0] = Q[0][2]; e2[1] = Q[1][2]; e2[2] = Q[2][2];
                
                while   np.dot(e0, e1) +  np.dot(e0, e2) + np.dot(e1, e2) > 10e-30:
                        e2 = np.cross(e0, e1);
                        e0 = np.cross(e1, e2);
                Qa = np.array([[e0[0], e1[0], e2[0] ],[ e0[1], e1[1], e2[1]],[ e0[2], e1[2], e2[2]]])
                return Qa 
        
        ##############################################################################################
        #Function 3 - This function rotates the axis calculated in the previous steps
        ##############################################################################################
        
        def orient(Q, al, bt):
                arr1 =np.array( [[np.cos(al), -np.sin(al), 0  ],[np.sin(al), np.cos(al), 0],[0, 0, 1]]);
                arr2 = np.array([[1, 0, 0],[0, np.cos(bt), np.sin(bt) ],[0, -np.sin(bt), np.cos(bt) ]]);
                out = np.dot(Q, arr1, arr2)
                return out
        
        ##################################################################################################
        #Function 4 - This function calculates quarternions and interpolates these quarternions
        #################################################################################################
        
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
        
            qa = pyq.Quaternion(matrix=Qa)
            qb = pyq.Quaternion(matrix=Qb)
        
        
            val = np.zeros(8)
        
            quat_i = pyq.Quaternion(0,1,0,0)
            quat_j = pyq.Quaternion(0,0,1,0)
            quat_k = pyq.Quaternion(0,0,0,1)
        
            quat_array = [qa, -qa, qa*quat_i, -qa*quat_i, qa*quat_j, -qa*quat_j, qa*quat_k, -qa*quat_k] 
        
            cnt = 0
            qb_arr = np.array([qb[0], qb[1], qb[2], qb[3]])
        
            for qt in quat_array:
                #val[cnt] = qt.dot(qb)
                qt_arr = np.array([qt[0], qt[1], qt[2], qt[3]])
                val[cnt] = np.dot(qt_arr, qb_arr)
                cnt = cnt + 1
        
            qt = quat_array[val.argmax(axis=0)] 
        
            if(t < 0):
                t = 0.0 
        
            #qm = slerp(t, qt, qb)
            qm = pyq.Quaternion.slerp(qt, qb, t)
            qm = qm.normalised
            Qm_M = qm.rotation_matrix
        
            return Qm_M
        
        
        
        #################################################################################################
        ######Generating results and using the functions ##########
        #################################################################################################
        # DOF map  
        parameters["form_compiler"]["quadrature_degree"] = degree
        parameters["form_compiler"]["representation"] = 'quadrature'
        
        V = FunctionSpace(mesh, VectorElement("Quadrature", mesh.ufl_cell(), degree=degree, quad_scheme="default"))
        S = FunctionSpace(mesh, FiniteElement("Quadrature", mesh.ufl_cell(), degree=degree, quad_scheme="default"))
        
        dof_coordinates = V.tabulate_dof_coordinates()
        Sdof_coordinates = S.tabulate_dof_coordinates()
        n = V.dim()
        d = mesh.geometry().dim()
        dof_coordinates.resize((n, d))
        Sdof_coordinates.resize((n,d))
        
        x_my_first, x_my_last = V.sub(0).dofmap().ownership_range()
        x_dofs = np.arange(0, x_my_last-x_my_first, d)
        y_dofs = np.arange(1, x_my_last-x_my_first, d)
        z_dofs = np.arange(2, x_my_last-x_my_first, d)
        
        
        jj = np.array([ 1, 1, 1] ); hh = np.array([1, 1, 1]); 
        
        
        ####### Solving the poisson equation ############ 
        
        # case 1)   epi -> u = 0 and rv/lv -> u = 1
        
        #if(MPI.rank(mpi_comm_world()) == 0): 
        if(MPI.comm_world.Get_rank() == 0): 
                print("Solve Poisson Eq. 1, epi -> u = 0 and rv/lv -> u = 1")
        par1 = [0, 1, 0]; epi, dd = lap(mesh, boundaries, par1, 0, 'phi_epi')
        
        #total_dd = comm.allreduce(len(dd), op=pyMPI.SUM)
        #num_of_nodes = total_dd/3
        #scalar_array = np.zeros(num_of_nodes)
        #scalar_dof = S.dofmap().dofs()
        
        total_dd = len(dd)
        S_my_first, S_my_last = S.dofmap().ownership_range()
        scalar_dof = filter(lambda dof: S.dofmap().local_to_global_index(dof) not in S.dofmap().local_to_global_unowned(),
                        range(S_my_last-S_my_first))
        
        
        # case 2)   lv -> u = 0 and rv/epi -> u = 1
        
        #if(MPI.rank(mpi_comm_world()) == 0): 
        if(MPI.comm_world.Get_rank() == 0): 
                print("Solve Poisson Eq. 2, lv -> u = 0 and rv/epi -> u = 1")
        par2 = [1, 0, 0]; lv, dd2 = lap(mesh, boundaries, par2, 0, 'phi_lv'); #dd2 = -1*dd2
        
        # case 3)   rv -> u = 0 and epi/lv -> u = 1
        
        #if(MPI.rank(mpi_comm_world()) == 0): 
        if(MPI.comm_world.Get_rank() == 0): 
                print("Solve Poisson Eq. 3, rv -> u = 0 and epi/lv -> u = 1")
        par3 = [0, 0, 1]; rv, dd3 = lap(mesh, boundaries, par3, 0, 'phi_rv')
        
        #case 4) from the top to the bottom
        
        #if(MPI.rank(mpi_comm_world()) == 0): 
        if(MPI.comm_world.Get_rank() == 0): 
                print("Solve Poisson Eq. 4, from the top to the bottom")
        par1 = [0, 1, 0]; b, dd4 = lap(mesh, boundaries, par1, 1, 'phi_ab')
        
        
        
        #  Start calculating the fiber orientation
        cnt = 0; c = 0;
        vector_array = np.zeros(V.dim())
        
        func_of_vector = Function(V); 
        func_of_scalar = Function(S);
        
        func_of_e0 = Function(V); 
        func_of_e1 = Function(V); 
        func_of_e2 = Function(V); 
        
        func_of_e0_vector = func_of_e0.vector()
        func_of_e1_vector = func_of_e1.vector() 
        func_of_e2_vector = func_of_e2.vector()  
        
        e0_fiber = func_of_e0_vector.get_local()
        e1_fiber = func_of_e1_vector.get_local()
        e2_fiber = func_of_e2_vector.get_local()
        
        vec_dd = np.zeros(3); vec_dd2 = np.zeros(3); vec_dd3 = np.zeros(3); vec_dd4 = np.zeros(3);
        
        Qepi = np.zeros((total_dd,3)); Qrv = np.zeros((total_dd,3)); Qlv = np.zeros((total_dd,3));
        Qendo = np.zeros((total_dd,3))
        Qfiber = np.zeros((total_dd,3))
        
        e0_epi = np.zeros(total_dd); e1_epi = np.zeros(total_dd); e2_epi = np.zeros(total_dd);
        e0_rv = np.zeros(total_dd); e1_rv = np.zeros(total_dd); e2_rv = np.zeros(total_dd);
        e0_lv = np.zeros(total_dd); e1_lv = np.zeros(total_dd); e2_lv = np.zeros(total_dd);
        check = np.zeros(total_dd); check2 = np.zeros(total_dd);
        ds =  np.zeros(int(total_dd/3));  al_s =  np.zeros(int(total_dd/3));   b_s = np.zeros(int(total_dd/3));
        al_w =  np.zeros(int(total_dd/3));   b_w = np.zeros(int(total_dd/3));
        
        e0_endo = np.zeros(total_dd); e1_endo = np.zeros(total_dd); e2_endo = np.zeros(total_dd);
        
        
        # Open MTV fiber and sheet files
        #if( MPI.rank(mpi_comm_world()) == 0):
        if(MPI.comm_world.Get_rank() == 0): 
        
                if(mtvfiberfilename):
                        mtvfiberfile = open(mtvfiberfilename, 'w')
                        print >>mtvfiberfile, int(total_dd/3), 10
                
        
                if(mtvsheetfilename):
                        mtvsheetfile = open(mtvsheetfilename, 'w')
                        print >>mtvsheetfile, int(total_dd/3)
        
        
        for x_dof, y_dof, z_dof, scl in zip(x_dofs, y_dofs, z_dofs, scalar_dof):
        
        
                if(abs(rv[scl]) < 1e-9 and abs(rv[scl] + lv[scl]) < 1e-9):
                        ds[scl] = 0.0
                else:
                        ds[scl] = rv[scl]/(lv[scl]+rv[scl])
        
                
                pt = Point(np.array([Sdof_coordinates[scl][0], Sdof_coordinates[scl][1], Sdof_coordinates[scl][2]]))
                meshid = mesh.bounding_box_tree().compute_first_collision(pt)
        
                al_endo = lv_fiber_angle[0];  b_endo = lv_sheet_angle[0];
                al_epi = lv_fiber_angle[1];  b_epi = lv_sheet_angle[1];
        
                if matid is not None:
                        if(matid.array()[meshid] == septummatid):
                                al_endo = septum_fiber_angle[0];  b_endo = septum_sheet_angle[0];
                                al_epi = septum_fiber_angle[1];  b_epi = septum_sheet_angle[1];                 
                        elif(matid.array()[meshid] == rvmatid):
                                al_endo = rv_fiber_angle[0];  b_endo = rv_sheet_angle[0];
                                al_epi = rv_fiber_angle[1];  b_epi = rv_sheet_angle[1];                 
                        
        
                al_s[scl] = al_endo*(1 - ds[scl]) - al_endo*ds[scl]; b_s[scl] = b_endo*(1 - ds[scl]) - b_endo*ds[scl];
                al_w[scl] = al_endo*(1 - epi[scl]) + al_epi*epi[scl]; b_w[scl] = b_endo*(1 - epi[scl]) + b_epi*epi[scl];
        
                vec_dd[0] = dd[x_dof]; vec_dd[1] = dd[y_dof]; vec_dd[2] = dd[z_dof];
                vec_dd2[0] = dd2[x_dof]; vec_dd2[1] = dd2[y_dof]; vec_dd2[2] = dd2[z_dof];
                vec_dd3[0] = dd3[x_dof]; vec_dd3[1] = dd3[y_dof]; vec_dd3[2] = dd3[z_dof];
                vec_dd4[0] = dd4[x_dof]; vec_dd4[1] = dd4[y_dof]; vec_dd4[2] = dd4[z_dof];
                
        
                Qlv[c:c+3][0:3]= axisf(vec_dd4, -1.0*vec_dd2);
                Qlv[c:c+3][0:3] = orto3(Qlv[c:c+3][0:3], cnt); 
                Qlv[c:c+3][0:3]= orient(Qlv[c:c+3][0:3], al_s[scl], b_s[scl])
        
        
                e0_lv[x_dof] =  Qlv[c][0]; e0_lv[y_dof] =  Qlv[c+1][0]; e0_lv[z_dof] =  Qlv[c+2][0];
                e1_lv[x_dof] =  Qlv[c][1]; e1_lv[y_dof] =  Qlv[c+1][1]; e1_lv[z_dof] =  Qlv[c+2][1];
                e2_lv[x_dof] =  Qlv[c][2]; e2_lv[y_dof] =  Qlv[c+1][2]; e2_lv[z_dof] =  Qlv[c+2][2];
        
                Qrv[c:c+3][0:3] = axisf(vec_dd4, -1.0*vec_dd3);
                Qrv[c:c+3][0:3] = orto3(Qrv[c:c+3][0:3], cnt);
                Qrv[c:c+3][0:3]= orient(Qrv[c:c+3][0:3], -al_s[scl], -b_s[scl]);
                
                e0_rv[x_dof] =  Qrv[c][0]; e0_rv[y_dof] =  Qrv[c+1][0]; e0_rv[z_dof] =  Qrv[c+2][0];
                e1_rv[x_dof] =  Qrv[c][1]; e1_rv[y_dof] =  Qrv[c+1][1]; e1_rv[z_dof] =  Qrv[c+2][1];
                e2_rv[x_dof] =  Qrv[c][2]; e2_rv[y_dof] =  Qrv[c+1][2]; e2_rv[z_dof] =  Qrv[c+2][2];
        
        
                #Qendo[c:c+3][0:3]  = bislerp(Qlv[c:c+3][0:3], Qrv[c:c+3][0:3], ds[scl])
                Qendo[c:c+3][0:3]  = my_bislerp(Qlv[c:c+3][0:3], Qrv[c:c+3][0:3], ds[scl])
                
                e0_endo[x_dof] =  Qendo[c][0]; e0_endo[y_dof] =  Qendo[c+1][0]; e0_endo[z_dof] =  Qendo[c+2][0];
                e1_endo[x_dof] =  Qendo[c][1]; e1_endo[y_dof] =  Qendo[c+1][1]; e1_endo[z_dof] =  Qendo[c+2][1];
                e2_endo[x_dof] =  Qendo[c][2]; e2_endo[y_dof] =  Qendo[c+1][2]; e2_endo[z_dof] =  Qendo[c+2][2];
        
        
                Qepi[c:c+3][0:3] = axisf(vec_dd4, vec_dd);
                Qepi[c:c+3][0:3] = orto3(Qepi[c:c+3][0:3], cnt); 
                Qepi[c:c+3][0:3] = orient(Qepi[c:c+3][0:3], al_w[scl], b_w[scl]);
        
                e0_epi[x_dof] =  Qepi[c][0]; e0_epi[y_dof] =  Qepi[c+1][0]; e0_epi[z_dof] =  Qepi[c+2][0];
                e1_epi[x_dof] =  Qepi[c][1]; e1_epi[y_dof] =  Qepi[c+1][1]; e1_epi[z_dof] =  Qepi[c+2][1];
                e2_epi[x_dof] =  Qepi[c][2]; e2_epi[y_dof] =  Qepi[c+1][2]; e2_epi[z_dof] =  Qepi[c+2][2];
                
                #Qfiber[c:c+3][0:3]  = bislerp(Qendo[c:c+3][0:3], Qepi[c:c+3][0:3], epi[scl])
                Qfiber[c:c+3][0:3]  = my_bislerp(Qendo[c:c+3][0:3], Qepi[c:c+3][0:3], epi[scl])
        
                e0_fiber[x_dof] =  Qfiber[c][0]; e0_fiber[y_dof] =  Qfiber[c+1][0]; e0_fiber[z_dof] =  Qfiber[c+2][0];
                e1_fiber[x_dof] =  Qfiber[c][1]; e1_fiber[y_dof] =  Qfiber[c+1][1]; e1_fiber[z_dof] =  Qfiber[c+2][1];
                e2_fiber[x_dof] =  Qfiber[c][2]; e2_fiber[y_dof] =  Qfiber[c+1][2]; e2_fiber[z_dof] =  Qfiber[c+2][2];
        
        
                cnt = cnt + 1;
                c = c + 3;
        
                if(isrotatept):
                        points = [(-Sdof_coordinates[scl][2]+maxz)/10.0, Sdof_coordinates[scl][1]/10.0, Sdof_coordinates[scl][0]/10.0]
                        fvectors =  [-1.0*e0_fiber[z_dof], e0_fiber[y_dof], e0_fiber[x_dof]]
                        svectors =  [-1.0*e2_fiber[z_dof], e2_fiber[y_dof], e2_fiber[x_dof]]
        
                else:
                        points = [Sdof_coordinates[scl][0], Sdof_coordinates[scl][1], Sdof_coordinates[scl][2]]
                        fvectors =  [e0_fiber[x_dof], e0_fiber[y_dof], e0_fiber[z_dof]]
                        svectors =  [e2_fiber[x_dof], e2_fiber[y_dof], e2_fiber[z_dof]]
        
        
                #print >>mtvfiberfile, points[0], points[1], points[2], fvectors[0], fvectors[1], fvectors[2]
                #print >>mtvsheetfile, points[0], points[1], points[2], svectors[0], svectors[1], svectors[2]
        
        
        vtkoutfile = outdirectory+outfilename+"_e0_fiber"
        func_of_e0_vector.set_local(e0_fiber)
        func_of_e0_vector.apply("insert")
        vtk_py.convertQuadDataToVTK(mesh, V, func_of_e0, vtkoutfile)
        
        vtkoutfile = outdirectory+outfilename+"_e1_fiber"
        func_of_e1_vector.set_local(e1_fiber)
        func_of_e1_vector.apply("insert")
        vtk_py.convertQuadDataToVTK(mesh, V, func_of_e1, vtkoutfile)
        
        vtkoutfile = outdirectory+outfilename+"_e2_fiber"
        func_of_e2_vector.set_local(e2_fiber)
        func_of_e2_vector.apply("insert")
        vtk_py.convertQuadDataToVTK(mesh, V, func_of_e2, vtkoutfile)
        
        #if(MPI.rank(mpi_comm_world()) == 0): 
        if(MPI.comm_world.Get_rank() == 0): 
                print(isrotatept)
                print(outdirectory + outfilename)
                print("*******************************************************")
        
        if(isreturn):
        
                return func_of_e0, func_of_e1, func_of_e2

if (__name__ == "__main__"):

        parser = argparse.ArgumentParser()
        parser.add_argument('--xml_folder', type=str, required=True)
        parser.add_argument('--xml_meshfilename', type=str, required=True)
        parser.add_argument('--xml_facetfilename', type=str, required=True)
        parser.add_argument('--mtv_grid_directory', type=str, required=True)
        parser.add_argument('--mtv_basename', type=str, required=True)
        parser.add_argument('--isrotatept', type=bool, required=True)
        parser.add_argument('--al_endo', type=float, required=True)
        parser.add_argument('--al_epi', type=float, required=True)
        parser.add_argument('--b_endo', type=float, required=True)
        parser.add_argument('--b_epi', type=float, required=True)
        args = parser.parse_args()
        
        print("************* Entering SetBiVFiber.py *****************")
        
        xmlmeshfilename = os.path.join(args.xml_folder, args.xml_meshfilename)
        xmlfacetfilename = os.path.join(args.xml_folder, args.xml_facetfilename)
        outdirectory = args.mtv_grid_directory
        outfilename = args.mtv_basename
        isrotatept = args.isrotatept
        al_endo = args.al_endo
        al_epi = args.al_epi
        b_endo = args.b_endo
        b_epi = args.b_epi

        mesh = Mesh(xmlmeshfilename)
        boundaries = MeshFunction("size_t", mesh, xmlfacetfilename)
        mtvfiberfilename = outdirectory+outfilename+"_fiber_rotated.axis"
        mtvsheetfilename = outdirectory+outfilename+"_sheet_rotated.axis"
        
        SetBiVFiber(mesh, boundaries, mtvfiberfilename, mtvsheetfilename, al_endo, al_epi, b_endo, b_epi, isrotatept, outfilename, outdirectory)

