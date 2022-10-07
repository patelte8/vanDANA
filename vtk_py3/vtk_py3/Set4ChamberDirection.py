########################################################################
import argparse
import numpy as np
from dolfin import *
import math
########################################################################


def Set4ChamberDirection(mesh, outfilename, apexC, apexR, basalC, basalN, outdirectory="./"):


        def BCbase(x, on_boundary):
        	basalNorm = np.linalg.norm(basalN)
        	basalNN = 1.0/basalNorm*np.array(basalN)
        	
        	return (x[0] - basalC[0])*basalNN[0] + (x[1] - basalC[1])*basalNN[1]  < DOLFIN_EPS  and  on_boundary
        
        def BCapex(x, on_boundary):
        	return (x[0] - apexC[0])**2 + (x[1] - apexC[1])**2 < apexR[0]**2  and  on_boundary


        V = FunctionSpace(mesh, 'Lagrange', 1)
        u0 = Constant(1.0)
        u1 = Constant(0.0)
        
        bc =[DirichletBC(V, u1, BCapex), DirichletBC(V, u0, BCbase)] 
        
        # Define variational problem
        u = TrialFunction(V)
        v = TestFunction(V)
        f = Constant(0)
        a = inner(nabla_grad(u), nabla_grad(v))*dx
        L = f*v*dx
        
        # Compute solution
        u = Function(V)
        solve(a == L, u, bc)
        u_a = u.vector().array()
        
        # Compute gradient
        V_g = VectorFunctionSpace(mesh, 'Lagrange', 1)
        v = TestFunction(V_g)
        w = TrialFunction(V_g)
        
        a = inner(w, v)*dx
        L = inner(grad(u), v)*dx
        grad_u = Function(V_g)
        solve(a == L, grad_u)
        #normalize_grad_u = project(grad_u/sqrt(dot(grad_u, grad_u)), V_g)
        grad_u.rename('matdir', 'matdir')
        
        #plot(normalize_grad_u, interactive=True)
        pvdoutfile = outdirectory+outfilename+"_matdir.pvd"
        file15 = File(pvdoutfile)
        file15 << grad_u
        
        #plot(normalize_grad_u, interactive=True)
        pvdoutfile2 = outdirectory+outfilename+"_dirsoln.pvd"
        file16 = File(pvdoutfile2)
        file16 << u 


