from dolfin import sym, nabla_grad, Identity, inv, det, ln
from .functions import *

# Define fluid stress tensor
def sigma(Re, u, p):
    
    return (2/Re)*sym(nabla_grad(u)) - p*Identity(len(u))

# Define solid stress : incompressible
def stress_inc(D_R, ps_R, Sm): 
    
    B = D_R
    return -1*ps_R*inv(F(B)) + Sm*(F(B).T - inv(F(B)))

# Define solid stress matrix : compressible
def stress_c(D_R, Ld, Sm): 
    
    B = D_R
    return Ld*ln(J(F(B)))*inv(F(B)) + Sm*(F(B).T - inv(F(B)))    

# Deformation gradient
def F(x):
    
    a = x.geometric_dimension()
    I = Identity(a)
    F = I + nabla_grad(x)
    return F

# Jacobian
def J(F):
    
    return det(F)    