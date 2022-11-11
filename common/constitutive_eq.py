from dolfin import sym, nabla_grad, Identity, inv, det, ln
from .functions import *

# Define fluid stress tensor
def sigma(Re, u, p):
    
    return (2/Re)*sym(nabla_grad(u)) - p*Identity(len(u))

# Define solid stress : incompressible
def stress_inc(D_R, ps_R, Sm): 
    
    B = F(D_R)
    return -1*(ps_R + Sm)*inv(B) + Sm*B.T 

# Define solid stress matrix : compressible
def stress_c(D_R, Ld, Sm): 
    
    B = F(D_R)
    return (Ld*ln(J(B)) - Sm)*inv(B) + Sm*B.T    

# Deformation gradient
def F(x):
    
    F = Identity(x.geometric_dimension()) + nabla_grad(x)
    return F

# Jacobian
def J(F):
    
    return det(F)    
