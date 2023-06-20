from dolfin import sym, nabla_grad, Identity, inv, det, ln
from .functions import *

# Deformation gradient
def F(x):
    
    F = Identity(x.geometric_dimension()) + nabla_grad(x)
    return F

# Jacobian
def J(F):
    
    return det(F)   

# Cauchy's strain tensor
def E(B, dim):
    
    E = B + B.T
    E *= 0.5
    E -= Identity(dim)
    return E

# Define fluid stress tensor
# ----------------------------------------------------------

def sigma(Re, u, p):
    
    return (2/Re)*sym(nabla_grad(u)) - p*Identity(len(u))

# Linearly elastic material
# ----------------------------------------------------------

# Define solid stress : incompressible
def stress_lr_elastic_inc(D_R, ps_R, Sm): 
    
    B = F(D_R)
    E = E(B, D_R.geometric_dimension())

    return (-1*ps_R + 2*Sm*E)*inv(B)  

# Define solid stress : compressible
def stress_lr_elastic_c(D_R, nw, Sm): 
    
    B = F(D_R)
    E = E(B, D_R.geometric_dimension())
    lame1 = (2*Sm*nw)/(1-(2*nw))

    return (lame1*tr(E) + 2*Sm*E)*inv(B)

# Neohookean material
# ----------------------------------------------------------

# Define solid stress : incompressible
def stress_inc(D_R, ps_R, Sm): 
    
    B = F(D_R)
    return -1*(ps_R + Sm)*inv(B) + Sm*B.T 

# Define solid stress : compressible
def stress_c(D_R, Ld, Sm): 
    
    B = F(D_R)
    return (Ld*ln(J(B)) - Sm)*inv(B) + Sm*B.T    
