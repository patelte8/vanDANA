from dolfin import sqrt, dot, Identity, \
                   nabla_grad, inner, outer, DOLFIN_EPS, conditional
from ufl import Max, ne

# SUPG/PSPG
def tau(alpha, vel, h, Num, dt):

    vnorm = dot(vel, vel)
    quant = Num*h*h 
    tau = alpha/sqrt((4.0/(dt*dt)) + (4*((vnorm)/(h*h))) + (144.0/(quant*quant)))
    return tau

def Pop(u, w):

    return dot(u, nabla_grad(w))

# Crosswind
def tau_cw(C_cw, fx, h, Num, Rs):

    vnorm = sqrt(inner(nabla_grad(fx), nabla_grad(fx)))
    res_norm = sqrt(dot(Rs, Rs))
    VX = C_cw - ((2*vnorm)/(h*Num*(res_norm + DOLFIN_EPS)))
    tau_cw = 0.5*Max(0, VX)*h*(res_norm/(vnorm + DOLFIN_EPS))
    return tau_cw

def D(u): 

    umag = dot(u, u)
    TX = conditional(ne(umag, 0), Identity(len(u)) - (outer(u, u)/(umag)), 0.0*Identity(len(u)))    # Projection orthogonal to streamline direction
    return TX      

def Pop_CW(u, fx):

    return dot(D(u), nabla_grad(fx))

# LSIC
def tau_lsic(Re):

    return 2/(3*Re)