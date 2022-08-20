from dolfin import sym, nabla_grad, Identity

# Define fluid stress tensor
def sigma(Re, u, p):
    
    return (2/Re)*sym(nabla_grad(u)) - p*Identity(len(u))