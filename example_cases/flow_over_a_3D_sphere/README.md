The drag and lift coefficients are defined as follows :

- Cd = assemble(dot(traction, self.nx)*ds(4))/(0.5*(PI)/4)
- CL = assemble(dot(traction, self.ny)*ds(4))/(0.5*(PI)/4)

Nusselt number :

- area_sphere = assemble(1*self.ds(4))
- Nu = assemble(dot(nabla_grad(T_[0]), n)*ds(4))/self.area_sphere