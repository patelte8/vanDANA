The commands to mesh the pipe are as follows:

- gmsh -3 pipe.geo -format msh2
- gmsh -refine pipe.msh -format msh2
- dolfin-convert pipe.msh pipe.xml
- python3 formpi_gmsh.py