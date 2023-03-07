# vanDANA

vanDANA is a highly efficient FEM Immersed Boundary (IB) based Flow-thermal FSI solver utilizing the [FEniCS](https://fenicsproject.org/) library (version 2019.2.0). The FSI solver is based on the [Distibuted Langrange Multiplier based Fictitious Domain method](https://www.sciencedirect.com/science/article/pii/S0021999105000148) and is extended to deal with [heat transfer](https://www.sciencedirect.com/science/article/pii/S0021999106000167). The interpolation of variables is conducted using the smeared [delta-functions](https://www.sciencedirect.com/science/article/pii/S0021999109004136). Additionally, the flow solver has the option of using various stabilization schemes : SUPG, PSPG, LSIC, Crosswind and Backflow; and the structure can be set as either incompressible/compressible.

<div align="center">
    <img src="/src/vanDANA.png" width="520px"> 
</div>
