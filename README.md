# vanDANA

vanDANA is a highly efficient FEM Immersed Boundary (IB) based Flow-thermal FSI solver utilizing the [FEniCS](https://fenicsproject.org/) library (version 2019.2.0). 

<div align="center">
    <img src="/src/vanDANA.png" width="520px"> 
</div>

The FSI solver is based on the [Distibuted Langrange Multiplier based Fictitious Domain method](https://www.sciencedirect.com/science/article/pii/S0021999105000148) and is extended to deal with [heat transfer](https://www.sciencedirect.com/science/article/pii/S0021999106000167). The interpolation of variables is conducted using the smeared [delta-functions](https://www.sciencedirect.com/science/article/pii/S0021999109004136). Additionally, the flow solver is incompressible and has the option of choosing from various stabilization schemes : SUPG, PSPG and Crosswind; and the structure can be set as either incompressible/compressible.

<p align="center">
    <img src="/src/turek_benchmark.gif" width="600" height="310"/>
</p>
<p align="center">
  The classical Turek benchmark (FSI2)
</p>

## License

vanDANA (all versions) is licensed under the GNU GPL and is Copyright (2022) by the authors.

## Organization of the code

- [common](./common) module
- [utilities](./utilities) module
- [user_inputs](./user_inputs) module

## Installation and Running the code

A singularity build file is provided that will install necessary libraries & setup the environment to run the code.

1. Install singularity by following the instruction in [here](https://docs.sylabs.io/guides/3.6/admin-guide/installation.html).

2. Build a singularity container using the [build file](./src/fenics_2019_dev) by
```
sudo singularity build <container_name>.img fenics_2019_dev
```

3. Once the container is built, you can launch the Singularity container by
```
 singularity run <container_name>.img
```

4. The code can be run within the singularity container by simply running :
```
 mpirun.mpich -np <# processors> python3 vanDANA.py
```

## Documentation

For an introduction and tutorial to vanDANA code, please refer to the [documentation](https://deepnote.com/@research-2834/vanDANA-User-Manual-dcbd70e8-f8a8-4cc9-84ba-10cbae5aa5a5).


## Authors

- Tejas Patel (contact : patelte8@gmail.com, https://patelte8.github.io/cfd_gallery/)
- Kai Thin

## Contribute

Please report bugs and other issues through the issue tracker at: https://github.com/patelte8/vanDANA/issues. We welcome suggestions and if you wish to contribute/improve vanDANA solver, contact the author via mail.
