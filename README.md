# vanDANA

vanDANA is a high-performance FEM based Flow-thermal CFD solver utilizing the [FEniCS](https://fenicsproject.org/) library. 

<div align="center">
    <img src="/src/vanDANA.png" width="520px"> 
</div>

Future efforrs are oriented towards developing a FSI solver coupled with Heat transfer based on the Distibuted Langrange Multiplier based Fictiious Domain method. 

### Installation and Running the code

A singularity build file is provided that will install necessary libraries & setup the environment to run the code.

1. Install singularity by following the instruction in [here](https://docs.sylabs.io/guides/3.6/admin-guide/installation.html).

2. Build a singularity container using the [build file](./src/fenics_2019_dev) with
```
sudo singularity build <container_name>.img fenics_2019_dev
```

3. Once the container is built, you can launch the Singularity container by
```
 singularity run <container_name>.img
```

4. The code can be run within the singularity container by simply running :
```
python3 vanDANA.py
```
or in parallel
```
mpirun.mpich -np <# processors> python3 vanDANA.py
```

### Organization of the code

- [common](./common) module
- [utilities](./utilities) module
- [user_inputs](./user_inputs) module

### Authors

* Tejas Patel 
  
  (contact : patelte8@gmail.com)



