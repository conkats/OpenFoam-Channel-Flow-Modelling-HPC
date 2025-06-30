# OpenFoam-CFD modelling
## Channel-Flow-SGE-UoM-HPC
- The main objective of this repository is to gain some experience in developing code, implementing methods and using an open-source CFD tool, `OpenFOAM` v9 and v2.1.1 distributed by OpenFOAM Foundation. The simulations concern a simple channel flow and templates are generated to run the turbulent flow simulations and post-process them on the HPC clusters of the University of Manchester, configured with sun-grid-enginee (`SGE`) and `SLURM`. 

- OpenFoam is written in C++, a high-level object-oriented programming language. For the purposes of this tutorial, the fully turbulent channel flow of Retau=395 will be simulated as 3-D since it forms the basis of the Large Eddy Simulation tutorial provided by OpenFoam v9 and also an 1-D computation is used to compared the results with published DNS data. The models include the standard k-epsilon model, RSM and an implementation of a pressure gradient solver.


### Reference:
Abe, H., Kawamura, H. and  Matsuo, Y., "Direct numerical simulation of a fully developed turbulent channel flow with respect to Reynolds number dependence," ASME J. Fluids Eng., vol. 123, pp. 382-393, 2001.
