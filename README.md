# dendriticGrowth 

### About
dendriticGrowth is a C++ code for C++ code developed to model complex dendrite observed in the solidification of the pure metal and binary alloy using the phase-field models and finite-element method. The master branch has a finite-element code which models heat diffusion and the diffusion of the phase-field order parameter in the evolution of a single equiaxed dendrite of a pure melt. The other branches has a code for different types of dendrites such as single equiaxed, single columnar and multi-columnar dendrites.  
The salient features of the computational implementation are: high-fidelity MPI based parallel implementation, support for direct and iterative solvers with preconditioning, fully implicit time-stepping schemes. 

### Installation:

dendriticGrowth code builds on top of the deal.II library.

1) Install CMake, PETSc, Trilinos, SLEPc, p4est, and deal.II (version 9.3.0 recommended)<br>

2) Clone the neuronalActionPotential GitHub repository <br>
```
$ git clone git clone https://github.com/cmmg/dendriticGrowth.git
$ cd dendriticGrowth
$ git checkout master
$ cmake .
$ make -j nprocs
  ```
[here nprocs denotes the number of processors]

### Visualization:

  Output of the primary fields is in the vtk format (parallel:*.pvtu, serial:*.vtu files). These can be visualized with the following open source applications:
  1. VisIt (https://visit-dav.github.io/visit-website/releases-as-tables/)
  2. Paraview (http://www.paraview.org/download/)


License
-------
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.
