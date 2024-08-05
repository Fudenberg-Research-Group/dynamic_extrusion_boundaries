# loop extrusion with dynamic boundaries

### description
This github page is designed for polymer simulation of chromatin loop extrusion by Cohesin with dynamic model of CTCF. For this simulaitons, polychrom library (https://github.com/open2c/polychrom) and OpenMM
https://github.com/openmm/openmm are used.
### usage
loops are implemented as bonds between monomers, which can be updated to connect more distal monomers, indicating loop extrusion. To indicate the monomers involve in the loop, first pairs of monomer indeces is determined by performing simulations on a one-dimensional lattice, which provide input data for coarse-grained molecular dynamic simulations. The one dimensional simulation is a python based code in 'script/simulations'. 
To use this repository, we first xxx
##### requirement
-Python
-Polychrom
-Openmm

### example