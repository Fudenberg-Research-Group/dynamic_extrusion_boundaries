# loop extrusion with dynamic boundaries

### description
This github page is designed for polymer simulation of chromatin loop extrusion by Cohesin with dynamic model of CTCF.
### usage
loops are implemented as bonds between monomers, which can be updated to connect more distal monomers, indicating chromatin loop extrusion. To indicate the monomers involve in the loop, first pairs of monomer indeces is determined by performing simulations on a one-dimensional lattice, which provide 'bonds' input data for coarse-grained molecular dynamic simulations. The one dimensional simulation is a python based code in 'script/simulations'.  
#### requirement
- Python
- Polychrom (https://github.com/open2c/polychrom)
- Openmm (https://github.com/openmm/openmm)

### example