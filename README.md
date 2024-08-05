# loop extrusion with dynamic boundaries

### description
This github page is designed for polymer simulation of chromatin loop extrusion by Cohesin with dynamic model of CTCF.

#### requirement
- Python
- Polychrom (https://github.com/open2c/polychrom)
- Openmm (https://github.com/openmm/openmm)

### usage
#### running simulaitons 
loops are implemented as bonds between monomers, which can be updated to connect more distal monomers, indicating chromatin loop extrusion. To indicate the monomers involve in the loop, first pairs of monomer indeces is determined by performing simulations on a one-dimensional lattice, which provide 'bonds' input data for coarse-grained molecular dynamic simulations. The one dimensional simulation is a python based code in 'script/simulations'. The parameters of loop extrusion, including the lifetime, velocity, and density of extruders (Cohesins), as well as MD parameters, can be modified in ... file. 
#### processing simulation data
The trajectory data from 1d or md simulations can be processed to provide virtual chip-seq or contact maps. The python codes for these process is provided in 'script/processing'. 

#### analysis
From processed data, the analysis for quantifying features such as FRiP, TADs, peaks, and vermicelli is provided at 'analysis' as jupyter notebooks.  


