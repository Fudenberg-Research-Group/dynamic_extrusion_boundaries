# Chromatin loop extrusion with dynamic boundaries
![Alt Text](./output/cartoons/output.gif)
%![Alt Text](./output/cartoons/output_mult_seq_bar.gif)
### description
This GitHub repository contains tools for simulating chromatin loop extrusion using a dynamic model of CTCF to better understand how binding of CTCF contributes to chromatin organization and regulation.

### structure of the repository
The structure of this repository follows as below:
- inputs : this contains simulations and experimental data.
- output : this contain files after processing and analyzing the input data.
- analysis: this contains the notebooks and codes for analyzing simulations and experimenal data.
- script: this contains files required for performing simulations.
  
### requirement
- Python\Python libraries
- Polychrom: A toolkit for polymer simulations. (https://github.com/open2c/polychrom)
- Openmm: A library for molecular simulations. (https://github.com/openmm/openmm)
- Open2c analysis packages. (https://github.com/open2c)


### usage
#### running simulaitons 
1. One-Dimensional Lattice Simulation:
Before running the full molecular dynamics simulations, implement loops on a one-dimensional lattice to determine pairs of monomer indices involved in loop extrusion. This step helps define the harmonic 'bonds' input data for the coarse-grained molecular dynamics simulations. Parameters such as the lifetime, velocity, and density of extruders (Cohesins) can be adjusted in the configuration file located at script/simulations/config/simulation_params.py. For CTCF barriers, you can determine dynamic or static barriers, and adjust their parameters such as binding time, unbound time, and/or the occupancy. You can find the Python script for this simulation in the script/simulations directory.

2. Run the coarse-grained molecular dynamics simulations to model loop extrusion. 

#### processing simulation data
After running the simulations, process the simulated trajectory data to generate virtual ChIP-seq profiles, contact maps, and images. The scripts for data processing are available in the script/processing directory. For any processing code the instruction is added to the relevant python code.
##### example 
To generate contact maps: 
'''
python script/processing/generate_maps.py
'''

#### analysis
Once the data is processed, perform analyses to quantify observable features such as:

- FRiP (Fraction of Reads in Peaks)
- TADs (Topologically Associating Domains)
- Dots (loops between barriers)
- Vermicelli: Analysis of accumulation of extruders on axial structures.
- P(s): contact frequency vs genomic distance
  
The analysis scripts are provided as Jupyter notebooks in the analysis directory.
Each notebook includes detailed instructions and examples to guide you through the analysis process.




