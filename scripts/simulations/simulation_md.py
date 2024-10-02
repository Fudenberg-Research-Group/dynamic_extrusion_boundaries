#salloc --partition=debug --gres=gpu --mem-per-cpu=2GB --cpus-per-task=8

import numpy as np
import os
import h5py 
import time
import sys

from lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary

import funcs


# for MD simulations
import polychrom
from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
from polychrom.lib.extrusion import  bondUpdater
import warnings
import ast





filename = sys.argv[-1]

print('this is file name %s'%filename)

params = [ast.literal_eval(i) for i in filename.split('folder_')[1].split('_')[1::2]]
face, back, Clife, Cof, life, slife, birth, pause, sep, site, monomer, replica, steps, vel = params

paramdict={
            'CTCF_facestall':[face],
            'CTCF_backstall':[back],
            'CTCF_lifetime':[Clife],
            'CTCF_offtime':[Cof],
            'LEF_lifetime':[life],
            'LEF_stalled_lifetime':[slife],
            'LEF_birth':[birth],
            'LEF_pause':[pause],
            'LEF_separation':sep,
            'sites_per_monomer':site,
            'monomers_per_replica':monomer,
            'number_of_replica':replica,
            'steps':steps,
            'velocity_multiplier':vel
            }

paramdict_keys={
                'CTCF_facestall':'face',
                'CTCF_backstall':'back',
                'CTCF_lifetime':'Clife',
                'CTCF_offtime':'Cof',
                'LEF_lifetime':'life',
                'LEF_stalled_lifetime':'slife',
                'LEF_birth':'birth',
                'LEF_pause':'pause',
                'LEF_separation':'sep',
                'sites_per_monomer':'site',
                'monomers_per_replica':'monomer',
                'number_of_replica':'replica',
                'steps':'steps',
                'velocity_multiplier':'vel'
                }

file_name = funcs.paramdict_to_filename(paramdict, paramdict_keys)
folder_name = '/sims/'+'folder_' + file_name.split('file_')[1]
folder = os.getcwd() + folder_name
if os.path.exists(folder):
    print("already exist")
else:
    os.mkdir(folder)

monomers_per_replica = paramdict['monomers_per_replica'] 
sites_per_monomer = paramdict['sites_per_monomer']
sites_per_replica = monomers_per_replica * sites_per_monomer

typedict = {'A': 0}
monomer_types = np.zeros(monomers_per_replica, dtype=int)
site_types = np.repeat(monomer_types, sites_per_monomer)

# Create some CTCF boundary sites
TAD_size = 50 # in monomers
CTCF_left_positions = np.arange(0, sites_per_replica, TAD_size * sites_per_monomer)
CTCF_right_positions = np.arange(1, sites_per_replica, TAD_size * sites_per_monomer)

########### 1d simulation parameters for lattice ###########
Trajn = 7500 # trajectory length in monomer 
trajectory_length = Trajn * paramdict['sites_per_monomer'] #trajectory length in lattice land
num_dummy_steps = trajectory_length // 5 #dummy steps in lattice land
blocksteps = 5 
bins = np.linspace(0, trajectory_length, blocksteps, dtype=int)
N = (paramdict['monomers_per_replica']*paramdict['number_of_replica'])
LEFNum = N // paramdict['LEF_separation']



translocator = funcs.make_translocator(LEFTranslocatorDynamicBoundary, 
                                 site_types,
                                 CTCF_left_positions,
                                 CTCF_right_positions, 
                                 **paramdict)

with h5py.File(folder+"/LEFPositions.h5", mode='w') as myfile:
    dset = myfile.create_dataset("positions", 
                                 shape=(trajectory_length, LEFNum, 2), #edited
                                 dtype=np.int32, 
                                 compression="gzip")
    
    translocator.steps(0)
    
    for st, end in zip(bins[:-1], bins[1:]):
        cur = []
        for i in range(st, end):
            translocator.step()        
            cur.append(translocator.LEFs.copy())
        cur = np.array(cur)
        dset[st:end] = cur
    myfile.attrs["N"] = N * paramdict['sites_per_monomer']
    myfile.attrs["LEFNum"] = LEFNum

### Molecular dynamics simulaiton ###
myfile = h5py.File(folder + "/LEFPositions.h5", mode='r')
sites_per_monomer = paramdict['sites_per_monomer']
N = myfile.attrs["N"] // sites_per_monomer
print(N)
LEFNum = myfile.attrs["LEFNum"]
LEFpositions = myfile["positions"][::sites_per_monomer]// sites_per_monomer
Nframes = LEFpositions.shape[0]

# Md simulation characteristics
stiff = 1
dens = 0.2
box = (N / dens) ** 0.33  # density = 0.1.

smcStepsPerBlock = 1  # now doing 1 SMC step per block 
# initialize positions
data = grow_cubic(N, int(box) - 2)  # creates a compact conformation 
block = 0  # starting block 
steps= paramdict['steps']



# new parameters because some things changed 
saveEveryBlocks = 10   # save every 10 blocks (saving every block is now too much almost)
restartSimulationEveryBlocks = 100

# parameters for smc bonds
smcBondWiggleDist = 0.2
smcBondDist = 0.5

# assertions for easy managing code below 
assert (Nframes % restartSimulationEveryBlocks) == 0 
assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0

savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
simInitsTotal  = (Nframes) // restartSimulationEveryBlocks 


tstp = 70 # timestep for integrator in fs
tmst = 0.01 # thermostat for integrator

milker = polychrom.lib.extrusion.bondUpdater(LEFpositions)

reporter = HDF5Reporter(folder=folder, max_data_length=100, overwrite=True, blocks_only=False)

for iteration in range(simInitsTotal):
    # simulation parameters are defined below 
    a = Simulation(
            platform="cuda",
            integrator='langevin',  timestep=tstp, collision_rate=tmst,
            error_tol=0.01,  
            GPU="0",
            N = len(data),
            reporters=[reporter],
            PBCbox=[box, box, box],
            precision="mixed")  # timestep not necessary for variableLangevin
    ############################## New code ##############################
    a.set_data(data)  # loads a polymer, puts a center of mass at zero

    a.add_force(
        forcekits.polymer_chains(
            a,
            chains=[(0, None, 0)],

                # By default the library assumes you have one polymer chain
                # If you want to make it a ring, or more than one chain, use self.setChains
                # self.setChains([(0,50,1),(50,None,0)]) will set a 50-monomer ring and a chain from monomer 50 to the end

            bond_force_func=forces.harmonic_bonds,
            bond_force_kwargs={
                'bondLength':1.0,
                'bondWiggleDistance':0.1, # Bond distance will fluctuate +- 0.05 on average
             },

            angle_force_func=forces.angle_force,
            angle_force_kwargs={
                'k':1.5
                # K is more or less arbitrary, k=4 corresponds to presistence length of 4,
                # k=1.5 is recommended to make polymer realistically flexible; k=8 is very stiff
            },

            nonbonded_force_func=forces.polynomial_repulsive,
            nonbonded_force_kwargs={
                'trunc':1.5, # this will let chains cross sometimes
                'radiusMult':1.05, # this is from old code
                #'trunc':10.0, # this will resolve chain crossings and will not let chain cross anymore
            },
            except_bonds=True,
    ))
    # ------------ initializing milker; adding bonds ---------
    # copied from addBond
    kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
    bondDist = smcBondDist * a.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)

    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(bondForce=a.force_dict['harmonic_bonds'],
                blocks=restartSimulationEveryBlocks)
    #print(milker.allBonds[0])
    for t,l in enumerate(milker.allBonds):
        for b in l:
            if (b[0] == 11296) or (b[1] == 11296):
                print(t,b)
    # If your simulation does not start, consider using energy minimization below
    if iteration==0:
        a.local_energy_minimization() 
    else:
        a._apply_forces()

    for i in range(restartSimulationEveryBlocks):        
       # print("restart#",i)
        if i % saveEveryBlocks == (saveEveryBlocks - 1):  
            a.do_block(steps=steps)
        else:
            a.integrator.step(steps)  # do steps without getting the positions from the GPU (faster)
        if i < restartSimulationEveryBlocks - 1: 
            curBonds, pastBonds = milker.step(a.context)  # this updates bonds. You can do something with bonds here
    data = a.get_data()  # save data and step, and delete the simulation
    del a

    reporter.blocks_only = True  # Write output hdf5-files only for blocks

    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

reporter.dump_data()

myfile.close()








