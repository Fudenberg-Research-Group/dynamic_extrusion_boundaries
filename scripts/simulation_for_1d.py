#salloc --partition=debug --gres=gpu --mem-per-cpu=2GB --cpus-per-task=8

import numpy as np
import matplotlib.pylab as plt
import os
import h5py 
import time
import sys

from lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary

import funcs






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
#make a random config
#TAD_size = 50 # in monomers
#CTCF_left_positions = np.arange(0, sites_per_replica, TAD_size * sites_per_monomer)
#CTCF_right_positions = np.arange(1, sites_per_replica, TAD_size * sites_per_monomer)
ctcfrightlist = [284, 302, 867, 1005, 2185, 2526, 3760, 3945, 4530, 4986, 5570, 6041, 6183, 6621, 6752, 8084, 9752]
ctcfleftlist = [557, 2130, 2608, 2608, 2787, 2899, 3259, 3327, 3641, 3646, 4300, 4518, 5172, 5783, 7112, 7940, 8905]
CTCF_right_positions = np.array(ctcfrightlist)
CTCF_left_positions = np.array(ctcfleftlist)

########### 1d simulation parameters for lattice ###########
Trajn = 100 # trajectory length in monomer 
trajectory_length = Trajn * paramdict['sites_per_monomer'] #trajectory length in lattice land
pause_multiplier = 1/(1-pause)
trajectory_length = trajectory_length * pause_multiplier
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








