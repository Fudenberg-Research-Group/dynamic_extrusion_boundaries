import numpy as np
import os
import h5py 
import time
import sys

from lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary #import from ./engines

import funcs #import from ./utils




filename = sys.argv[-1]

print('file name: %s'%filename)

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

# Create CTCF boundary sites
CTCF_sites_right = np.array([284, 302, 867, 1005, 2185, 2526, 3760, 3945, 4530, 4986, 5570, 6041, 6183, 6621, 6752, 8084, 9752])
CTCF_sites_left = np.array([557, 2130, 2608, 2608, 2787, 2899, 3259, 3327, 3641, 3646, 4300, 4518, 5172, 5783, 7112, 7940, 8905])

########### 1d simulation parameters for lattice ###########
Trajn = 100 # trajectory length in monomer 
trajectory_length = Trajn * paramdict['sites_per_monomer'] #trajectory length in lattice land
pause_multiplier = 1/(1-pause)
trajectory_length = trajectory_length * pause_multiplier
num_dummy_steps = int(trajectory_length // 5) #dummy steps in lattice land for loops equilibration
blocksteps = 5 
bins = np.linspace(0, trajectory_length, blocksteps, dtype=int)
N = (paramdict['monomers_per_replica']*paramdict['number_of_replica'])
LEFNum = N // paramdict['LEF_separation']



translocator = funcs.make_translocator(LEFTranslocatorDynamicBoundary, 
                                 site_types,
                                 CTCF_sites_left,
                                 CTCF_sites_right, 
                                 **paramdict)

with h5py.File(folder+"/LEFPositions.h5", mode = 'w') as myfile:
    # creating data set for LEF positions
    dset_loop_positions = myfile.create_dataset("positions", 
                                 shape = (trajectory_length, LEFNum, 2), #edited
                                 dtype = np.int32, 
                                 compression = "gzip")

    # creating data sets for boundary elements possible sites
    dset_ctcf_sites_right = myfile.create_dataset("CTCF_sites_right",
                                                 shape = (len(CTCF_sites_right)), 
                                                 compression = "gzip", 
                                                 data=CTCF_sites_right.copy())


    dset_ctcf_sites_left = myfile.create_dataset("CTCF_sites_left",
                                                shape = len(CTCF_sites_left), 
                                                compression="gzip",
                                                data=CTCF_sites_left.copy())

    # creating data sets for boundary elements positions
    dset_ctcf_positions_right = myfile.create_dataset("CTCF_positions_right",
                                      shape = (trajectory_length, len(CTCF_sites_right), 1), 
                                     dtype = np.bool, 
                                     compression = "gzip")
    dset_ctcf_positions_left = myfile.create_dataset("CTCF_positions_left",
                                     shape = (trajectory_length, len(CTCF_sites_left), 1), 
                                     dtype = np.bool, 
                                     compression = "gzip")

    translocator.steps(num_dummy_steps)
    
    for st, end in zip(bins[:-1], bins[1:]):
        loop_positions = []
        ctcf_right_cur= []
        ctcf_left_cur = []
        for i in range(st, end):
            translocator.step()        

            loop_positions.append(translocator.LEFs.copy())

            ctcf_positions_right = (translocator.stallProbRight)[CTCF_sites_right]*1
            ctcf_positions_left = (translocator.stallProbLeft)[CTCF_sites_left]*1
            ctcf_right_cur.append(ctcf_positions_right.reshape(len(ctcf_positions_right),1))
            ctcf_left_cur.append(ctcf_positions_left.reshape(len(ctcf_positions_left),1))

        loop_positions = np.array(loop_positions)
        ctcf_right_cur = np.array(ctcf_right_cur)
        ctcf_left_cur = np.array(ctcf_left_cur)
        
        dset_loop_positions[st:end] = loop_positions

        dset_ctcf_positions_right[st:end] = ctcf_right_cur
        dset_ctcf_positions_left[st:end] = ctcf_left_cur

    myfile.attrs["N"] = N * paramdict['sites_per_monomer']
    myfile.attrs["LEFNum"] = LEFNum








