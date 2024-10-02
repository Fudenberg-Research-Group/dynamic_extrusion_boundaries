import glob
import os
import time
import numpy as np
import ast

import polychrom

from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
import polychrom.contactmaps

import cooltools
import cooltools.lib.plotting
import sys


path_dict = {}

directory = '/scratch1/rahmanin/dynamical_boundary_simulations/simulation_md/simulation_with_pause_mod/simulations_with_hsteps/'
for name  in glob.glob(directory+'sims/folder_*'):
    path_dict[name.split('/sims/')[1][:]]= name
path_dict = dict(sorted(path_dict.items()))


monomer_per_replica = 1000

mapN = 3 * monomer_per_replica #number of monomer for averaging
total_monomers = 10000
mapstarts = (np.arange(0,8000 , monomer_per_replica))
min_time = 0
freq = 1
refresh = True
if refresh == True:
    map_dict_eq = {}
i = 1
for name in list(path_dict.keys())[:]:
    i += 1
    if os.path.isfile('../maps/maps_pause_mod_hsteps/%s.npz' % name)==True: continue
    params = [ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
    face, back, clife, cof, life, slife, birth, pause, sep, site, mon, rep, step, vel = params
    URIs = polychrom.hdf5_format.list_URIs(path_dict[name])
    URIs_eq = np.array(URIs)[np.array([int(i.split("::")[-1]) for i in URIs]) > min_time][::freq]
    mrc = polychrom.contactmaps.monomerResolutionContactMapSubchains(
        URIs_eq,
        mapstarts,
        mapN,
        cutoff = 2.3,
        n = 8)
    map_dict_eq[name] = mrc
    np.savez_compressed('../maps/maps_pause_mod_hsteps/%s.npz' % name, mrc)

    # To save in cooler format 
    cool_uri = ('inputs/maps/maps_points_cools_eq/%s.mcool'%name)
    mrc_new = cms.coolify(mrc,
            cool_uri,
            chrom_dict = {},
            binsize=2500,
            chunksize = 10000000)
    
    clr = cooler.Cooler(cool_uri+'.2500.cool')
    # to convert format to a new resolution
    base_uri =cool_uri+'.2500.cool'
    output_uri = cool_uri+'.10000.cool'
    factor = 4
    chunksize = 10000000
    clr_10 = cooler.coarsen_cooler(base_uri, output_uri, factor, chunksize)





