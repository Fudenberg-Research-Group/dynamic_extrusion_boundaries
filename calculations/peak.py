import pickle
import os
import time
import numpy as np
import ast

import polychrom

from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
import polychrom.contactmaps

import cooltools
import cooltools.lib.plotting

import pandas as pd
import warnings
import h5py 
import glob
import sys





#directory to the simulations
path_dict = {}
#directory = '/scratch1/rahmanin/dynamical_boundary_simulations/simulation_md/simulation_with_pause_mod/sims_a/'
directory = '/project/fudenber_735/polychrom/Dynamic_boundary_model/3d_sims_dynamic_boundary/layout_random/simulation_with_pause_mod/'
for name  in glob.glob(directory+'sims/folder_*'):
    path_dict[name.split('/sims/')[1][:]]= name
path_dict = dict(sorted(path_dict.items()))


#List of the position of boundary elements on the main diagonal. 
ctcfrightlist = [314, 579, 1195, 3717, 3772, 3921, 4451, 5193, 5723, 6302, 6574, 6779, 7000, 9232, 9310, 9861]
ctcfleftlist = [495, 865, 1404, 2164, 3143, 3615, 3971, 4069, 4480, 4938, 5300, 5587, 6401, 7725, 8764, 9619]
lst_a = np.sort((ctcfrightlist + ctcfleftlist))
lst_a = [int(lst_a[i]/10 ) for i in range(len(lst_a))]
new_lst = [i+1000 for i in lst_a]
lst_t = lst_a + new_lst
lst = np.sort(lst_t)
lst = np.unique(lst)

#making bins for distances
mindist = 5
band_edges = np.append([], 10+(4 * 1.3 ** np.arange(1,25)))
band_edges_array = band_edges.astype(int)


band_edge_list=list(band_edges_array)

for name in list(path_dict.keys()):

    params=[ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
    face, back, clife, cof, life, slife, birth, pause, sep, site, mon, rep, step, vel = params
    f=open('../data/peakscore/peakscore_%s.csv' % name,'w')
    f.write('chr_distance,peak_score\n')
    
    #importing saved maps as matrices 
    data = np.load('../maps/maps_pause_mod/%s.npz' % name)
    mrc = data['arr_0']
    mrc  = mrc.astype(float)
    mrc /= np.median(np.diag(mrc,2))
    mat = np.log10(mrc)
    pile = chrmap.get_offdiagonal_pileup_binlist(mrc,lst,band_edge_list,30)
    
    
    for i in range(len(band_edge_list)-1):
        c+=1
        pseudo_count=(1/np.sum(pile[i][1]))
        peakscore=chrsco.peak_score(pile[i][1],6,10,pseudo_count=pseudo_count)
        chr_dist=pile[i][0]*(2.5)
        f.write('%s,%s\n'%(chr_dist, peakscore))
    f.close()








