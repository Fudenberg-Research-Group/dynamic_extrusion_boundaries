import glob
import os
import time
import numpy as np
import ast


import cooltools
import cooltools.lib.plotting

import pandas as pd
import warnings
import h5py 
import sys


path_dict = {}

directory='./inputs/sims/'

for fname  in glob.glob(directory+'folder*'):
    path_dict[fname.split('sims/')[1][:]]= fname
path_dict = dict(sorted(path_dict.items()))


CTCF_sites_right = [314, 579, 1195, 3717, 3772, 3921, 4451, 5193, 5723, 6302, 6574, 6779, 7000, 9232, 9310, 9861]
CTCF_sites_left = [495, 865, 1404, 2164, 3143, 3615, 3971, 4069, 4480, 4938, 5300, 5587, 6401, 7725, 8764, 9619]
lst = np.array(CTCF_sites_right + CTCF_sites_left)


### list of boundary elements on all replications

rep = 10 
mon = 1000
site = 10
lst_t = []
for i in range(rep):
    lst_t += list(np.array(lst)+ i * mon * site)


def FRiP(num_sites_t, lef_positions, peak_positions ):
   """
   Calculate the Fraction of Reads in Peaks (FRiP) score.

   Args:
       num_sites_t (int): total number of genomic sites (lattice sites)
       lef_positions (np.array): positions of loop extrusion factors legs.
       peak_positions (np.array): Indices of positions corresponding to peaks

    Returns:
       float: The FRiP score, which is the fraction of LEF positions that fall within the peak regions.
   
   """
    
    hist, edges = np.histogram(  lef_positions  , np.arange(num_sites_t + 1) )
    return np.sum(hist[peak_positions] )/len(lef_positions)

def peak_positions(boundary_list, window_sizes=[1]):
    """
    Calculate peak positions based on a boundary_list within window_sizes.

    Args:
        boundary_list (list): List of boundary values.
        window_sizes (list, optional): List of window sizes. Defaults to [1].

    Returns:
        np.ndarray: Array containing peak positions.
    """
    peak_monomers = np.array([])

    for i in window_sizes:
        inds_to_add = [boundary + i for boundary in boundary_list]
        peak_monomers = np.hstack((peak_monomers, inds_to_add))

    return peak_monomers.astype(int)



window_size = 1
min_time = -500000
file = open('./data/fripscore.csv', 'w')
file.write('lifetime, velocity, clife, cof,sep, fripscore\n')
i = 1
for name in list(path_dict.keys())[:]:
    params=[ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
    face, back, clife, cof, life, slife, birth, pause, sep, site, mon, rep, step, vel = params

    i += 1
    mapN=mon*site
    lefs = h5py.File(path_dict[name]+'/LEFPositions.h5','r')["positions"]

    lef_lefts = lefs[min_time:,:,0].flatten()
    lef_rights = lefs[min_time:,:,1].flatten()
    lef_positions = np.hstack((lef_lefts,lef_rights))


    peak_monomers = peak_positions(lst_t,window_sizes=np.arange(-window_size,(window_size)+1) )
    frip = FRiP(mapN * rep, lef_positions, peak_monomers)
    score = FRiP(mapN * rep, lef_positions, peak_monomers)
    file.write('%s,%s,%s,%s,%s,%s\n'%(life, vel, clife, cof, sep, score))
file.close()
