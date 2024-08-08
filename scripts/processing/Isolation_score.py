### Import necessary libraries
import os
import ast
import glob
import h5py
import numpy as np
import warnings
import cooler

### Import specific functions for analysis
from scipy.optimize import fsolve
from scipy.ndimage import gaussian_filter1d
from polychrom.polymer_analyses import slope_contact_scaling

### Import chromoscore library functions for peak score calculation
import chromoscores.maputils as chrmap
import chromoscores.scorefunctions as chrsco
import chromoscores.snipping as chrsn

### Define the directory containing simulation maps
directory = '/home1/rahmanin/start/polychrom/projects/Dynamic_boundary_elements/analysis/maps/maps_points_cools_eq/'

### Create a dictionary to store paths to simulation files
path_dict = {}
for name in glob.glob(directory + 'folder_*'):
    path_dict[name.split('_eq/')[1].split('.mcool')[0]] = name
path_dict = dict(sorted(path_dict.items()))
print('The number of simulations is %s' % len(path_dict))

### Define boundary element positions on the main diagonal
ctcfrightlist = [314, 579, 1195, 3717, 3772, 3921, 4451, 5193, 5723, 6302, 6574, 6779, 7000, 9232, 9310, 9861]
ctcfleftlist = [495, 865, 1404, 2164, 3143, 3615, 3971, 4069, 4480, 4938, 5300, 5587, 6401, 7725, 8764, 9619]

### Calculate boundary element list for 10 kb resolution
lst = np.sort(ctcfrightlist + ctcfleftlist)
lst = [int(lst[i] / (10 * 4)) for i in range(len(lst))]
lst = np.array(lst)

### Create list of boundary elements for all replications
rep = 3  # Number of replications
mon = 1000 // 4  # Length of the monomer
lst_t = []
for i in range(rep):
    lst_t += list(np.array(lst) + i * mon)
lst_t = np.array(lst_t)

### calculating pile-up isolation score
delta_val = 1
diag_offset_val = 5
max_distance_val = 30


file = open('../../Outputs/isolation_score_traiangles_coolres10.csv','w')

file.write('life,velocity,clife,cof,isoscore\n')


for name in list(path_dict.keys())[:]:
    print(name)
    params=[ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
    face, back,clife, cof, life, slife, birth, pause, sep, site, mon, rep, step, vel = params

    #importing saved maps as matrices 
    cool_uri = (directory+'%s.mcool'%name)
    region = 'chr_sim:0-7500000' 
    cname = (cooler.Cooler(cool_uri+'.10000.cool').matrix(balance=False).fetch(region))
    mrcn = cname
    mrc  = mrcn.astype(float)
    mat = np.log10(mrc)
    mrc_exp = chrmap.get_observed_over_expected(mrc)
    pile = chrmap.get_diagonal_pileup(mrc_exp, lst_t[6:91],60)   #A function to pile up snippets around investigated features (here boundary elements)

    clr = cooler.Cooler(cool_uri+'.10000.cool')
    pseudocount = 0
    score = chrsco.isolation_score(pile,delta=delta_val,diag_offset=diag_offset_val,max_dist=max_distance_val,snippet_shapes='triangle',pseudo_count=pseudocount)
    file.write('%s,%s,%s,%s,%s\n'%(life, vel,clife,cof, score))
    
file.close()