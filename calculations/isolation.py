import numpy as np
import warnings
import matplotlib.pyplot as plt
import pandas as pd
import glob
import ast
import chromoscores.maputils as chrmap
import chromoscores.snipping as chrsnip
import chromoscores.scorefunctions as chrscores
import seaborn as sns



#directory to the simulations
path_dict = {}
directory = '/scratch1/rahmanin/dynamical_boundary_simulations/simulation_md/simulation_with_pause_mod/sims_a/'
#directory = '/project/fudenber_735/polychrom/Dynamic_boundary_model/3d_sims_dynamic_boundary/layout_random/simulation_with_pause_mod/'
for name  in glob.glob(directory+'sims/folder_*'):
    path_dict[name.split('/sims/')[1][:]]= name
path_dict = dict(sorted(path_dict.items()))


#List of the position of boundary elements on the main diagonal. 
ctcfrightlist = [314, 579, 1195, 3717, 3772, 3921, 4451, 5193, 5723, 6302, 6574, 6779, 7000, 9232, 9310, 9861]
ctcfleftlist = [495, 865, 1404, 2164, 3143, 3615, 3971, 4069, 4480, 4938, 5300, 5587, 6401, 7725, 8764, 9619]
lst = np.sort((ctcfrightlist + ctcfleftlist))
lst = [int(lst[i]/10 ) for i in range(len(lst))]
lst = np.array(lst)
### list of boundary elements on all replications
rep = 3 
mon = 1000
site = 10
lst_t = []
for i in range(rep):
    lst_t += list(np.array(lst)+i*mon)
lst_t = np.array(lst_t)

### calculating the isolation score
delta_val = 1
diag_offset_val = 15
max_distance_val = 45


file=open('data/isoscore_for_dynamic_random_pause_mod_%s_%s_%s_traiangle.csv'%(delta_val,diag_offset_val,max_distance_val),'w')

file.write('life,velocity,clife,cof,sep,isoscore\n')
for name in list(path_dict.keys()):
    print(name)
    params=[ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
    face, back, clife, cof, life, slife, birth, pause, sep, site, mon, rep, step, vel = params
    
    #importing saved maps as matrices 
    data=np.load('../maps/maps_pause_mod/%s.npz' % name)
    mrcn=data['arr_0']
    mrc  = mrcn.astype(float)
    mrc /= np.median(np.diag(mrc,2))
    mat = np.log10(mrc)
    mrc_exp=chrmap.get_observed_over_expected(mrc)
    pile=chrmap.get_diagonal_pileup(mrc_exp, lst_t[1:75],96)   #A function to pile up snippets around investigated features (here boundary elements)

    #pseudocount=1/np.nansum(pile)

    pseudocount = np.nanmean(pile)/1000
    score=chrscores.isolation_score(pile,delta=delta_val,diag_offset=diag_offset_val,max_dist=max_distance_val,snippet_shapes='triangle',pseudo_count=pseudocount)
    print(score)
    file.write('%s,%s,%s,%s,%s,%s\n'%(life, vel, clife,cof, sep, score))

file.close()



