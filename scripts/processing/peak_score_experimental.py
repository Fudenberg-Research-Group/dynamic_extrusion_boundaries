import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 96
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# import libraries for biological data analysis
from coolpuppy import coolpup
from coolpuppy.lib import numutils
from coolpuppy.lib.puputils import divide_pups
from coolpuppy import plotpup
import cooler
import bioframe
import cooltools
from cooltools import expected_cis, expected_trans
from cooltools.lib import plotting

ctcf = bioframe.read_table('/project/fudenber_735/yxiao977/frip_sm_data/ChIP_fastqs_maps/Justice_2020_mm10/SRR10099910/SRR10099910.q30.mm10.sort_peaks.narrowPeak', schema='bed')

ctcf['mid']=(ctcf.end+ctcf.start)/2
ctcf['quartile_score']=pd.qcut(ctcf['score'],4, labels=False) + 1
ctcf_quart=ctcf[ctcf['quartile_score']==4]

direct='/project/fudenber_735/motifs/mm10/jaspar/MA0139.1.tsv.gz'
motif=bioframe.read_table(direct)
motif.head()
motif=motif.rename(columns={0: 'chrom', 1: 'start',2:'end',3:'name',4:'score',5:'pval',6:'strand'})


peaks_motifs = bioframe.overlap(ctcf,motif,how='inner')
peaks_motifs = peaks_motifs[(peaks_motifs['chrom']!= 'chrX')&(peaks_motifs['chrom']!= 'chrY')]

peaks_motifs=peaks_motifs.rename(columns={'strand':'strand_','strand_':'strand'})
peaks_motifs=peaks_motifs[['chrom','start','end','name','score','mid','quartile_score','strand']]

sites = peaks_motifs
peaks_motifs_bin = bioframe.cluster(peaks_motifs, min_dist=10000)\
    .drop_duplicates('cluster')\
    .reset_index(drop=True)

peaks_motifs_cluster = bioframe.cluster(peaks_motifs, min_dist=10000)#\

# Group by 'chuster_start' and 'cluster_end' and get the index of rows with the max 'score' within each group
idx = peaks_motifs_cluster.groupby(['cluster_start', 'cluster_end'])['score'].idxmax()

# Filter the dataframe to keep only the rows with the max 'score' in each group
peaks_motifs_bin_s = peaks_motifs_cluster.loc[idx].reset_index(drop=True)
print('length of bin s is %s'%len(peaks_motifs_bin_s))


res = 10_000
bonev_file = '/project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/HiC_ES.mm10.mapq_30.1000.mcool'
bonev_cooler = cooler.Cooler(bonev_file+'::resolutions/'+str(res))
view_df_bonev = cooltools.lib.make_cooler_view(bonev_cooler)[:19]

prefix_dir_liu = '/project/fudenber_735/GEO/liu_deWit_GSE181848/'
cooler_prefix_liu = prefix_dir_liu
sample_dict_liu ={
    'Wapl-0h':'GSM5512837_HiC.01_WAPL_0h.mcool',
    'Wapl-6h':'GSM5512838_HiC.02_WAPL_6h.mcool'
}
sample='Wapl-6h'
mcool_path_liu =  cooler_prefix_liu + sample_dict_liu[sample]
liu_cooler = cooler.Cooler(mcool_path_liu+'::resolutions/'+str(res))
view_df_liu = cooltools.lib.make_cooler_view(liu_cooler)[:19]


# data from Liu 2021
prefix_dir_liu_new = '/project/fudenber_735/GEO/liu_deWit_2021_GSE135180/distiller-mm10/results/coolers_library/'
cooler_prefix_liu_new = prefix_dir_liu_new
sample_dict_liu_new ={
    'Wapl-0h_new':'liu_wapl0h.mm10.mapq_30.1000.mcool',
    'Wapl-24h_new':'liu_wapl24h.mm10.mapq_30.1000.mcool'
}
sample='Wapl-24h_new'
mcool_path_liu_new =  cooler_prefix_liu_new + sample_dict_liu_new[sample]
liu_new_cooler = cooler.Cooler(mcool_path_liu_new+'::resolutions/'+str(res))
view_df_liu_new = cooltools.lib.make_cooler_view(liu_new_cooler)[:19]


from chromoscores.snipping import tad_snippet_sectors

from chromoscores.snipping import *
from chromoscores.snipping import peak_snipping
from chromoscores.scorefunctions import peak_score


beginning = 100_000
mindist = 90_000
band_edges = np.append([ beginning, 110_000], (100_000 + 15_000 * 1.3 ** np.arange(1,34)))
band_edge_list=list(band_edges[:])
print(band_edge_list[:31])
print(len(band_edge_list))

sites=peaks_motifs_bin_s


print(band_edge_list)

sample_dict = ['WT']
for sample in sample_dict:
    f=open('dots_vs_distance_peaks_motifs_%s_Bonev_bin_s_Justice.csv'%sample,'w')
    f.write('orientation,dist,n,peak_score\n')
    clr =  bonev_cooler
    pup = coolpup.pileup(clr, 
                         sites, 
                         features_format='bed', view_df=view_df_bonev[:1],
                         flip_negative_strand=True,
                         by_distance=np.array(band_edges),
                         by_strand=True, mindist=mindist, maxdist=10_000_000,
                        flank=50_000, min_diag=2,
                        nproc=19
                        )
    for i in range(len(pup['data'])):
        pseudocount = np.mean(np.isfinite(pup['data'][i]))/1000
        score=peak_score(pup['data'][i],3,4,pseudocount)
        if pup['separation'][i]=='all':
            continue
        
        dist=np.mean(pup['distance_band'][i])
        orientation=pup['orientation'][i]
        n = pup['n'][i]
        f.write('%s,%s,%s,%s\n'%(orientation,dist,n,score))
    fg = plotpup.plot(pup,  rows='orientation',cols='separation',
                  row_order=['-+', '--', '++', '+-'],score=False, 
                  cmap='fall', 
                  scale='log', sym=False,
                  height=3)
    
    plt.show()
    f.close()



sample_dict_liu ={
    'Wapl-0h':'GSM5512837_HiC.01_WAPL_0h.mcool',
    'Wapl-6h':'GSM5512838_HiC.02_WAPL_6h.mcool'
}

mindist = 100_000
for sample in sample_dict_liu:
    f=open('dots_vs_distance_peaks_motifs_%s_Liu_bin_s_pseudo1000.csv'%sample,'w')
    f.write('orientation,dist,n,peak_score\n')
    mcool_path_liu =  cooler_prefix_liu + sample_dict_liu[sample]
    liu_cooler = cooler.Cooler(mcool_path_liu+'::resolutions/'+str(res))
    view_df_liu = cooltools.lib.make_cooler_view(liu_cooler)[:19]
    clr =  liu_cooler
    pup = coolpup.pileup(clr, 
                         sites, 
                         features_format='bed', view_df=view_df_liu[:19],
                         flip_negative_strand=True,
                         by_distance=np.array(band_edges),
                         by_strand=True, mindist=mindist, maxdist=10_000_000,
                        flank=50_000, min_diag=2,
                        nproc=19
                        )
    for i in range(len(pup['data'])):
        pseudocount = np.mean(np.isfinite(pup['data'][i]))/1000
        score=peak_score(pup['data'][i],3,4,pseudocount)
        if pup['separation'][i]=='all':
            continue
        
        dist=np.mean(pup['distance_band'][i])
        orientation=pup['orientation'][i]
        n = pup['n'][i]
        f.write('%s,%s,%s,%s\n'%(orientation,dist,n,score))
    fg = plotpup.plot(pup,  rows='orientation',cols='separation',
                  row_order=['-+', '--', '++', '+-'],score=False, 
                  cmap='fall', 
                  scale='log', sym=False,
                  height=3)
    
    plt.show()
    f.close()
















