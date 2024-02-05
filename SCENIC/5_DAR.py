#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:42:15 2023

@author: nrq2
"""
import os
import scenicplus
import scanpy
import pycisTopic
## 1. Initialize cisTopic object
from pycisTopic.cistopic_class import *
from pycisTopic.utils import *
from pycisTopic.lda_models import CistopicLDAModel
from pycisTopic.cistopic_class import *
from pycisTopic.topic_binarization import *
from pycisTopic.clust_vis import *
work_dir="/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/"

cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))



## Differentially Accessible Regions (DARs)
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

os.mkdir(work_dir+'DARs')
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj,
                                           min_disp = 0.05,
                                           min_mean = 0.0125,
                                           max_mean = 3,
                                           max_disp = np.inf,
                                           n_bins=20,
                                           n_top_features=None,
                                           plot=True,
                                           save= work_dir + 'DARs/HVR_plot.pdf')
## There is a total of 41011 variable regions
len(variable_regions)

markers_dict= find_diff_features(cistopic_obj,
                      imputed_acc_obj,
                      variable='CellType',
                      var_features=variable_regions,
                      contrasts=None,
                      adjpval_thr=0.05,
                      log2fc_thr=np.log2(2),
                      n_cpu=5,
                      _temp_dir='/data/rachel/tempModels_retina',
                      split_pattern = '-')


# Save
with open(work_dir + 'DARs/Imputed_accessibility.pkl', 'wb') as f:
  pickle.dump(imputed_acc_obj, f)
with open(work_dir + 'DARs/DARs.pkl', 'wb') as f:
  pickle.dump(markers_dict, f)


## plot region accessibility into the cell-topic UMAP
from pycisTopic.clust_vis import *
plot_imputed_features(cistopic_obj,
                    reduction_name='UMAP',
                    imputed_data=imputed_acc_obj,
                    features=[markers_dict[x].index.tolist()[0] for x in ['T3']],
                    scale=False,
                    num_columns=4,
                    save= work_dir + 'DARs/example_best_DARs.pdf')



  
## Gene activity
import pybiomart as pbm
import pyranges as pr

dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')


annot = dataset.query(attributes=['chromosome_name', 'start_position', 'end_position', 'strand', 'external_gene_name', 'transcription_start_site', 'transcript_biotype'])
annot['Chromosome/scaffold name'] = 'chr' + annot['Chromosome/scaffold name'].astype(str)
annot.columns=['Chromosome', 'Start', 'End', 'Strand', 'Gene','Transcription_Start_Site', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']
annot.Strand[annot.Strand == 1] = '+'
annot.Strand[annot.Strand == -1] = '-'
pr_annotation = pr.PyRanges(annot.dropna(axis = 0))


# Get chromosome sizes
import pandas as pd
import requests
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
chromsizes=pr.PyRanges(chromsizes)


from pycisTopic.gene_activity import *
gene_act, weigths = get_gene_activity(imputed_acc_obj, # Region-cell probabilities
                pr_annotation, # Gene annotation
                chromsizes, # Chromosome size
                use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
                upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                                      #these bp will be taken (1kbp here)
                downstream=[1000,100000], # Search space downstream
                distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
                decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
                extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                                      #this weight)
                extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
                gene_size_weight=False, # Whether to add a weights based on the length of the gene
                gene_size_scale_factor='median', # Dividend to calculate the gene size weigth. Default is the median value of all genes
                                      #in the genome
                remove_promoters=False, # Whether to remove promoters when computing gene activity scores
                average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                                      #activity score
                scale_factor=1, # Value to multiply for the final gene activity matrix
                extend_tss=[10,10], # Space to consider a promoter
                gini_weight = True, # Whether to add a gini index weigth. The more unique the region is, the higher this weight will be
                return_weights= True, # Whether to return the final weights
                project='Gene_activity') # Project name for the gene activity object

os.mkdir(work_dir+'DAGs')
plot_imputed_features(cistopic_obj,
                    reduction_name='UMAP',
                    imputed_data=gene_act,
                    features=['RXRG', 'CRX', 'NRL', 'ONECUT2'],
                    scale=True,
                    num_columns=4,
                    save= work_dir + 'DAGs/example_best_DAGs.pdf')


#How many DAGs do we have per cell type?
x = [print(x + ': '+ str(len(markers_dict[x]))) for x in markers_dict.keys()]

# Save
with open(work_dir + 'DAGs/Gene_activity.pkl', 'wb') as f:
  pickle.dump(gene_act, f)
with open(work_dir + 'DAGs/DAGs.pkl', 'wb') as f:
  pickle.dump(markers_dict, f)


