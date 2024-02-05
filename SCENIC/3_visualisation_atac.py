#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 09:27:48 2023

@author: nrq2
"""

import scenicplus
import scanpy
import pycisTopic
## 1. Initialize cisTopic object
from pycisTopic.cistopic_class import *
from pycisTopic.utils import *
from pycisTopic.lda_models import CistopicLDAModel
import pickle
from pycisTopic.cistopic_class import *
from pycisTopic.clust_vis import *



work_dir="/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/"
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
## make cell names match with RNAseq
ct = cistopic_obj.cell_data.CellType
ct = ct.replace('Gabaergic amacrine cells', 'ACs')
ct = ct.replace('glycinergic amacrine cells ', 'ACs')
ct = ct.replace('starburst amacrine cells', 'ACs')
cistopic_obj.cell_data.CellType =ct

models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/mix_mm_models_500_iter_LDA.pkl'), 'rb'))

model=evaluate_models(models,
                     select_model=None,
                     return_model=True,
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= work_dir + 'scATAC/models/model_selection.pdf')


# Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)




# umap
run_umap(cistopic_obj, target  = 'cell', scale=True)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

os.mkdir(work_dir +'/visualisation')

cistopic_obj.cell_data.index = cistopic_obj.cell_names

plot_metadata(cistopic_obj,
                 reduction_name='UMAP',
                 variables=['CellType'], # Labels from RNA and new clusters
                 target='cell', num_columns=3,
                 text_size=10,
                 dot_size=5,
                 figsize=(15,5),
                 save= work_dir + 'visualisation/dimensionality_reduction_label.pdf')



plot_topic(cistopic_obj,
            reduction_name = 'UMAP',
            target = 'cell',
            num_columns=5,
            save= work_dir + 'visualisation/dimensionality_reduction_topic_contr.pdf')



cell_topic_heatmap(cistopic_obj,
                     variables = ['CellType'],
                     scale = False,
                     legend_loc_x = 1.05,
                     legend_loc_y = -1.2,
                     legend_dist_y = -1,
                     figsize=(10,10),
                     save = work_dir + 'visualisation/heatmap_topic_contr.pdf')



pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))



