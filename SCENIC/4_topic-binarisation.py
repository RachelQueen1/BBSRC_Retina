#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 10:30:42 2023

@author: nrq2
"""
print("load libraries")
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

work_dir="/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/"

os.mkdir(work_dir+'topic_binarization')

cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
#binarize the topic-region distributions
region_bin_topics = binarize_topics(cistopic_obj, method='otsu', ntop=3000, plot=True, num_columns=5, save= work_dir + 'topic_binarization/otsu.pdf')

# binarize the cell-topic distribions
binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li', plot=True, num_columns=5, nbins=60)


from pycisTopic.topic_qc import *

topic_qc_metrics = compute_topic_metrics(cistopic_obj)


fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)


fig=plt.figure(figsize=(40, 43))


i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
fig.savefig(work_dir + 'topic_binarization/Topic_qc.pdf', bbox_inches='tight')
plt.show()

# ## change early late to RPC
# s=cistopic_obj.cell_data.CellType
# s=s.replace('early RPCs', 'RPCs')
# s=s.replace('late RPCs', 'RPCs')
# cistopic_obj.cell_data.CellType=s


topic_annot = topic_annotation(cistopic_obj, annot_var='CellType', binarized_cell_topic=binarized_cell_topic, general_topic_thr = 0.2)

topic_qc_metrics = pd.concat([topic_annot[['CellType']], topic_qc_metrics], axis=1)
topic_qc_metrics


# Save
with open(work_dir + 'topic_binarization/Topic_qc_metrics_annot.pkl', 'wb') as f:
  pickle.dump(topic_qc_metrics, f)
with open(work_dir + 'topic_binarization/binarized_cell_topic.pkl', 'wb') as f:
  pickle.dump(binarized_cell_topic, f)
with open(work_dir + 'topic_binarization/binarized_topic_region.pkl', 'wb') as f:
  pickle.dump(region_bin_topics, f)
  
  

pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

