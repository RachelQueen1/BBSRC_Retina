# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 14:46:54 2023

@author: nrq2

"""


## ssh rachel@nostromo.ncl.ac.uk
## spyder-scenic-ad3
# cd /data/rachel/Linda_Lako/Retina/AD3/ATAC_analysis/SCENIC/
# conda activate scenicplus
# python -m spyder_kernels.console
# /home/rachel/.local/share/jupyter/runtime/kernel-2436331.json
# https://docs.spyder-ide.org/current/_images/console-connect-remote-step2.gif

import scenicplus
import scanpy
import pycisTopic
## 1. Initialize cisTopic object
from pycisTopic.cistopic_class import *
from pycisTopic.utils import *
from pycisTopic.lda_models import CistopicLDAModel
adata = scanpy.read_h5ad("AD3_all_cells_RNA.h5ad")


# Load count matrix
matrix_path="feathers/model_to_pycisTopic/count_matrix.feather"
fragment_matrix = pd.read_feather(matrix_path)

regions=pd.read_csv("regions.txt")
region_names = regions["x"].tolist()
fragment_matrix.index = region_names

cisTopic_obj = pycisTopic.cistopic_class.create_cistopic_object(fragment_matrix)

# Also add the cell annotation, cell_data should be a pandas df with cells as rows (cell names as index) and variables as columns
cell_data = pd.read_feather("feathers/model_to_pycisTopic/celldata.feather")
cell_data.index = fragment_matrix.columns

cisTopic_obj.cell_data.index = fragment_matrix.columns
cisTopic_obj.add_cell_data(cell_data)

cisTopic_obj.cell_names = fragment_matrix.columns


## 2. Add model
## check that the data can be plotted after loading
## https://pycistopic.readthedocs.io/en/latest/Toy_melanoma-RTD.html

from pycisTopic.cistopic_class import *
sys.stderr = open(os.devnull, "w")  # silence stderr
work_dir="/data/rachel/Linda_Lako/Retina/AD3/ATAC_analysis/SCENIC/"
models=run_cgs_models(cisTopic_obj,
                    n_topics=[16,17,18,19,20,21,22,23,24],
                    n_cpu=6,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None, _temp_dir = work_dir)


if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
    open(os.path.join(work_dir, 'scATAC/models/mix_mm_models_500_iter_LDA.pkl'), 'wb'))


numTopics = 20
model = evaluate_models(models,
                     select_model = numTopics,
                     return_model = True,
                     metrics = ['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics = False)


cistopic_obj.add_LDA_model(model)



from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cisTopic_obj, method='otsu')

region_bin_topics = binarize_topics(cisTopic_obj, target="region",
                                    method='otsu', 
                                    ntop=3000, 
                                    plot=True, 
                                    num_columns=5, 
                                    save= 'topic_binarization/otsu.pdf')

pickle.dump(cisTopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))


# from pycisTopic.utils import *
# path_to_cisTopic_model_matrices="feathers/model_to_pycisTopic/"
# cell_topic = pd.read_feather(path_to_cisTopic_model_matrices + "cell_topic.feather")
# topic_region = pd.read_feather(path_to_cisTopic_model_matrices + "topic_region.feather")
# topic_region.index = ["Topic" + str(x) for x in range(1, topic_region.shape[0] + 1)]
# topic_region = topic_region.T
# parameters = None
# metrics = None
# metrics = None
# coherence = None
# marg_topic = None
# topic_ass = None
# model = pycisTopic.lda_models.CistopicLDAModel(metrics,  coherence, marg_topic, topic_ass, cell_topic, topic_region,parameters)
# cisTopic_obj.add_LDA_model(model)

# from pycisTopic.clust_vis import *
# run_umap(cisTopic_obj, target  = 'cell', scale=True)

# plot_metadata(cisTopic_obj, reduction_name = "UMAP", variables=['CellType'])
# ## create a Pycistarget motif enrichment dictionary
# # https://pycistarget.readthedocs.io/en/stable/pycistarget_scenic%2B_wrapper.html

# By default, we use binarized topics and DARs

# cisTopic_obj.selected_model.topic_ass["Regions_in_binarized_topic"]
# Traceback (most recent call last):
#   Input In [19] in <cell line: 1>
#     cisTopic_obj.selected_model.topic_ass["Regions_in_binarized_topic"]
# TypeError: 'NoneType' object is not subscriptable
