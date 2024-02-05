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
work_dir="/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/"


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

pickle.dump(cisTopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))



