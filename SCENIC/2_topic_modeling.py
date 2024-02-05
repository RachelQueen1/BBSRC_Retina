#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 09:53:30 2023

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

work_dir="/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/"
if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))


print("load data")
import pickle
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))


print("load data")
#sys.stderr = open(os.devnull, "w")  # silence stderr

print("run models..")
models=run_cgs_models(cistopic_obj,
                    n_topics=[2,5,10,15,30,45],
                    n_cpu=6,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = '/data/rachel/tempModels')
#sys.stderr = sys.__stderr__  # unsilence stderr

print("write models..")


pickle.dump(models,
    open(os.path.join(work_dir, 'scATAC/models/mix_mm_models_500_iter_LDA.pkl'), 'wb'))


cistopic_obj.add_LDA_model(model)


### run umap
run_umap(cistopic_obj, target = 'cell', scale = True)




pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))



