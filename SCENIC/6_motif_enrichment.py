#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:08:21 2023

@author: nrq2
"""
# https://pycistarget.readthedocs.io/en/latest/pycistarget_scenic%2B_wrapper.html
#%matplotlib inline
import pycistarget
pycistarget.__version__
import os

import pyranges as pr
import os
from scenicplus.wrappers.run_pycistarget import *

# Load region binarized topics
import pickle
work_dir = '/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/'
outDir = '/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/'
infile = open(outDir+'topic_binarization/binarized_topic_region.pkl', 'rb')
binarized_topic_region = pickle.load(infile)
infile.close()
# Load DARs

infile = open(outDir+'DARs/DARs.pkl', 'rb')
DARs_dict = pickle.load(infile)
infile.close()
# Format region sets
import re
import pyranges as pr
from pycistarget.utils import *
region_sets = {}
region_sets['Topics'] = {key: pr.PyRanges(region_names_to_coordinates(binarized_topic_region[key].index.tolist())) for key in binarized_topic_region.keys()}
region_sets['DARs'] = {re.sub('[^A-Za-z0-9]+', '_', key): pr.PyRanges(region_names_to_coordinates(DARs_dict[key].index.tolist())) for key in DARs_dict.keys()}


## precomputed database from https://resources.aertslab.org/cistarget/
DBpath='/data/rachel/Linda_Lako/Retina/cistargetDB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
annoPath='/data/rachel/Linda_Lako/Retina/cistargetDB/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl'

savePath=os.path.join(work_dir, "SCREEN_cluster_V10_V2")
#os.mkdir(savePath)


demOuts=os.path.join(work_dir, "cluster_SCREEN.regions_vs_motifs.scores.feather")
run_pycistarget(region_sets,
                  ctx_db_path = DBpath,
                  species = 'homo_sapiens',
                  save_path = savePath,
                  dem_db_path = DBpath,
                  run_without_promoters = True,
                  biomart_host = 'http://www.ensembl.org',
                  promoter_space = 500,
                  ctx_auc_threshold = 0.005,
                  ctx_nes_threshold = 3.0,
                  ctx_rank_threshold = 0.05,
                  dem_log2fc_thr = 0.5,
                  dem_motif_hit_thr = 3.0,
                  dem_max_bg_regions = 500,
                  path_to_motif_annotations = annoPath,
                  annotation_version = 'v10nr_clust',
                  annotation = ['Direct_annot', 'Orthology_annot'],
                  n_cpu = 1,
                  _temp_dir = '/data/rachel/')
















