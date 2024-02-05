#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:48:59 2023

@author: nrq2
"""

#%matplotlib inline
import scenicplus
scenicplus.__version__
import pickle
import scanpy
# Load functions
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *
from pycistarget.motif_enrichment_dem import *
work_dir = '/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/'

# Create SCENIC+ object

## ATAC - cisTopic object
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

## Precomputed imputed data
imputed_acc_obj = pickle.load(open(os.path.join(work_dir, 'DARs/Imputed_accessibility.pkl'), 'rb'))

## RNA - anndata
rna_anndata = scanpy.read_h5ad("retina_all_cells_RNA.h5ad")


## check intersection
rna=set(rna_anndata.obs.CellType)
atac=set(cistopic_obj.cell_data.CellType)
rna.intersection(atac)


## make cell names match with ATACseq
ct = rna_anndata.obs.CellType
#ct = ct.replace('BCs', 'BPs')
#ct = ct.replace('Cones', 'cones')
#ct = ct.replace('MCs', 'MG')
#ct = ct.replace('Putative cillary margin', 'pCM')
#ct = ct.replace('Rods', 'rods')
rna_anndata.obs.CellType = ct






##  motif enrichment results
motif_enrichment_dict=pickle.load(open(os.path.join(work_dir, 'SCREEN_cluster_V10_V2/menr.pkl'), 'rb'))


## create a scenic plus object
scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = rna_anndata,
        cisTopic_obj = cistopic_obj,
        imputed_acc_obj = imputed_acc_obj,
        menr = motif_enrichment_dict,
        multi_ome_mode = False,
        key_to_group_by = 'CellType',
        ACC_prefix = 'ACC_',
        GEX_prefix = 'GEX_',
        bc_transform_func = lambda x: x,
        normalize_imputed_acc = False)


pickle.dump(scplus_obj,
            open(os.path.join(work_dir, 'scATAC/scplus_obj.pkl'), 'wb'))





