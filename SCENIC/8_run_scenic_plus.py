#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 11:26:40 2023

@author: nrq2
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import scenicplus
scenicplus.__version__
import pickle
import scanpy
# Load functions
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *
from pycistarget.motif_enrichment_dem import *
from scenicplus.wrappers.run_scenicplus import run_scenicplus

work_dir = '/data/rachel/Linda_Lako/Retina/AD3/ATAC_analysis/SCENIC/'

## load scenic plus object
scplus_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/scplus_obj.pkl'), 'rb'))

# =============================================================================
# 
# =============================================================================
## check with which biomart host our gene names match
# ensembl_version_dict = {'105': 'http://www.ensembl.org',
#                         '104': 'http://may2021.archive.ensembl.org/',
#                         '103': 'http://feb2021.archive.ensembl.org/',
#                         '102': 'http://nov2020.archive.ensembl.org/',
#                         '101': 'http://aug2020.archive.ensembl.org/',
#                         '100': 'http://apr2020.archive.ensembl.org/',
#                         '99': 'http://jan2020.archive.ensembl.org/',
#                         '98': 'http://sep2019.archive.ensembl.org/',
#                         '97': 'http://jul2019.archive.ensembl.org/',
#                         '96': 'http://apr2019.archive.ensembl.org/',
#                         '95': 'http://jan2019.archive.ensembl.org/',
#                         '94': 'http://oct2018.archive.ensembl.org/',
#                         '93': 'http://jul2018.archive.ensembl.org/',
#                         '92': 'http://apr2018.archive.ensembl.org/',
#                         '91': 'http://dec2017.archive.ensembl.org/',
#                         '90': 'http://aug2017.archive.ensembl.org/',
#                         '89': 'http://may2017.archive.ensembl.org/',
#                         '88': 'http://mar2017.archive.ensembl.org/',
#                         '87': 'http://dec2016.archive.ensembl.org/',
#                         '86': 'http://oct2016.archive.ensembl.org/',
#                         '80': 'http://may2015.archive.ensembl.org/',
#                         '77': 'http://oct2014.archive.ensembl.org/',
#                         '75': 'http://feb2014.archive.ensembl.org/',
#                         '54': 'http://may2009.archive.ensembl.org/'}
# import pybiomart as pbm
# def test_ensembl_host(scplus_obj, host, species):
#     dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
#     annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
#     annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
#     annot['Chromosome'] = annot['Chromosome'].astype('str')
#     filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
#     annot = annot[~filter]
#     annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
#     gene_names_release = set(annot['Gene'].tolist())
#     ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
#     print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
#     return ov
# n_overlap = {}
# for version in ensembl_version_dict.keys():
#     print(f'host: {version}')
#     try:
#         n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
#     except:
#         print('Host not reachable')
# v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
# print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

## set host
biomart_host = "http://jul2018.archive.ensembl.org"

## SCENIC+ workflow
## from here
## get list of of TF from nostromo
## https://scenicplus.readthedocs.io/en/latest/mix_melanoma_cell_lines.html#inferring-enhancer-driven-Gene-Regulatory-Networks-(eGRNs)-using-SCENIC+
## run scenic plus
from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    sys.stderr = open(os.devnull, "w")  # silence stderr
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['CellType'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = 'utoronto_human_tfs_v_1.01.txt',
        save_path = os.path.join(work_dir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = False,
        calculate_DEGs_DARs = False,
        export_to_loom_file = False,
        export_to_UCSC_file = False,
        n_cpu = 12,
        _temp_dir = '/data/rachel/temp')
    sys.stderr = sys.__stderr__  # unsilence stderr
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus_results/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)
