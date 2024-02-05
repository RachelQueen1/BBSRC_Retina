#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:30:56 2023
https://scenicplus.readthedocs.io/en/latest/Scenicplus_step_by_step-RTD.html#6.-Exploring-SCENIC+-results
@author: nrq2
"""

## open eGRNs
import dill
import scanpy
from scenicplus.utils import format_egrns
from scenicplus.eregulon_enrichment import *
from scenicplus.utils import format_egrns
from scenicplus.dimensionality_reduction import *
import pandas
from scenicplus.networks import *
import networkx as nx
from scenicplus.diff_features import *
from itertools import compress
"""
A. Generate eRegulon metadata
As a first step, we will format the eGRN metadata. This will integrate the 
results from the inference of enhancer-to-gene and TF-to-gene relationships 
and the eRegulon construction in one pandas dataframe that can be used 
for further exploration.
"""

## open eGRNs
work_dir = '/data/rachel/Linda_Lako/Retina/ATAC_analysis/SCENIC/'
scenic_plus_path=os.path.join(work_dir, 'scenicplus/scplus_obj.pkl')
scplus_obj= pickle.load(open(scenic_plus_path, 'rb'))



gene_ranking_path=os.path.join(work_dir, 'scenicplus/gene_ranking.pkl')
gene_ranking= pickle.load(open(gene_ranking_path, 'rb'))
# infile = open(scenic_plus_path, 'rb')
# scplus_obj = dill.load(infile)
# infile.close()


  
#look at eRegulaon data:
eregulons = scplus_obj.uns['eRegulon_metadata']






"""
D. Identification of high quality regulons
We will select a subset of regulons based on the correlation between the 
region based and the gene based regulons. We will only use the extended 
eRegulon if there is not a direct eRegulon available.
"""

# Correlation between region based regulons and gene based regulons
df1 = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
df2 = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()
df1.columns = [x.split('_(')[0] for x in df1.columns]
df2.columns = [x.split('_(')[0] for x in df2.columns]
correlations = df1.corrwith(df2, axis = 0)
correlations = correlations[abs(correlations) > 0.6]
# Kepp only R2G +
keep = [x for x in correlations.index if '+_+' in x] + [x for x in correlations.index if '-_+' in x]
# Keep extended if not direct
extended = [x for x in keep if 'extended' in x]
direct = [x for x in keep if not 'extended' in x]
keep_extended = [x for x in extended if not x.replace('extended_', '') in direct]
keep = direct + keep_extended
# Keep regulons with more than 10 genes
keep_gene = [x for x in scplus_obj.uns['eRegulon_AUC']['Gene_based'].columns if x.split('_(')[0] in keep]
keep_gene = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 10)]
keep_all = [x.split('_(')[0] for x in keep_gene]
keep_region = [x for x in scplus_obj.uns['eRegulon_AUC']['Region_based'].columns if x.split('_(')[0] in keep]
scplus_obj.uns['selected_eRegulons'] = {}
scplus_obj.uns['selected_eRegulons']['Gene_based'] = keep_gene
scplus_obj.uns['selected_eRegulons']['Region_based'] = keep_region

len(keep_gene)

"""
E. Overlap between eRegulons
In addition, to assess which eRegulons tend to be enriched in the same group 
of cells we can generate a correlation plot as well.
"""
# from scenicplus.plotting.correlation_plot import *
# correlation_heatmap(scplus_obj,
#                     auc_key = 'eRegulon_AUC',
#                     signature_keys = ['Gene_based'],
#                     selected_regulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
#                     fcluster_threshold = 0.1,
#                     fontsize = 3)
# """
# We can also check the overlap between eRegulons:
# """
# jaccard_heatmap(scplus_obj,
#                     gene_or_region_based = 'Gene_based',
#                     signature_key = 'eRegulon_signatures',
#                     selected_regulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
#                     fcluster_threshold = 0.1,
#                     fontsize = 3,
#                     method='intersect')





"""
eGRN dimensionality reduction
"""
plot_metadata(scplus_obj,
                 reduction_name='eRegulons_UMAP',
                 variables=['CellType'],
                 num_columns=1,
                 text_size=10,
                 dot_size=5)


"""
Find clusters based on eRegulon data
"""
find_clusters(scplus_obj,
              signature_keys=['Gene_based', 'Region_based'],
              k = 10,
              res = [0.6, 1.2, 1.5],
              prefix = 'SCENIC+_',
              scale = True)


plot_metadata(scplus_obj,
                 reduction_name='eRegulons_UMAP',
                 variables=['CellType', 'SCENIC+_leiden_10_1.5'],
                 num_columns=2,
                 text_size=10,
                 dot_size=5)

"""
G. eRegulon specificity scores
"""
from scenicplus.RSS import *
regulon_specificity_scores(scplus_obj,
                         'CellType',
                         signature_keys=['Gene_based'],
                         selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'],
                         out_key_suffix='_gene_based',
                         scale=False)
plot_rss(scplus_obj, 'CellType_gene_based', num_columns=4, top_n=10)



## rename cones to photoreceptor precursors to match current annotation
# celltype = scplus_obj.uns['RSS']['CellType_gene_based']
# celltype = celltype.rename(index={'cones': 'photoreceptors'})
# scplus_obj.uns['RSS']['CellType_gene_based'] = celltype

# scplus_obj.metadata_cell = scplus_obj.metadata_cell.replace('cones', 'photoreceptors')


from scenicplus.plotting.dotplot import *
from plotnine import ggplot, aes, geom_line, ggsave

## heatmap of cell type specific eRegulons
plt_heatmap=(heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['RSS']['CellType_gene_based'],
        color_matrix = scplus_obj.to_df('EXP'),
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'CellType',
        subset_eRegulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
        figsize = (15, 7),
        orientation = 'horizontal',
        split_repressor_activator=True))

ggsave(plt_heatmap, 'scenic_heatmap.tiff', format='tiff', dpi=600)



# ### add DEGs and DARs
# get_differential_features(scplus_obj, 'CellType', use_hvg = False, contrast_type = ['DEGs', 'DARs'])
# #get_differential_features(scplus_obj, 'ACC_Seurat_cell_type', use_hvg = True, contrast_type = ['DARs', 'DEGs'])

# from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
# apply_std_filtering_to_eRegulons(scplus_obj)

# from scenicplus.networks import *

# from pycisTopic.diff_features import find_highly_variable_features
# hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Region']))], n_top_features=300, plot = False)
# rna_anndata = scanpy.read_h5ad("AD3_all_cells_RNA.h5ad")
# scanpy.pp.highly_variable_genes(rna_anndata, n_top_genes = 1000, flavor = "seurat")
# hvg = rna_anndata.var.index.tolist()

# ### DE genes for each cell type
# df = pd.read_excel("Table S4.xlsx",  sheet_name='integrated WT2-ROs', skiprows=2)
# #df = df[df['avg_log2FC']>0.7]


# RPCs =df.loc[ df['cell fate'].str.contains("RPC"), "gene"].unique().tolist()

# PRs_types = ["BCs", "PR precursors", "Rods", "Cones"]
# intersection_df = df[df['cell fate'].isin(PRs_types)]  
# PRs= intersection_df["gene"].unique().tolist()                                            

# #### PLOT NETWORKS
# cat_cols = {'CRX': 'Orange', 'NEUROD1': 'Purple', 'ATOH7': 'Red', 'HMGA2': 'Orange', 'FOS': 'Blue', 'FOSB': 'Blue'}


# import matplotlib
# matplotlib.use("Agg")

# ### create function for plotting
# import matplotlib.pyplot as plt
# from matplotlib.patches import Patch

# def plot_networkx_save(G, pos, saveName):
#     nx.draw_networkx_nodes(G, pos, node_color=nx.get_node_attributes(G,'color').values(),
#                            node_size=list(nx.get_node_attributes(G,'size').values()),
#                            node_shape = 'D')
#     nx.draw_networkx_edges(G, pos, edge_color = nx.get_edge_attributes(G,'color').values(),
#                            width = list(nx.get_edge_attributes(G,'width').values()))
#     fontsize_d = {y:x['size'] for x,y in zip(list(nx.get_node_attributes(G,'font').values()),list(nx.get_node_attributes(G,'label').values())) if x['size'] != 0.0}
#     fontcolor_d = {y:x['color'] for x,y in zip(list(nx.get_node_attributes(G,'font').values()),list(nx.get_node_attributes(G,'label').values())) if x['size'] != 0.0}
#     for node, (x, y) in pos.items():
#         if node in fontsize_d.keys():
#             plt.text(x, y, node, fontsize=fontsize_d[node], color=fontcolor_d[node],  ha='center', va='center')
#     ax = plt.gca()
#     ax.margins(0.11)
#     plt.tight_layout()
#     plt.axis("off")
    
   
#     # this line is added in because the files were saved as blank
#     plt.gcf()
#     plt.savefig(saveName, dpi=600) 
#     plt.show()
    




# def visualize_network(scplus_obj, eRegs, hvr, hvg, cat_cols, saveName, filterRegulons=True):
    
    
#     if filterRegulons:
#         nx_tables = create_nx_tables(scplus_obj,
#                                  eRegulon_metadata_key='eRegulon_metadata_filtered',
#                                  subset_eRegulons=eRegs,
#                                  subset_regions=hvr,
#                                  subset_genes=hvg,
#                                  add_differential_gene_expression=True,
#                                  add_differential_region_accessibility=True,
#                                  differential_variable=['CellType'])
#         # 2a: Filter R2G Importance Table
#         R2Gimportance_table = nx_tables['Edge']['R2G']
#         nx_tables['Edge']['R2G'] = R2Gimportance_table[R2Gimportance_table['R2G_importance'] > 0.07]
#         genes_plot = R2Gimportance_table[R2Gimportance_table['R2G_importance'] > 0.08]["Gene"].tolist()
#     else:
#         nx_tables = create_nx_tables(scplus_obj,
#                                   eRegulon_metadata_key='eRegulon_metadata',
#                                   subset_eRegulons=eRegs,
#                                   subset_regions=hvr,
#                                   subset_genes=None,
#                                   add_differential_gene_expression=True,
#                                   add_differential_region_accessibility=True,
#                                   differential_variable=['CellType'])
#         # 2b: No Filter genes
#         genes_plot = nx_tables['Edge']['R2G']["Gene"].tolist()
    
#     # 3: Create NetworkX Graph
#     G_c, pos_c, edge_tables_c, node_tables_c = create_nx_graph(nx_tables,
#                                                                use_edge_tables=['TF2R', 'R2G'],
#                                                                color_edge_by={'TF2R': {'variable': 'TF', 'category_color': cat_cols},
#                                                                               'R2G': {'variable': 'R2G_rho',
#                                                                                       'continuous_color': 'viridis',
#                                                                                       'v_min': -1, 'v_max': 1}},
#                                                                transparency_edge_by={'R2G': {'variable': 'R2G_importance',
#                                                                                               'min_alpha': 0.1,
#                                                                                               'v_min': 0}},
#                                                                width_edge_by={'R2G': {'variable': 'R2G_importance',
#                                                                                       'max_size': 1.5, 'min_size': 1}},
#                                                                color_node_by={'TF': {'variable': 'TF',
#                                                                                     'category_color': cat_cols},
#                                                                               'Region': {'variable': 'CellType_Log2FC_pCM',
#                                                                                          'continuous_color': 'viridis'},
#                                                                               'Gene': {'variable': 'fixed_color',
#                                                                                        'fixed_color': '#F5F5F5'}
#                                                                               },
#                                                                transparency_node_by={'R2G': {'variable': 'R2G_importance',
#                                                                                                'min_alpha': 0.1,
#                                                                                                'v_min': 0}},
#                                                                size_node_by={'TF': {'variable': 'fixed_size', 
#                                                                                     'fixed_size': 30},
#                                                                              'Gene': {'variable': 'fixed_size',
#                                                                                       'fixed_size': 1},
#                                                                              'Region': {'variable': 'fixed_size',
#                                                                                         'fixed_size': 20}},
#                                                                shape_node_by={'TF': {'variable': 'fixed_shape',
#                                                                                       'fixed_shape': 'ellipse'},
#                                                                              'Gene': {'variable': 'fixed_shape',
#                                                                                       'fixed_shape': 'ellipse'},
#                                                                              'Region': {'variable': 'fixed_shape',
#                                                                                         'fixed_shape': 'diamond'}},
#                                                                label_size_by={'TF': {'variable': 'fixed_label_size',
#                                                                                     'fixed_label_size': 30.0},
#                                                                               'Gene': {'variable': 'fixed_label_size',
#                                                                                        'fixed_label_size': 10.0},
#                                                                               'Region': {'variable': 'fixed_label_size',
#                                                                                          'fixed_label_size': 0}},
#                                                                layout='kamada_kawai_layout',
#                                                                scale_position_by=250)

#     # Step 5: Save the Graph
#     plt.figure(figsize=(10, 10))
#     plot_networkx_save(G_c, pos_c, saveName)



# visualize_network(scplus_obj, ['CRX'], hvr, hvg, cat_cols, "CRX_network.tiff")
# visualize_network(scplus_obj, ['ATOH7'], hvr, hvg, cat_cols, "ATOH7_network.tiff")
# visualize_network(scplus_obj, ['NEUROD1'], hvr, hvg, cat_cols, "NEUROD1_network.tiff")



# visualize_network(scplus_obj, ['FOS'], hvr, hvg, cat_cols, "FOS_network.tiff", filterRegulons=False)
# visualize_network(scplus_obj, ['FOSB'], hvr, hvg, cat_cols, "FOSB_network.tiff", filterRegulons=False)
# visualize_network(scplus_obj, ['HMGA2'], hvr, hvg, cat_cols, "HMGA2_network.tiff", filterRegulons=False)
