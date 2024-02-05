# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:46:34 2024

@author: nrq2
"""

import matplotlib
matplotlib.use("Agg")
from scenicplus.networks import *
import networkx as nx

### create function for plotting
import matplotlib.pyplot as plt

def plot_networkx_save(G, pos, saveName):
    nx.draw_networkx_nodes(G, pos, node_color=nx.get_node_attributes(G,'color').values(),
                           node_size=list(nx.get_node_attributes(G,'size').values()),
                           node_shape = 'D')
    nx.draw_networkx_edges(G, pos, edge_color = nx.get_edge_attributes(G,'color').values(),
                           width = list(nx.get_edge_attributes(G,'width').values()))
    fontsize_d = {y:x['size'] for x,y in zip(list(nx.get_node_attributes(G,'font').values()),list(nx.get_node_attributes(G,'label').values())) if x['size'] != 0.0}
    fontcolor_d = {y:x['color'] for x,y in zip(list(nx.get_node_attributes(G,'font').values()),list(nx.get_node_attributes(G,'label').values())) if x['size'] != 0.0}
    for node, (x, y) in pos.items():
        if node in fontsize_d.keys():
            plt.text(x, y, node, fontsize=fontsize_d[node], color=fontcolor_d[node],  ha='center', va='center')
    ax = plt.gca()
    ax.margins(0.11)
    plt.tight_layout()
    plt.axis("off")
    plt.show()
    # this line is added in because the files were saved as blank
    plt.gcf()
    plt.savefig(saveName, dpi=600) 
    




def visualize_network(scplus_obj, eRegs, hvr, hvg, celltype, cat_cols, saveName):
    # Step 1: Create NetworkX Tables
    nx_tables = create_nx_tables(scplus_obj,
                                 eRegulon_metadata_key='eRegulon_metadata_filtered',
                                 subset_eRegulons=eRegs,
                                 subset_regions=hvr,
                                 subset_genes=hvg,
                                 add_differential_gene_expression=True,
                                 add_differential_region_accessibility=True,
                                 differential_variable=['CellType'])

    # Step 2: Filter R2G Importance Table
    R2Gimportance_table = nx_tables['Edge']['R2G']
    nx_tables['Edge']['R2G'] = R2Gimportance_table[R2Gimportance_table['R2G_importance'] > 0.07]
    genes_plot = R2Gimportance_table[R2Gimportance_table['R2G_importance'] > 0.08]["Gene"].tolist()

    # Step 3: Extract Other Importance Tables
    TF2Rimportance_table = nx_tables['Edge']['TF2R']
    TF2Gimportance_table = nx_tables['Edge']['TF2G']

    # Step 4: Create NetworkX Graph
    region="CellType_Log2FC_" + celltype
    G_c, pos_c, edge_tables_c, node_tables_c = create_nx_graph(nx_tables,
                                                               use_edge_tables=['TF2R', 'R2G'],
                                                               color_edge_by={'TF2R': {'variable': 'TF', 'category_color': cat_cols},
                                                                              'R2G': {'variable': 'R2G_rho',
                                                                                      'continuous_color': 'viridis',
                                                                                      'v_min': -1, 'v_max': 1}},
                                                               transparency_edge_by={'R2G': {'variable': 'R2G_importance',
                                                                                              'min_alpha': 0.1,
                                                                                              'v_min': 0.1}},
                                                               width_edge_by={'R2G': {'variable': 'R2G_importance',
                                                                                      'max_size': 1.5, 'min_size': 1}},
                                                               color_node_by={'TF': {'variable': 'TF',
                                                                                    'category_color': cat_cols},
                                                                              'Region': {'variable': region,
                                                                                         'continuous_color': 'viridis'},
                                                                              'Gene': {'variable': 'fixed_color',
                                                                                       'fixed_color': '#FFFFFF'}
                                                                              },
                                                               transparency_node_by={'R2G': {'variable': 'R2G_importance',
                                                                                               'min_alpha': 0.1,
                                                                                               'v_min': 0}},
                                                               size_node_by={'TF': {'variable': 'fixed_size', 'fixed_size': 20},
                                                                             'Gene': {'variable': 'fixed_size',
                                                                                      'fixed_size': 0.001},
                                                                             'Region': {'variable': 'fixed_size',
                                                                                        'fixed_size': 1}},
                                                               shape_node_by={'TF': {'variable': 'fixed_shape',
                                                                                      'fixed_shape': 'ellipse'},
                                                                             'Gene': {'variable': 'fixed_shape',
                                                                                      'fixed_shape': 'ellipse'},
                                                                             'Region': {'variable': 'fixed_shape',
                                                                                        'fixed_shape': 'circle'}},
                                                               label_size_by={'TF': {'variable': 'fixed_label_size',
                                                                                    'fixed_label_size': 10.0},
                                                                              'Gene': {'variable': 'fixed_label_size',
                                                                                       'fixed_label_size': 10.0},
                                                                              'Region': {'variable': 'fixed_label_size',
                                                                                         'fixed_label_size': 0}},
                                                               layout='kamada_kawai_layout',
                                                               scale_position_by=250)

    # Step 5: Visualize the Graph
    plt.figure(figsize=(10, 10))
    plot_networkx_save(G_c, pos_c, saveName)
    plt.figure(figsize=(10, 10))
    plt.show()
    return(G_c, pos_c, edge_tables_c, node_tables_c)




def visualize_network_no_filter(scplus_obj, eRegs, hvr, hvg, celltype, cat_cols, saveName):
    # Step 1: Create NetworkX Tables
    nx_tables = create_nx_tables(scplus_obj,
                                 eRegulon_metadata_key='eRegulon_metadata',
                                 subset_eRegulons=eRegs,
                                 subset_regions=hvr,
                                 subset_genes=None,
                                 add_differential_gene_expression=True,
                                 add_differential_region_accessibility=True,
                                 differential_variable=['CellType'])

    # Step 2: Filter R2G Importance Table
    R2Gimportance_table = nx_tables['Edge']['R2G']
    nx_tables['Edge']['R2G'] = R2Gimportance_table[R2Gimportance_table['R2G_importance'] > 0.07]
    genes_plot = R2Gimportance_table[R2Gimportance_table['R2G_importance'] > 0.08]["Gene"].tolist()

    # Step 3: Extract Other Importance Tables
    TF2Rimportance_table = nx_tables['Edge']['TF2R']
    TF2Gimportance_table = nx_tables['Edge']['TF2G']
    
    region="CellType_Log2FC_" + celltype

    # Step 4: Create NetworkX Graph
    G_c, pos_c, edge_tables_c, node_tables_c = create_nx_graph(nx_tables,
                                                               use_edge_tables=['TF2R', 'R2G'],
                                                               color_edge_by={#'TF2R': {'variable': 'TF', 'category_color': cat_cols},
                                                                              'R2G': {'variable': 'R2G_rho',
                                                                                      'continuous_color': 'viridis',
                                                                                      'v_min': -1, 'v_max': 1}},
                                                               transparency_edge_by={'R2G': {'variable': 'R2G_importance',
                                                                                              'min_alpha': 0.1,
                                                                                              'v_min': 0.1}},
                                                               width_edge_by={'R2G': {'variable': 'R2G_importance',
                                                                                      'max_size': 1.5, 'min_size': 1}},
                                                               color_node_by={#'TF': {'variable': 'TF',
                                                                              #      'category_color': cat_cols},
                                                                              # 'Region': {'variable': 'CellType_Log2FC_pCM',
                                                                              #            'continuous_color': 'viridis'},
                                                                              'Gene': {'variable': 'fixed_color',
                                                                                       'fixed_color': '#FFFFFF'}
                                                                              },
                                                               transparency_node_by={'R2G': {'variable': 'R2G_importance',
                                                                                               'min_alpha': 0.1,
                                                                                               'v_min': 0}},
                                                               size_node_by={'TF': {'variable': 'fixed_size', 'fixed_size': 20},
                                                                             'Gene': {'variable': 'fixed_size',
                                                                                      'fixed_size': 0.001},
                                                                             'Region': {'variable': 'fixed_size',
                                                                                        'fixed_size': 1}},
                                                               shape_node_by={'TF': {'variable': 'fixed_shape',
                                                                                      'fixed_shape': 'ellipse'},
                                                                             'Gene': {'variable': 'fixed_shape',
                                                                                      'fixed_shape': 'ellipse'},
                                                                             'Region': {'variable': 'fixed_shape',
                                                                                        'fixed_shape': 'circle'}},
                                                               label_size_by={'TF': {'variable': 'fixed_label_size',
                                                                                    'fixed_label_size': 10.0},
                                                                              'Gene': {'variable': 'fixed_label_size',
                                                                                       'fixed_label_size': 10.0},
                                                                              'Region': {'variable': 'fixed_label_size',
                                                                                         'fixed_label_size': 0}},
                                                               layout='kamada_kawai_layout',
                                                               scale_position_by=250)
    # Step 5: Visualize the Graph
    plt.figure(figsize=(10, 10))
    plot_networkx_save(G_c, pos_c, saveName)
    
    
