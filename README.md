# "Spatiotemporal single cell analyses reveal a transient population of retinal progenitor cells in the ciliary margin of developing human retina"

This repository contains the analysis code needed to reproduce the results in the above paper.

## Data

The raw and processed data can be downloaded from GEO under the
following accession number: GSE234971.

## Software

Software for the analysis used:

Seurat (Version 4.3.0) Cran

DoubletFinder (Version 2.0.3) Bioconductor

Harmony (Version 0.1.1) Bioconductor

Monocle 3 (Version 1.3.1) Bioconductor

Signac (Version 1.6) Bioconductor

ComplexHeatmap (Version 2.14) Bioconductor

Chromvar (Version 1.20) Bioconductor

Spaniel (Version 1.12) Bioconductor

SCENIC+ (Version 1.0.1) <https://github.com/aertslab>

pycisTopic (Version 1.0.3) <https://github.com/aertslab>

pycistarget (Version 1.0.3) <https://github.com/aertslab>

scanpy (Version 1.9.5) <https://github.com/scverse/scanpy> RRID:SCR_018139

cell2Location (Version 0.13) <https://github.com/BayraktarLab/cell2location>

## Installation of Analysis Packages

Cran packages can be intalled with the command.

```{r}
install.packages(PACKAGENAME) 
```

Bioconductor packages can be installed

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(PACKAGENAME)
```

Installation instructions of the other packages can be found in the github links.

The scripts presented here don't require installation to run (no extra installation time).

## scRNAseq Integrated Analysis

1.  developmental_retina/QC.Rmd & Eye/QC.Rmd

    These script performs the filtering of retinal and eye samples and was used to generate the cell numbers and QC metrics in **Table S1**

2.  developmental_retina/Cluster_Analysis.Rmd

    & Eye/QC.Rmd

    These scripts were used to perform clustering analysis of individual samples and generate marker lists (**Table S1**). The script also

    and individual clustering analysis

3.  Eye/Extract_Retina_Clusters.Rmd

    This script was used to extract retinal clusters from whole eye samples

4.  developmental_retina_and_eye_integrated/filter_hb_mt_integrate_find_markers_Retina.R

    This script was used to integrate the datasets together and generate the UMAP in Figure 1A. Clusters were identified and marker lists generated for the integrated analysis (**Table S2** - integrated UMAP ) Retinal marker genes (**Table S4**) were used to assign cell identity.

    ## scRNAseq Pseudotime Analysis

5.  developmental_retina_and_eye_integrated/Pseudotime/subset_trajectory.R

    Clusters were subset from the integrated datasets and reintegrated using this script. UMAP plots are are shown in **Figure 1C** and **Figure S2**.

6.  developmental_retina_and_eye_integrated/Pseudotime/convert_seurat_to_cds.R

    script used to convert filtered seurat object into cds object (the format used by Monocle 3)

7.  developmental_retina_and_eye_integrated/Pseudotime/Pseudotime_Template.Rmd a markdown template used to run trajectory analysis and create plots in **Figure 1E-1G** (NOTE:

    developmental_retina_and_eye_integrated/Pseudotime/StartNode contains the start node used for the pseudotime analysis and developmental_retina_and_eye_integrated/Pseudotime/annotations contains the cell type annotations).

8.  Marker lists were generated from this analysis (**Table S2** - RPCs_T1_T2_T3, RPC_T1_T2_RGCs-HCs_ACs, RPC_T1_T2_HCs_ACs, RPCs_T1_T2_T3).

## ATAC analysis

1.  ATAC_analysis/QC.Rmd - indivdual script used to import scATAC samples, QC and filter data. an overview of the QC metrics can be found in (**Table S5).** The script also clusters each individual dataset.
2.  ATAC_analysis/R_scripts/gene_activity_individual.R - gene activity scores were added to the objects using this script.
3.  Annotations from scRNAseq were transfered from the scRNAseq datasets using this script. This was used to guide the cluster annotation (ATAC_analysis/annotation).
4.  ATAC_analysis/R_scripts/Select_Retina_Cells.R - retinal cells were selected from the individual samples (shown in colour in **Figure S6**)
5.  ATAC_analysis/R_scripts/import_shared_peaks.R script to import individual samples after cellranger reanalysis with a shared peak dataset.
6.  ATAC_analysis/R_scripts/shared_peaks_with_13PCW_sample.R - script to merge data and cluster. Clustered data is shown in **Figure S6C.** The script also generates lists of Differentially accessible peaks (**Table S5**). The script also adds motif information to the object using ChromVAR.
7.   ATAC_analysis/R_scripts/chromvar_motif_enrichment.R script to identify enriched motifs (**Table S7)**
8.  ATAC_analysis/R_scripts/Footprinting_top_motifs.R script to generate footprint profiles for selected TFs.
9.  ATAC_analysis/ATACseq_figures.Rmd - code used to plot the data shown in **Figures 4 and 5.**

## Spatial Analysis

1.  Spatial/individual_results.rmd

2.  Spatial/fetal_retina/hiRes.Rmd

## Network Analysis

Numbered scripts in the SCENIC directory were run to create the SCENIC networks shown in Figure

## Analysis Time

The run time for the scripts ranged from 10 minutes to 24 hours on a server for the full datasets with 80 cores.
