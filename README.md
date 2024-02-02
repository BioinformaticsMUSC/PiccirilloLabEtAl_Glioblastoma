# PiccirilloEtAl_Glioblastoma
The Tumor Microenvironment Analysis of Human Glioblastoma and the Sub-Ventricular Zone - scRNAseq Data Analysis Code

![](InputData/TOP_PANEL.pdf)
==========================
This repository contains analysis code for the single cell RNA-seq project carried out by researchers at the [Sara G.M. Piccirillo, UNM](https://braintumortranslationallab.com) and [Berto Lab, MUSC](https://bertolab.org/)

Glioblastoma (GBM) clinical management is challenging due to its heterogeneous nature, invasive potential, and poor response to radio- and chemo-therapy. As a result, GBM inevitably recurs and only a minority of patients survive 5 years post-diagnosis. We were the first to show that in the majority of GBM patients, the sub-ventricular zone (SVZ) of the lateral ventricles is a reservoir of cancer stem-like cells (CSCs) that show distinct patterns of treatment resistance when compared to matched CSCs from the tumor mass and contribute to the seeding of the recurrent tumor. However, the sampling and the characterization of the SVZ in GBM patients pose several challenges, as this area is extremely small and needs to be objectively identified during tumor surgical resection.

To overcome this challenge, by using our fluorescence-guided multiple sampling scheme, we built a single-nucleus RNA-sequencing-based microenvironment landscape of the SVZ in GBM patients (T_SVZ). By performing a systematic comparison with tumor mass (T_Mass) samples isolated from the same patients and using two histologically normal SVZ (N_SVZ) samples through a number of computational tools and experimental methods, we identified novel therapeutic vulnerabilities in the T_SVZ of GBM patients.
## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: 

## Access the data with an app:

Here a webapp to analyze the data:

[Glioblastoma-SVZs](https://bioinformatics-musc.shinyapps.io/sara_piccirillo_glioblastoma/)

## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`InputData`](InputData/) | Input data of the initial processing and quality check. ||
| [`Initial_Clustering`](Initial_Clustering/) | Code and Output data of the initial clustering and integration. | 01_Seurat_Clustering_Nsvz.R \ 02_Seurat_Clustering_TM.R \ 03_Seurat_Clustering_Tsvz.R \ 04_Doubletting_Nsvz.R \ 05_Doubletting_TM.R \ 06_Doubletting_Tsvz.R \ 07_CellCycleScore.R \ 08_FindMarkers.R \ 09_InferCNV_Tsvz.R \ 10_InferCNV_TM \ 11_EnrichmentPlot.R |
| [`Final_Labelling`](output_reclust/) | Code and Output data of the Final Labelling | 12_Relable_Nsvz.R \ 13_Relable_TM.R \ 14_Relable_Tsvz.R |
| [`DGE`](output_dge/) | Output data of the DGE analyses. |  |
| [`ShinyApp`](ShinyApp/) | Output of the ShinyApp. | |

