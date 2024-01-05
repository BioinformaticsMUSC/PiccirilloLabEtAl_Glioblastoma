suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridis)
  library(scGate)
  library(Seurat)
  library(Nebulosa)
  library(RColorBrewer)
  library(ggpubr)
})

setwd("/users/SuganyaSubramanian/Sara_Final")

load("output_Nsvz/08_Nsvz_withHarmonyWithCellCycle.RData")
load("output_Tsvz/08_Tsvz_withHarmonyWithCellCycle.RData")
load("output_TumorMass/08_TM_withHarmonyWithCellCycle.RData")

seu_NSVZ <- seuOject_Nsvz_nodoub_withHarmony
seu_TSVZ <- seuOject_Tsvz_nodoub_withHarmony
seu_TM <- seuOject_TM_nodoub_withHarmony

#markers
all_markers_clustID <- presto::wilcoxauc(seu_NSVZ, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_Nsvz/NSVZ_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

all_markers_clustID <- presto::wilcoxauc(seu_TSVZ, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_Tsvz/TSVZ_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

all_markers_clustID <- presto::wilcoxauc(seu_TM, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_TumorMass/TM_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")