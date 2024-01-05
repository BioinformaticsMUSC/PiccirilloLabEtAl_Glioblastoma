suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(speckle)
  library(magrittr)
  library(broom)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  library(RColorBrewer)
  library(MetBrewer)
  library(scales)
  library(ComplexHeatmap)
  library(randomcoloR)
  library(png)
  library(patchwork)
  library(randomcoloR)
  library(Nebulosa)
})
load("output_TumorMass/08_TM_withHarmonyWithCellCycle.RData")

seu_TM <- seuOject_TM_nodoub_withHarmony

seu_TM <- FindClusters(seu_TM, resolution = c(0.7,0.8,0.9,1.0,1.1,1.2,1.3), 
                    algorithm = 1, n.iter = 1000,
                    graph.name = "SCT_snn")

Idents(seu_TM) <- "SCT_snn_res.1.1"

seu_TM$seurat_clusters <- seu_TM$SCT_snn_res.1.1

dir.create("Relabeled")

Idents(seu_TM) <- "seurat_clusters"

Labels <- read.table("Relabeled/TM_Cell_labels.txt",header=T,sep="\t")

#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Cell)

##Step 3: Rename cluster with new Labels
seu_TM@active.ident <- plyr::mapvalues(x = seu_TM@active.ident, 
                                         from = current.cluster.ids, 
                                         to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
seu_TM@meta.data$Cell <- seu_TM@active.ident

seu_TM@meta.data <- seu_TM@meta.data %>%
  mutate(Cell_Class = case_when(grepl("MDM", Cell) ~ "MDM",
                                grepl("Microglia", Cell) ~ "Microglia",
                                grepl("Oligo", Cell) ~ "Oligodendrocytes",
                                grepl("Endo", Cell) ~ "Endothelial",
                                grepl("GBMcc", Cell) ~ "CancerCell",
                                grepl("GBMac", Cell) ~ "GBMac",
                                grepl("GBMmes", Cell) ~ "GBMmes",
                                grepl("GBMopc", Cell) ~ "GBMopc",
                                grepl("GBMnpc", Cell) ~ "GBMnpc",
                                grepl("Neuron", Cell) ~ "Neurons")) %>%
  mutate(Tumor_Class = case_when(grepl("GBM", Cell) ~ "tumor",.default = "normal"))

save(seu_TM,file = "Relabeled/02_TM_Res11_labelled.RData")

pdf("Relabeled/TM_01_UMAP_Labelled.pdf", width = 6, height = 6)
DimPlot(object = seu_TM, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

color <- distinctColorPalette(10)
pdf("Relabeled/TM_01_UMAP_Cell_Class.pdf", width = 8, height = 6)
DimPlot(object = seu_TM, reduction = "umap",
        group.by = "Cell_Class", label = FALSE, 
        pt.size = 0.5, cols = color) #+ theme(legend.position="none")
dev.off()

pdf("Relabeled/TM_02_UMAP_Labelled_Groupby.pdf", width = 12, height = 6)
p1 <- DimPlot(object = seu_TM, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seu_TM, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/TM_02_UMAP_Labelled_TumorClass.pdf", width = 12, height = 6)
p1 <- DimPlot(object = seu_TM, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seu_TM, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Tumor_Class")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/TM_03_UMAP_Labelled_withclusters.pdf", width = 12, height = 6)
Idents(seu_TM) <- "seurat_clusters"
p1 <- DimPlot(object = seu_TM, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(seu_TM) <- "Cell"
p2 <- DimPlot(object = seu_TM, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/TM_04_UMAP_Labelled_Splitby_Genotype.pdf", width = 20, height = 15)
DimPlot(object = seu_TM, reduction = "umap", label = TRUE,ncol = 4, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

order <- c("Oligo_1","Oligo_2",
           "Microglia_1","MDM_1","MDM_2",
           "Endo_1","Neuron_1","Neuron_2",
           "GBMnpc_1","GBMnpc_2",
           "GBMac_1","GBMac_2","GBMac_3","GBMac_4","GBMac_5",
           "GBMmes_1","GBMmes_2",
           "GBMopc_1","GBMopc_2",
           "GBMcc_1","GBMcc_2","GBMcc_3")

Idents(seu_TM) <- factor(Idents(seu_TM), levels = order) # reorder the factor based on cell type

ct_markers <- c("MBP","MEF2C","PTPRC","SPARCL1","SYT1","IL1RAPL2",
                "GFAP","AQP4","CD44",
                "VIM","DCN","SOX6","SOX5",
                "EGFR","CNTN5","GLIS3","VEGFA","TOP2A")

pdf("Relabeled/TM_06_DotPlot_Celltype.pdf", width = 15, height = 6)
DotPlot(seu_TM, features = ct_markers) + rotate_x_text(angle = 45)
dev.off()

#Density Plot
p1 <- plot_density(seu_TM, ct_markers,reduction="umap")
p2 <- p1 + plot_layout(ncol = 4)
pdf("Relabeled/TM_07_Nebulosa_DensityPlots.pdf", width = 15, height = 6)
p2
dev.off()

#Find Markers by cell lables

all_markers_clustID <- presto::wilcoxauc(seu_TM, 'Cell', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "Relabeled/TM_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

color <- c("#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

#reorder cells then plot
seu_TM$Cell <- factor(seu_TM$Cell, levels = order)

prop_cell <-  seu_TM@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(Cell) %>%
  dplyr::group_by(Cell,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(Cell) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(Cell) %>%
  ggbarplot("Cell", "percent",
            fill = "Genotype", color = "Genotype", #palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="right") +
  xlab("")

ggsave("Relabeled/TM_09_CellProportion_By_Sample.pdf", plot = prop_cell, width = 12, height = 6, units = "in", dpi = 150)

cc_order <- c("Microglia","Neurons","MDM","Oligodendrocytes","Endothelial",
              "GBMac","GBMopc","GBMmes","GBMnpc","CancerCell")
seu_TM$Cell_Class <- factor(seu_TM$Cell_Class, levels = cc_order)

prop_cell <-  seu_TM@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(Cell_Class) %>%
  dplyr::group_by(Cell_Class,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(Cell_Class) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(Cell_Class) %>%
  ggbarplot("Cell_Class", "percent",
            fill = "Genotype", color = "Genotype", #palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="right") +
  xlab("")

ggsave("Relabeled/TM_09_CellProportion_By_cellclass.pdf", plot = prop_cell, width = 12, height = 6, units = "in", dpi = 150)

seu_TM_slim <- DietSeurat(seu_TM, 
                          counts = TRUE, 
                          data = TRUE, 
                          scale.data = FALSE,
                          assays="RNA",
                          dimreducs = c("pca","umap"))
save(seu_TM_slim, file = "Relabeled/02_TM_Res11_labelled_slim.RData")

modules <- read.csv('All_Cell_State_Markers.csv')
mod_list = list()
for (col in colnames(modules)) {
  mod_list[col] <- c(modules[col])
}

p3 <- SCpubr::do_EnrichmentHeatmap(sample = seu_TM_slim,
                                   input_gene_list = mod_list,
                                   group.by = "Cell",
                                   viridis_direction = -1,
                                   use_viridis = FALSE,
                                   enforce_symmetry = TRUE)

pdf("Relabeled/TM_10_EnrichHeatMap_Cell.pdf", width = 15, height = 10)  
p3
dev.off()






