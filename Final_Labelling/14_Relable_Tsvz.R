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

load("output_Tsvz/08_Tsvz_withHarmonyWithCellCycle.RData")
seu_TSVZ <- seuOject_Tsvz_nodoub_withHarmony

#taking res3
seu_TSVZ$seurat_clusters <- seu_TSVZ$SCT_snn_res.0.3

Idents(seu_TSVZ) <- "seurat_clusters"

Labels <- read.table("Relabeled/Tsvz_Cell_labels.txt",header=T,sep="\t")

#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Cell)

##Step 3: Rename cluster with new Labels
seu_TSVZ@active.ident <- plyr::mapvalues(x = seu_TSVZ@active.ident, 
                                         from = current.cluster.ids, 
                                         to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
seu_TSVZ@meta.data$Cell <- seu_TSVZ@active.ident

##Add inferCNV class
seu_TSVZ@meta.data <- seu_TSVZ@meta.data %>%
  mutate(Cell_Class = case_when(grepl("MDM", Cell) ~ "MDM",
                                grepl("Microglia", Cell) ~ "Microglia",
                                grepl("Oligo", Cell) ~ "Oligodendrocytes",
                                grepl("Endo", Cell) ~ "Endothelial",
                                grepl("GBMcc", Cell) ~ "CancerCell",
                                grepl("GBMac", Cell) ~ "GBMac",
                                grepl("GBMmes", Cell) ~ "GBMmes",
                                grepl("GBMopc", Cell) ~ "GBMopc",
                                grepl("GBMnpc", Cell) ~ "GBMnpc",
                                grepl("NPC", Cell) ~ "NPC")) %>%
  mutate(Tumor_Class = case_when(grepl("GBM", Cell) ~ "tumor",.default = "normal"))


dir.create("Relabeled")
save(seu_TSVZ,file = "Relabeled/03_Tsvz_Res03_labelled.RData")

pdf("Relabeled/Tsvz_01_UMAP_Labelled.pdf", width = 6, height = 6)
DimPlot(object = seu_TSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

color <- distinctColorPalette(10)
pdf("Relabeled/Tsvz_01_UMAP_Cell_Class.pdf", width = 8, height = 6)
DimPlot(object = seu_TSVZ, reduction = "umap",
        group.by = "Cell_Class", label = FALSE, 
        pt.size = 0.5, cols = color) #+ theme(legend.position="none")
dev.off()

pdf("Relabeled/Tsvz_02_UMAP_Labelled_Groupby.pdf", width = 12, height = 6)
p1 <- DimPlot(object = seu_TSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seu_TSVZ, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/Tsvz_02_UMAP_Labelled_TumorClass.pdf", width = 12, height = 6)
p1 <- DimPlot(object = seu_TSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seu_TSVZ, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Tumor_Class")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/Tsvz_03_UMAP_Labelled_withclusters.pdf", width = 12, height = 6)
Idents(seu_TSVZ) <- "seurat_clusters"
p1 <- DimPlot(object = seu_TSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(seu_TSVZ) <- "Cell"
p2 <- DimPlot(object = seu_TSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/Tsvz_04_UMAP_Labelled_Splitby_Genotype.pdf", width = 20, height = 15)
DimPlot(object = seu_TSVZ, reduction = "umap", label = TRUE,ncol = 4, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

order <- c("Oligo_1","NPC_1",
           "Microglia_1","MDM_1",
           "Endo_1",
           "GBMopc_1","GBMopc_2","GBMopc_3",
           "GBMmes_1","GBMmes_2",
           "GBMnpc_1","GBMnpc_2",
           "GBMcc_1","GBMcc_2","GBMcc_3")

Idents(seu_TSVZ) <- factor(Idents(seu_TSVZ), levels = order) # reorder the factor based on cell type

ct_markers <- c("MBP","GFAP","MEF2C","PTPRC","SPARCL1",
                "OLIG1","VIM","CD44","SOX2",
                "TRIO","GALNT17")

pdf("Relabeled/Tsvz_06_DotPlot_Celltype.pdf", width = 15, height = 6)
DotPlot(seu_TSVZ, features = ct_markers) + rotate_x_text(angle = 45)
dev.off()

#Density plot
p1 <- plot_density(seu_TSVZ, ct_markers,reduction="umap")
p2 <- p1 + plot_layout(ncol = 4)
pdf("Relabeled/Tsvz_07_Nebulosa_DensityPlots.pdf", width = 15, height = 6)
p2
dev.off()

all_markers_clustID <- presto::wilcoxauc(seu_TSVZ, 'Cell', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "Relabeled/Tsvz_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

color <- c("#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

#reorder cells then plot
seu_TSVZ$Cell <- factor(seu_TSVZ$Cell, levels = order)

prop_cell <-  seu_TSVZ@meta.data %>%
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

ggsave("Relabeled/Tsvz_09_CellProportion_By_Sample.pdf", plot = prop_cell, width = 12, height = 6, units = "in", dpi = 150)

cc_order <- c("Microglia","MDM","NPC","Oligodendrocytes","Endothelial",
              "GBMopc","GBMmes","GBMnpc","CancerCell")

seu_TSVZ$Cell_Class <- factor(seu_TSVZ$Cell_Class, levels = cc_order)

prop_cell <-  seu_TSVZ@meta.data %>%
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

ggsave("Relabeled/TSVZ_09_CellProportion_By_cellclass.pdf", plot = prop_cell, width = 12, height = 6, units = "in", dpi = 150)


seu_TSVZ_slim <- DietSeurat(seu_TSVZ, 
                            counts = TRUE, 
                            data = TRUE, 
                            scale.data = FALSE,
                            assays="RNA",
                            dimreducs = c("pca","umap"))
save(seu_TSVZ_slim, file = "Relabeled/03_Tsvz_Res03_labelled_slim.RData")

modules <- read.csv('/users/SuganyaSubramanian/Sara_FinalPtn/All_Cell_State_Markers.csv')
mod_list = list()
for (col in colnames(modules)) {
  mod_list[col] <- c(modules[col])
}

p2 <- SCpubr::do_EnrichmentHeatmap(sample = seu_TSVZ_slim,
                                   input_gene_list = mod_list,
                                   group.by = "Cell",
                                   viridis_direction = -1,
                                   use_viridis = FALSE,
                                   enforce_symmetry = TRUE)

pdf("Relabeled/Tsvz_10_EnrichHeatMap_Cell.pdf", width = 15, height = 10)  
p2
dev.off()