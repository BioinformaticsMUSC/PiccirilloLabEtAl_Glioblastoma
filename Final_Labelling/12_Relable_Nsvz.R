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

load("output_Nsvz/08_Nsvz_withHarmonyWithCellCycle.RData")

seu_NSVZ <- seuOject_Nsvz_nodoub_withHarmony

current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12")
new.cluster.ids <- c("Oligo_1","Oligo_1","OPC_1","Astro_1","Microglia_1","Oligo_1","OPC_2",
                     "Neurons_1","Oligo_1","NPC_1","OPC_3","Endo_1","Epen_1")

##Step 3: Rename cluster with new Labels
seu_NSVZ@active.ident <- plyr::mapvalues(x = seu_NSVZ@active.ident, 
                                         from = current.cluster.ids, 
                                         to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
seu_NSVZ@meta.data$Cell <- seu_NSVZ@active.ident

##Save new Seurat data
##Add inferCNV class
seu_NSVZ$Tumor_Class <- "normal"

seu_NSVZ@meta.data <- seu_NSVZ@meta.data %>%
  mutate(Cell_Class = case_when(grepl("Oligo", Cell) ~ "Oligodendrocytes",
                                grepl("OPC", Cell) ~ "OPC",
                                grepl("Astro", Cell) ~ "Astrocytes",
                                grepl("Microglia", Cell) ~ "Microglia",
                                grepl("Neuron", Cell) ~ "Neurons",
                                grepl("NPC", Cell) ~ "NPC",
                                grepl("Endo", Cell) ~ "Endothelial",
                                grepl("Epen", Cell) ~ "Ependymal"))


dir.create("Relabeled")

save(seu_NSVZ,file = "Relabeled/01_NSvz_Res03_labelled.RData")

pdf("Relabeled/NSvz_01_UMAP_Labelled.pdf", width = 6, height = 6)
DimPlot(object = seu_NSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

color <- distinctColorPalette(10)
pdf("Relabeled/NSvz_01_UMAP_Cell_Class.pdf", width = 7, height = 6)
DimPlot(object = seu_NSVZ, reduction = "umap",
        group.by = "Cell_Class", label = FALSE, 
        pt.size = 0.5, cols = color) #+ theme(legend.position="none")
dev.off()

pdf("Relabeled/NSvz_02_UMAP_Labelled_Groupby.pdf", width = 12, height = 6)
p1 <- DimPlot(object = seu_NSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seu_NSVZ, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/NSvz_02_UMAP_Labelled_TumorClass.pdf", width = 12, height = 6)
p1 <- DimPlot(object = seu_NSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seu_NSVZ, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Tumor_Class")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/NSvz_03_UMAP_Labelled_withclusters.pdf", width = 12, height = 6)
Idents(seu_NSVZ) <- "seurat_clusters"
p1 <- DimPlot(object = seu_NSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(seu_NSVZ) <- "Cell"
p2 <- DimPlot(object = seu_NSVZ, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("Relabeled/NSvz_04_UMAP_Labelled_Splitby_Genotype.pdf", width = 12, height = 6)
DimPlot(object = seu_NSVZ, reduction = "umap", label = TRUE,ncol = 2, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

order <- c("Astro_1","Microglia_1",
           "Oligo_1","OPC_1","OPC_2","OPC_3",
           "Neurons_1","NPC_1",
           "Endo_1","Epen_1")

Idents(seu_NSVZ) <- factor(Idents(seu_NSVZ), levels = order) # reorder the factor based on cell type

ct_markers <- c("GFAP", #Astrocytes
                "MEF2C", #Microglia
                "MBP", #Oligodendrocytes
                "SOX5",#OPC
                "SYT1",#Neurons"DLX2","MEG3",
                "CDK4",#NPC" SOX2","NES",
                "CLDN5",#Endo #"CLDN5", "IFI27",
                "S100B")#MSC

pdf("Relabeled/Nsvz_06_DotPlot_Celltype.pdf", width = 15, height = 6)
DotPlot(seu_NSVZ, features = ct_markers) + rotate_x_text(angle = 45)
dev.off()

### Density Plot
p1 <- plot_density(seu_NSVZ, ct_markers,reduction="umap")
p2 <- p1 + plot_layout(ncol = 4)
pdf("Relabeled/Nsvz_07_Nebulosa_DensityPlots.pdf", width = 15, height = 6)
p2
dev.off()

#Find Markers by cell lables

all_markers_clustID <- presto::wilcoxauc(seu_NSVZ, 'Cell', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "Relabeled/Nsvz_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

color <- c("#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

#reorder cells then plot
seu_NSVZ$Cell <- factor(seu_NSVZ$Cell, levels = order)

prop_cell <-  seu_NSVZ@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(Cell) %>%
  dplyr::group_by(Cell,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(Cell) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(Cell) %>%
  ggbarplot("Cell", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="right") +
  xlab("")

ggsave("Relabeled/NSvz_09_CellProportion_By_Sample.pdf", plot = prop_cell, width = 12, height = 6, units = "in", dpi = 150)
