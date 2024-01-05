suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  library(randomcoloR)
  library(harmony)
})

source("/home/subramanian/R/x86_64-pc-linux-gnu-library/4.1/Utils.R")
setwd("/zfs/musc3/Sara/01Final")

load("output_TumorMass/04_seuObject_Res0.5.RData")
seuOject_Integrated_TM <- seuObject_integrated

#Convert into Single Cell Experiment object
sce <- as.SingleCellExperiment(seuOject_Integrated_TM)

# Doubleting by genotype
sce <- scDblFinder(sce,
                   samples="Genotype", 
                   BPPARAM=MulticoreParam(3),
                   nfeatures = 3000,
                   dims = 30,
                   dbr.sd = 1)

# add in the new column "Doublets" to make if it doublet or singlet
seuOject_Integrated_TM@meta.data$Doublets <- sce$scDblFinder.class

save(seuOject_Integrated_TM, file = "output_TumorMass/05_SeuratObj_withDoublets.RData")


df <- seuOject_Integrated_TM@meta.data %>% as.data.frame()
tmp <- table(df$Doublets,df$Genotype) %>% as.data.frame()

pdf("output_TumorMass/08_GeneExp_Singlet-Doublet_Genotype.pdf",width = 25,height = 6)
ggplot(tmp,aes(x=Var2, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  #scale_fill_manual(values=c("orange2", "green3")) +
  geom_text(aes(label=Freq), size=3.5) +
  ggtitle("Singlet/Doublet by Genotype") +
  xlab("Sample") +
  ylab("Gene Expression") +
  theme(legend.position="none")
dev.off()

#Subset only singlets (remove Doublets)
seuOject_TM_nodoub <- subset(seuOject_Integrated_TM, subset = Doublets == "singlet")

seuOject_TM_nodoub <- processing_seurat_sctransform(seuOject_TM_nodoub, 
                                                      vars_to_regress = c("nCount_RNA","pMito","pRibo"), 
                                                      npcs = 30, 
                                                      res = c(0.1,0.2,0.3,0.4,0.5,0.6))

DefaultAssay(seuOject_TM_nodoub) <- "RNA"
seuOject_TM_nodoub <- NormalizeData(object = seuOject_TM_nodoub, 
                                      normalization.method = "LogNormalize", 
                                      scale.factor = 10000)


save(seuOject_TM_nodoub, file = "output_TumorMass/06_SeuratObj_SCT_30pcs_06res_NoDoublet_Reclustered.RData")

seuOject_TM_nodoub_slim <- DietSeurat(seuOject_TM_nodoub, 
                                        counts = TRUE, 
                                        data = TRUE, 
                                        scale.data = FALSE,
                                        assays="RNA",
                                        dimreducs = c("pca","umap"))

save(seuOject_TM_nodoub_slim, file = "output_TumorMass/06_SeuratObj_SCT_30pcs_06res_NoDoublet_Reclustered_slim.RData")

load("output_TumorMass/06_SeuratObj_SCT_30pcs_06res_NoDoublet_Reclustered_slim.RData")
pdf("output_TumorMass/09_UMAP_TM_NoDoublets_Reclustered_Genotype.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuOject_TM_nodoub_slim, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuOject_TM_nodoub_slim, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

prop_cell <-  seuOject_TM_nodoub_slim@meta.data %>%
  as.data.frame() %>%
  arrange(seurat_clusters) %>%
  group_by(seurat_clusters,Genotype) %>% 
  summarise(cnt = n()) %>% 
  group_by(seurat_clusters) %>%  
  mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  group_by(seurat_clusters) %>%
  ggbarplot("seurat_clusters", "percent",
            fill = "Genotype", color = "Genotype", #palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_TumorMass/10_CellProportion_By_Genotype_TM.pdf", plot = prop_cell, width = 25, height = 15, units = "in", dpi = 150)


pdf("output_TumorMass/11_UMAP_Nodoublets_SplitBy_Genotype.pdf", width = 15, height = 30)
DimPlot(seuOject_TM_nodoub_slim,split.by = "Genotype",ncol = 3) + theme(legend.position="none")
dev.off()

library(harmony)
#seuOject_TM_nodoub_slim object looses PCA reduction
seuOject_TM_nodoub_withHarmony <- seuOject_TM_nodoub

seuOject_TM_nodoub_withHarmony <- RunHarmony(seuOject_TM_nodoub_withHarmony, assay.use="RNA", group.by.vars = "Set")
seuOject_TM_nodoub_withHarmony <- RunUMAP(seuOject_TM_nodoub_withHarmony, reduction = "harmony", dims = 1:30)
seuOject_TM_nodoub_withHarmony <- FindNeighbors(seuOject_TM_nodoub_withHarmony, reduction = "harmony", dims = 1:30) 
seuOject_TM_nodoub_withHarmony <- FindClusters(seuOject_TM_nodoub_withHarmony, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), 
                                               algorithm = 1, n.iter = 1000,
                                               graph.name = "SCT_snn")

Idents(object = seuOject_TM_nodoub_withHarmony) <- "SCT_snn_res.0.3"

pdf("output_TumorMass/12_UMAP_TM_NoDoublets_Reclustered_afterHarmonyWithSet_Res0.3.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuOject_TM_nodoub_withHarmony, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuOject_TM_nodoub_withHarmony, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

Idents(object = seuOject_TM_nodoub_withHarmony) <- "SCT_snn_res.0.1"

pdf("output_TumorMass/12_UMAP_TM_NoDoublets_Reclustered_afterHarmonyWithSet_Res0.1.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuOject_TM_nodoub_withHarmony, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuOject_TM_nodoub_withHarmony, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

Idents(object = seuOject_TM_nodoub_withHarmony) <- "SCT_snn_res.0.2"

pdf("output_TumorMass/12_UMAP_TM_NoDoublets_Reclustered_afterHarmonyWithSet_Res0.2.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuOject_TM_nodoub_withHarmony, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuOject_TM_nodoub_withHarmony, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

seuOject_TM_nodoub_withHarmony$seurat_clusters <- seuOject_TM_nodoub_withHarmony$SCT_snn_res.0.3
Idents(seuOject_TM_nodoub_withHarmony) <- "seurat_clusters"

save(seuOject_TM_nodoub_withHarmony, file = "output_TumorMass/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")

prop_cell <-  seuOject_TM_nodoub_withHarmony@meta.data %>%
  as.data.frame() %>%
  arrange(seurat_clusters) %>%
  group_by(seurat_clusters,Genotype) %>% 
  summarise(cnt = n()) %>% 
  group_by(seurat_clusters) %>%  
  mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  group_by(seurat_clusters) %>%
  ggbarplot("seurat_clusters", "percent",
            fill = "Genotype", color = "Genotype", #palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_TumorMass/10_CellProportion_By_Genotype_TM_withHarmonyBySet.pdf", plot = prop_cell, width = 25, height = 15, units = "in", dpi = 150)
