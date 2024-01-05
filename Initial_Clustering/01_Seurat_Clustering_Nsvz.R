suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(sctransform)
  library(hdf5r)
  library(ggrastr)
  library(clustree)
  library(cowplot)
})

setwd("/Users/SuganyaSubramanian/Sara_Final")

pmBrE <- Read10X_h5("cellbender_out_Nsvz/Set1_E_cellbender_filtered.h5", use.names = TRUE, unique.features = TRUE)
MIS2D <- Read10X_h5("cellbender_out_Nsvz/Set3_MIS2D_filtered.h5", use.names = TRUE, unique.features = TRUE)

#Filter for protein coding genes
load("/Users/SuganyaSubramanian/HgProteinCodingGenesSara.rda")
pmBrE <- pmBrE[rownames(pmBrE)%in%ptn_genes,]
MIS2D <- MIS2D[rownames(MIS2D)%in%ptn_genes,]


pmBrE_obj <- CreateSeuratObject(counts = pmBrE,min.features = 100)
pmBrE_obj@meta.data$Region <- "E"
pmBrE_obj@meta.data$PatientID <- "pmBrain"
pmBrE_obj@meta.data$Sex <- "Female"
pmBrE_obj@meta.data$Age <- "30"
pmBrE_obj@meta.data$Set <- "Set1"
pmBrE_obj@meta.data$Label <- "Putative-normal-sub-ventricular-zone"
pmBrE_obj@meta.data$Abbreviation <- "N-SVZ"

MIS2D_obj <- CreateSeuratObject(counts = MIS2D,min.features = 100)
MIS2D_obj@meta.data$Region <- "D"
MIS2D_obj@meta.data$PatientID <- "MIS2"
MIS2D_obj@meta.data$Sex <- "Female"
MIS2D_obj@meta.data$Age <- "73"
MIS2D_obj@meta.data$Set <- "Set3"
MIS2D_obj@meta.data$Label <- "Putative-normal-sub-ventricular-zone"
MIS2D_obj@meta.data$Abbreviation <- "N-SVZ"


##### Merge all ########
seuObject <- merge(pmBrE_obj,
                   y = c(MIS2D_obj),
                   add.cell.ids = c("pmBrE","MIS2D"),
                   project = "Sara_final_Nsvz")

########## Mito column ##########
seuObject[["pMito"]] <- PercentageFeatureSet(seuObject, pattern = "^MT-")

seuObject[["pRibo"]] <- PercentageFeatureSet(seuObject, pattern = "^RP[SL]")

seuObject@meta.data <- seuObject@meta.data %>%
  rownames_to_column("Cell") %>%
  mutate(Genotype = sapply(X = strsplit(colnames(seuObject), split = "_"), FUN = "[", 1)) %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
  column_to_rownames("Cell")

dir.create("output_Nsvz")

save(seuObject,file="output_Nsvz/01_SeuObj_Nsvz_Unfilt.RData")

#QC plot 1
pdf("output_Nsvz/01_Quality_Control_plots_Genotype.pdf", width=6,height=5)
feats <- c("nUMI", "nGene", "pMito")
VlnPlot(seuObject, group.by = "Genotype", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
dev.off()

### Filter Mito ######
mito_filtered <- seuObject@assays$RNA@counts[-grep("^MT-",rownames(seuObject@assays$RNA@counts)),]

#Initialize the Seurat object with the raw (non-normalized data).
seuObject_final <- CreateSeuratObject(counts = mito_filtered, project = "Sara_set2")

## Add pMito info from meta data for all cells before filtering
metaAll <- as.data.frame(seuObject@meta.data)
seuObject_final <- AddMetaData(object = seuObject_final, metadata = as.data.frame(seuObject@meta.data))
seuObject_final@meta.data$nCount_RNA <- NULL
seuObject_final@meta.data$nFeature_RNA <- NULL
seuObject_filt <- subset(x = seuObject_final, subset = nUMI < 25000 & pMito < 5)

save(seuObject_filt,file="output_Nsvz/02_SeuObj_Nsvz_MtFiltered_u25000_m5.RData")

df <- seuObject_filt@meta.data %>% as.data.frame()
tmp <- table(df$Genotype) %>% as.data.frame()

pdf("output_Nsvz/04-GeneExp_barplot_Genotype.pdf",width = 5,height = 7)
ggplot(tmp,aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") +
  #scale_fill_manual(values=c("orange2", "green3")) +
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  ggtitle("Gene Expression by Genotype") +
  xlab("Sample") +
  ylab("Gene Expression") +
  theme(legend.position="none")
dev.off()

######### STEP 10 Data Integration ###############

#load("output_Nsvz/02_SeuObj_Nsvz_MtFiltered_u25000_m5.RData")
seuObject_split <- SplitObject(seuObject_filt, split.by = "Genotype")

seuObject_split <- seuObject_split[c("pmBrE","MIS2D")]

for (i in 1:length(seuObject_split)) {
  seuObject_split[[i]] <- SCTransform(seuObject_split[[i]],
                                      vars.to.regress = c("nUMI","pMito","pRibo"),
                                      verbose = FALSE)
}

integ_features <- SelectIntegrationFeatures(object.list = seuObject_split,
                                            nfeatures = 4000)

seuObject_split <- PrepSCTIntegration(object.list = seuObject_split,
                                      anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = seuObject_split,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)

seuObject_integrated <- IntegrateData(
  anchorset = integ_anchors,
  new.assay.name = "integrated",
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 100,
  sd.weight = 1,
  eps = 0.5,
  verbose = TRUE
)

DefaultAssay(seuObject_integrated) <- "integrated"

seuObject_integrated <- RunPCA(object = seuObject_integrated,
                               features=NULL,
                               weight.by.var = TRUE,
                               ndims.print = 1:5,
                               nfeatures.print = 30,
                               npcs = 30,
                               reduction.name = "pca")

seuObject_integrated <- FindNeighbors(object = seuObject_integrated,
                                      reduction = "pca",
                                      dims = 1:30,
                                      nn.eps = 0.5)


seuObject_integrated <- FindClusters(object = seuObject_integrated,
                                     resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2),
                                     algorithm = 1,
                                     n.iter = 1000)

save(seuObject_integrated, file = "output_Nsvz/04_Integrated_allRes.RData")

pdf("output_Nsvz/05_Clustree-Data_Integrated.pdf", width = 15, height = 6)
clustree(seuObject_integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Select the RNA counts slot to be the default assay
DefaultAssay(seuObject_integrated) <- "RNA"

Idents(object = seuObject_integrated) <- "integrated_snn_res.0.5"

seuObject_integrated <- RunUMAP(object = seuObject_integrated,
                                reduction = "pca",
                                dims = 1:30)

seuObject_integrated <- NormalizeData(object = seuObject_integrated,
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

pdf("output_Nsvz/06_UMAP-Data_Integrated_Genotype.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_Nsvz/07_UMAP-Data_Integrated_splitby_Genotype.pdf", width = 10, height = 5)
DimPlot(object = seuObject_integrated, reduction = "umap", ncol = 3, label = FALSE, pt.size = 0.5, split.by = "Genotype")
dev.off()

save(seuObject_integrated, file = "output_Nsvz/04_seuObject_Res0.5.RData")