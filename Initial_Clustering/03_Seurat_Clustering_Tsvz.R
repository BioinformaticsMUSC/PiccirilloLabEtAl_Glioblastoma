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

setwd("/zfs/musc3/Sara/01Final")

G10B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_02GBM10B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G12B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_05GBM12B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G16B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_08GBM16B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G17B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_11GBM17B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G20B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_17GBM20B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G22B <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_GBM22B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G7B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_19GBM7B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G8B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_22GBM8B_filtered.h5", use.names = TRUE, unique.features = TRUE)
G9B <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_25GBM9B_filtered.h5", use.names = TRUE, unique.features = TRUE)
MIS1B <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_MIS1B_filtered.h5", use.names = TRUE, unique.features = TRUE)
MIS2B <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_MIS2B_filtered.h5", use.names = TRUE, unique.features = TRUE)
MIS3B <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_MIS3B_filtered.h5", use.names = TRUE, unique.features = TRUE)

load("HgProteinCodingGenesSara.rda")

G10B <- G10B[rownames(G10B)%in%ptn_genes,]
G12B <- G12B[rownames(G12B)%in%ptn_genes,]
G16B <- G16B[rownames(G16B)%in%ptn_genes,]
G17B <- G17B[rownames(G17B)%in%ptn_genes,]
G20B <- G20B[rownames(G20B)%in%ptn_genes,]
G22B <- G22B[rownames(G22B)%in%ptn_genes,]
G7B <- G7B[rownames(G7B)%in%ptn_genes,]
G8B <- G8B[rownames(G8B)%in%ptn_genes,]
G9B <- G9B[rownames(G9B)%in%ptn_genes,]
MIS1B <- MIS1B[rownames(MIS1B)%in%ptn_genes,]
MIS2B <- MIS2B[rownames(MIS2B)%in%ptn_genes,]
MIS3B <- MIS3B[rownames(MIS3B)%in%ptn_genes,]


G10B_obj <- CreateSeuratObject(counts = G10B,min.features = 100)
G10B_obj@meta.data$Region <- "B"
G10B_obj@meta.data$PatientID <- "GBM10"
G10B_obj@meta.data$Sex <- "Male"
G10B_obj@meta.data$Age <- "64"
G10B_obj@meta.data$Set <- "Set2"
G10B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G10B_obj@meta.data$Abbreviation <- "T-SVZ"

G12B_obj <- CreateSeuratObject(counts = G12B,min.features = 100)
G12B_obj@meta.data$Region <- "B"
G12B_obj@meta.data$PatientID <- "GBM12"
G12B_obj@meta.data$Sex <- "Female"
G12B_obj@meta.data$Age <- "66"
G12B_obj@meta.data$Set <- "Set2"
G12B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G12B_obj@meta.data$Abbreviation <- "T-SVZ"

G16B_obj <- CreateSeuratObject(counts = G16B,min.features = 100)
G16B_obj@meta.data$Region <- "B"
G16B_obj@meta.data$PatientID <- "GBM16"
G16B_obj@meta.data$Sex <- "Male"
G16B_obj@meta.data$Age <- "77"
G16B_obj@meta.data$Set <- "Set2"
G16B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G16B_obj@meta.data$Abbreviation <- "T-SVZ"

G17B_obj <- CreateSeuratObject(counts = G17B,min.features = 100)
G17B_obj@meta.data$Region <- "B"
G17B_obj@meta.data$PatientID <- "GBM17"
G17B_obj@meta.data$Sex <- "Female"
G17B_obj@meta.data$Age <- "56"
G17B_obj@meta.data$Set <- "Set2"
G17B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G17B_obj@meta.data$Abbreviation <- "T-SVZ"

G20B_obj <- CreateSeuratObject(counts = G20B,min.features = 100)
G20B_obj@meta.data$Region <- "B"
G20B_obj@meta.data$PatientID <- "GBM20"
G20B_obj@meta.data$Sex <- "Male"
G20B_obj@meta.data$Age <- "71"
G20B_obj@meta.data$Set <- "Set2"
G20B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G20B_obj@meta.data$Abbreviation <- "T-SVZ"

G22B_obj <- CreateSeuratObject(counts = G22B,min.features = 100)
G22B_obj@meta.data$Region <- "B"
G22B_obj@meta.data$PatientID <- "GBM22"
G22B_obj@meta.data$Sex <- "Female"
G22B_obj@meta.data$Age <- "80"
G22B_obj@meta.data$Set <- "Set3"
G22B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G22B_obj@meta.data$Abbreviation <- "T-SVZ"

G7B_obj <- CreateSeuratObject(counts = G7B,min.features = 100)
G7B_obj@meta.data$Region <- "B"
G7B_obj@meta.data$PatientID <- "GBM7"
G7B_obj@meta.data$Sex <- "Male"
G7B_obj@meta.data$Age <- "77"
G7B_obj@meta.data$Set <- "Set2"
G7B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G7B_obj@meta.data$Abbreviation <- "T-SVZ"

G8B_obj <- CreateSeuratObject(counts = G8B,min.features = 100)
G8B_obj@meta.data$Region <- "B"
G8B_obj@meta.data$PatientID <- "GBM8"
G8B_obj@meta.data$Sex <- "Male"
G8B_obj@meta.data$Age <- "64"
G8B_obj@meta.data$Set <- "Set2"
G8B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G8B_obj@meta.data$Abbreviation <- "T-SVZ"

G9B_obj <- CreateSeuratObject(counts = G9B,min.features = 100)
G9B_obj@meta.data$Region <- "B"
G9B_obj@meta.data$PatientID <- "GBM9"
G9B_obj@meta.data$Sex <- "Female"
G9B_obj@meta.data$Age <- "71"
G9B_obj@meta.data$Set <- "Set2"
G9B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
G9B_obj@meta.data$Abbreviation <- "T-SVZ"

MIS1B_obj <- CreateSeuratObject(counts = MIS1B,min.features = 100)
MIS1B_obj@meta.data$Region <- "B"
MIS1B_obj@meta.data$PatientID <- "MIS1"
MIS1B_obj@meta.data$Sex <- "Male"
MIS1B_obj@meta.data$Age <- "21"
MIS1B_obj@meta.data$Set <- "Set3"
MIS1B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
MIS1B_obj@meta.data$Abbreviation <- "T-SVZ"

MIS2B_obj <- CreateSeuratObject(counts = MIS2B,min.features = 100)
MIS2B_obj@meta.data$Region <- "B"
MIS2B_obj@meta.data$PatientID <- "MIS2"
MIS2B_obj@meta.data$Sex <- "Female"
MIS2B_obj@meta.data$Age <- "73"
MIS2B_obj@meta.data$Set <- "Set3"
MIS2B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
MIS2B_obj@meta.data$Abbreviation <- "T-SVZ"

MIS3B_obj <- CreateSeuratObject(counts = MIS3B,min.features = 100)
MIS3B_obj@meta.data$Region <- "B"
MIS3B_obj@meta.data$PatientID <- "MIS3"
MIS3B_obj@meta.data$Sex <- "Female"
MIS3B_obj@meta.data$Age <- "51"
MIS3B_obj@meta.data$Set <- "Set3"
MIS3B_obj@meta.data$Label <- "Tumor_Sub-Ventricular_Zone"
MIS3B_obj@meta.data$Abbreviation <- "T-SVZ"

##### Merge all ########
seuObject <- merge(G10B_obj,
                   y = c(G12B_obj,G16B_obj,G17B_obj,G20B_obj,
                         G22B_obj,G7B_obj,G8B_obj,G9B_obj,MIS1B_obj,MIS2B_obj,MIS3B_obj),
                   add.cell.ids = c("GBM10B","GBM12B","GBM16B","GBM17B",
                                    "GBM20B","GBM22B","GBM7B","GBM8B","GBM9B","MIS1B",
                                    "MIS2B","MIS3B"),
                   project = "Sara_final_Tsvz")

########## Mito column ##########
seuObject[["pMito"]] <- PercentageFeatureSet(seuObject, pattern = "^MT-")

seuObject[["pRibo"]] <- PercentageFeatureSet(seuObject, pattern = "^RP[SL]")

seuObject@meta.data <- seuObject@meta.data %>%
  rownames_to_column("Cell") %>%
  mutate(Genotype = sapply(X = strsplit(colnames(seuObject), split = "_"), FUN = "[", 1)) %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
  column_to_rownames("Cell")

dir.create("output_Tsvz")

save(seuObject,file="output_Tsvz/01_SeuObj_Tsvz_Unfilt.RData")

#QC plot 1
pdf("output_Tsvz/01_Quality_Control_plots_Genotype.pdf", width=25,height=5)
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

save(seuObject_filt,file="output_Tsvz/02_SeuObj_Tsvz_MtFiltered_u25000_m5.RData")

df <- seuObject_filt@meta.data %>% as.data.frame()
tmp <- table(df$Genotype) %>% as.data.frame()

pdf("output_Tsvz/04-GeneExp_barplot_Genotype.pdf",width = 25,height = 7)
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
seuObject_split <- SplitObject(seuObject_filt, split.by = "Genotype")

seuObject_split <- seuObject_split[c("GBM10B","GBM12B","GBM16B","GBM17B",
                                     "GBM20B","GBM22B","GBM7B","GBM8B","GBM9B","MIS1B",
                                     "MIS2B","MIS3B")]

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

save(seuObject_integrated, file = "output_Tsvz/04_Integrated_allRes.RData")

pdf("output_Tsvz/05_Clustree-Data_Integrated.pdf", width = 15, height = 6)
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

pdf("output_Tsvz/06_UMAP-Data_Integrated_Genotype.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_Tsvz/07_UMAP-Data_Integrated_splitby_Genotype.pdf", width = 15, height = 30)
DimPlot(object = seuObject_integrated, reduction = "umap", ncol = 3, label = FALSE, pt.size = 0.5, split.by = "Genotype")
dev.off()
save(seuObject_integrated, file = "output_Tsvz/04_seuObject_Res0.5.RData")