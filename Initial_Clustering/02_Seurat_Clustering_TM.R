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

G10A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_01GBM10A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G12A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_04GBM12A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G14A <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_GBM14A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G16A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_07GBM16A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G17A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_10GBM17A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G18A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_13GBM18A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G20A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_16GBM20A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G22A <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_GBM22A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G7A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_19GBM7A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G8A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_21GBM8A_filtered.h5", use.names = TRUE, unique.features = TRUE)
G9A <- Read10X_h5("/zfs/musc3/Sara/set2_cellbender_out/Set2_24GBM9A_filtered.h5", use.names = TRUE, unique.features = TRUE)
MIS1A <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_MIS1A_filtered.h5", use.names = TRUE, unique.features = TRUE)
MIS2A <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_MIS2A_filtered.h5", use.names = TRUE, unique.features = TRUE)
MIS3A <- Read10X_h5("/zfs/musc3/Sara/set3_cellbender_out/Set3_MIS3A_filtered.h5", use.names = TRUE, unique.features = TRUE)

load("HgProteinCodingGenesSara.rda")

G10A <- G10A[rownames(G10A)%in%ptn_genes,]
G12A <- G12A[rownames(G12A)%in%ptn_genes,]
G14A <- G14A[rownames(G14A)%in%ptn_genes,]
G17A <- G17A[rownames(G17A)%in%ptn_genes,]
G18A <- G18A[rownames(G18A)%in%ptn_genes,]
G20A <- G20A[rownames(G20A)%in%ptn_genes,]
G22A <- G22A[rownames(G22A)%in%ptn_genes,]
G7A <- G7A[rownames(G7A)%in%ptn_genes,]
G8A <- G8A[rownames(G8A)%in%ptn_genes,]
G9A <- G9A[rownames(G9A)%in%ptn_genes,]
MIS1A <- MIS1A[rownames(MIS1A)%in%ptn_genes,]
MIS2A <- MIS2A[rownames(MIS2A)%in%ptn_genes,]
MIS3A <- MIS3A[rownames(MIS3A)%in%ptn_genes,]

G10A_obj <- CreateSeuratObject(counts = G10A,min.features = 100)
G10A_obj@meta.data$Region <- "A"
G10A_obj@meta.data$PatientID <- "GBM10"
G10A_obj@meta.data$Sex <- "Male"
G10A_obj@meta.data$Age <- "64"
G10A_obj@meta.data$Set <- "Set2"
G10A_obj@meta.data$Label <- "Tumor_Mass"
G10A_obj@meta.data$Abbreviation <- "T"

G12A_obj <- CreateSeuratObject(counts = G12A,min.features = 100)
G12A_obj@meta.data$Region <- "A"
G12A_obj@meta.data$PatientID <- "GBM12"
G12A_obj@meta.data$Sex <- "Female"
G12A_obj@meta.data$Age <- "66"
G12A_obj@meta.data$Set <- "Set2"
G12A_obj@meta.data$Label <- "Tumor_Mass"
G12A_obj@meta.data$Abbreviation <- "T"

G14A_obj <- CreateSeuratObject(counts = G14A,min.features = 100)
G14A_obj@meta.data$Region <- "A"
G14A_obj@meta.data$PatientID <- "GBM14"
G14A_obj@meta.data$Sex <- "Male"
G14A_obj@meta.data$Age <- "61"
G14A_obj@meta.data$Set <- "Set3"
G14A_obj@meta.data$Label <- "Tumor_Mass"
G14A_obj@meta.data$Abbreviation <- "T"

G16A_obj <- CreateSeuratObject(counts = G16A,min.features = 100)
G16A_obj@meta.data$Region <- "A"
G16A_obj@meta.data$PatientID <- "GBM16"
G16A_obj@meta.data$Sex <- "Male"
G16A_obj@meta.data$Age <- "77"
G16A_obj@meta.data$Set <- "Set2"
G16A_obj@meta.data$Label <- "Tumor_Mass"
G16A_obj@meta.data$Abbreviation <- "T"

G17A_obj <- CreateSeuratObject(counts = G17A,min.features = 100)
G17A_obj@meta.data$Region <- "A"
G17A_obj@meta.data$PatientID <- "GBM17"
G17A_obj@meta.data$Sex <- "Female"
G17A_obj@meta.data$Age <- "56"
G17A_obj@meta.data$Set <- "Set2"
G17A_obj@meta.data$Label <- "Tumor_Mass"
G17A_obj@meta.data$Abbreviation <- "T"

G18A_obj <- CreateSeuratObject(counts = G18A,min.features = 100)
G18A_obj@meta.data$Region <- "A"
G18A_obj@meta.data$PatientID <- "GBM18"
G18A_obj@meta.data$Sex <- "Female"
G18A_obj@meta.data$Age <- "78"
G18A_obj@meta.data$Set <- "Set2"
G18A_obj@meta.data$Label <- "Tumor_Mass"
G18A_obj@meta.data$Abbreviation <- "T"

G20A_obj <- CreateSeuratObject(counts = G20A,min.features = 100)
G20A_obj@meta.data$Region <- "A"
G20A_obj@meta.data$PatientID <- "GBM20"
G20A_obj@meta.data$Sex <- "Male"
G20A_obj@meta.data$Age <- "71"
G20A_obj@meta.data$Set <- "Set2"
G20A_obj@meta.data$Label <- "Tumor_Mass"
G20A_obj@meta.data$Abbreviation <- "T"

G22A_obj <- CreateSeuratObject(counts = G22A,min.features = 100)
G22A_obj@meta.data$Region <- "A"
G22A_obj@meta.data$PatientID <- "GBM22"
G22A_obj@meta.data$Sex <- "Female"
G22A_obj@meta.data$Age <- "80"
G22A_obj@meta.data$Set <- "Set3"
G22A_obj@meta.data$Label <- "Tumor_Mass"
G22A_obj@meta.data$Abbreviation <- "T"

G7A_obj <- CreateSeuratObject(counts = G7A,min.features = 100)
G7A_obj@meta.data$Region <- "A"
G7A_obj@meta.data$PatientID <- "GBM7"
G7A_obj@meta.data$Sex <- "Male"
G7A_obj@meta.data$Age <- "77"
G7A_obj@meta.data$Set <- "Set2"
G7A_obj@meta.data$Label <- "Tumor_Mass"
G7A_obj@meta.data$Abbreviation <- "T"

G8A_obj <- CreateSeuratObject(counts = G8A,min.features = 100)
G8A_obj@meta.data$Region <- "A"
G8A_obj@meta.data$PatientID <- "GBM8"
G8A_obj@meta.data$Sex <- "Male"
G8A_obj@meta.data$Age <- "64"
G8A_obj@meta.data$Set <- "Set2"
G8A_obj@meta.data$Label <- "Tumor_Mass"
G8A_obj@meta.data$Abbreviation <- "T"

G9A_obj <- CreateSeuratObject(counts = G9A,min.features = 100)
G9A_obj@meta.data$Region <- "A"
G9A_obj@meta.data$PatientID <- "GBM9"
G9A_obj@meta.data$Sex <- "Female"
G9A_obj@meta.data$Age <- "71"
G9A_obj@meta.data$Set <- "Set2"
G9A_obj@meta.data$Label <- "Tumor_Mass"
G9A_obj@meta.data$Abbreviation <- "T"

MIS1A_obj <- CreateSeuratObject(counts = MIS1A,min.features = 100)
MIS1A_obj@meta.data$Region <- "A"
MIS1A_obj@meta.data$PatientID <- "MIS1"
MIS1A_obj@meta.data$Sex <- "Male"
MIS1A_obj@meta.data$Age <- "21"
MIS1A_obj@meta.data$Set <- "Set3"
MIS1A_obj@meta.data$Label <- "Tumor_Mass"
MIS1A_obj@meta.data$Abbreviation <- "T"

MIS2A_obj <- CreateSeuratObject(counts = MIS2A,min.features = 100)
MIS2A_obj@meta.data$Region <- "A"
MIS2A_obj@meta.data$PatientID <- "MIS2"
MIS2A_obj@meta.data$Sex <- "Female"
MIS2A_obj@meta.data$Age <- "73"
MIS2A_obj@meta.data$Set <- "Set3"
MIS2A_obj@meta.data$Label <- "Tumor_Mass"
MIS2A_obj@meta.data$Abbreviation <- "T"

MIS3A_obj <- CreateSeuratObject(counts = MIS3A,min.features = 100)
MIS3A_obj@meta.data$Region <- "A"
MIS3A_obj@meta.data$PatientID <- "MIS3"
MIS3A_obj@meta.data$Sex <- "Female"
MIS3A_obj@meta.data$Age <- "51"
MIS3A_obj@meta.data$Set <- "Set3"
MIS3A_obj@meta.data$Label <- "Tumor_Mass"
MIS3A_obj@meta.data$Abbreviation <- "T"

##### Merge all ########
seuObject <- merge(G10A_obj,
                   y = c(G12A_obj,G14A_obj,G16A_obj,G17A_obj,G18A_obj,G20A_obj,
                         G22A_obj,G7A_obj,G8A_obj,G9A_obj,MIS1A_obj,MIS2A_obj,MIS3A_obj),
                   add.cell.ids = c("GBM10A","GBM12A","GBM14A","GBM16A","GBM17A","GBM18A",
                                    "GBM20A","GBM22A","GBM7A","GBM8A","GBM9A","MIS1A",
                                    "MIS2A","MIS3A"),
                   project = "Sara_final_TM")

########## Mito column ##########
seuObject[["pMito"]] <- PercentageFeatureSet(seuObject, pattern = "^MT-")

seuObject[["pRibo"]] <- PercentageFeatureSet(seuObject, pattern = "^RP[SL]")

seuObject@meta.data <- seuObject@meta.data %>%
  rownames_to_column("Cell") %>%
  mutate(Genotype = sapply(X = strsplit(colnames(seuObject), split = "_"), FUN = "[", 1)) %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
  column_to_rownames("Cell")

dir.create("output_TumorMass")

save(seuObject,file="output_TumorMass/01_SeuObj_TumorMass_Unfilt.RData")

#QC plot 1
pdf("output_TumorMass/01_Quality_Control_plots_Genotype.pdf", width=25,height=5)
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

save(seuObject_filt,file="output_TumorMass/02_SeuObj_TumorMass_MtFiltered_u25000_m5.RData")

df <- seuObject_filt@meta.data %>% as.data.frame()
tmp <- table(df$Genotype) %>% as.data.frame()

pdf("output_TumorMass/04-GeneExp_barplot_Genotype.pdf",width = 25,height = 7)
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

seuObject_split <- seuObject_split[c("GBM10A","GBM12A","GBM14A","GBM16A","GBM17A",
                                     "GBM18A","GBM20A","GBM22A","GBM7A",
                                     "GBM8A","GBM9A","MIS1A","MIS2A","MIS3A")]

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

save(seuObject_integrated, file = "output_TumorMass/04_Integrated_allRes.RData")

pdf("output_TumorMass/05_Clustree-Data_Integrated.pdf", width = 15, height = 6)
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

pdf("output_TumorMass/06_UMAP-Data_Integrated_Genotype.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_TumorMass/07_UMAP-Data_Integrated_splitby_Genotype.pdf", width = 15, height = 30)
DimPlot(object = seuObject_integrated, reduction = "umap", ncol = 3, label = FALSE, pt.size = 0.5, split.by = "Genotype")
dev.off()

save(seuObject_integrated, file = "output_TumorMass/04_seuObject_Res0.5.RData")