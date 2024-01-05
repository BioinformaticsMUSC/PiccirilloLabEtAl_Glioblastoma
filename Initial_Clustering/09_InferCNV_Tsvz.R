suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(sctransform)
  library(hdf5r) #install.packages("hdf5r")
  library(ggrastr)
  library(clustree)
  library(cowplot)
  #library(copykat)
  #library(SCEVAN)
  library(infercnv)
})


setwd("/zfs/musc3/Sara/01Final")
dir.create("New_inferCNV")

#nsvz
load("output_Nsvz/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")

#tsvz
load("output_Tsvz/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res02.RData")

seuOject_Nsvz_nodoub_withHarmony@meta.data <- seuOject_Nsvz_nodoub_withHarmony@meta.data %>%
  unite("ForAnno", c(seurat_clusters,Abbreviation),remove = F)

seuOject_Tsvz_nodoub_withHarmony@meta.data <- seuOject_Tsvz_nodoub_withHarmony@meta.data %>%
  unite("ForAnno", c(seurat_clusters,Abbreviation),remove = F)

#################################################
#####  Normal svz vs tumor svz  subset 1 ########
#################################################

### Obejct was large to run so subsetting by few clusters

sub1 <- subset(seuOject_Tsvz_nodoub_withHarmony, subset = (seurat_clusters %in% c("0","1","2","3","4","5")))

merge_Nsvz_sub1Tsvz <- merge(seuOject_Nsvz_nodoub_withHarmony,
                           y = c(sub1))

annoTsvz <- merge_Nsvz_sub1Tsvz@meta.data %>% dplyr::select(c("ForAnno"))
table(merge_Nsvz_sub1Tsvz$ForAnno)

write.table(annoTsvz, file = "New_inferCNV/Annotation_NSVZ_vs_TSVZsub1.txt", sep="\t",quote=F,
            row.names = TRUE, col.names = FALSE)

count_mxt <- as.matrix(merge_Nsvz_sub1Tsvz@assays$RNA@counts)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_mxt,
                                    annotations_file="New_inferCNV/Annotation_NSVZ_vs_TSVZsub1.txt",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c("0_N-SVZ","1_N-SVZ","10_N-SVZ","11_N-SVZ","12_N-SVZ",
                                                      "2_N-SVZ","3_N-SVZ","4_N-SVZ","5_N-SVZ",
                                                      "6_N-SVZ","7_N-SVZ","8_N-SVZ","9_N-SVZ"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="New_inferCNV/NSVZ_vs_TSVZsub1", 
                             cluster_by_groups=TRUE, 
                             plot_steps=F,
                             denoise=TRUE,
                             output_format="pdf",
                             #resume_mode = FALSE,
                             HMM=TRUE)

#################################################
#####  Normal svz vs tumor svz  subset 2 ########
#################################################

sub2 <- subset(seuOject_Tsvz_nodoub_withHarmony, subset = (seurat_clusters %in% c("6","7","8","9","10","11","12")))

merge_Nsvz_sub2Tsvz <- merge(seuOject_Nsvz_nodoub_withHarmony,
                           y = c(sub2))

annoTsvz <- merge_Nsvz_sub2Tsvz@meta.data %>% dplyr::select(c("ForAnno"))

write.table(annoTsvz, file = "New_inferCNV/Annotation_NSVZ_vs_TSVZsub2.txt", sep="\t",quote=F,
            row.names = TRUE, col.names = FALSE)

count_mxt <- as.matrix(merge_Nsvz_sub2Tsvz@assays$RNA@counts)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_mxt,
                                    annotations_file="New_inferCNV/Annotation_NSVZ_vs_TSVZsub2.txt",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c("0_N-SVZ","1_N-SVZ","10_N-SVZ","11_N-SVZ","12_N-SVZ",
                                                      "2_N-SVZ","3_N-SVZ","4_N-SVZ","5_N-SVZ",
                                                      "6_N-SVZ","7_N-SVZ","8_N-SVZ","9_N-SVZ"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="New_inferCNV/NSVZ_vs_TSVZsub2", 
                             cluster_by_groups=TRUE, 
                             plot_steps=F,
                             denoise=TRUE,
                             output_format="pdf",
                             #resume_mode = FALSE,
                             HMM=TRUE)
