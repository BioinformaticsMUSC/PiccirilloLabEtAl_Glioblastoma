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

#nsvz
load("output_Nsvz/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")

#tm
load("output_TumorMass/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")

seuOject_Nsvz_nodoub_withHarmony@meta.data <- seuOject_Nsvz_nodoub_withHarmony@meta.data %>%
  unite("ForAnno", c(seurat_clusters,Abbreviation),remove = F)

seuOject_TM_nodoub_withHarmony@meta.data <- seuOject_TM_nodoub_withHarmony@meta.data %>%
  unite("ForAnno", c(seurat_clusters,Abbreviation),remove = F)

#################################################
#####  Normal svz vs Tumor Mass  subset 1 #######
#################################################

sub1 <- subset(seuOject_TM_nodoub_withHarmony, subset = (seurat_clusters %in% c("0","1","2","3")))

merge_Nsvz_sub1TM <- merge(seuOject_Nsvz_nodoub_withHarmony,
                           y = c(sub1))

annoTM <- merge_Nsvz_sub1TM@meta.data %>% dplyr::select(c("ForAnno"))

write.table(annoTM, file = "New_inferCNV/Annotation_NSVZ_vs_sub1TM.txt", sep="\t",quote=F,
            row.names = TRUE, col.names = FALSE)

count_mxtTM <- as.matrix(merge_Nsvz_sub1TM@assays$RNA@counts)

infercnv_obj1TM = CreateInfercnvObject(raw_counts_matrix=count_mxtTM,
                                       annotations_file="New_inferCNV/Annotation_NSVZ_vs_sub1TM.txt",
                                       delim="\t",
                                       gene_order_file="gencode_v19_gene_pos.txt",
                                       ref_group_names=c("0_N-SVZ","1_N-SVZ","10_N-SVZ","11_N-SVZ","12_N-SVZ",
                                                         "2_N-SVZ","3_N-SVZ","4_N-SVZ","5_N-SVZ",
                                                         "6_N-SVZ","7_N-SVZ","8_N-SVZ","9_N-SVZ"))

infercnv_obj1TM = infercnv::run(infercnv_obj1TM,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir="New_inferCNV/Nsvz_sub1TM", 
                                cluster_by_groups=TRUE, 
                                plot_steps=F,
                                denoise=TRUE,
                                resume_mode = FALSE,
                                HMM=TRUE)

#################################################
#####  Normal svz vs Tumor Mass  subset 2 #######
#################################################
sub2 <- subset(seuOject_TM_nodoub_withHarmony, subset = (seurat_clusters %in% c("4","5","6","7","8","9")))

merge_Nsvz_sub2TM <- merge(seuOject_Nsvz_nodoub_withHarmony,
                           y = c(sub2))

annoTM <- merge_Nsvz_sub2TM@meta.data %>% dplyr::select(c("ForAnno"))

write.table(annoTM, file = "New_inferCNV/Annotation_NSVZ_vs_sub2TM.txt", sep="\t",quote=F,
            row.names = TRUE, col.names = FALSE)

count_mxtTM <- as.matrix(merge_Nsvz_sub2TM@assays$RNA@counts)

infercnv_obj2TM = CreateInfercnvObject(raw_counts_matrix=count_mxtTM,
                                       annotations_file="New_inferCNV/Annotation_NSVZ_vs_sub2TM.txt",
                                       delim="\t",
                                       gene_order_file="gencode_v19_gene_pos.txt",
                                       ref_group_names=c("0_N-SVZ","1_N-SVZ","10_N-SVZ","11_N-SVZ","12_N-SVZ",
                                                         "2_N-SVZ","3_N-SVZ","4_N-SVZ","5_N-SVZ",
                                                         "6_N-SVZ","7_N-SVZ","8_N-SVZ","9_N-SVZ"))

infercnv_obj2TM = infercnv::run(infercnv_obj2TM,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir="New_inferCNV/Nsvz_sub2TM", 
                                cluster_by_groups=TRUE, 
                                plot_steps=F,
                                denoise=TRUE,
                                resume_mode = FALSE,
                                HMM=TRUE)

#################################################
#####  Normal svz vs Tumor Mass  subset 3 #######
#################################################
sub3 <- subset(seuOject_TM_nodoub_withHarmony, 
  subset = (seurat_clusters %in% c("10","11","12","13","14","15","16","17","18","19","20","21")))

merge_Nsvz_sub3TM <- merge(seuOject_Nsvz_nodoub_withHarmony,
                           y = c(sub3))

annoTM <- merge_Nsvz_sub3TM@meta.data %>% dplyr::select(c("ForAnno"))

write.table(annoTM, file = "New_inferCNV/Annotation_NSVZ_vs_sub3TM.txt", sep="\t",quote=F,
            row.names = TRUE, col.names = FALSE)

count_mxtTM <- as.matrix(merge_Nsvz_sub3TM@assays$RNA@counts)

infercnv_obj3TM = CreateInfercnvObject(raw_counts_matrix=count_mxtTM,
                                       annotations_file="New_inferCNV/Annotation_NSVZ_vs_sub3TM.txt",
                                       delim="\t",
                                       gene_order_file="gencode_v19_gene_pos.txt",
                                       ref_group_names=c("0_N-SVZ","1_N-SVZ","10_N-SVZ",
                                                         "2_N-SVZ","3_N-SVZ","4_N-SVZ","5_N-SVZ",
                                                         "6_N-SVZ","7_N-SVZ","8_N-SVZ","9_N-SVZ"))

infercnv_objTM = infercnv::run(infercnv_obj3TM,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="New_inferCNV/Nsvz_sub3TM", 
                               cluster_by_groups=TRUE, 
                               plot_steps=F,
                               denoise=TRUE,
                               resume_mode = FALSE,
                               HMM=TRUE)