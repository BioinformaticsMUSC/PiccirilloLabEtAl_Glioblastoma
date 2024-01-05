setwd("/users/SuganyaSubramanian/Sara_Final")

#nsvz
load("output_Nsvz/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")

#tm
load("output_TumorMass/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")

#tsvz
load("output_Tsvz/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res02.RData")

## All module genes combined ####
modules <- read.csv('All_Cell_State_Markers.csv')

mod_list = list()
for (col in colnames(modules)) {
  mod_list[col] <- c(modules[col])
}

p1 <- SCpubr::do_EnrichmentHeatmap(sample = seuOject_Nsvz_nodoub_withHarmony,
                                   input_gene_list = mod_list,
                                   group.by = "seurat_clusters",
                                   viridis_direction = -1,
                                   use_viridis = FALSE,
                                   enforce_symmetry = TRUE)

pdf("output_Nsvz/15_EnrichHeatMap_clusters.pdf", width = 15, height = 10)  
p1
dev.off()

p2 <- SCpubr::do_EnrichmentHeatmap(sample = seuOject_Tsvz_nodoub_withHarmony,
                                   input_gene_list = mod_list,
                                   group.by = "seurat_clusters",
                                   viridis_direction = -1,
                                   use_viridis = FALSE,
                                   enforce_symmetry = TRUE)

pdf("output_Tsvz/15_EnrichHeatMap_clusters.pdf", width = 15, height = 10)  
p2
dev.off()

p3 <- SCpubr::do_EnrichmentHeatmap(sample = seuOject_TM_nodoub_withHarmony,
                                   input_gene_list = mod_list,
                                   group.by = "seurat_clusters",
                                   viridis_direction = -1,
                                   use_viridis = FALSE,
                                   enforce_symmetry = TRUE)

pdf("output_TumorMass/15_EnrichHeatMap_clusters.pdf", width = 15, height = 10)  
p3
dev.off()