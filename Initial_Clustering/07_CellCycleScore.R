library(RCurl)
library(AnnotationHub)
library(ensembldb)

setwd("/users/SuganyaSubramanian/Sara_FinalPtn")
#nsvz
load("output_Nsvz/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")
#tsvz
load("output_Tsvz/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res02.RData")
#tm
load("output_TumorMass/07_SeuratObj_NoDoublet_Reclustered_withHarmonyBySet_Res03.RData")

seuOject_Nsvz_nodoub_withHarmony <- CellCycleScoring(seuOject_Nsvz_nodoub_withHarmony,
                             g2m.features = g2m_genes,
                             s.features = s_genes)

seuOject_Tsvz_nodoub_withHarmony <- CellCycleScoring(seuOject_Tsvz_nodoub_withHarmony,
                                                     g2m.features = g2m_genes,
                                                     s.features = s_genes)

seuOject_TM_nodoub_withHarmony <- CellCycleScoring(seuOject_TM_nodoub_withHarmony,
                                                     g2m.features = g2m_genes,
                                                     s.features = s_genes)

save(seuOject_Nsvz_nodoub_withHarmony, file ="output_Nsvz/08_Nsvz_withHarmonyWithCellCycle.RData")
save(seuOject_Tsvz_nodoub_withHarmony, file ="output_Tsvz/08_Tsvz_withHarmonyWithCellCycle.RData")
save(seuOject_TM_nodoub_withHarmony, file ="output_TumorMass/08_TM_withHarmonyWithCellCycle.RData")

pdf("output_Nsvz/14_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuOject_Nsvz_nodoub_withHarmony,
        reduction = "umap",
        group.by= "Phase")
dev.off()

pdf("output_Tsvz/14_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuOject_Tsvz_nodoub_withHarmony,
        reduction = "umap",
        group.by= "Phase")
dev.off()

pdf("output_TM/14_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuOject_TM_nodoub_withHarmony,
        reduction = "umap",
        group.by= "Phase")
dev.off()








