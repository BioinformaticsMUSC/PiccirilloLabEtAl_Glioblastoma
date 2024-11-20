library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(MoMAColors)

setwd("/Users/bryanwgranger/biocm/projects/sara_gbm/main_analysis/scenic/new/big_heatmap")

tm_rss <- read.csv("RSS_TM_full.csv", row.names = 1) |> t() |> as.data.frame()
rownames(tm_rss) <- str_split_i(rownames(tm_rss), pattern = "_", 1)
tsvz_rss <- read.csv("RSS_TSVZ_full.csv", row.names = 1) |> t() |> as.data.frame()
rownames(tsvz_rss) <- str_split_i(rownames(tsvz_rss), pattern = "_", 1)

tm_auc <- read.csv("TM_full_AUC.csv", row.names = 1)
colnames(tm_auc) <- str_split_i(colnames(tm_auc), pattern = "_", 1)
tsvz_auc <- read.csv("TSVZ_full_AUC.csv", row.names = 1)
colnames(tsvz_auc) <- str_split_i(colnames(tsvz_auc), pattern = "_", 1)

tm_meta <- read.csv("TM_full_meta.csv", row.names = 1) |>
  dplyr::select(Cell_Class)

tsvz_meta <- read.csv("TSVZ_full_meta.csv", row.names = 1) |>
  dplyr::select(Cell_Class)

#####################
get_top_TF <- function (rss_data, n = 10) {
  top_tf <- list()
  for (cell in colnames(rss_data)) {
    top_regs <- rss_data |>
      dplyr::arrange(-!!as.name(cell)) |>
      dplyr::select(!!as.name(cell)) |>
      head(n = n) |> 
      rownames() |>
      as.vector() |>
      unlist()
    top_tf[[cell]] <- top_regs
    
  }
  return (top_tf)
}
tm_top <- get_top_TF(tm_rss, n=5)
tsvz_top <- get_top_TF(tsvz_rss, n=5)

tm_auc_filt <- tm_auc[,colnames(tm_auc) %in% unlist(tm_top)]
tsvz_auc_filt <- tsvz_auc[,colnames(tsvz_auc) %in% unlist(tsvz_top)]
#rm(tm_auc, tsvz_auc)

tm_auc_filt <- cbind(tm_auc_filt, tm_meta)
tsvz_auc_filt <- cbind(tsvz_auc_filt, tsvz_meta)

tm_data <- tm_auc_filt |>
  dplyr::group_by(Cell_Class) |>
  dplyr::summarise_if(is.numeric, mean) |>
  tibble::column_to_rownames("Cell_Class") |>
  t() |>
  as.data.frame()

tm_data_Z <- tm_data |>
  tibble::rownames_to_column("TF") |>
  reshape2::melt()

# |>
#   dplyr::filter(Var1 %in% unlist(tm_top)) 

names(tm_top)
cell_order = c("CancerCell", "GBMac", "GBMmes", "GBMnpc", "GBMopc",
               "Endothelial", "MDM", "Microglia", "Neurons", "Oligodendrocytes")

new_tm_top = list()
for (c in cell_order) {
  new_tm_top[[c]] = tm_top[[c]]
}

new_tm_top

colnames(tm_data_Z) <- c("TF", "Cell_Class", "Value")
tm_data_Z$TF <- factor(tm_data_Z$TF, levels = unique(unlist(new_tm_top)))
tm_data_Z$Cell_Class <- factor(tm_data_Z$Cell_Class, levels = cell_order)

#add RSS data
tm_rss <- tm_rss |>
  tibble::rownames_to_column("TF")


tm_rss_melted <- reshape2::melt(tm_rss)

colnames(tm_rss_melted) <- c("TF", "Cell_Class", "RSS_value")
all_data_tm <- dplyr::left_join(tm_data_Z, tm_rss_melted, by = c("TF", "Cell_Class"))
all_data_tm$TF <- factor(all_data_tm$TF, levels = unique(unlist(new_tm_top)))
all_data_tm$Cell_Class <- factor(all_data_tm$Cell_Class, levels = cell_order)

tumor_celltypes = c("CancerCell", "GBMac", "GBMmes", "GBMnpc", "GBMopc")
norm_celltypes = c("Endothelial", "MDM", "Microglia", "Neurons", "Oligodendrocytes")

max_tm = round(max(tm_data_Z$Value),2)

hm1_tum <- tm_data_Z |>
  dplyr::filter(Cell_Class %in% tumor_celltypes) |>
#all_data_tm |>
  ggplot(mapping = aes(
           x = Cell_Class,
           y = TF,
           fill = Value
         )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color = "grey") +
  #scale_fill_continuous(limits=c(0,0.5)) +
  #viridis::scale_fill_viridis(option = "H", limits=c(0,0.43)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.6, hjust=1),
        legend.position = "none") +
  scale_fill_gradient2(low = "white", high="red", limits = c(0,max_tm)) +
  #scale_fill_gradientn(colors=moma.colors("Fritsch")) +
  # theme(axis.text.x = element_text(angle=90, vjust=0.6, hjust=1),
  #       legend.position = "none") +
  ggtitle("TM Regulons AUC Heatmap")
hm1_tum

ggsave("NEW_tm_hm1_RED_tumor.pdf", height=9, width = 5)

##################
tsvz_data <- tsvz_auc_filt |>
  dplyr::group_by(Cell_Class) |>
  dplyr::summarise_if(is.numeric, mean) |>
  tibble::column_to_rownames("Cell_Class") |>
  t() |>
  as.data.frame()

tsvz_data_Z <- tsvz_data |>
  tibble::rownames_to_column("TF") |>
  reshape2::melt() 


tsvz_order <- c("CancerCell", "GBMmes", "GBMnpc", "GBMopc",
                "Endothelial", "MDM", "Microglia", "NPC", "Oligodendrocytes")
new_tsvz_top <- list()
for (c in tsvz_order) {
  new_tsvz_top[[c]] <- tsvz_top[[c]]
}

colnames(tsvz_data_Z) <- c("TF", "Cell_Class", "Value")
tsvz_data_Z$TF <- factor(tsvz_data_Z$TF, levels = unique(unlist(new_tsvz_top)))
tsvz_data_Z$Cell_Class <- factor(tsvz_data_Z$Cell_Class, levels = tsvz_order)

#add RSS data
tsvz_rss <- tsvz_rss |>
  tibble::rownames_to_column("TF")

tsvz_rss_melted <- reshape2::melt(tsvz_rss)

colnames(tsvz_rss_melted) <- c("TF", "Cell_Class", "RSS_value")
all_data_tsvz <- dplyr::left_join(tsvz_data_Z, tsvz_rss_melted, by = c("TF", "Cell_Class"))
all_data_tsvz$TF <- factor(all_data_tsvz$TF, levels = unique(unlist(new_tsvz_top)))
all_data_tsvz$Cell_Class <- factor(all_data_tsvz$Cell_Class, levels = tsvz_order)

hm2 <- tsvz_data_Z |>
#all_data_tsvz |>
  ggplot(mapping = aes(
    x = Cell_Class,
    y = TF,
    fill = Value
  )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color="grey") +
  #viridis::scale_fill_viridis(option = "A", limits=c(0,0.43)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_gradient2(low = "white", high="red") + 
  theme(axis.text.x = element_text(angle=90, vjust=0.6, hjust=1)) +
  ggtitle("TSVZ Regulon AUC Heatmap")
hm2
ggsave("NEW_tsvz_hm2_no_RED.pdf", height=9, width = 5)

hm1 + hm2
ggsave("NEW_both_hm.pdf", width = 10, height = 9)


#### SCALING
min_max_scale <- function (x) {(x - min(x)) / (max(x) - min(x))}

tm_data_Z$Value_scaled <- min_max_scale(tm_data_Z$Value)
tsvz_data_Z$Value_scaled <- min_max_scale(tsvz_data_Z$Value)

hm1_scaled <- tm_data_Z |>
  #all_data_tm |>
  ggplot(mapping = aes(
    x = Cell_Class,
    y = TF,
    fill = Value_scaled
  )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color = "grey") +
  #viridis::scale_fill_viridis(option = "A") +
  #scale_fill_continuous(limits=c(0,0.5)) +
  viridis::scale_fill_viridis(option = "A", limits=c(0,1)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.6, hjust=1)) +
  ggtitle("TM Regulons AUC Heatmap")

hm2_scaled <- tsvz_data_Z |>
  #all_data_tsvz |>
  ggplot(mapping = aes(
    x = Cell_Class,
    y = TF,
    fill = Value_scaled
  )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color="grey") +
  viridis::scale_fill_viridis(option = "A", limits=c(0,1)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.6, hjust=1)) +
  ggtitle("TSVZ Regulon AUC Heatmap")
#ggsave("NEW_tsvz_hm2_scaled.pdf", height=9, width = 5)

hm1_scaled + hm2_scaled
ggsave("NEW_both_scaled.pdf")


##### SPLIT BETWEEN CELL REGIONS
max_tsvz = round(max(tsvz_data_Z$Value),2)
max_tm = round(max(tm_data_Z$Value),2)

max_metric = max(max_tm, max_tsvz)
### T MASS
tm_tumor_celltypes = c("CancerCell", "GBMac", "GBMmes", "GBMnpc", "GBMopc")
tm_norm_celltypes = c("Endothelial", "MDM", "Microglia", "Neurons", "Oligodendrocytes")

tm_tumor_TFs <- new_tm_top[tm_tumor_celltypes]
tm_normal_TFs <- new_tm_top[tm_norm_celltypes]

hm1_tum <- tm_data_Z |>
  dplyr::filter(Cell_Class %in% tm_tumor_celltypes) |>
  dplyr::filter(TF %in% unlist(tm_tumor_TFs)) |>
  #all_data_tm |>
  ggplot(mapping = aes(
    x = Cell_Class,
    y = TF,
    fill = Value
  )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color = "grey") +
  #scale_fill_continuous(limits=c(0,0.5)) +
  #viridis::scale_fill_viridis(option = "H", limits=c(0,0.43)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.6, hjust=1),
        legend.position = "none") +
  scale_fill_gradient2(low = "white", high="red", limits = c(0,max_metric)) +
  ggtitle("T_Mass Tumor Cells AUC Heatmap")
hm1_tum


hm1_norm <-  tm_data_Z |>
  dplyr::filter(Cell_Class %in% tm_norm_celltypes) |>
  dplyr::filter(TF %in% unlist(tm_normal_TFs)) |>
  #all_data_tm |>
  ggplot(mapping = aes(
    x = Cell_Class,
    y = TF,
    fill = Value
  )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color = "grey") +
  #scale_fill_continuous(limits=c(0,0.5)) +
  #viridis::scale_fill_viridis(option = "H", limits=c(0,0.43)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.6, hjust=1),
        legend.position = "none") +
  scale_fill_gradient2(low = "white", high="red", limits = c(0,max_metric)) +
  #scale_fill_gradientn(colors=moma.colors("Fritsch")) +
  # theme(axis.text.x = element_text(angle=90, vjust=0.6, hjust=1),
  #       legend.position = "none") +
  ggtitle("T_Mass Normal Cells AUC Heatmap")
hm1_norm

hm1_tum + hm1_norm
ggsave("NEW_heatmap_red_TM_SPLIT.pdf", height = 7, width = 9)

### T SVZ

tsvz_tumor_celltypes = c("CancerCell", "GBMmes", "GBMnpc", "GBMopc")
tsvz_norm_celltypes = c("Endothelial", "MDM", "Microglia", "NPC", "Oligodendrocytes")

tsvz_tumor_TFs <- new_tsvz_top[tsvz_tumor_celltypes]
tsvz_normal_TFs <- new_tsvz_top[tsvz_norm_celltypes]



hm1_ts <- tsvz_data_Z |>
  dplyr::filter(Cell_Class %in% tsvz_tumor_celltypes) |>
  dplyr::filter(TF %in% unlist(tsvz_tumor_TFs)) |>
  #all_data_tm |>
  ggplot(mapping = aes(
    x = Cell_Class,
    y = TF,
    fill = Value
  )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color = "grey") +
  #scale_fill_continuous(limits=c(0,0.5)) +
  #viridis::scale_fill_viridis(option = "H", limits=c(0,0.43)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.6, hjust=1)) +
  scale_fill_gradient2(low = "white", high="red", limits = c(0,max_metric)) +
  #scale_fill_gradientn(colors=moma.colors("Fritsch")) +
  # theme(axis.text.x = element_text(angle=90, vjust=0.6, hjust=1),
  #       legend.position = "none") +
  ggtitle("T_SVZ Tumor Cells AUC Heatmap")
hm1_ts


hm1_ts_norm <-  tsvz_data_Z |>
  dplyr::filter(Cell_Class %in% tsvz_norm_celltypes) |>
  dplyr::filter(TF %in% unlist(tsvz_normal_TFs)) |>
  #all_data_tm |>
  ggplot(mapping = aes(
    x = Cell_Class,
    y = TF,
    fill = Value
  )) +
  geom_tile() +
  #geom_point(aes(size=RSS_value), color = "grey") +
  #scale_fill_continuous(limits=c(0,0.5)) +
  #viridis::scale_fill_viridis(option = "H", limits=c(0,0.43)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.6, hjust=1)) +
  scale_fill_gradient2(low = "white", high="red", limits = c(0,max_tsvz)) +
  #scale_fill_gradientn(colors=moma.colors("Fritsch")) +
  # theme(axis.text.x = element_text(angle=90, vjust=0.6, hjust=1),
  #       legend.position = "none") +
  ggtitle("T_SVZ Normal Cells AUC Heatmap")
hm1_ts_norm


###COMBINE
#TUMOR
hm1_tum + hm1_ts
ggsave("NEW_heatmap_red_TUMORCELLS_SPLIT.pdf", height = 7, width = 9)

##NORMAL
hm1_norm + hm1_ts_norm
ggsave("NEW_heatmap_red_NORMALCELLS_SPLIT.pdf", height = 7, width = 9)
