#paired.r
#make correlation data

#cell2location cells
library(psych)
library(ggplot2)
library(patchwork)
library(stringr)
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(MoMAColors)
library(ggrepel)
setwd("/Users/bryanwgranger/biocm/projects/sara_gbm/spatial/cell2location")

nsvz_cells = c("Astrocytes", "Endothelial", "Ependymal", "Microglia", "NPC", "Neurons",      
               "OPC", "Oligodendrocytes")

tm_cells = c("CancerCell", "GBMac", "GBMmes", "GBMnpc", "GBMopc",  "Endothelial",      
             "MDM", "Microglia", "Neurons", "Oligodendrocytes")

tsvz_cells = c("CancerCell",  "GBMmes", "GBMnpc", "GBMopc", "Endothelial",
               "MDM", "Microglia", "NPC", "Oligodendrocytes")

all_cells = unique(append(append(nsvz_cells, tm_cells), tsvz_cells))
all_cells2 <- unique(append(append(tm_cells, tsvz_cells), nsvz_cells))
cell_order = c("Astrocytes",
               "Endothelial",
               "Ependymal", 
               "Neurons",
               "NPC",
               "Oligodendrocytes", 
               "OPC",
               "CancerCell", 
               "GBMac", 
               "GBMmes", 
               "GBMnpc", 
               "GBMopc",
               "MDM")

c2l_files <- list.files(path = ".", pattern = "c2l_cells.csv", recursive = T, full.names = T)

big_df = data.frame()
sample_size = list()
for (f in c2l_files) {
  data = read.csv(f, row.names=1)
  sample_list = unlist(str_split(f, pattern = "/"))[length(unlist(str_split(f, pattern = "/")))] |>
    str_split(pattern = "_") |>
    unlist() 
  sample = sample_list[1]
  print(sample)
  sample_size[[sample]] = nrow(data)
  
  if (str_detect(sample, "A")) {
    cell_list = tm_cells
  } else if (str_detect(sample, "B")) {
    cell_list = tsvz_cells
  } else {
    cell_list = nsvz_cells
  }
  
  t = cor(data[, cell_list]) |> as.data.frame() |>
    select(Microglia)
  colnames(t) <- sample
  t <- t(t) |> as.data.frame()
  big_df = dplyr::bind_rows(big_df, t)
}

big_df_melt <- big_df |>
  tibble::rownames_to_column("sample") |>
  dplyr::select(-Microglia) |>
  melt()

big_df_melt$variable <- factor(big_df_melt$variable, levels = cell_order)
big_df_melt$sample <- factor(big_df_melt$sample, levels = c("GBM4A", "GBM7A", "GBM8A", "GBM9A",
                                                            "GBM4B", "GBM7B", "GBM8B", "GBM9B",
                                                            "NSVZ"))

big_df_melt <- big_df_melt |>
  dplyr::mutate(patient = stringr::str_sub(sample, start = 1, end = -2)) |>
  dplyr::mutate(region = dplyr::case_when(stringr::str_ends(sample, "A") ~ "TM",
                                          stringr::str_ends(sample, "B") ~ "TSVZ",
                                          .default = "NSVZ"))
common_cells = c("CancerCell",  "GBMmes", "GBMnpc", "GBMopc",  "Endothelial",      
                 "MDM",  "Oligodendrocytes")
t_samples <- c("GBM4", "GBM7", "GBM8", "GBM9")


res_list4 <- list()
for (c in common_cells){
  for (pat in t_samples) {
    id = str_glue("{pat}_{c}")
    xz <- big_df_melt |>
      dplyr::filter(variable == c, region == "TM", patient == pat)
    
    xy <- big_df_melt |>
      dplyr::filter(variable == c, region == "TSVZ", patient == pat)
    print(xy)
    print(xz)
    print(sample_size[[str_glue("{pat}A")]])
    print(sample_size[[str_glue("{pat}B")]])
    res_list4[[id]] <- r.test(r12 = xy$value,
                                r34 = xz$value,
                                n=sample_size[[str_glue("{pat}B")]],
                                n2=sample_size[[str_glue("{pat}A")]],
                                twotailed = TRUE)
  }
}

res_df4 <- data.frame()
for (n in names(res_list4)){
  results <- res_list4[[n]]
  tdf <- data.frame(p=results$p, z=results$z, cell=n)
  print(tdf)
  res_df4 <- rbind(res_df4, tdf)
  rownames(res_df4) <- NULL
}
write.csv(res_df4, "paired_R_test.csv")

res_df4 <- res_df4 |>
  tidyr::separate_wider_delim(cell, delim="_", names = c("patient", "cell_class"))

cor_data <- big_df_melt |>
  dplyr::filter(variable %in% common_cells, sample != "NSVZ")
  
cor_data <- cor_data |>
  tidyr::pivot_wider(id_cols = c("patient", "variable"),
                     names_from = "region",
                     values_from = "value") |>
  dplyr::mutate(diff = TSVZ - TM) |>
  dplyr::mutate(patient_cell = str_c(patient, variable, sep = "_"))

res_df4$pt <- res_df4$p + 1e-30
res_df4 <- res_df4 |>
  dplyr::mutate(patient_cell = str_c(patient, cell_class, sep = "_"))

res_df4 <- res_df4 |>
  dplyr::left_join(cor_data, by = "patient_cell")

min_val = res_df4 |> dplyr::filter(pt < 0.05) |>
  dplyr::select(pt) |>
  log10() |>
  max()
min_val = -1 * min_val

res_df4$cell_class <- factor(res_df4$cell_class, levels = c("CancerCell", "GBMmes", "GBMnpc", "GBMopc",  "Endothelial",      
                                                            "MDM",   "Oligodendrocytes"))
res_df4 |>
  ggplot(aes(
    y=forcats::fct_rev(cell_class),
    x=patient.x,
    fill=diff
  )) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid = "white", high="red", name = "Difference in correlation,\nTSVZ - TM") +
  geom_point(data = res_df4 |> dplyr::filter(pt < 0.05), aes(size=-log10(pt))) +
  geom_hline(yintercept=3.5, color="black") +
  #geom_text(data = res_df4 |> dplyr::filter(pt > 0.05), label = "n.s.") +
  scale_size_continuous(name = "-log10(p.value)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  #coord_flip() +
  ylab("Patient") +
  xlab("Cell class")
ggsave("Paired_R_updated5.pdf")

res_df4 |>
  ggplot(aes(
    x=diff,
    y=-log10(pt),
    label=stringr::str_replace(patient_cell, "_", " ")
  )) +
  geom_point(aes(color=cell_class), size=3) +
  geom_text_repel(size=3, min.segment.length = 0.2, force = 20) +
  geom_hline(yintercept=1.3, color='gray', linetype="dashed") +
  theme_bw() +
  ylab("-nlog10(p.value)") +
  xlab("Difference in Pearson's correlation, T_SVZ - T_Mass") +
  labs(color="Cell Class") +
  xlim(-0.75, 0.75)
ggsave("Paired_R_scatter.pdf", height = 6)
