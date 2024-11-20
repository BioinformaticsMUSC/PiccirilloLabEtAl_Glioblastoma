#cellrank - biocm_python 1.0.6
#https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/kernels/200_rna_velocity.html

import numpy as np
import pandas as pd
import cellrank as cr
import scanpy as sc
import scvelo as scv
import os
import matplotlib.pyplot as plt

os.chdir("/zfs/musc3/Sara/cellrank")

region = "TSVZ_tumor"
#adata = sc.read("/zfs/musc3/Sara/seurat_objects/newest/full/tm_adata.h5ad")
#adata = sc.read("/zfs/musc3/Sara/seurat_objects/newest/full/tsvz_adata.h5ad")
adata = sc.read("/zfs/musc3/Sara/Set2_velo/individual_velos/tsvz/tsvz_merged_all_genes.h5ad")
#adata = sc.read("/zfs/musc3/Sara/Set2_velo/individual_velos/tm/tm_merged_all_genes.h5ad")

#adata = adata[adata.obs['Tumor_Class'] == "normal"]

sc.pp.normalize_total(adata)

# scv.pp.filter_and_normalize(
#     adata, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False
# )

scv.pp.filter_and_normalize(
    adata, min_shared_counts=20, n_top_genes=5000, subset_highly_variable=True,
    retain_genes = ["ZEB1", "RBPJ", "NRF1", "SOX5"]
    #retain_genes = ["SOX2", "TCF4", "NFIA", "TCF7L2"]
)

sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

scv.tl.recover_dynamics(adata, n_jobs=24) 

scv.tl.velocity(adata, mode="dynamical")
scv.tl.latent_time(adata)

sc.external.pp.magic(adata)

## CELLRANK

vk = cr.kernels.VelocityKernel(adata)

vk.compute_transition_matrix()

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()

combined_kernel = 0.8 * vk + 0.2 * ck
print(combined_kernel)

vk.write_to_adata()


adata.write("/zfs/musc3/Sara/Set2_velo/individual_velos/tm/tm_FULL_merged_5k_hvg_cellrank.h5ad")
#adata.write("/zfs/musc3/Sara/Set2_velo/individual_velos/tsvz/tsvz_FULL_merged_5k_hvg_cellrank.h5ad")

###advanced
g2 = cr.estimators.GPCCA(vk)
print(g2)
g2.compute_schur()
g2.plot_spectrum(real_only=True)
plt.tight_layout()
plt.savefig(f"{region}_schur_decomp.pdf")

N_STATES = 6

g2.compute_macrostates(n_states=N_STATES, cluster_key="Cell_Class")
g2.plot_macrostates(which="all", legend_loc="right", s=100)
plt.tight_layout()
plt.savefig(f"{region}_macrostates_all_n{N_STATES}.pdf")

g2.predict_terminal_states()
g2.plot_macrostates(which="terminal", legend_loc="right", s=100)
plt.tight_layout()
plt.savefig(f"{region}_terminal_n{N_STATES}.pdf")

g2.predict_initial_states(allow_overlap=True)
g2.plot_macrostates(which="initial", s=100)
plt.tight_layout()
plt.savefig(f"{region}_initial_n{N_STATES}.pdf")


g2.compute_fate_probabilities(use_petsc=True, preconditioner="ilu") #tol=1e-7 for TM
g2.plot_fate_probabilities(same_plot=False, ncols=3)
plt.tight_layout()
plt.savefig(f"{region}_fate_probs.pdf")

g2.plot_fate_probabilities(same_plot=True, legend_loc="right margin")
plt.tight_layout()
plt.savefig(f"{region}_fate_probs_all.pdf")

cr.pl.circular_projection(adata, keys=["Cell_Class"], legend_loc="right", ncols=1)
plt.tight_layout()
plt.savefig(f"{region}_circle_proj_newest.pdf")

#### DRIVER GENES ######

######### HEATMAP
model = cr.models.GAM(adata, n_knots=6)

### all macrostates
for m in g2.terminal_states.cat.categories:
	drivers = g2.compute_lineage_drivers(lineages=m)
	cr.pl.heatmap(adata, model=model, lineages=m, cluster_key="Cell_Class", show_fate_probabilities=True, data_key="X", genes=drivers.head(40).index, time_key="latent_time", figsize=(12, 10), show_all_genes=True, weight_threshold=(1e-3, 1e-3), save=f"{region}_heatmap_{m}.pdf")


for m in g2.terminal_states.cat.categories:
    drivers = g2.compute_lineage_drivers(lineages=m)
    adata.obs[f"fate_probabilities_{m}"] = g2.fate_probabilities[m].X.flatten()
    sc.pl.embedding(
        adata,
        basis="umap",
        color=[f"fate_probabilities_{m}"] + list(drivers.index[:8]),
        color_map="viridis",
        s=20,
        ncols=3,
        vmax="p96",
        save=f"{m}_drivers.pdf"
    )

drivers = g2.compute_lineage_drivers(lineages="GBMmes_1")
drivers2 = g2.compute_lineage_drivers(lineages="GBMmes_2")

genes_oi = {"GBMmes_1": list(drivers.index[:10]),
            "GBMmes_2": list(drivers2.index[:10])}

adata.var["mean expression"] = adata.X.mean(axis=0)  

comb = [c for c in combinations(g2.terminal_states.cat.categories, r=2)]
for c in comb:
    print(c[0], "vs", c[1])
    drivers = g2.compute_lineage_drivers(lineages=c[0])
    drivers2 = g2.compute_lineage_drivers(lineages=c[1])
    genes_oi = {c[0]: list(drivers.index[:10]),
                c[1]: list(drivers2.index[:10])}
    both_drivers = g2.compute_lineage_drivers(lineages=[c[0], c[1]])
    g2.plot_lineage_drivers_correlation(
        lineage_x=c[0],
        lineage_y=c[1],
        adjust_text=True,
        gene_sets=genes_oi,
        color="mean expression",
        legend_loc="none",
        figsize=(5, 5),
        dpi=150,
        fontsize=9,
        size=50,
    )
    plt.savefig(f"CORR_{c[0]}_{c[1]}.pdf")

adata.write("/zfs/musc3/Sara/Set2_velo/individual_velos/tsvz/tsvz_merged_5k_hvg_cellrank.h5ad")
#######
for m in g2.terminal_states.cat.categories:
    g.plot_lineage_drivers_correlation(
        lineage_x="GBMac",
        lineage_y="Alpha",
        adjust_text=True,
        gene_sets=genes_oi,
        color="mean expression",
        legend_loc="none",
        figsize=(5, 5),
        dpi=150,
        fontsize=9,
        size=50,
    )

all_drivers = pd.DataFrame()
for m in g2.terminal_states.cat.categories:
    drivers = g2.compute_lineage_drivers(lineages=m)
    all_drivers = pd.concat([all_drivers, drivers], axis=1)


cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["SOX2", "TCF4", "NFIA", "TCF7L2"],
    same_plot=True,
    ncols=2,
    time_key="latent_time",
    hide_cells=True,
    weight_threshold=(1e-3, 1e-3),
)
plt.tight_layout()
plt.savefig("tm_gene_trends.pdf")

##### scenic regulons
 df = pd.read_csv("/zfs/musc3/Sara/scenic/new/normal/heatmap_data_4regions.csv", index_col = [0])
 tm_reg = df.sort_values(by="TM_tumor", ascending=False).index[:40]
 tm_reg_genes = [t.replace("_(+)", "") for t in tm_reg]

df = pd.read_csv("/zfs/musc3/Sara/scenic/new/normal/heatmap_data_4regions.csv", index_col = [0])
ts_reg = df.sort_values(by="TSVZ_tumor", ascending=False).index[:40]
ts_reg_genes = [t.replace("_(+)", "") for t in ts_reg]

same_genes = np.intersect1d(ts_reg_genes, adata.var_names)
 ### all macrostates
for m in g2.terminal_states.cat.categories:
    drivers = g2.compute_lineage_drivers(lineages=m)
    cr.pl.heatmap(adata, model=model, lineages=m, cluster_key="Cell_Class", show_fate_probabilities=True, data_key="X", genes=same_genes, time_key="latent_time", figsize=(12, 10), show_all_genes=True, weight_threshold=(1e-3, 1e-3), save=f"scenic_regs_{region}_heatmap_{m}.pdf")


