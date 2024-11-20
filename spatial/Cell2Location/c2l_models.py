import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import squidpy as sq
import pandas as pd

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import os


base_dir = "/mnt/Sara/spatial/4_patients/analysis/"

samples = ["GBM7A", "GBM7B", "GBM8A",
			"GBM8B", "GBM9A", "GBM9B"]

adata_tm = sc.read("/mnt/Sara/seurat_objects/newest/full/tm_adata.h5ad")

adata_tsvz = sc.read("/mnt/Sara/seurat_objects/newest/full/tsvz_adata.h5ad")

cell_count_df = pd.read_csv("/mnt/Sara/spatial/4_patients/analysis/gbm_cell_count.csv", index_col=[0])
print(cell_count_df)

# create paths and names to results folders for reference regression and cell2location models
results_folder = './results_v2'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

for s in samples:
	if not os.path.exists(os.path.join(base_dir, s, "cell2location")):
		os.mkdir(os.path.join(base_dir, s, "cell2location"))
	os.chdir(os.path.join(base_dir, s, "cell2location"))

	print("loading adatas...")
	adata_vis = sc.read_visium(
            path = f"/mnt/Sara/spatial/4_patients/run/{s}/outs")
	adata_vis.var_names_make_unique()

	if "A" in s:
		adata_ref = adata_tm.copy()

	else:
		adata_ref = adata_tsvz.copy()

	if s in adata_ref.obs['Genotype'].values:
		adata_ref = adata_ref[adata_ref.obs['Genotype'] == s]

	print("filtering genes...")
	#filter genes to use
	
	selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

	# filter the object
	adata_ref = adata_ref[:, selected].copy()

	cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='Genotype',
                        # cell type, covariate used for constructing signatures
                        labels_key='Cell_Class',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        #categorical_covariate_keys=['Method']
                       )

	
	mod = RegressionModel(adata_ref)

	# view anndata_setup as a sanity check
	mod.view_anndata_setup()

	mod.train(max_epochs=500, use_gpu=False)
	mod.plot_history()
	plt.show()
	plt.savefig("01_mod_history.pdf")
	plt.clf()

	#export the estimated cell abundance (summary of the posterior distribution).
	adata_ref = mod.export_posterior(
		adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
	)
	print("saving reference model...")
	# Save model
	mod.save(f"{ref_run_name}", overwrite=True)
	# Save anndata object with results
	adata_file = f"{ref_run_name}/sc.h5ad"
	adata_ref.write(adata_file)

	if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
		inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()
	else:
		inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()

	inf_aver.columns = adata_ref.uns['mod']['factor_names']
	inf_aver.iloc[0:5, 0:5]

	# find shared genes and subset both anndata and reference signatures
	intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
	adata_vis = adata_vis[:, intersect].copy()
	inf_aver = inf_aver.loc[intersect, :].copy()

	cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

	N_cells = cell_count_df.loc[s, "mean_cell_count"]
	print(f"Training model using {int(N_cells)} cells per location")

	mod = cell2location.models.Cell2location(
	    adata_vis, cell_state_df=inf_aver,
	    # the expected average cell abundance: tissue-dependent
	    # hyper-prior which can be estimated from paired histology:
	    N_cells_per_location=int(N_cells),
	    # hyperparameter controlling normalisation of
	    # within-experiment variation in RNA detection:
	    detection_alpha=20
	)


	mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=False,
         )

	mod.plot_history(1000)
	plt.legend(labels=['full data training'])
	plt.savefig("02_full_data_training.pdf")
	plt.clf()

	adata_vis = mod.export_posterior(
	    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
	)

	# Save model
	mod.save(f"{run_name}", overwrite=True)
	# Save anndata object with results
	adata_file = f"{run_name}/sp.h5ad"
	adata_vis.write(adata_file)
