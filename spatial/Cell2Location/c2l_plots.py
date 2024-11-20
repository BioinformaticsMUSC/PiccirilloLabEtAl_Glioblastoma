import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import squidpy as sq
import pandas as pd
import cell2location
from cell2location.plt import plot_spatial
import os


base_dir = "/zfs/musc3/Sara/spatial/4_patients/analysis/"

samples = ["GBM7A", "GBM7B", "GBM8A",
			"GBM8B", "GBM9A", "GBM9B"]

samples = ["GBM7A"]

samples = ["GBM4A", "GBM4B"]

# adata_tm = sc.read("/mnt/Sara/seurat_objects/newest/full/tm_adata.h5ad")

# adata_tsvz = sc.read("/mnt/Sara/seurat_objects/newest/full/tsvz_adata.h5ad")

results_folder = './results_v2'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

###################
for s in samples: 
	if not os.path.exists(os.path.join(base_dir, s, "cell2location")):
		os.mkdir(os.path.join(base_dir, s, "cell2location"))
	os.chdir(os.path.join(base_dir, s, "cell2location"))
	os.getcwd()
	#load adata_ref
	adata_ref_file = f"{ref_run_name}/sc.h5ad"
	adata_ref = sc.read_h5ad(adata_ref_file)
	adata_ref
	#load adata_vis
	adata_vis_file = f"{run_name}/sp.h5ad"
	adata_vis = sc.read_h5ad(adata_vis_file)
	adata_vis
	#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
	adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
	sc.pl.spatial(adata_vis, cmap='viridis',
		# show first 8 cell types
		color=adata_ref.obs['Cell_Class'].unique().tolist(),
		ncols=5, size=1.3,
		img_key='hires',
		# limit color scale at 99.2% quantile of cell abundance
		vmin=0, vmax='p99.2',
		save=f"_{s}_cell_abundance.pdf"
	)
	# select up to 6 clusters
	clust_labels = ['Microglia', 'MDM', 'CancerCell','GBMmes', 'GBMnpc', 'GBMopc']
	clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
	with mpl.rc_context({'figure.figsize': (15, 15)}):
	    fig = plot_spatial(
	        adata=adata_vis,
	        # labels to show on a plot
	        color=clust_col, labels=clust_labels,
	        show_img=True,
	        # 'fast' (white background) or 'dark_background'
	        style='fast',
	        # limit color scale at 99.2% quantile of cell abundance
	        max_color_quantile=0.992,
	        # size of locations (adjust depending on figure size)
	        circle_diameter=6,
	        colorbar_position='right'
	  	)
	plt.savefig(f"03_{s}_spatial_plot.pdf")
	sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf', n_neighbors = 15)
	# Cluster spots into regions using scanpy
	sc.tl.leiden(adata_vis, resolution=1.1)
	# add region as categorical variable
	adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")
	sc.pl.spatial(adata_vis, color="region_cluster", save="region_cluster.pdf")
	#rank genes here
	adata_vis.raw = adata_vis
	sc.pp.normalize_total(adata_vis)
	sc.pp.log1p(adata_vis)
	sc.tl.rank_genes_groups(adata_vis, groupby="region_cluster", method="wilcoxon")
	df = pd.DataFrame(adata_vis.uns['rank_genes_groups']['names'])
	df.head(20).to_csv(f"{s}_new_markers.csv")
	sq.pl.spatial_scatter(adata_vis, color=["IL1RAP", "IL1B", "WNT5A", "FZD3"], 
		ncols=2,
		save='marker_genes.pdf')
	adata_vis.write(adata_vis_file)



###############
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

	results_folder = './results'

	# create paths and names to results folders for reference regression and cell2location models
	ref_run_name = f'{results_folder}/reference_signatures'
	run_name = f'{results_folder}/cell2location_map'

	print("filtering genes...")
	#filter genes to use
	from cell2location.utils.filtering import filter_genes
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

	from cell2location.models import RegressionModel
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

	mod = cell2location.models.Cell2location(
	    adata_vis, cell_state_df=inf_aver,
	    # the expected average cell abundance: tissue-dependent
	    # hyper-prior which can be estimated from paired histology:
	    N_cells_per_location=20,
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


s = "GBM7A"
for s in ["GBM7A"]: 
	if not os.path.exists(os.path.join(base_dir, s, "cell2location")):
		os.mkdir(os.path.join(base_dir, s, "cell2location"))
	os.chdir(os.path.join(base_dir, s, "cell2location"))
	os.getcwd()
	#load adata_ref
	adata_ref_file = f"{ref_run_name}/sc.h5ad"
	adata_ref = sc.read_h5ad(adata_ref_file)
	adata_ref
	#load adata_vis
	adata_vis_file = f"{run_name}/sp.h5ad"
	adata_vis = sc.read_h5ad(adata_vis_file)
	adata_vis
	#clust_options = ['Oligodendrocytes', 'NPC', "GBMac", "Endothelial"]
	#clust_labels = [c for c in clust_options if c in adata_ref.obs['Cell_Class'].cat.categories]
	clust_labels = ['Microglia', 'CancerCell']
	clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
	with mpl.rc_context({'figure.figsize': (15, 15)}):
	    fig = plot_spatial(
				adata=adata_vis,
				# labels to show on a plot
				color=clust_col, labels=clust_labels,
				show_img=False,
				# 'fast' (white background) or 'dark_background'
				style='fast',
				# limit color scale at 99.2% quantile of cell abundance
				max_color_quantile=0.992,
				# size of locations (adjust depending on figure size)
				circle_diameter=6,
				colorbar_position='right',
				reorder_cmap=[0,2,1,3,4,5,6])
	plt.savefig(f"{s}_spatial_test1.pdf")

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
			adata=adata_vis,
			# labels to show on a plot
			color=clust_col, labels=clust_labels,
			show_img=False,
			# 'fast' (white background) or 'dark_background'
			style='fast',
			# limit color scale at 99.2% quantile of cell abundance
			max_color_quantile=0.992,
			# size of locations (adjust depending on figure size)
			circle_diameter=6,
			colorbar_position='right',
			reorder_cmap=[1,3,0,4,2,5,6])

plt.savefig(f"{s}_spatial_test1.pdf")

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
			adata=adata_vis,
			# labels to show on a plot
			color=clust_col, labels=clust_labels,
			show_img=True,
			img_alpha=0.2,
			# 'fast' (white background) or 'dark_background'
			style='fast',
			# limit color scale at 99.2% quantile of cell abundance
			max_color_quantile=0.992,
			# size of locations (adjust depending on figure size)
			circle_diameter=6,
			colorbar_position='right',
			reorder_cmap=[1,3,0,4,2,5,6],
			plt_axis=False)

plt.savefig(f"{s}_spatial_test1.pdf")

def plot_pair(adata, feature1, feature2, sample):
	clust_labels = [feature1, feature2]
	clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
	plt.clf()
	with mpl.rc_context({'figure.figsize': (15, 15)}):
	    fig = plot_spatial(
				adata=adata,
				# labels to show on a plot
				color=clust_col, labels=clust_labels,
				show_img=True,
				img_alpha=0.2,
				# 'fast' (white background) or 'dark_background'
				style='fast',
				# limit color scale at 99.2% quantile of cell abundance
				max_color_quantile=0.992,
				# size of locations (adjust depending on figure size)
				circle_diameter=6,
				colorbar_position='right',
				reorder_cmap=[1,3,0,4,2,5,6],
				plt_axis=False)
	plt.savefig(f"sp_{sample}_{feature1}_{feature2}.pdf")

for s in samples: 
	if not os.path.exists(os.path.join(base_dir, s, "cell2location")):
		os.mkdir(os.path.join(base_dir, s, "cell2location"))
	os.chdir(os.path.join(base_dir, s, "cell2location"))
	os.getcwd()
	#load adata_ref
	adata_ref_file = f"{ref_run_name}/sc.h5ad"
	adata_ref = sc.read_h5ad(adata_ref_file)
	adata_ref
	#load adata_vis
	adata_vis_file = f"{run_name}/sp.h5ad"
	adata_vis = sc.read_h5ad(adata_vis_file)
	adata_vis
	tc = ['CancerCell','GBMac', 'GBMmes', 'GBMnpc', 'GBMopc', 'GBMac']
	for t in tc:
		if t in adata_ref.obs['Cell_Class'].cat.categories:
			plot_pair(adata_vis, "Microglia", t, s)

s="NSVZ"

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
sc.pl.spatial(adata_vis, cmap='viridis',
	# show first 8 cell types
	color=adata_ref.obs['Cell_Class'].unique().tolist(),
	ncols=5, size=1.3,
	img_key='hires',
	# limit color scale at 99.2% quantile of cell abundance
	vmin=0, vmax='p99.2',
	save=f"_{s}_cell_abundance.pdf"
)
	# select up to 6 clusters
clust_labels = ['Microglia', 'Neurons', 'Endothelial','Astrocytes', 'Oligodendrocytes', 'Ependymal']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
  	)

plt.savefig(f"03_{s}_spatial_plot.pdf")
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf', n_neighbors = 15)
# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=1.1)
# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")
sc.pl.spatial(adata_vis, color="region_cluster", save="region_cluster.pdf")
#rank genes here
adata_vis.raw = adata_vis
sc.pp.normalize_total(adata_vis)
sc.pp.log1p(adata_vis)
sc.tl.rank_genes_groups(adata_vis, groupby="region_cluster", method="wilcoxon")
df = pd.DataFrame(adata_vis.uns['rank_genes_groups']['names'])
df.head(20).to_csv(f"{s}_new_markers.csv")
sq.pl.spatial_scatter(adata_vis, color=["IL1RAP", "IL1B", "WNT5A", "FZD3"], 
	ncols=2,
	save='marker_genes.pdf')
adata_vis.write(adata_vis_file)



###############
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

	results_folder = './results'

	# create paths and names to results folders for reference regression and cell2location models
	ref_run_name = f'{results_folder}/reference_signatures'
	run_name = f'{results_folder}/cell2location_map'

	print("filtering genes...")
	#filter genes to use
	from cell2location.utils.filtering import filter_genes
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

	from cell2location.models import RegressionModel
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

	mod = cell2location.models.Cell2location(
	    adata_vis, cell_state_df=inf_aver,
	    # the expected average cell abundance: tissue-dependent
	    # hyper-prior which can be estimated from paired histology:
	    N_cells_per_location=20,
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



#### create cell csvs
for s in samples: 
	if not os.path.exists(os.path.join(base_dir, s, "cell2location")):
		os.mkdir(os.path.join(base_dir, s, "cell2location"))
	os.chdir(os.path.join(base_dir, s, "cell2location"))
	os.getcwd()
	#load adata_ref
	adata_ref_file = f"{ref_run_name}/sc.h5ad"
	adata_ref = sc.read_h5ad(adata_ref_file)
	adata_ref
	#load adata_vis
	adata_vis_file = f"{run_name}/sp.h5ad"
	adata_vis = sc.read_h5ad(adata_vis_file)
	adata_vis
	adata_vis.obs.to_csv(f"{s}_c2l_cells.csv")
