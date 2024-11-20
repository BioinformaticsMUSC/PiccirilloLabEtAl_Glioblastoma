#cell2location?
import ncem as nc
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import os

from ncem.interpretation import InterpreterDeconvolution
from ncem.train import TrainModelLinearDeconvolution
from ncem.data import get_data_custom, customLoaderDeconvolution

base_dir = "/zfs/musc3/Sara/spatial/4_patients/analysis/"

samples = ["GBM4A", "GBM4B", "GBM7A", "GBM7B", "GBM8A",
			"GBM8B", "GBM9A", "GBM9B"]

for s in samples:
	datadir = f"/zfs/musc3/Sara/spatial/4_patients/analysis/{s}/cell2location/"
	os.chdir(datadir)
	print(f"working on {s}")

	adata = sc.read(datadir + f'cell2location_ncem_{s}.h5ad')
	adata.obs['library_id'] = s

	ncem_ip = InterpreterDeconvolution()

	ncem_ip.data = customLoaderDeconvolution(
	    adata=adata, patient=None, library_id='library_id', radius=None
	)

	get_data_custom(interpreter=ncem_ip, deconvolution=True)

	ncem_ip.get_sender_receiver_effects()

	ncem_ip.cv_idx = s + "v2"

	type_coupling = ncem_ip.type_coupling_analysis_circular(
	    edge_attr='magnitude', figsize=(9,8), text_space=1.28, de_genes_threshold=150, save=s,
	    edge_width_scale = 1.5
	)

	type_coupling.to_csv(f"ncem_{s}_type_coupling2.csv")
