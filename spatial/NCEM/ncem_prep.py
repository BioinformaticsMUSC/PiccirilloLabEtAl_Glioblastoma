

cell_types = ['CancerCell', 'Endothelial', 'GBMac', 'GBMmes', 'GBMnpc', 'GBMopc', 'MDM', 'Microglia', 'Neurons', 'Oligodendrocytes']
prop = adata_vis.obs[cell_types]

cell_expression = []
node_types = []
proportions = []
spatial = []
for i, ct in enumerate(cell_types):
    proportions.append(prop)
    cell_expression.append(adata_vis.layers[ct].toarray())
    nt = np.zeros((prop.shape[0], len(cell_types)))
    nt[:, i] = 1
    node_types.append(nt)
    spatial.append(adata_vis.obsm['spatial'])
    
proportions = pd.DataFrame(np.concatenate(proportions), columns=cell_types)
cell_expression = pd.DataFrame(np.concatenate(cell_expression), columns=adata_vis.var_names)
node_types = pd.DataFrame(np.concatenate(node_types), columns=cell_types)
spatial = pd.DataFrame(np.concatenate(spatial))

