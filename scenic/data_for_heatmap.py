import os
import pandas as pd
import loompy as lp

# set a working directory
wdir = '${SINGULARITY_BIND}/Sara/scenic/new/normal'
os.chdir(wdir)

#open loom files, get AUC data
#tsvz normal
f_final_loom = 'final_scenic_tsvz_normal.loom'
lf = lp.connect( f_final_loom, mode='r', validate=False )
auc_mtx_tsvz_normal = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

#tm normal
f_final_loom = 'final_scenic_tm_normal.loom'
lf = lp.connect( f_final_loom, mode='r', validate=False )
auc_mtx_tm_normal = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

#tsvz tumor
f_final_loom = '../final_scenic_tsvz_tumor.loom'
lf = lp.connect( f_final_loom, mode='r', validate=False )
auc_mtx_tsvz_tumor = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

#tm tumor
f_final_loom = '../final_scenic_tm_tumor.loom'
lf = lp.connect( f_final_loom, mode='r', validate=False )
auc_mtx_tm_tumor = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

#AUC mean heatmap
#both auc_mtxs loaded separately
tm_normal_avg = auc_mtx_tm_normal.mean()
tsvz_normal_avg = auc_mtx_tsvz_normal.mean()
tm_tumor_avg = auc_mtx_tm_tumor.mean()
tsvz_tumor_avg = auc_mtx_tsvz_tumor.mean()

data = pd.concat([tm_tumor_avg, tsvz_tumor_avg, tm_normal_avg, tsvz_normal_avg], axis=1)
data.columns = ['TM_tumor', 'TSVZ_tumor', 'TM_normal', 'TSVZ_normal']

# z scores
for col in data.columns:
    avg = data[col].mean()
    sd = data[col].std()
    data[f"{col}_z"] = data[col].apply(lambda x: (x - avg) / sd)

tm_subset = data.sort_values(by="TM_tumor", ascending=False)[["TM_tumor_z", "TM_normal_z"]][:30]
tsvz_subset = data.sort_values(by="TSVZ_tumor", ascending=False)[["TSVZ_tumor_z", "TSVZ_normal_z"]][:30]

heat_data = pd.concat([tm_subset, tsvz_subset], axis=1)
heat_data = heat_data[['TM_tumor_z', "TSVZ_tumor_z", "TM_normal_z", "TSVZ_normal_z"]]
heat_data.columns = ['TM_tumor', 'TSVZ_tumor', 'TM_normal', 'TSVZ_normal']

heat_data = heat_data[['TM_tumor', "TM_normal", "TSVZ_tumor", "TSVZ_normal"]]
heat_data.to_csv("auc_heatmap_data_processed.csv")