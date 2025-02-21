# %%
##############################################################################
# Examine ESL beta values in different human tissues using the WGBS cell atlas
# published by Loyfer et al.

# Requires a directory containing .beta files for all samples in the cell atlas,
# which can be downloaded from the Gene Expression Omnibus database (ID: GSE186458)
##############################################################################
import numpy as np
import pandas as pd
import os
import pyreadr
from import_data import stable_um, stable_m

# %%
################################################
# LOAD METHYLATION CELL ATLAS FROM LOYFER ET AL.
################################################
cpg_nums = pd.read_csv('data/GSE186458/hg38.ilmn2CpG.tsv', sep='\t', names=['probe','cpg','array'])
unmeth_nums = pd.merge(pd.DataFrame({'probe': stable_um}), cpg_nums).dropna()
meth_nums =  pd.merge(pd.DataFrame({'probe': stable_m}), cpg_nums).dropna()

atlas_dir = 'data/GSE186458/'
atlas_samples = [x for x in os.listdir(atlas_dir) if 'GSM' in x]
atlas_samples = sorted(atlas_samples, key=lambda x: x.split('_')[1])

# %%
### Create dataframe with ESL beta values for all samples
for i, probe_set in enumerate([unmeth_nums, meth_nums]):
    probes = probe_set
    all_samples_beta = []

    for file in atlas_samples:
        content = np.fromfile(atlas_dir+file, dtype=np.uint8).reshape((-1, 2))
        sample_i = probes.copy()
        sample_i['meth'] = content[sample_i['cpg'].astype('int32'),0]
        sample_i['depth'] = content[sample_i['cpg'].astype('int32'),1]
        
        id = file.split('_')[1].split('.')[0]
        sample_i[id] = sample_i['meth'] / sample_i['depth']
        all_samples_beta.append(sample_i[id])

    all_samples_beta = pd.concat(all_samples_beta, axis=1)
    all_samples_beta.index = probes['probe']

    ### Group biological replicates of the same cell type
    samples = pd.DataFrame(np.array(all_samples_beta.columns), columns=['sample'])
    samples['cell type'] = ['-'.join(col.split('-')[:-1]) for col in samples['sample']]

    ### Keep youngest sample from each cell type
    youngest = []
    metadata = pd.read_table('data/GSE186458/GSE186458_series_matrix_fmt.txt')
    metadata['id'] = [x.split('-')[-1] for x in metadata['sample']]
    for celltype in pd.unique(samples['cell type']):
        this_type = samples.loc[samples['cell type'] == celltype, 'sample']
        this_type_id = [x.split('-')[-1] for x in this_type]
        idx = np.argmin([metadata.loc[metadata['id'] == x, 'age'].iloc[0] for x in this_type_id])
        youngest.append(this_type.iloc[idx])

    ### Save data for plotting
    youngest_beta = all_samples_beta[youngest]
    youngest_beta.columns = ['-'.join(x.split('-')[:-1]) for x in youngest_beta.columns]
    youngest_beta.to_csv(f"plots/source_data/Fig_1c_methylation_atlas_beta_{['U','M'][i]}.csv")
