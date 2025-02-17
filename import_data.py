# %%
##############################################################
# Auxiliary data that is used throughout the various analyses.
# Other scripts will import this data from here when needed.

# Requires:
# - humanmethylation450_15017482_v1-2.csv (Illumina manifest)
# - pan-cancer_beta_binarized.pkl
# - pan-cancer_beta.pkl
# Download these files from Zenodo: <URL>
##############################################################
import pandas as pd
import numpy as np
import pickle

# %%
### Load 450K manifest
manifest = pd.read_csv('data/humanmethylation450_15017482_v1-2.csv',
                      skiprows=7, skipfooter=916, engine='python')
manifest['Relation_to_UCSC_CpG_Island'].fillna('Open Sea', inplace=True)
manifest['Relation_to_UCSC_CpG_Island'].replace('S_Shore', 'Shore', inplace=True)
manifest['Relation_to_UCSC_CpG_Island'].replace('N_Shore', 'Shore', inplace=True)
manifest['Relation_to_UCSC_CpG_Island'].replace('S_Shelf', 'Shelf', inplace=True)
manifest['Relation_to_UCSC_CpG_Island'].replace('N_Shelf', 'Shelf', inplace=True)

### Only keep "valid" probes (no XY probes, no SNP loci, no epigenetic clock CpGs, etc.)
valid_probes = pd.read_csv('data/valid_probes.csv', header=None).iloc[:,0]
manifest = manifest[manifest['IlmnID'].isin(valid_probes)]
manifest.set_index('IlmnID', inplace=True)

# %%
### Load stable sites
with open('data/stable_sites.pkl', 'rb') as fin:
    stable_sites = pickle.load(fin)
stable_um = stable_sites['ALL'][0]
stable_m = stable_sites['ALL'][1]
non_esl = pd.Series(list(set(manifest.index) - set(stable_um) - set(stable_m)))

# %%
### Calculate recurrence from pan-cancer binarized matrix
with open('data/pan-cancer_beta_binarized.pkl', 'rb') as fin:
    pan_cancer_binarized = pickle.load(fin)
recurrence_counts = np.sum(pan_cancer_binarized, axis=1)
recurrence = recurrence_counts / pan_cancer_binarized.shape[1]
recurrence_above_5 = recurrence[recurrence > 0.05]
recurrence_below_5 = recurrence[recurrence < 0.05]

### Original pan-cancer beta value matrix
with open('data/pan-cancer_beta.pkl', 'rb') as fin:
    pan_cancer_beta = pickle.load(fin)
