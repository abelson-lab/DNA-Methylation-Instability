# %%
#############################################################################
# Create the matrix of ESL beta values across the pan-blood cancer cohort.
#
# Requires:
# - TCGA-LAML_beta_values.RDS
# - GSE159907_Dx_beta_values.RDS
# - GSE124413_AML_beta_values.RDS
# - GSE153347_beta_values.RDS
# - GSE147667_beta_values.RDS
# - GSE155333_beta_values.RDS
# - GSE49031_T-ALL_beta_values.RDS
# - GSE49031_BCP-ALL_beta_values.RDS
# - GSE143411_timepoint1_beta_values.RDS
# - GSE136724_beta_values.RDS
# - FL_beta_values.RDS
# - GSE104770_beta_values.RDS
# - GSE197696_beta_values.RDS
# - GSE113545_beta_values.RDS
# - MDS_beta_values.RDS
# - CML_beta_values.RDS
# Download these files from Zenodo: <URL>

# RESTRICTED ACCESS:
# - EGA_CLL_490_beta_values.RDS
#     - original data at European Genome-Phenome Archive: EGAD00010001975
# - EGA_CLL_40_beta_values.RDS
#     - original data at European Genome-Phenome Archive: EGAD00010000254
#############################################################################
import numpy as np
import pandas as pd
import pyreadr
import copy
import pickle
from import_data import stable_sites, stable_um, stable_m

# %%
########################################
# LOAD ALL HEMATOLOGICAL CANCER DATASETS
########################################
# Acute myeloid leukemia (n = 997)
tcga_laml = pyreadr.read_r('data/TCGA-LAML/beta_values.RDS')[None]
with open('data/GSE159907/Dx_beta_values.pkl', 'rb') as fin: 
    gse159907 = pickle.load(fin)
gse124413 = pyreadr.read_r('data/GSE124413/AML_beta_values.RDS')[None]
with open('data/GSE153347/unique_pts_beta_values.pkl', 'rb') as fin: 
    gse153347 = pickle.load(fin)
aml = [tcga_laml, gse159907, gse124413, gse153347]
aml = pd.concat(aml, axis=1)

# T-cell acute lymphoblastic leukemia (n = 353)
gse147667 = pyreadr.read_r('data/GSE147667/beta_values.RDS')[None]
gse155333 = pyreadr.read_r('data/GSE155333/beta_values.RDS')[None]
gse49031_tall = pyreadr.read_r('data/GSE49031/T-ALL_beta_values.RDS')[None]
tall = pd.concat([gse147667, gse155333, gse49031_tall], axis=1)

# B-cell precursor acute lymphoblastic leukemia (n = 663)
bcpall = pyreadr.read_r('data/GSE49031/BCP-ALL_beta_values.RDS')[None]

# Chronic lymphocytic leukemia (n = 612)
gse143411 = pyreadr.read_r('data/GSE143411/timepoint1_beta_values.RDS')[None]
gse136724 = pyreadr.read_r('data/GSE136724/beta_values.RDS')[None]
egaCLL490 = pyreadr.read_r('data/EGAD00010001975/beta_values.RDS')[None]
with open('data/EGAD00010000254/Dx_beta_values.pkl', 'rb') as fin:
    egaCLL40 = pickle.load(fin)
cll = pd.concat([gse143411, gse136724, egaCLL490, egaCLL40], axis=1)

# Follicular lymphoma (n = 246)
fl = pyreadr.read_r('data/FL/beta_values.RDS')[None]
# Diffuse Large B-cell Lymphoma (n = 96)
dlbcl1 = pyreadr.read_r('data/GSE255869/beta_values.RDS')[None]
# Primary plasma cell leukemia (n = 14)
ppcl = pyreadr.read_r('data/GSE104770/beta_values.RDS')[None]
# Myeloid / natural killer precursor leukemia (n = 7)
mnkpl = pyreadr.read_r('data/GSE197696/beta_values.RDS')[None]
# Mixed phenotype acute leukemia (n = 31)
mpal = pyreadr.read_r('data/GSE113545/beta_values.RDS')[None]

# %%
############################################################
# COMBINE DATASETS, SUBSET TO ESLs, CREATE BINARIZED VERSION
############################################################
all_cancers = pd.concat([aml, tall, bcpall, cll, fl, dlbcl1, ppcl, mnkpl, mpal], axis=1)                                                                       # mds, cml
all_cancers = all_cancers.loc[stable_um]
all_cancers_binarized = (all_cancers > stable_sites['ALL'][2]).astype(int)

# %%
### Save data
# with open('data/general/pan-cancer-3019_beta.pkl', 'wb') as fout:
#     pickle.dump(all_cancers, fout, pickle.HIGHEST_PROTOCOL)
# with open('data/general/pan-cancer-3019_beta_binarized.pkl', 'wb') as fout:
#     pickle.dump(all_cancers_binarized, fout, pickle.HIGHEST_PROTOCOL)
# all_cancers_binarized.to_csv('supplementary/SupplementaryTable2.csv')

# %%
###### ANNOTATION OF COLUMNS ######
dfs = [tcga_laml, gse159907, gse124413, gse153347, 
       gse147667, gse155333, gse49031_tall,
       bcpall,
       gse143411, gse136724, egaCLL490, egaCLL40,
       fl, dlbcl1, ppcl, mnkpl, mpal]

names = ['AML']*4 + \
        ['T-ALL']*3 + \
        ['BCP-ALL'] + \
        ['CLL']*4 + \
        ['FL', 'DLBCL', 'pPCL', 'MNKPL', 'MPAL']

ids = ['TCGA-LAML', 'GSE159907', 'GSE124413', 'GSE153347',
       'GSE147667', 'GSE155333', 'GSE49031',
       'GSE49031',
       'GSE143411', 'GSE136724', 'EGAD00010001975', 'EGAD00010000254',
       'Shelton et al., 2024', 'GSE255869', 'GSE104770', 'GSE197696', 'GSE113545']

column_anno = pd.DataFrame(columns = ['column_name', 'cancer', 'dataset'])
for df, name, id in zip(dfs, names, ids):
    cols = df.columns
    cancer = [name] * df.shape[1]
    dataset = [id] * df.shape[1]
    this_df = pd.DataFrame({'column_name': cols,
                            'cancer': cancer,
                            'dataset': dataset})
    column_anno = pd.concat([column_anno, this_df], axis=0)

column_anno['index'] = range(1, column_anno.shape[0]+1)
column_anno = column_anno[['index'] + list(column_anno.columns[:-1])]
# column_anno.to_csv('supplementary/SupplementaryTable2_column_annotations.csv', index=False)
