# %%
############################################################
# Generate list of unmethylated and methylated ESLs.
#
# Requires:
# - GSE105018_beta_values.RDS
# - GSE124413_healthy_BM_beta_values.RDS
# - GSE141682_beta_values.RDS
# - GSE132203_beta_values.RDS
# - GSE103541_beta_values.RDS
# Download these files from Zenodo: <URL>
############################################################
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyreadr as pr
import pickle
from import_data import manifest, valid_probes

# %%
####################################################################
# LOAD DISCOVERY COHORT, HEALTHY CONTROLS, PURIFIED BLOOD CELLS DATA
####################################################################
# Discovery cohort (GSE105018)
with open('data/GSE105018/beta_values.pkl', 'rb') as fin:
    gse105018_beta = pickle.load(fin)
gse105018_normal_blood_counts = pd.read_csv('data/GSE105018/GSE105018_normal_blood_counts.csv', header=None).iloc[:,0]
gse105018_beta = gse105018_beta.loc[:, gse105018_normal_blood_counts]
sd = np.std(gse105018_beta, axis=1) # get SD of beta values for each CpG
mean = np.mean(gse105018_beta, axis=1) # get mean beta value for each CpG

# Healthy control datasets (Unknown, Chinese, and African descent individuals)
# Take the mean beta values of each cohort and combine the three vectors.
# This will be used to identify the 99.9th percentile of beta values in healthy individuals.
unknown = np.mean(pr.read_r('data/GSE124413/healthy_BM_beta_values.RDS')[None], axis=1)
chinese = np.mean(pr.read_r('data/GSE141682/beta_values.RDS')[None], axis=1)
african = np.mean(pr.read_r('data/GSE132203/beta_values.RDS')[None], axis=1)
healthy = pd.concat([unknown, chinese, african], axis=1)

# Purified blood cells dataset - to be used in filtering
cell_bvals = [pr.read_r(file)[None] for file in [
    'data/GSE103541/Bcells/beta_values.RDS',
    'data/GSE103541/CD4Tcell/beta_values.RDS',
    'data/GSE103541/CD8Tcell/beta_values.RDS',
    'data/GSE103541/Granulocytes/beta_values.RDS',
    'data/GSE103541/Monocytes/beta_values.RDS']
]

# %%
################################################################
# DETERMINE THE MOST STABLE (LOWEST VARIANCE) SITES ON THE ARRAY
################################################################
stable = {} # a dictionary to store SMS sets and beta value thresholds
name = 'ALL'
# get the top 10% most stable sites (i.e., the bottom 10% in terms of standard deviation)
sd_cutoff = np.quantile(sd.loc[valid_probes], 0.10)
stable_sites = sd.loc[valid_probes]
stable_sites = stable_sites[stable_sites < sd_cutoff]
selection = mean.loc[stable_sites.index]

# split them into stable unmethylated and stable methylated
stable_um = pd.Series(selection[selection < 0.5].index)
stable_m = pd.Series(selection[selection > 0.5].index)
stable[name] = [stable_um, stable_m]

# %%
##################################################################################
# FURTHER FILTERING: REMOVE CPGS THAT ARE OUTLIERS IN ANY PURIFIED BLOOD CELL TYPE
##################################################################################
for group in stable.keys():
    for num in [0,1]: # unmethylated, methylated
        sites = stable[group][num]
        mean_bvals = [np.mean(df.loc[sites,:], axis=1) for df in cell_bvals]

        # save the unfiltered data for plotting
        purified_cells_bvals = pd.DataFrame(index=stable[group][num])
        purified_cells_bvals['B'] = mean_bvals[0]
        purified_cells_bvals['CD4T'] = mean_bvals[1]
        purified_cells_bvals['CD8T'] = mean_bvals[2]
        purified_cells_bvals['Granulocyte'] = mean_bvals[3]
        purified_cells_bvals['Monocyte'] = mean_bvals[4]
        purified_cells_bvals.to_csv(f"plots/source_data/SuppFig_2_purified_cells_beta_{['U','M'][num]}_unfiltered.csv")

        # remove outliers
        all_outliers = pd.Series()
        for cell in mean_bvals:
            q1 = cell.quantile(0.25)
            q3 = cell.quantile(0.75)
            iqr = q3 - q1
            lower = q1 - 1.5*iqr
            upper = q3 + 1.5*iqr
            outliers = cell[(cell < lower) | (cell > upper)]
            all_outliers = pd.concat([all_outliers, pd.Series(outliers.index)])            
        all_outliers.drop_duplicates(inplace=True)
        stable[group][num] = pd.Series(list(set(sites) - set(all_outliers)))

        # save the filtered data for plotting
        purified_cells_bvals = pd.DataFrame(index=stable[group][num])
        purified_cells_bvals['B'] = mean_bvals[0]
        purified_cells_bvals['CD4T'] = mean_bvals[1]
        purified_cells_bvals['CD8T'] = mean_bvals[2]
        purified_cells_bvals['Granulocyte'] = mean_bvals[3]
        purified_cells_bvals['Monocyte'] = mean_bvals[4]
        purified_cells_bvals.to_csv(f"plots/source_data/Fig_1a_purified_cells_beta_{['U','M'][num]}_filtered.csv")

### Add destabilization thresholds (99.9th percentile of healthy cohorts) for each category
for group in stable.keys():
    for num in [0,1]:
        sites = stable[group][num]
        bvals = healthy.loc[sites,:].values.flatten()
        limit = np.quantile(bvals, 0.999) if num == 0 else np.quantile(bvals, 0.001)
        stable[group].append(limit)

# %%
###########
# SAVE DATA
###########
# Stable sites for use in all subsequent analyses
# with open('data/general/stable_sites.pkl', 'wb') as fout:
#     pickle.dump(stable, fout, pickle.HIGHEST_PROTOCOL)

stable['ALL'][0].to_csv('plots/source_data/ESL_U.csv', header=None, index=False)
stable['ALL'][1].to_csv('plots/source_data/ESL_M.csv', header=None, index=False)

# Mean beta values at ESLs in the discovery cohort
mean.loc[pd.concat(stable['ALL'][:2])].to_csv('plots/source_data/Fig_1b_GSE105018_mean_ESL_beta.csv', header=['Beta'])
