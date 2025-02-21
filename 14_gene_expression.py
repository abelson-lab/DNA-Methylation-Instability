# %%
########################################################################
# Conduct methylation vs. expression, age vs. expression, and promoter
# hypermethylation analyses.
#
# Download the required input files from Zenodo:
########################################################################
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pyreadr as pr
import pickle
from scipy import stats as st
from scipy.spatial import distance
from scipy.cluster import hierarchy
from import_data import stable_um, recurrence, manifest
from joblib import Parallel, delayed

# %%
### LOAD PANCANCER METHYLATION VS. EXPRESSION RESULTS
with open('data/gene_expression/pan-cancer_meth_vs_expression.pkl', 'rb') as fin:
    methExp = pickle.load(fin)
mani = manifest.loc[methExp.index].copy()

# %%
### LOAD HEALTHY AGING AND AML RNAseq DATA
healthySamples = pd.read_csv('data/gene_expression/GTEx_blood_sample_annotations.csv')
ages = pd.read_csv('data/gene_expression/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t')
healthySamples['SUBJID'] = ['-'.join(healthySamples.at[i, 'SAMPID'].split('-')[:2]) for i in healthySamples.index]
healthyAge = pd.merge(healthySamples, ages, on='SUBJID')
healthyTpms = pd.read_csv('data/gene_expression/GTEx_blood_sample_TPMs.csv')
healthyTpms.set_index('Description', inplace=True)

amlTpms = pd.read_csv('data/gene_expression/TCGA-LAML_TPMs.csv')
amlTpms.set_index('gene_name', inplace=True)

groups = sorted(list(pd.unique(healthyAge['AGE'])))
labels_n = [x+'\nn='+str(len(healthyAge[healthyAge['AGE'] == x])) for x in groups] \
    + [f'AML\nn={len(amlTpms.columns)-1}']

# %%
## GET LIST OF ALL GENES PRESENT IN HEALTHY AND AML DATASETS
healthy_genes = pd.Series(healthyTpms.index).drop_duplicates()
aml_genes = pd.Series(amlTpms.index).drop_duplicates()
tpm_genes = healthy_genes[healthy_genes.isin(aml_genes)]
tpm_genes

# %%
### GET ALL CGI-GENE PAIRS
wrk = pd.merge(methExp, manifest[['UCSC_CpG_Islands_Name']],
                                 left_index=True, right_index=True)
wrk = wrk[wrk['gene'].isin(tpm_genes)]
wrk.dropna(inplace=True)
wrk

# %%
gene_cgi = wrk.drop_duplicates(subset=['gene', 'UCSC_CpG_Islands_Name']).copy()
gene_cgi['cpgs'] = pd.NA
gene_cgi['mean_r'] = pd.NA
# get list of ESLs corresponding to each cgi-gene pair
for i in gene_cgi.index[:]:
    cgi = gene_cgi.at[i, 'UCSC_CpG_Islands_Name']
    cgi_esls = mani.loc[mani['UCSC_CpG_Islands_Name'] == cgi, 'Name']
    gene_cgi.at[i, 'cpgs'] = list(cgi_esls.index)

    gene_cgi.at[i, 'mean_r'] = np.mean(methExp.loc[cgi_esls.index, 'r'])
    gene_cgi.at[i, 'max_p'] = np.max(methExp.loc[cgi_esls.index, 'p'])

gene_cgi.drop(labels=['r', 'p'], axis=1, inplace=True)
gene_cgi

# %%
gene_cgi['ageExpR'] = pd.NA
gene_cgi['ageExpP'] = pd.NA

### ADD AGE VS. EXPRESSION CORRELATION
for i in gene_cgi.index:
    tpmdata = []
    ages = []
    gene = gene_cgi.at[i, 'gene']
    for agegroup in groups:
        try:
            thisgroup = healthyAge.loc[healthyAge['AGE'] == agegroup, 'SAMPID']
            tpmdata.extend(list(healthyTpms.loc[gene, thisgroup]))
            ages.extend( len(thisgroup) * [int(agegroup[:2])] )
        except:
             print(i, gene)
    try:
        r, p = st.pearsonr(ages, tpmdata)
        gene_cgi.at[i, 'ageExpR'] = r
        gene_cgi.at[i, 'ageExpP'] = p
    except:
        pass # both will remain NA
gene_cgi

# %%
### ADD MEAN RECURRENCE ###
gene_cgi['mean_recurrence'] = pd.NA
for i in gene_cgi.index:
    cpgs = gene_cgi.at[i, 'cpgs']
    gene_cgi.at[i, 'mean_recurrence'] = np.mean(recurrence.loc[cpgs])
gene_cgi

# %%
## ADD HG38 CGI COORDINATES
hg38cgis = pd.read_table('data/gene_expression/all_cgis_450k_toHg38.bed',
                         header=None, names=['chr','start','end','id'], index_col='id')

gene_cgi['chr'] = pd.NA
gene_cgi['start'] = pd.NA
gene_cgi['end'] = pd.NA
for i in gene_cgi.index:
    chrm, start, end = hg38cgis.loc[gene_cgi.at[i, 'UCSC_CpG_Islands_Name']]
    gene_cgi.at[i, 'chr'] = chrm
    gene_cgi.at[i, 'start'] = start
    gene_cgi.at[i, 'end'] = end

# %%
##### CALCULATE LOG2FC HYPERMETHYLATION IN AML WGBS ####
# chroms = ['chr'+str(i) for i in range(1,23)] # + ['chrX', 'chrY']
# with open('data/gene_expression/aml_and_control_WGBS.pkl', 'rb') as fin:
#     allSamples = pickle.load(fin)

# # iterate through CGIs and calculate p value, logFC for each
# def dmr_cgis(chrm):
#     results = hg38cgis.copy()
#     results['-log10p'] = np.nan
#     results['log2FC'] = np.nan
#     thisSamples = allSamples[allSamples['chr'] == chrm]
#     results = results[results['chr'] == chrm]
#     for i in results.index:
#         chrom, start, end = results.loc[i, ['chr','start','end']]
#         start, end = int(start), int(end)
#         selection = thisSamples[(thisSamples['chr'] == chrom) &
#                             (thisSamples['start'] >= start) &
#                             (thisSamples['end'] <= end)]
        
#         healthy = selection.iloc[:,3:9].values.flatten()
#         healthy = healthy[~np.isnan(healthy)]
#         aml = selection.iloc[:,9:15].values.flatten()
#         aml = aml[~np.isnan(aml)]
#         try:
#             result = st.mannwhitneyu(healthy, aml)
#             results.at[i, '-log10p'] = -np.log10(result[1])
#             results.at[i, 'log2FC'] = np.log2(np.mean(aml) / np.mean(healthy))
#         except:
#             results.at[i, '-log10p'] = pd.NA
#             results.at[i, 'log2FC'] = pd.NA
#     print(chrm)
#     return results

### run in parallel
# dfs = Parallel(n_jobs=10)(delayed(dmr_cgis)(c) for c in chroms[:])
# results = pd.concat(dfs, axis=0)
# results.to_csv('data/gene_expression/WGBS_results_MWU.csv')

# %%
### ADD RESULTS FROM AML WGBS HYPERMETHYLATION ANALYSIS ###
results = pd.read_csv('data/gene_expression/WGBS_results_MWU.csv', index_col='id')
gene_cgi['-log10p'] = pd.NA
gene_cgi['log2FC'] = pd.NA
for i in gene_cgi.index:
    cgi = gene_cgi.at[i, 'UCSC_CpG_Islands_Name']
    score, fc = results.loc[cgi, ['-log10p', 'log2FC']]
    gene_cgi.at[i, '-log10p'] = score
    gene_cgi.at[i, 'log2FC'] = fc
gene_cgi

# %%
### WRITE DATAFRAME TO SUPPLEMENTARY TABLE 3
supp3 = gene_cgi.rename(columns={'mean_r': 'meth_vs_exp_R_mean',
                         'max_p': 'meth_vs_exp_maxPvalue',
                         'ageExpR': 'age_vs_exp_R',
                         'ageExpP': 'age_vs_exp_Pvalue',
                         'log2FC': 'WGBS_log2FC',
                         '-log10p': 'WGBS_-log10(Pvalue)',
                         'chr': 'CGI_chr',
                         'start': 'CGI_start',
                         'end': 'CGI_end',
                         })

supp3 = supp3[['gene', 'CGI_chr', 'CGI_start', 'CGI_end',
               'mean_recurrence', 'meth_vs_exp_R_mean', 'meth_vs_exp_maxPvalue',
               'age_vs_exp_R', 'age_vs_exp_Pvalue',
               'WGBS_log2FC', 'WGBS_-log10(Pvalue)',               
               'UCSC_CpG_Islands_Name', 'cpgs'
               ]]

# bonferroni for ageExpP value
supp3['age_vs_exp_Pvalue'] = supp3['age_vs_exp_Pvalue'] * gene_cgi['gene'].drop_duplicates().shape[0]

supp3.to_csv('data/gene_expression/SupplementaryTable3_FULL.csv', index=False)

# %%
### 3D SCATTERPLOT ###
df = pd.read_csv('data/gene_expression/SupplementaryTable3_FULL.csv')
df = df[df['mean_recurrence'] > 0.05]
df.sort_values(by='mean_recurrence', ascending=False, inplace=True)
df.drop_duplicates(subset='gene', inplace=True)
genes234 = pd.read_csv('data/gene_expression/SupplementaryTable3.csv').iloc[:,0].drop_duplicates()
top_genes = df[df['gene'].isin(genes234)]
all_other_genes = df[~df['gene'].isin(genes234)]

### PLOT 3 VARIABLES ON 3D AXES
fig = plt.figure() #figsize=(5,5)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(all_other_genes['age_vs_exp_R'], all_other_genes['meth_vs_exp_R_mean'], all_other_genes['mean_recurrence'],
           c='lightblue', marker='.',)
ax.scatter(top_genes['age_vs_exp_R'], top_genes['meth_vs_exp_R_mean'], top_genes['mean_recurrence'],
           c='red', marker='.', zorder=99999)

ax.tick_params(labelsize=10)
ax.set_xlabel(ax.get_xlabel(), fontsize=11)
ax.set_ylabel(ax.get_ylabel(), fontsize=11)
ax.set_zlabel(ax.get_zlabel(), fontsize=11)
ax.set_xlim(-0.41,0.29)
ax.set_zlim(0,0.82)

## LABEL A POINT
pt_idx = df[df['gene'] == 'PRDM5'].index[0]
prdm5_x = df.at[pt_idx, 'age_vs_exp_R']
prdm5_y = df.at[pt_idx, 'meth_vs_exp_R_mean']
prdm5_z = df.at[pt_idx, 'mean_recurrence']

ax.text(x=prdm5_x,
        y=prdm5_y,
        z=prdm5_z,
        s='PRDM5', color='black', zorder=999,
        rotation=45,
        ha='left', fontsize=12)

ax.view_init(elev=30, azim=225)
ax.set_xlabel('Age vs. Expression R')
ax.set_ylabel('Methylation vs. Expression R')
ax.set_zlabel('Recurrence')

fig.savefig('plots/figures/Fig_6c_gene_expression_3D_scatter.png')
source_data = pd.concat([top_genes, all_other_genes])[[
    'gene', 'mean_recurrence', 'meth_vs_exp_R_mean',
    'meth_vs_exp_maxPvalue', 'age_vs_exp_R', 'age_vs_exp_Pvalue']]
source_data.to_csv('plots/source_data/Fig_6c_gene_expression_3D_scatter.csv', index=False)
