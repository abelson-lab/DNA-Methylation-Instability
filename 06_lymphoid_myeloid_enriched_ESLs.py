# %%
#######################################################
# Determine lymphoid-enriched and myeloid-enriched ESLs,
# subsequently used for scATAC-seq analysis
#######################################################
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pyreadr as pr
import pickle
from scipy import stats as st
from scipy.spatial import distance
from scipy.cluster import hierarchy
from import_data import pan_cancer_binarized, pan_cancer_beta

# %%
########################################
# Select AML, T-ALL, and BCP-ALL cohorts
########################################
aml = pan_cancer_binarized.iloc[:, :997]
tall = pan_cancer_binarized.iloc[:, 997 : 997+353]
bcp = pan_cancer_binarized.iloc[:, 997+353 : 997+353+663]

aml_beta = pan_cancer_beta.iloc[:, :997]
tall_beta = pan_cancer_beta.iloc[:, 997 : 997+353]
bcp_beta = pan_cancer_beta.iloc[:, 997+353 : 997+353+663]

# %%
##############################################################
# FUNCTIONS TO CALCULATE ESL ENRICHMENT FOR CANCER COMPARISONS
##############################################################
def esl_enrichment(a, b, names):
    '''Use Fisher's test to calculate odds ratio for each ESL between two cancers'''
    a_name, b_name = names
    results = pd.DataFrame(columns=[a_name, b_name, 'Difference', 'OR', 'p'])

    for i in a.index:
        a_yes = np.sum(a.loc[i])
        a_no = a.shape[1] - a_yes
        b_yes = np.sum(b.loc[i])
        b_no = b.shape[1] - b_yes
        a_proportion = a_yes / a.shape[1]
        b_proportion = b_yes / b.shape[1]
        diff = a_proportion - b_proportion        
        odds_ratio, p_value = st.fisher_exact([[a_yes, a_no], [b_yes, b_no]])
        results.loc[i] = a_proportion, b_proportion, diff, odds_ratio, p_value    
    
    results['p'] = results['p'] * a.shape[0] # bonferroni 
    return results

def get_enriched_depleted_cpgs(results, min_proportion):
    # select all sites with OR > 1 and min_proportion in column A
    or_values = results.loc[(results['OR'] > 1) & (results.iloc[:,0] > min_proportion) , 'OR']
    threshold = np.quantile(or_values, 0.5)
    enriched = pd.Series(results[(results['OR'] > threshold) & (results.iloc[:,0] > min_proportion)].index)
    # select all sites with OR < 1 and min_proportion in column B
    or_values = results.loc[(results['OR'] < 1) & (results.iloc[:,1] > min_proportion) , 'OR']
    threshold = np.quantile(or_values, 0.5)
    depleted = pd.Series(results[(results['OR'] < threshold) & (results.iloc[:,1] > min_proportion)].index)
    return enriched, depleted

# %%
#####################################################
# IDENTIFY LYPHOID-ENRICHED AND MYELOID-ENRICHED ESLs
#####################################################
# Identify sites enriched in lymphoid cancers vs myeloid cancers
lymphoid_vs_myeloid = esl_enrichment(pd.concat([tall, bcp], axis=1), aml, ['Lymphoid', 'Myeloid'])
lymphoid_enriched = get_enriched_depleted_cpgs(lymphoid_vs_myeloid, 0.5)[0]

# To determine myeloid-enriched, compare AML to T-ALL and BCP-ALL separately instead of as one "lymphoid" group
# Due to the cohort sizes of T-ALL (353) and BCP-ALL (663), a site may appear as myeloid-enriched relative to 
# the combined lymphoid group, when actually T-ALL is enriched for that site but BCP-ALL is not.
aml_vs_tall = esl_enrichment(aml, tall, ['AML', 'T-ALL'])
aml_vs_bcp = esl_enrichment(aml, bcp, ['AML', 'BCP-ALL'])

# Minimum proportion threshold is set to 10% since AML has much lower destabilization than ALL
# Take the intersection of the two lists to get a final myeloid-enriched set
aml_vs_tall_enriched = get_enriched_depleted_cpgs(aml_vs_tall, 0.1)[0]
aml_vs_bcp_enriched = get_enriched_depleted_cpgs(aml_vs_bcp, 0.1)[0]
myeloid_enriched = aml_vs_tall_enriched[aml_vs_tall_enriched.isin(aml_vs_bcp_enriched)]

# %%
##############################################
# SAVE CPG SETS FOR USE IN scATAC-seq ANALYSIS
##############################################
hg38cpgs = pd.read_table('data/scATAC/hg38/450k_cpgs_hg38.bed', header=None,
                         names=['chr','start','end','cpg'])
hg19cpgs = pd.read_table('data/scATAC/hg19/450k_cpgs_hg19.bed', header=None,
                         names=['chr','start','end','cpg'])

# for sites, name in zip([lymphoid_enriched, myeloid_enriched], ['LymphoidEnriched', 'MyeloidEnriched']):
#     f = f'data/scATAC/{name}-{sites.shape[0]}'
#     hg38cpgs[hg38cpgs['cpg'].isin(sites)].to_csv(f'{f}-hg38.bed', index=False, header=None, sep='\t')
#     hg19cpgs[hg19cpgs['cpg'].isin(sites)].to_csv(f'{f}-hg19.bed', index=False, header=None, sep='\t')

# %%
#####################################
# HEATMAP OF LINEAGE-ASSOCIATED SITES
#####################################
colors = ['#b2df8a', '#6a3d9a', '#ff7f00']
colColors = []
for i, cancer in enumerate([aml, tall, bcp]): colColors.extend([colors[i]]*cancer.shape[1])

# sites, fname = lymphoid_enriched, 'SuppFig_4a_lymphoid_enriched_sites' 
sites, fname = myeloid_enriched, 'SuppFig_4b_myeloid_enriched_sites'

df = pd.concat([aml, tall, bcp], axis=1).loc[sites]

# plot heatmap
cpgHM = sns.clustermap(
                    df.values,
                    col_colors=colColors,
                    col_cluster=False,
                    row_cluster=False,
                    cbar_pos=None,
                    # cmap='coolwarm'
                    )
cpgHM.ax_row_dendrogram.set_visible(False)
cpgHM.ax_col_dendrogram.set_visible(False)
cpgHM.ax_heatmap.axes.set_ylabel('CpG sites', fontsize=17)
cpgHM.ax_heatmap.axes.set_xlabel('Patient Samples', fontsize=17)

countby = 100 if df.shape[0] > 500 else 10
cpgHM.ax_heatmap.axes.set_yticks(range(0, df.shape[0], countby))
cpgHM.ax_heatmap.axes.set_yticklabels(range(0, df.shape[0], countby),
                                      fontsize=15)

cpgHM.ax_heatmap.axes.tick_params(labelbottom = False, bottom=False,
                                labelright=True, right=True)
# save figure and source data
df.to_csv(f'plots/source_data/{fname}.csv')
cpgHM.savefig(f'plots/figures/{fname}.png')
