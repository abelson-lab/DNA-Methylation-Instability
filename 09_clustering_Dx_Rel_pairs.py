# %%
#########################################################################
# Cluster diagnosis and relapse pairs based on the most destabilized ESLs
# found at one time point.
#
# Download the required input files from Zenodo:
#########################################################################
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Patch
import scipy.stats as st
import pyreadr as pr
import pickle
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns
import random
from import_data import stable_um, recurrence_below_5

mpl.rcParams['legend.handlelength'] = 1 # 2
mpl.rcParams['legend.handleheight'] = 1 # 0.7

# %%
################################################
# LOAD ALL DX/REL PAIRED SAMPLES FOR EACH CANCER
#################################################
## AML
aml = pr.read_r('data/GSE159907/beta_values.RDS')[None]
aml.columns = [x.split('_')[0] for x in aml.columns]
with open('data/GSE159907/dx2rel_BEATAML.pkl', 'rb') as fin:
    aml_dx2rel = pickle.load(fin)

## BCP-ALL
bcp = pd.concat([pr.read_r(file)[None] for file in [
    'data/GSE49031/BCP-ALL_beta_values.RDS',
    'data/GSE49031/relapse1_beta_values.RDS']], axis=1) 
bcp.columns = [x.split('_')[0] for x in bcp.columns]
with open('data/GSE49031/dx2rel_BCP-ALL.pkl', 'rb') as fin:
    bcp_dx2rel = pickle.load(fin)

## CLL
cll = pr.read_r('data/EGAD00010000254/beta_values.RDS')[None]
with open('data/EGAD00010000254/dx2rel_CLL.pkl', 'rb') as fin:
    cll_dx2rel = pickle.load(fin)

### Create reciprocal relapse: diagnosis pairings
aml_rel2dx = {v:k for k,v in aml_dx2rel.items()}
bcp_rel2dx = {v:k for k,v in bcp_dx2rel.items()}
cll_rel2dx = {v:k for k,v in cll_dx2rel.items()}

# %%
###########################################################################
# CLUSTER ALL SAMPLES BASED ON UNIQUE SMSs IDENTIFIED IN EACH DX/REL SAMPLE
###########################################################################
def get_unique_sites(cancer, dx2rel):
    dx2cpg = {}
    for dx in dx2rel.keys():
        dx_bvals = cancer.loc[recurrence_below_5.index, dx].sort_values(ascending=False)
        dx2cpg[dx] = dx_bvals
    return dx2cpg

def cluster_dx_rel_pairs(cancer, dx2rel, num_sites, label):
    dx2cpg = get_unique_sites(cancer, dx2rel)
    cluster_cpgs = pd.concat([pd.Series(x[:num_sites].index) for x in dx2cpg.values()])
    cluster_cpgs.drop_duplicates(inplace=True)
    cluster_df = cancer.loc[cluster_cpgs, list(dx2rel.keys()) + list(dx2rel.values()) ]
    
    # cluster
    row_linkage = hierarchy.linkage(
        distance.pdist(cluster_df, metric='correlation'), method='ward')
    col_linkage = hierarchy.linkage(
        distance.pdist(cluster_df.T, metric='correlation'), method='ward')

    # check how many were correctly paired using linkage matrix
    link_mtx = pd.DataFrame(col_linkage, columns=['clustA','clustB','dist','new_n'])
    link_mtx = link_mtx[link_mtx['new_n'] == 2]
    link_mtx['diff'] = abs(link_mtx['clustA'] - link_mtx['clustB'])
    correct = link_mtx[link_mtx['diff'] == len(dx2rel)]
    # assign correct/incorrect colors
    correct_idx = list(map(int, pd.concat([correct['clustA'], correct['clustB']])))
    pair_colors = []
    for i in range(cluster_df.shape[1]):
        if i in correct_idx:
            pair_colors.append('green')
        else:
            pair_colors.append('orange')

    # plot
    cpgHM = sns.clustermap(cluster_df,
                        col_linkage=col_linkage,
                        row_linkage=row_linkage,
                        
                        # cbar_pos=(0.94,0.12,0.02,0.3), # BCP
                        # cbar_pos=(0.92,0.12,0.02,0.3), # AML
                        cbar_pos=(0.92,0.176,0.02,0.3), # CLL
                        
                        cmap='coolwarm',
                        col_colors=pair_colors,
                        # figsize=(13,10),
                        )
    cpgHM.ax_row_dendrogram.set_visible(False)

    cpgHM.ax_heatmap.axes.set_ylabel('Non-recurrent ESLs', fontsize=20)
    cpgHM.ax_heatmap.axes.yaxis.set_label_position('left')
    cpgHM.ax_heatmap.axes.set_xlabel('Patient Samples', fontsize=20)
    cpgHM.ax_heatmap.axes.tick_params(labelbottom = False, bottom=False,
                                    labelright=False, right=False)

    # nicer color bar
    cpgHM.ax_cbar.axes.set_ylabel('β-value', fontsize=17)
    cpgHM.ax_cbar.axes.tick_params(right=False)
    cpgHM.ax_cbar.axes.set_yticklabels(
        cpgHM.ax_cbar.axes.get_yticklabels(), fontsize=17)

    cpgHM.ax_col_colors.axes.legend(loc=(1.02, -1.4), frameon=False,
        handles=[Patch(facecolor='orange', edgecolor='black', label='Incorrect'),
                 Patch(facecolor='green', edgecolor='black', label='Correct')],
        fontsize=17)

    # add borders to color annotation squares
    for i, col in enumerate(cluster_df.columns):
        square = patches.Rectangle((i, 0), 1, 1, fill=None,
                                   edgecolor='black', linewidth=1)
        cpgHM.ax_col_colors.axes.add_patch(square)
    
    return cpgHM, cluster_df.iloc[cpgHM.dendrogram_row.reordered_ind,
                                  cpgHM.dendrogram_col.reordered_ind]

def random_cluster_dx_rel_pairs(cancer, dx2rel, num_sites, label, iter=10):
    dx2cpg = get_unique_sites(cancer, dx2rel)
    scores = []
    for i in range(iter):
        rand_idx = random.sample(range(len(list(dx2cpg.values())[0])), num_sites)
        cluster_cpgs = pd.concat([pd.Series(x.iloc[rand_idx].index) for x in dx2cpg.values()])
        cluster_cpgs.drop_duplicates(inplace=True)
        cluster_df = cancer.loc[cluster_cpgs, list(dx2rel.keys()) + list(dx2rel.values()) ]

        # cluster
        row_linkage = hierarchy.linkage(
            distance.pdist(cluster_df, metric='correlation'), method='ward')
        col_linkage = hierarchy.linkage(
            distance.pdist(cluster_df.T, metric='correlation'), method='ward')

        # check how many were correctly paired using linkage matrix
        link_mtx = pd.DataFrame(col_linkage, columns=['clustA','clustB','dist','new_n'])
        link_mtx = link_mtx[link_mtx['new_n'] == 2]
        link_mtx['diff'] = abs(link_mtx['clustA'] - link_mtx['clustB'])
        correct = link_mtx[link_mtx['diff'] == len(dx2rel)]
        scores.append(correct.shape[0])
    
    return scores

# %%
##################
# PERMUTATION TEST
##################
# random.seed(420)
# bcp_scores = pd.Series(random_cluster_dx_rel_pairs(bcp, bcp_rel2dx, 15, 'BCP-ALL', 1000))
# aml_scores = pd.Series(random_cluster_dx_rel_pairs(aml, aml_rel2dx, 15, 'AML', 1000))
# cll_scores = pd.Series(random_cluster_dx_rel_pairs(cll, cll_rel2dx, 15, 'CLL', 1000))

# %%
# P VALUES
# print((bcp_scores[bcp_scores >= 24].shape[0] + 1) / 1001)
# print((aml_scores[aml_scores >= 10].shape[0] + 1) / 1001)
# print((cll_scores[cll_scores >= 31].shape[0] + 1) / 1001)

# %%
### Pair samples based on most destabilized ESLs found at diagnosis
n = 15
aml_heatmap, aml_data = cluster_dx_rel_pairs(aml, aml_dx2rel, n, 'AML')
bcp_heatmap, bcp_data = cluster_dx_rel_pairs(bcp, bcp_dx2rel, n, 'BCP-ALL')
cll_heatmap, cll_data = cluster_dx_rel_pairs(cll, cll_dx2rel, n, 'CLL')

aml_heatmap.savefig('plots/figures/SuppFig_6b_AML_DxRel_clustering.png')
bcp_heatmap.savefig('plots/figures/SuppFig_6a_BCP-ALL_DxRel_clustering.png')
cll_heatmap.savefig('plots/figures/SuppFig_6c_CLL_DxRel_clustering.png')

aml_data.to_csv('plots/source_data/SuppFig_6b_AML_DxRel_clustering.csv')
bcp_data.to_csv('plots/source_data/SuppFig_6a_BCP-ALL_DxRel_clustering.csv')
cll_data.to_csv('plots/source_data/SuppFig_6c_CLL_DxRel_clustering.csv')

# %%
### Pair samples based on most destabilized ESLs found at relapse
n = 15
aml_heatmap, aml_data = cluster_dx_rel_pairs(aml, aml_rel2dx, n, 'AML')
bcp_heatmap, bcp_data = cluster_dx_rel_pairs(bcp, bcp_rel2dx, n, 'BCP-ALL')
cll_heatmap, cll_data = cluster_dx_rel_pairs(cll, cll_rel2dx, n, 'CLL')

aml_heatmap.savefig('plots/figures/ExtFig_3a_AML_DxRel_clustering.png')
bcp_heatmap.savefig('plots/figures/Fig_4a_BCP-ALL_DxRel_clustering.png')
cll_heatmap.savefig('plots/figures/ExtFig_3b_CLL_DxRel_clustering.png')

aml_data.to_csv('plots/source_data/ExtFig_3a_AML_DxRel_clustering.csv')
bcp_data.to_csv('plots/source_data/Fig_4a_BCP-ALL_DxRel_clustering.csv')
cll_data.to_csv('plots/source_data/ExtFig_3b_CLL_DxRel_clustering.csv')

# %%
###########################################################################
# CHECK REMISSION SAMPLES OF 2 PATIENTS TO CONFIRM THAT METHYLATION DROPPED
###########################################################################
all_bcp_samples = pr.read_r('data/GSE49031/patients_257_464_beta_values.RDS')[None]
all_bcp_samples.columns = [x.split('_')[0] for x in all_bcp_samples.columns]

# %%
# Patient 257 Dx, Rem, and Rel time points
patient_257 = 'GSM1203645', 'GSM1204183', 'GSM1203646'
# Patient 464 Dx, Rem, and Rel time points
patient_464 = 'GSM1203866', 'GSM1204207', 'GSM1203867'

def compare_dx_rem_rel(patient, num_sites, label):
    # dx_unique_sites = get_unique_sites(bcp, bcp_dx2rel)[patient[0]][:num_sites]    
    # patient_sites = all_bcp_samples.loc[dx_unique_sites.index, patient]
    
    fig, ax = plt.subplots(figsize=(6.5,3))
    rel_unique_sites = get_unique_sites(bcp, bcp_rel2dx)[patient[2]][:num_sites]
    patient_sites = all_bcp_samples.loc[rel_unique_sites.index, patient]
    for cpg in patient_sites.index:
        line = ax.plot(['Diagnosis', 'Remission', 'Relapse'], patient_sites.loc[cpg],
                       linestyle='solid', marker='o', color='black')

    ax.set_xlim(-0.2, 2.2)
    ax.set_ylabel('β-value', fontsize=15)
    ax.tick_params(labelsize=13)
    ax.set_title(label, fontsize=14)

    return fig, patient_sites

patient_257_fig, patient_257_data = compare_dx_rem_rel(patient_257, 15, 'Patient 257')
patient_464_fig, patient_464_data = compare_dx_rem_rel(patient_464, 15, 'Patient 464')

patient_257_fig.savefig('plots/figures/Fig_4b_patient_257.png')
patient_464_fig.savefig('plots/figures/Fig_4c_patient_464.png')
patient_257_data.to_csv('plots/source_data/Fig_4b_patient_257.csv')
patient_464_data.to_csv('plots/source_data/Fig_4c_patient_464.csv')
