# %%
##############################################################################
# Examine the correlation between recurrence and beta values in the pan-cancer 
# cohort and in specific cancer types.
##############################################################################
import numpy as np
import pandas as pd
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.stats as st
from import_data import pan_cancer_beta, pan_cancer_binarized, recurrence

# %%
##############################
# Select specific cancer types
##############################
def get_cancer_specific_data(start, end):
    '''generate cancer_beta, cancer_binarized, and recurrence from a subset of the pan-cancer data'''
    cancer_beta = pan_cancer_beta.iloc[:, start:end]
    cancer_binarized = pan_cancer_binarized.iloc[:, start:end]
    recurrence_counts = np.sum(cancer_binarized, axis=1)
    recurrence_prop = recurrence_counts / cancer_binarized.shape[1]
    return cancer_beta, cancer_binarized, recurrence_prop # use cohort-specific recurrence

aml = get_cancer_specific_data(0, 997)
tall = get_cancer_specific_data(997, 997+353)
bcp = get_cancer_specific_data(997+353, 997+353+663)
cll = get_cancer_specific_data(997+353+663, 997+353+663+612)
fl = get_cancer_specific_data(997+353+663+612, 997+353+663+612+246)
pan_cancer = pan_cancer_beta, pan_cancer_binarized, recurrence

def get_mean_destabilized_beta(cancer_beta, cancer_binarized, recurrence_values, r):
    cpgs = recurrence_values[recurrence_values == r].index
    all_cpg_means = []
    for cpg in cpgs:
        patient_destabilized = cancer_binarized.loc[cpg]
        patient_destabilized = patient_destabilized[patient_destabilized == 1].index
        cpg_mean_beta = np.mean([cancer_beta.at[cpg, patient] for patient in patient_destabilized])
        all_cpg_means.append(cpg_mean_beta)
    grouped_mean_beta = np.mean(all_cpg_means)
    return grouped_mean_beta

def get_data(cancer_beta, cancer_binarized, recurrence_values, label):
    data = pd.DataFrame({'recurrence': pd.unique(recurrence_values),
                        'mean_beta': [get_mean_destabilized_beta(cancer_beta, cancer_binarized, recurrence_values, r) \
                                    for r in pd.unique(recurrence_values)]})
    data.dropna(inplace=True)
    return data

def plot(cancer_beta, cancer_binarized, recurrence_values, label, data, gs=None):
    if gs == None:
        fig, (ax, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [5, 1], 'hspace':0.15})
    else:
        ax, ax2 = [plt.subplot(x) for x in gs]

    # plot beta values against recurrence
    ax.scatter(data['recurrence'], data['mean_beta'], marker='.', s=8, color='blue')
    ax.get_xaxis().set_visible(False)
    ax.set_ylabel('Mean β-value', fontsize=14)
    # of destabilized sites

    # get correlation
    slope, intercept, r, p = st.linregress(data['recurrence'], data['mean_beta'])[:4]
    print(r, p)
    ax.set_title(f'{label}: r = {round(r, 2)}, n = {cancer_binarized.shape[1]}',
                    fontsize=14)
    # plot regression line
    x_regress = np.linspace(0, np.max(data['recurrence']), 1000)
    y_regress = slope*x_regress + intercept
    ax.plot(x_regress, y_regress, color='green')
    ax.tick_params(labelsize=12)

    # plot # ESLs against recurrence
    log10_counts = np.log10([recurrence_values[recurrence_values == x].shape[0] for x in data['recurrence']])
    ax2.bar(data['recurrence'], log10_counts, width=0.003, color='red')
    ax2.set_ylabel('log\n(# ESLs)', fontsize=11)
    ax2.set_xlabel('Recurrence', fontsize=14)
    ax2.tick_params(labelsize=12)

    if gs == None:
        fig.tight_layout()

# %%
### generate data to plot
cancer_datas = []
for label, data in zip(['Pan-cancer', 'AML', 'T-ALL', 'BCP-ALL', 'CLL', 'FL'],
                       [pan_cancer, aml, tall, bcp, cll, fl]):
    cancer_datas.append(get_data(data[0], data[1], data[2], label))

# %%
### plot into 2x3 grid
fig = plt.figure(figsize=(18, 9))
outer = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.2)
gs = []
for i in range(6):
    gs.append(gridspec.GridSpecFromSubplotSpec(2, 1,
                                               height_ratios=[5,1],
                                               subplot_spec = outer[i],
                                               hspace=0.05))

for gs_i, label, data, preload in zip(gs,
                          ['Pan-cancer', 'AML', 'T-ALL', 'BCP-ALL', 'CLL', 'FL'],
                          [pan_cancer, aml, tall, bcp, cll, fl],
                          cancer_datas):
    plot(data[0], data[1], data[2], label, data=preload, gs=gs_i)
