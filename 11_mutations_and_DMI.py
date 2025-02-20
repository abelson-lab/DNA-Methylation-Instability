# %%
########################################################################
# Examine DMI levels and mutation VAFs in AML patients who were sampled
# at diagnosis and post-induction therapy

# Requires:
# - AML_PB_beta_values.RDS
# Download these files from Zenodo: <URL>
########################################################################
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import pyreadr as pr
from import_data import stable_um, recurrence_above_5

# %%
############################################### 
# LOAD METHYLATION, MUTATION, AND CLINICAL DATA
###############################################
beta_values = pr.read_r('data/longitudinal_AML/PB_beta_values.RDS')[None]
clinical_data = pd.read_table('data/longitudinal_AML/PB_sample_info.txt')
muts = pd.read_csv('data/longitudinal_AML/PB_all_mutations.csv')
muts = muts[muts['Reported at diagnosis'] == 'YES']
mut_pts = list(muts['altID'].drop_duplicates()) + ['MH_6745']

patients = {}
for sample in pd.unique(clinical_data['LTB patient ID']):
    patients[sample] = clinical_data[clinical_data['LTB patient ID'] == sample]

# %%
##############################################
# PLOT DMI AND MUTATION VAFs ACROSS TIMEPOINTS
##############################################
def get_dmi_and_vafs(patient):
    dmi_values = []
    for i in patient.index[:3]:
        sample = patient.loc[i, 'sample_ID']
        # only use sites above 5% recurrence to calculate DMI
        dmi = np.std(beta_values.loc[recurrence_above_5.index, sample])
        dmi_values.append(dmi)
    
    patient_id = patient['LTB patient ID'].iloc[0]
    result_df = pd.DataFrame()
    result_df['Patient'] = [patient_id] * 3
    result_df['Time Point'] = ['Diagnosis', 'T2', 'T3']
    result_df['DMI'] = dmi_values
    
    mutations = muts[muts['altID'] == patient_id]
    for pos in pd.unique(mutations['Position']):
        vafs = []
        mutation_i = mutations[mutations['Position'] == pos]
        mut_name = mutation_i['Gene ID'].iloc[0] + '_' + str(pos)
        for time in ['Dx', 'T1', 'T2']:
            if time in list(mutation_i['Detected at time point']):
                vaf = mutation_i.loc[mutation_i['Detected at time point'] == time, 'VAF'].iloc[0]
                vafs.append(vaf)
            else:
                vafs.append(pd.NA)
        result_df[mut_name] = vafs
    
    return result_df

# %%
### plotting
# fname, pts_to_plot = 'Fig_4e_mutations_and_DMI', ['2522_86899', 'MH_6726', '2522_85368', 'MH_6745']
fname, pts_to_plot = 'ExtFig_4_mutations_and_DMI', ['2522_86891', '2522_86210', 'MH_6729', 'MH_6747', '2522_85291', '2522_87221']

fig = plt.figure(figsize=(18, 10))
outer = gridspec.GridSpec(2, 3, figure=fig, hspace=0.2, wspace=0.25)

gs = []
for i in range(6):
    gs.append(gridspec.GridSpecFromSubplotSpec(2, 1,
                                               subplot_spec = outer[i],
                                               hspace=0.05))
for num, pt in enumerate(pts_to_plot):
    cell1 = gs[num][0]
    cell2 = gs[num][1]
    ax1 = plt.subplot(cell1)
    ax2 = plt.subplot(cell2)

    df = get_dmi_and_vafs(patients[pt])
    ax1.plot(df['Time Point'], df['DMI'], color='blue', marker='o')
    ax1.tick_params(labelsize=14)
    ax1.xaxis.set_visible(False)
    ax1.set_ylabel('DMI', fontsize=16)
    

    mut_vafs = df.iloc[:, 3:]
    handles = []
    if mut_vafs.shape[1]:
        nd = np.max(mut_vafs) * -0.06
        ax2.text(-0.07, nd*(-0.3), 'Detection Limit', fontsize=13, ha='left', va='bottom')

        for col in mut_vafs:
            vafs = mut_vafs[col].fillna(nd)
            l = ax2.plot(df['Time Point'], vafs, linestyle='solid', marker='o')
            label = col.split('_')[0]
            handles.append(
                Line2D([0], [0], marker='o', color=l[0]._color, markerfacecolor=l[0]._color, markersize=5, label=label)
            )
        ax2.axhline(y=0, linestyle='dashed', color='black')
        ax2.legend(handles=handles, frameon=False, fontsize=14)
    else:
        ax2.text(0.5, 0.5, 'No mutations detected', transform=ax2.transAxes,
                 fontsize=18, ha='center', va='center')

    ax2.set_ylabel('VAF', fontsize=16)
    ax2.tick_params(labelsize=14)
    ax2.set_xlim(-0.1, 2.1)
    ax2.set_xticks([0,1,2])
    ax2.set_xticklabels(['Diagnosis','T2','T3'], fontsize=15)

    renamed_pts = {'2522_86899': 5,
                   'MH_6726': 6,
                   '2522_85368': 7,
                   'MH_6745': 8,
                   '2522_86891': 9,
                   '2522_86210': 10,
                   'MH_6729': 11,
                   'MH_6747': 12,
                   '2522_85291': 13,
                   '2522_87221': 14}

    ax1.set_title(f'Patient {renamed_pts[pt]}', fontsize=16)

fig.savefig(f'plots/figures/{fname}.png')
# save source data
with open(fname, 'w') as f:
    pass  # Creates an empty file
for pt in pts_to_plot:
    df = get_dmi_and_vafs(patients[pt])
    df.to_csv(f'plots/source_data/{fname}.csv', mode='a', index=False)
