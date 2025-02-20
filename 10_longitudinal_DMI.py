# %%
########################################################################
# Examine DMI levels in AML patients who were sampled at different times
# throughout their disease (diagnosis, remissions, relapses)

# Requires:
# - AML_BM_beta_values.RDS
# Download these files from Zenodo: <URL>
########################################################################
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import pyreadr as pr
from import_data import stable_um, recurrence_above_5

# %%
################################################# 
# LOAD LONGITUDINAL AML SAMPLES AND CLINICAL DATA
#################################################
beta_values = pr.read_r('data/longitudinal_AML/BM_beta_values.RDS')[None]
clinical_data = pd.read_csv('data/longitudinal_AML/BM_clinical_data.csv')

patients = {0: clinical_data.iloc[:4,:].copy(),
            1: clinical_data.iloc[8:12,:].copy(),
            2: clinical_data.iloc[4:8,:].copy(),
            3: clinical_data.iloc[12:16,:].copy()}

# %%
############################
# PLOT DMI ACROSS TIMEPOINTS
############################
def longitudinal_dmi(patient, beta_values):
    dmi_values = []
    for i in patient.index:
        sample = patient.loc[i, 'barcode']
        # only use sites above 5% recurrence to calculate DMI
        dmi = np.std(beta_values.loc[recurrence_above_5.index, sample])
        dmi_values.append(dmi)
    days_since_dx = list(patient['time from dx'])
    disease_status = list(patient['Diagnosis_remission_relapse'])
    return days_since_dx, dmi_values, disease_status

## save source data
data = pd.DataFrame(columns=['patient', 'days_since_Dx', 'DMI', 'status'])

### plotting
fig, axs = plt.subplots(1, 4, figsize=(22, 3))
for i, ax in enumerate(axs.flatten()):
    days_since_dx, dmi_values, disease_status = longitudinal_dmi(patients[i], beta_values)
    for j in range(len(dmi_values)):
        data.loc[data.shape[0]] = (f'Patient {i+1}', days_since_dx[j], dmi_values[j], disease_status[j])

    ax.plot(days_since_dx, dmi_values, color='black')
    red_idx = [j for j in range(4) if disease_status[j] != 'remission']
    blue_idx = [j for j in range(4) if disease_status[j] == 'remission']
    for idx, color in zip([red_idx, blue_idx], ['red','blue']):
        ax.scatter([days_since_dx[j] for j in idx],
                [dmi_values[j] for j in idx],
                color=color, marker='o', s=100, zorder=2)

    handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=12, label='Diagnosis/\nRelapse'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=12, label='Remission')]
    if i == 0:
        ax.set_ylabel('DNA Methylation Instability', fontsize=12)
        ax.legend(handles=handles, frameon=False, fontsize=12)
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Days since diagnosis', fontsize=13)
    ax.set_title(f'Patient {i+1}', fontsize=14)

fig.tight_layout()
fig.savefig('plots/figures/Fig_4d_longitudinal_DMI.png')
data.to_csv('plots/source_data/Fig_4d_longitudinal_DMI.csv', index=False)
