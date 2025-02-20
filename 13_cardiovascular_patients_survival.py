# %%
########################################################################
# Conduct Kaplan-Meier and Cox proportional hazards regression analyses
# on the Framingham Heart Study cohort and the cardiogenic shock cohort

# Requires:
# - CS_beta_values.RDS
# Download these files from Zenodo: <URL>

# RESTRICTED ACCESS:
# - FHS_beta_values.RDS
#     - original data at database of Genotypes and Phenotypes (dbGaP):
#           phs000007.v32.p13
########################################################################
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyreadr as pr
import scipy.stats as st
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.compare import compare_survival
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from import_data import stable_um, recurrence_above_5

# %%
###########################################################
# LOAD FRAMINGHAM COHORT AND CARDIOGENIC SHOCK PATIENT DATA
###########################################################
### Framingham cohort
with open('data/Framingham/beta_values.pkl', 'rb') as fin:
    fhs_beta = pickle.load(fin)
fhs_clinical = pd.read_csv('data/Framingham/FHS_clinical.csv')
shared_id = list(set(fhs_beta.columns) & set(fhs_clinical['LABID']))
fhs_clinical = fhs_clinical.loc[fhs_clinical['LABID'].isin(shared_id), :]
fhs_clinical.set_index('LABID', inplace=True)
for col in ['died', 'cvd', 'chd', 'chf']:
    fhs_clinical[col] = fhs_clinical[col].astype(bool)

fhs_dmi = np.std(fhs_beta.loc[recurrence_above_5.index], axis=0); fhs_dmi.name = 'DMI'
fhs_clinical = pd.merge(fhs_clinical, fhs_dmi, left_index=True, right_index=True)
fhs_clinical['high_dmi'] = fhs_clinical['DMI'] > np.median(fhs_clinical['DMI'])

# only keep patients who did not experience any endpoint before sampling
fhs_clinical = fhs_clinical[(fhs_clinical['follow_up'] >= 0)
                            & (fhs_clinical['fu_cvd'] >= 0)
                            & (fhs_clinical['fu_chd'] >= 0)
                            & (fhs_clinical['fu_chf'] >= 0)]

# %%
### Cardiogenic shock patients
cs_beta = pr.read_r('data/CardiogenicShock/beta_values.RDS')[None]
cs_clinical = pd.read_csv('data/CardiogenicShock/SampleSheet.csv')
cs_clinical['id'] = list(cs_clinical['Sentrix_ID'].astype(int).astype(str) + '_' + cs_clinical['Sentrix_Position'])
cs_clinical = cs_clinical[cs_clinical['CS_HF'] == 'Cardiogenic shock']
cs_beta = cs_beta.loc[:,cs_clinical['id']]
cs_clinical.rename({'Seq_ID':'ARCHCVD'}, axis=1, inplace=True)
meta = pd.read_csv('data/CardiogenicShock/CS_clinical.csv')
meta = meta[['ARCHCVD','follow_up','death_1095']]
cs_clinical = pd.merge(cs_clinical, meta, on='ARCHCVD')

cs_clinical.set_index('id', inplace=True)
cs_clinical['died'] = cs_clinical['death_1095'].astype(bool)

cs_dmi = np.std(cs_beta.loc[recurrence_above_5.index], axis=0); cs_dmi.name = 'DMI'
cs_clinical = pd.merge(cs_clinical, cs_dmi, left_index=True, right_index=True)
cs_clinical['high_dmi'] = cs_clinical['DMI'] > np.median(cs_clinical['DMI'])

cs_clinical['DMI-CH'] = 0
for j in cs_clinical.index:
    if cs_clinical.at[j, 'high_dmi'] and cs_clinical.at[j, 'final_CH']:
        cs_clinical.at[j, 'DMI-CH'] = 4
    if cs_clinical.at[j, 'high_dmi'] and not(cs_clinical.at[j, 'final_CH']):
        cs_clinical.at[j, 'DMI-CH'] = 3
    if not(cs_clinical.at[j, 'high_dmi']) and cs_clinical.at[j, 'final_CH']:
        cs_clinical.at[j, 'DMI-CH'] = 2
    if not(cs_clinical.at[j, 'high_dmi']) and not(cs_clinical.at[j, 'final_CH']):
        cs_clinical.at[j, 'DMI-CH'] = 1

# %%
#######################################################
# PERFORM UNIVARIATE AND MULTIVARIATE SURVIVAL ANALYSIS
#######################################################
def kaplan_meier(clinical_data, outcome, follow_up, variable,
                 var_cats, labels, endpoint='Mortality',
                 colors=['#1f78b4', '#e31a1c']):
    mpl.rcParams['legend.handlelength'] = 2
    mpl.rcParams['legend.handleheight'] = 0.7 

    # require that the follow-up time is >0.
    df = clinical_data[clinical_data[follow_up] >= 0]
    
    # CoxPH with one binary variable (get stats for KM curve)
    cph = CoxPHFitter()
    selection = df[[follow_up, outcome, variable]]
    cph.fit(selection, duration_col=follow_up, event_col=outcome)
    hr = cph.summary['exp(coef)'].iloc[0]
    pval = cph.summary['p'].iloc[0]
    
    ### plotting
    fig, ax = plt.subplots(figsize=(6,4))
    ylab = 'Survival Probability'
    title = f'Endpoint: {endpoint}'
    textloc = 0.12
    
    # source data
    data = []
    for j, cat in enumerate(var_cats):
        df2 = df[df[variable] == cat]

        time, survival_prob, conf_int = kaplan_meier_estimator(
            df2[outcome], df2[follow_up], conf_type='log-log')
        time = time / 365
        time = np.concatenate((np.array([0]), time), axis=0)
        survival_prob = np.concatenate((np.array([1]), survival_prob), axis=0)

        if endpoint != 'Mortality':
            survival_prob = 1 - survival_prob
            conf_int[0] = 1 - conf_int[0]
            conf_int[1] = 1 - conf_int[1]
            textloc = 0.87
            ylab = f'Probability of {endpoint}'
        data.append(time)
        data.append(survival_prob)

        ax.set_ylabel(ylab, fontsize=14)
        ax.set_xlabel('Time (years)', fontsize=14)
        ax.tick_params(labelsize=12)
        ax.step(time, survival_prob, where='post', color=colors[j])
        
        if clinical_data.shape[0] != 64:
            ax.fill_between(time[1:], conf_int[0], conf_int[1],
                            alpha=0.25, color=colors[j], step='post')

    strpval = "{:0.3e}".format(pval) if pval < 0.0001 else round(pval, 4)
    hr_pval = f'HR: {round(hr, 2)}\np = {strpval}'

    text_fs = 12
    legend_fs = 12
    ### CARDIOGENIC SHOCK ###
    if clinical_data.shape[0] == 64:
        ax.set_xlim(-0.1,3)
        if variable == 'DMI-CH':
            pval = multivariate_logrank_test(
                cs_clinical['follow_up'], cs_clinical['DMI-CH'], cs_clinical['died']).p_value
            strpval = "{:0.3e}".format(pval) if pval < 0.0001 else round(pval, 4)
            hr_pval = f'p = {strpval}' #\nn = {df.shape[0]}
            ax.text(0.96, 0.9, hr_pval, horizontalalignment='right',
                verticalalignment='center', transform=ax.transAxes,
                fontsize=text_fs)
            mpl.rcParams['legend.handlelength'] = 1 # 2
            mpl.rcParams['legend.handleheight'] = 1 # 0.7
            ax.legend(handles=ax.lines, labels=labels,
                frameon=False, loc=(-0.1, -0.3), fontsize=11,
                ncol=4, columnspacing=1)
        else:
            ax.legend(handles=ax.lines, labels=labels,
                    frameon=False, fontsize=legend_fs)
            ax.text(0.05, textloc, hr_pval, horizontalalignment='left',
                verticalalignment='center', transform=ax.transAxes,
                fontsize=text_fs)
    ### FRAMINGHAM COHORT ###
    else:
        legend_loc = 'lower right' if endpoint in ['CVD','CHD'] else 'upper right'
        ax.legend(handles=ax.lines, labels=labels,
            frameon=False, loc=legend_loc, fontsize=legend_fs)
        ax.text(0.04, textloc, hr_pval, horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes,
            fontsize=text_fs)
        if outcome == 'cvd': ax.set_ylim(-0.0476*0.5,0.45)
        if outcome == 'chd': ax.set_ylim(-0.0476*0.15,0.15)
        if outcome == 'chf': ax.set_ylim(-0.0476*0.28,0.28)

    # return source data
    df = pd.DataFrame()
    for j, label in enumerate(labels):
        cat_df = pd.DataFrame({'variable': label,
                               'time': data[j*2],
                               'probability': data[1 + j*2]})
        df = pd.concat([df, cat_df])
    return fig, df

def coxph_multivariate(clinical_data, outcome, follow_up, covariates):
    summaries = []
    for i in range(len(covariates)):
        cph = CoxPHFitter()
        selection = clinical_data[[follow_up, outcome] + covariates[:i+1]]
        selection = selection[selection[follow_up] >= 0]
        cph.fit(selection, duration_col=follow_up, event_col=outcome)
        summaries.append(cph.summary)
    return summaries

def forest_plot(results):
    results = results.loc[results.index[::-1]]
    fig, (ax, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 0.7], 'wspace':0.1},
                                  figsize=(5,4))
    
    hr = results['exp(coef)']
    left_error = hr - results['exp(coef) lower 95%']
    right_error = results['exp(coef) upper 95%'] - hr

    ax.errorbar(hr, range(results.shape[0]), xerr=[left_error, right_error],
                fmt='o', capsize=5, color='black')
    ax.set_yticks(range(results.shape[0]))
    ax.set_yticklabels(results.index, fontsize=12)
    ax.tick_params(labelsize=12)
    ax.axvline(x=1, color='gray', linestyle='--')
    ax.set_xlabel('Hazard Ratio', fontsize=14)
    ax.set_xlim(0, ax.get_xlim()[1])

    # # results table
    ax2.set_ylim(ax.get_ylim())
    fs = 12
    for i, label in enumerate(results.index):
        x1, x2, x3 = 0, 3, 8.2
        headers_y = len(results)-1 + (len(results)-1)*0.07
        ax2.text(x1, headers_y, 'HR', fontsize=fs, weight='bold')
        ax2.text(x1, i, str(round(results.at[label, 'exp(coef)'], 2)),
                fontsize=fs, va='center')
        # CI
        ax2.text(x2, headers_y, '95% CI', fontsize=fs, weight='bold')
        ci = str(round(results.at[label, 'exp(coef) lower 95%'], 2)) + \
            '-' + str(round(results.at[label, 'exp(coef) upper 95%'], 2))
        ax2.text(x2, i, ci, fontsize=fs, va='center')
        # p value
        pval = results.at[label, 'p']
        strpval = "{:0.2e}".format(pval) if pval < 0.001 else round(pval, 3)
        ax2.text(x3, headers_y, 'P value', fontsize=fs, weight='bold')
        ax2.text(x3, i, strpval, fontsize=fs, va='center')

    ax2.set_xlim(0, 8)
    ax2.axis('off')
    return fig

# %%
### Cardiogenic shock
fig1, cardio_km_dmi = kaplan_meier(cs_clinical, 'died', 'follow_up', 'high_dmi',
             var_cats=[False, True], labels=['Low DMI', 'High DMI'])
fig2, cardio_km_ch = kaplan_meier(cs_clinical, 'died', 'follow_up', 'final_CH',
             var_cats=[False, True], labels=['CH-', 'CH+'])
fig3, cardio_km_dmi_ch = kaplan_meier(cs_clinical, 'died', 'follow_up', 'DMI-CH',
             var_cats=[4,3,2,1], labels=['DMI+, CH+', 'DMI+, CH-', 'DMI-, CH+', 'DMI-, CH-'],
             colors=['#e31a1c', '#ff7f00', '#33a02c', '#1f78b4'])

fig1.savefig('plots/figures/Fig_5g_cardiogenic_shock_DMI.png')
fig2.savefig('plots/figures/SuppFig_7_cardiogenic_shock_CH.png')
fig3.savefig('plots/figures/Fig_5h_cardiogenic_shock_DMI-CH.png')

cardio_km_dmi.to_csv('plots/source_data/Fig_5g_cardiogenic_shock_DMI.csv', index=False)
cardio_km_ch.to_csv('plots/source_data/SuppFig_7_cardiogenic_shock_CH.csv', index=False)
cardio_km_dmi_ch.to_csv('plots/source_data/Fig_5h_cardiogenic_shock_DMI-CH.csv', index=False)

# %%
s1 = coxph_multivariate(cs_clinical, 'died', 'follow_up',
                   covariates = ['high_dmi', 'age', 'final_CH'])
mortality = s1[-1]
mortality.index = ['DMI', 'Age', 'CH+']
forest_plot(mortality)

# %%
### Framingham
fig1, framingham_km_died = kaplan_meier(fhs_clinical, 'died', 'follow_up', 'high_dmi',
             var_cats=[False, True], labels=['Low DMI', 'High DMI'])
fig2, framingham_km_cvd = kaplan_meier(fhs_clinical, 'cvd', 'fu_cvd', 'high_dmi', endpoint='CVD',
             var_cats=[False, True], labels=['Low DMI', 'High DMI'])
fig3, framingham_km_chd = kaplan_meier(fhs_clinical, 'chd', 'fu_chd', 'high_dmi', endpoint='CHD',
             var_cats=[False, True], labels=['Low DMI', 'High DMI'])
fig4, framingham_km_chf = kaplan_meier(fhs_clinical, 'chf', 'fu_chf', 'high_dmi', endpoint='CHF',
             var_cats=[False, True], labels=['Low DMI', 'High DMI'])

fig1.savefig('plots/figures/Fig_5b_framingham_mortality.png')
fig2.savefig('plots/figures/Fig_5c_framingham_CVD.png')
fig3.savefig('plots/figures/Fig_5d_framingham_CHD.png')
fig4.savefig('plots/figures/Fig_5e_framingham_CHF.png')

framingham_km_died.to_csv('plots/source_data/Fig_5b_framingham_mortality.csv', index=False)
framingham_km_cvd.to_csv('plots/source_data/Fig_5c_framingham_CVD.csv', index=False)
framingham_km_chd.to_csv('plots/source_data/Fig_5d_framingham_CHD.csv', index=False)
framingham_km_chf.to_csv('plots/source_data/Fig_5e_framingham_CHF.csv', index=False)

# %%
s1 = coxph_multivariate(fhs_clinical, 'died', 'follow_up',
                   covariates = ['high_dmi', 'age8', 'sex'])
s2 = coxph_multivariate(fhs_clinical, 'cvd', 'fu_cvd',
                   covariates = ['high_dmi', 'age8', 'sex'])
s3 = coxph_multivariate(fhs_clinical, 'chd', 'fu_chd',
                   covariates = ['high_dmi', 'age8', 'sex'])
s4 = coxph_multivariate(fhs_clinical, 'chf', 'fu_chf',
                   covariates = ['high_dmi', 'age8', 'sex'])

age_sex_corrected = pd.concat([s1[-1].iloc[0], s2[-1].iloc[0],
                               s3[-1].iloc[0], s4[-1].iloc[0]], axis = 1)
age_sex_corrected.columns = ['Mortality', 'CVD', 'CHD', 'CHF']
age_sex_corrected = age_sex_corrected.T
fig = forest_plot(age_sex_corrected)
fig.savefig('plots/figures/Fig_5f_CoxPH_HRs.png')
