# %%
#############################################################################
# Calculate average beta value at ESLs in various cancer and control cohorts.
#
# Requires:
# - GSE124413_healthy_BM_beta_values.RDS
# - GSE141682_beta_values.RDS
# - GSE132203_beta_values.RDS
# - GSE124413_AML_beta_values.RDS
# - TCGA-LAML_beta_values.RDS
# - GSE49031_T-ALL_beta_values.RDS
# - GSE147667_beta_values.RDS
# - GSE49031_BCP-ALL_beta_values.RDS
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
# - EGA_CLL_165_beta_values.RDS
#     - original data at European Genome-Phenome Archive: EGAD00010000254
#############################################################################
import numpy as np
import pandas as pd
import pyreadr
import copy
import pickle
from import_data import stable_sites, stable_um, stable_m

# %%
########################################################
# LOAD HEALTHY CONTROL AND HEMATOLOGICAL CANCER DATASETS
########################################################
# Healthy controls
unknown = pyreadr.read_r('data/GSE124413/healthyBM/beta_values.RDS')[None]
chinese = pyreadr.read_r('data/GSE141682-Han/GSE141682_RAW/beta_values.RDS')[None]
african = pyreadr.read_r('data/GSE132203/beta_values.RDS')[None]

# Acute myeloid leukemia (pediatric)
aml1 = pyreadr.read_r('data/GSE124413/beta_values.RDS')[None]
# AML adult
aml2 = pyreadr.read_r('data/TCGA-LAML_idat/files/beta_values.RDS')[None]
# T-cell acute lymphoblastic leukemia (pediatric)
tall1 = pyreadr.read_r('data/GSE49031/T-ALL_beta_values.RDS')[None]
# T-ALL adult
tall2 = pyreadr.read_r('data/GSE147667/beta_values.RDS')[None]
# B-cell precursor acute lymphoblastic leukemia 
bcpall = pyreadr.read_r('data/GSE49031/BCP-ALL_beta_values.RDS')[None]
# Chronic lymphocytic leukemia
cll = pyreadr.read_r('data/EGA-CLL490/beta_values.RDS')[None]
cll_sex = pd.read_table('data/EGA-CLL490/CLL4_EGAD00010001975_predicted_sex.txt')
males = cll_sex.loc[cll_sex['sex'] == 'M', 'sentrix']
females = cll_sex.loc[cll_sex['sex'] == 'F', 'sentrix']
cll_m = cll.loc[:,males]
cll_f = cll.loc[:,females]
# Follicular lymphoma
fl = pyreadr.read_r('data/02_FL_Samples/fl/beta_values.RDS')[None]
# Primary plasma cell leukemia
ppcl = pyreadr.read_r('data/GSE104770/beta_values.RDS')[None]
# Myeloid / natural killer precursor leukemia
mnkpl = pyreadr.read_r('data/GSE197696/beta_values.RDS')[None]
# Mixed phenotype acute leukemia
mpal = pyreadr.read_r('data/GSE113545/beta_values.RDS')[None]
### MDS / CML
bvals = pyreadr.read_r('data/MDSCML/Methylation_idat/beta_values.RDS')[None]
anno = pd.read_csv('data/MDSCML/160629_JNU CSH_270sample sheet.csv', skiprows=7)
anno['id'] = anno['Sentrix_ID'].astype(str) +'_'+ anno['Sentrix_Position']
mds = anno[anno['Sample_Name'].str.contains('MDS')]
cml = anno[anno['Sample_Name'].str.contains('CML')]
# Myelodysplastic syndrome
mds = bvals.loc[:,mds['id'].iloc[:int(len(mds['id'])/2)]]
# Chronic myeloid leukemia
cml = bvals.loc[:,cml['id'].iloc[:int(len(cml['id'])/2)]]

not_chronic_phase = ['CML Dx-16', 'CML Dx-17', 'CML Dx-22', 'CML Dx-25', 'CML Dx-41',
                     'CML Dx-46', 'CML Dx-47', 'CML Dx-63', 'CML Dx-85']
not_chronic_phase = anno[anno['Sample_Name'].isin(not_chronic_phase)]['id']
cml = cml[list(set(cml.columns) - set(not_chronic_phase))]

# %%
datasets_og = {0: ['Control 1\nn=41', unknown, '#9ecae1', 'Control 1 (Unknown ethnicity)'],
            1: ['Control 2\nn=795', african, '#4292c6', 'Control 2 (African descent)'],
            2: ['Control 3\nn=42', chinese, '#084594', 'Control 3 (Chinese descent)'],
            3: ['CML\nn=69', cml, '#bbfc49', 'Chronic Myeloid Leukemia'],
            3.01: ['MDS\nn=57', mds, '#ffdf6b', 'Myelodysplastic Syndrome'],
            3.1: ['AML Ped.\nn=459', aml1, '#fc9272', 'Acute Myeloid Leukemia (Pediatric)'],
            3.2: ['AML Ad.\nn=194', aml2, '#de2d26', 'Acute Myeloid Leukemia (Adult)'],
            4: ['CLL (M)\nn=288', cll_m, '#756bb1', 'Chronic Lymphocytic Leukemia (Males)'],
            4.1: ['CLL (F)\nn=202', cll_f, '#9e9ac8', 'Chronic Lymphocytic Leukemia (Females)'],
            5: ['pPCL\nn=14', ppcl, '#2b8cbe', 'Primary Plasma Cell Leukemia'],
            6: ['FL\nn=246', fl, '#4eb3d3', 'B-cell Follicular Lymphoma'],
            7: ['MPAL\nn=31', mpal, '#7bccc4', 'Mixed Phenotype Acute Leukemia'],
            8: ['BCP-ALL\nn=663', bcpall, '#bae4b3', 'B-cell Precursor Acute Lymphoblastic Leukemia'],
            9: ['T-ALL Ped.\nn=101', tall1, '#74c476', 'T-cell Acute Lymphoblastic Leukemia (Pediatric)'],
            9.1: ['T-ALL Ad.\nn=143', tall2, '#31a354', 'T-cell Acute Lymphoblastic Leukemia (Adult)'],
            10: ['MNKPL\nn=7', mnkpl, '#006d2c', 'Myeloid/Natural Killer Precursor Acute Leukemia'],
           }

# %%
################################################
# Calculate mean ESL beta values in each dataset
################################################
def get_mean_beta(datasets_og, esl_group): # 0=unmethylated, 1=methylated
    datasets = copy.deepcopy(datasets_og)

    # convert datasets to mean beta value of all ESLs
    probes = stable_sites['ALL'][esl_group]
    limit = stable_sites['ALL'][esl_group+2]
    for i in datasets:
        datasets[i][1] = np.mean(datasets[i][1].loc[probes,:], axis=1)
    
    # save mean beta values
    group = []
    label = []
    for i in datasets:
        group.extend( [datasets[i][0]]*len(probes) )
        label.extend( [datasets[i][3]]*len(probes) )
    beta = pd.concat( [datasets[i][1] for i in datasets] )
    cpg = list(beta.index)
    df = pd.DataFrame({'group':group, 'cpg':cpg, 'beta':beta, 'label':label})
    df.to_csv(f"plotting/inputs/ESLs_in_cancer_{['U','M'][esl_group]}.csv", index=False)

get_mean_beta(datasets_og, 0) # unmethylated
get_mean_beta(datasets_og, 1) # methylated

# %%
##############################################
# GET MEDIAN # of DESTABILIZED ESLs PER COHORT
##############################################
def get_median_destabilized(datasets, esl_group):
    results = pd.DataFrame(columns=['median_destabilized'])
    for k, v in datasets.items():
        if esl_group == 0:
            df = v[1].loc[stable_um, :]
            df = df > stable_sites['ALL'][2]
        else:
            df = v[1].loc[stable_m, :]
            df = df < stable_sites['ALL'][3]
        df = np.sum(df, axis=0)
        results.loc[v[0]] = np.median(df)
    return results

get_median_destabilized(datasets_og, 0).to_csv(f"plotting/inputs/ESLs_in_cancer_median_destabilized_U.csv")
get_median_destabilized(datasets_og, 1).to_csv(f"plotting/inputs/ESLs_in_cancer_median_destabilized_M.csv")

# %%
### GET PERTURBED ESLs, GET BETA VALUE
# adult AML
aml_beta = aml2.loc[stable_um].copy()
aml_values = aml_beta[aml_beta > stable_sites['ALL'][2]].values.flatten()
aml_values = aml_values[~np.isnan(aml_values)]
np.mean(aml_values)

# adult T-ALL
tall_beta = tall2.loc[stable_um].copy()
tall_values = tall_beta[tall_beta > stable_sites['ALL'][2]].values.flatten()
tall_values = tall_values[~np.isnan(tall_values)]
np.mean(tall_values)
