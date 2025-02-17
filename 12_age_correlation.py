# %%
#############################################################################
# Calculate correlation between DMI and age in cohorts of healthy individuals

# Requires:
# - GSE87571_beta_values.RDS
# - GSE115278_beta_values.RDS
# - GSE40279_beta_values.RDS
# - GSE197676_beta_values.RDS
# Download these files from Zenodo: <URL>
########################################################################
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyreadr as pr
from scipy import stats as st
from import_data import stable_um, recurrence_above_5

# %%
#############################
# LOAD HEALTHY AGING DATASETS
#############################
gse87571 = pr.read_r('data/GSE87571/beta_values.RDS')[None]
ages87571 = pd.read_table('data/GSE87571/GSE87571_ages.txt')

gse115278 = pr.read_r('data/GSE115278/beta_values.RDS')[None]
ages115278 = pd.read_table('data/GSE115278/info_GSE115278-GPL16304.txt')
ages115278['age'] = [float(x.split(': ')[1]) for x in ages115278['age']]
ages115278 = ages115278[['geo','age']]

gse40279 = pr.read_r('data/GSE40279/beta_values.RDS')[None]
ages40279 = pd.read_table('data/GSE40279/info_GSE40279.txt')
ages40279.set_index('geo', inplace=True)
ages40279['age'] = [int(x.split(': ')[-1]) for x in ages40279['age']]

gse197676 = pr.read_r('data/GSE197676/beta_values.RDS')[None]
ages197676 = pd.read_csv('data/GSE197676/GSE197676_ages.csv')

for df in [gse87571, gse115278, gse197676]:
    df.columns = [x.split('_')[0] for x in df.columns]

ages87571.set_index('geo', inplace=True)
ages87571 = ages87571.dropna()
ages115278.set_index('geo', inplace=True)
ages197676.set_index('geo', inplace=True)

# %%
###################################################################
# STRATIFY BY AGE, CALCULATE AGE VS. DMI CORRELATION, AND PLOT DATA
###################################################################
def stratify_and_boxplot(df, age_col, dmi_col, bins, bin_labels, figsize=(6, 4)):
    """
    Parameters:
    - df (pd.DataFrame): The input DataFrame, indexed by IDs, with columns for age and DMI.
    - age_col (str): The name of the age column.
    - dmi_col (str): The name of the DMI column.
    - bins (list): The bin edges for stratification.
    - bin_labels (list): Labels for the age categories.
    - figsize (tuple): Size of the figure (width, height).
    - title (str): Title for the boxplot.
    """
    # Stratify the ages into discrete categories
    df['Age Category'] = pd.cut(df[age_col], bins=bins, labels=bin_labels, include_lowest=True)
    
    # Calculate means, counts, and standard errors for each category
    grouped = df.groupby('Age Category')[dmi_col]
    means = grouped.mean()
    counts = grouped.count()
    sems = grouped.apply(st.sem)  # Standard error of the mean

    # Calculate 95% confidence intervals
    ci_lower = means - 1.96*sems
    ci_upper = means + 1.96*sems
    
    # Map bin labels to bin midpoints
    bin_midpoints = means.index

    # Ensure alignment and remove NaNs
    # valid_mask = ~means.isna()  # Mask for non-NaN values
    valid_mask = counts > 1  # Mask for 
    midpoints = bin_midpoints[valid_mask]
    means = means[valid_mask]
    sems = sems[valid_mask]
    ci_lower = ci_lower[valid_mask]
    ci_upper = ci_upper[valid_mask]


    # Linear fit
    n = len(means)
    x = np.array(means.index)
    y = np.array(means.values)
    slope, yint, r, p = st.linregress(x, y)[:4]
    strpval = "{:0.3e}".format(p) if p < 0.0001 else round(p, 4)
    # Means of x and y
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    # Slope (b1) and intercept (b0)
    b1 = np.sum((x - x_mean) * (y - y_mean)) / np.sum((x - x_mean) ** 2)
    b0 = y_mean - b1 * x_mean
    # Predicted values
    y_pred = b0 + b1 * x
    # Residuals and residual standard error
    residuals = y - y_pred
    s_err = np.sqrt(np.sum(residuals**2) / (n - 2))
    # Critical t-value for the given confidence level
    confidence = 0.95
    t_value = t.ppf((1 + confidence) / 2, df=n - 2)
    # Confidence intervals
    x_var = np.sum((x - x_mean) ** 2)
    ci = t_value * s_err * np.sqrt(1 / n + (x - x_mean) ** 2 / x_var)
    line_ci_upper = y_pred + ci
    line_ci_lower = y_pred - ci


    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(means.index, means.values, marker='.', linestyle='-', color='black')
    
    # Plot the linear fit
    ax.plot(means.index, y_pred, linestyle='-', color='red')
    ax.fill_between(means.index, line_ci_upper, line_ci_lower, alpha=0.2)

    # plot the error bars for the means
    for i in range(means.shape[0]):
        ax.plot([list(means.index)[i], list(means.index)[i]],
                [ci_lower.iloc[i], ci_upper.iloc[i]],
                linewidth=1, color='black')

    ax.set_xlabel("Age", fontsize=14)
    ax.set_ylabel("DNA Methylation Instability", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.set_xlim(np.min(means.index)-1, np.max(means.index)+1)
    ax.text(0.05, 0.93, f'r = {round(r, 3)}\np = {strpval}',
            transform=ax.transAxes, va='top', ha='left',
            fontsize=13)

    return means

# %%
### plot age correlation for different healthy datasets
mygse = gse40279
myages = ages40279

# mygse = gse87571
# myages = ages87571

# mygse = gse115278
# myages = ages115278

# mygse = gse197676
# myages = ages197676

cpgs = list(set(mygse.index) & set(recurrence_above_5.index))
dmi = np.std(mygse.loc[cpgs], axis=0)
dmi.name = 'DMI'
df = pd.merge(myages, dmi, right_index=True, left_index=True)
labels = [i for i in range(0,95,5)]
mean_DMI = stratify_and_boxplot(df, 'age', 'DMI', labels + [labels[-1]+5], 
                        bin_labels=labels)
# source data
mean_DMI.to_csv('plots/source_data/Fig_5a')
# mean_DMI.to_csv('plots/source_data/ExtFig_5a')
# mean_DMI.to_csv('plots/source_data/ExtFig_5b')
# mean_DMI.to_csv('plots/source_data/ExtFig_5c')