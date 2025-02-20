# %%
import pandas as pd
import numpy as np
import re
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.transforms import Bbox
from matplotlib import patches
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.cluster import hierarchy

# %%
##############################################################
# LOAD ESLs AND JASPAR MOTIF HITS FOR ALL 234 CGIs, INTERSECT.
##############################################################
genes234 = pd.read_csv('data/gene_expression/SupplementaryTable3.csv').iloc[:,0].drop_duplicates()
esls = pd.read_table('data/gene_expression/topGenes234_CpGs.bed',
                     header=None, names=['chr','start','end','cpg', 'gene'])
motif_hits = pd.read_table('data/gene_expression/topGenes234_CGIs_Jaspar_motifs.bed', 
                           header=None, names=['chrom', 'chromStart', 'chromEnd', 'jasparID', 'score', 'strand', 'tf'])
motif_hits['tf'] = motif_hits['tf'].str.upper()

# %%
### for each ESL; for each motif overlap; add to total results df
all_overlaps = pd.DataFrame()
for i in esls.index[:]:
    c, s, e = esls.loc[i, ['chr', 'start', 'end']]
    overlap = motif_hits[(motif_hits['chrom'] == c)
                         & (motif_hits['chromStart'] <= s)
                         & (motif_hits['chromEnd'] > s)].copy()
    overlap['chr'] = c
    overlap['start'] = s
    overlap['end'] = e
    overlap['cpg'] = esls.at[i, 'cpg']
    overlap['gene'] = esls.at[i, 'gene']
    all_overlaps = pd.concat([all_overlaps, overlap], axis=0)

counts = all_overlaps['tf'].value_counts()
all_overlaps

# %%
### ONLY KEEP HITS WITH SCORE > 400 (i.e., P-value < 1e-4)
overlaps_400 = all_overlaps[all_overlaps['score'] >= 400].copy()
cpg_counts = overlaps_400['cpg'].value_counts()
tf_counts = overlaps_400['tf'].value_counts()
gene_counts = overlaps_400['gene'].value_counts()

# %%
### CREATE MOTIFS X GENE MATRIX
genes_motifs_mtx = pd.DataFrame(index=sorted(list(tf_counts.index)), columns=sorted(list(gene_counts.index)))
for tf in tf_counts.index:
    counts_per_gene = overlaps_400[overlaps_400['tf'] == tf]['gene'].value_counts()
    row = pd.Series(index=sorted(list(gene_counts.index)))
    for g in row.index:
        if g in list(counts_per_gene.index):
            row[g] = counts_per_gene[g]
        else:
            row[g] = 0
    genes_motifs_mtx.loc[tf] = list(row)
# binarized
genes_motifs_mtx_bin = (genes_motifs_mtx > 0).astype(int)

# %%
### get TF groupings for annotation bar
tfnames = pd.Series(genes_motifs_mtx_bin.index)
tf_groups = dict.fromkeys(['ATF','CDX', 'CREB','E2F','EGR','ELF','ELK','ERF','ETV',
                           'FOS', 'FOX', 'GLI', 'HES', 'HOX', 'JUN', 'KLF', 'MAF', 'NR',
                           'PAX', 'RAR', 'RFX', 'RXR', 'SP', 'STAT', 'TFAP2', 'TP', 'ZBTB',
                           'ZIC', 'ZNF', 'Other'])
for k in tf_groups: tf_groups[k] = []

for t in tfnames:
    for g in tf_groups:
        if t.startswith(g) and t not in ['NRF1', 'SPDEF', 'SPIC']:
            tf_groups[g].append(t)
            break
        if g == 'Other':
            tf_groups['Other'].append(t)

# reorder matrix rows
ordered_tfs = []
for k in tf_groups: ordered_tfs.extend(tf_groups[k])
genes_motifs_mtx_bin = genes_motifs_mtx_bin.loc[ordered_tfs, :]

yticks = []
vals = list(tf_groups.values())
val_lens = [len(x) for x in vals]
for i in range(len(tf_groups)):
    yticks.append(np.sum(val_lens[:i]) + 0.5*val_lens[i])

colors = [
    '#fb9a99',
    '#b2df8a',
    '#1f78b4',
    '#33a02c',
    '#a6cee3', 
    '#cab2d6',
    '#e31a1c',
    '#fdbf6f',
    '#6a3d9a',
    '#ff7f00',
    '#ffff99',
    '#b15928'
]
row_colors = []
for i, v in enumerate(vals):
    color = colors[i % 12]
    if i == 29:
        color = 'black'
    row_colors.extend( [color] * val_lens[i] )

# %%
##############
# PLOT HEATMAP
##############
mtx = genes_motifs_mtx_bin
cmap = ListedColormap(['lightgray', 'red'])

row_linkage = hierarchy.linkage(
    distance.pdist(mtx, metric='euclidean'), method='ward')
col_linkage = hierarchy.linkage(
    distance.pdist(mtx.T, metric='euclidean'), method='ward')

cpgHM = sns.clustermap(
                    mtx,
                    col_linkage=col_linkage,
                    row_cluster=False,
                    row_colors=row_colors,
                    cbar_pos=None,
                    cmap=cmap,
                    figsize=(12, 25)
                    )

cpgHM.ax_heatmap.axes.tick_params(labelbottom=False, bottom=False,
                                labelleft=False, left=False,
                                labelright=False, right=False,
                                labelsize=19)
cpgHM.ax_col_dendrogram.set_visible(False)
cpgHM.ax_row_dendrogram.set_visible(False)

# Add gridlines to the heatmap
heatmap_ax = cpgHM.ax_heatmap  # Get the heatmap Axes
# Customize gridlines to align with cells
heatmap_ax.set_xticks(np.arange(mtx.shape[1]) , minor=True)
heatmap_ax.set_yticks(np.arange(mtx.shape[0]) , minor=True)
heatmap_ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
heatmap_ax.tick_params(which="minor", top=False, bottom=False, left=False, right=False)

# add marker for RGS2
triangle_vertices = [[0, -3], [1.5, -1], [3, -3]]
triangle = patches.Polygon(triangle_vertices, closed=True, clip_on=False,
                           linewidth=2, edgecolor='red', facecolor='red')
heatmap_ax.add_patch(triangle)
heatmap_ax.text(1.5, -4, 'RGS2', size=18, ha='center', va='bottom')

# add labels for TF groups
cpgHM.ax_row_colors.set_yticks(yticks)
cpgHM.ax_row_colors.set_yticklabels(tf_groups.keys(),
                                    fontsize=18, rotation=0)
cpgHM.ax_row_colors.tick_params(left=False)
cpgHM.ax_row_colors.set_ylabel('TF Groups', fontsize=27)
cpgHM.ax_heatmap.axes.set_xlabel('Genes', fontsize=27)

cpgHM.savefig('plots/figures/Fig_6e_TF_heatmap.png')
mtx.iloc[:, cpgHM.dendrogram_col.reordered_ind].to_csv('plots/source_data/Fig_6e_TF_heatmap.csv')
