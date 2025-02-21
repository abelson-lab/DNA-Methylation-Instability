library(tidyr)
library(ggplot2)
library(ggpubr)
library(glue)
library(patchwork)
library(data.table)
library(pheatmap)
library(dplyr)

################################# SuppFig_1_GSE105018_estimated_cell_counts.png ###############################
### ESTIMATE CELL COUNTS IN DISCOVERY COHORT
df = read.csv('source_data/SuppFig_1_GSE105018_estimated_cell_counts.csv')
df = df[,2:(dim(df)[2] - 1)] %>%
  pivot_longer(cols = everything(), names_to = 'category', values_to = 'values')

bplot = ggplot(df, aes(x = category, y = values)) +
  theme_classic() +
  geom_boxplot(fill = '#33a02c') + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                         legend.position = 'none',
                         axis.title.x = element_text(size=12, face='bold'),
                         axis.title.y = element_text(size=12, face='bold')
                         ) + 
  ylab('Relative Proportion') + xlab('Cell Type') +
  scale_x_discrete(labels = c('B', 'CD4 T', 'CD8 T', 'Granulocyte', 'Monocyte', 'NK'))
bplot
ggsave('figures/SuppFig_1_GSE105018_estimated_cell_counts.png', bplot, width = 4, height = 5)



################################# Fig_1a_purified_cells_beta_filtered.png #####################################
################################# SuppFig_2_purified_cells_beta_unfiltered.png ################################
### ESLs IN PURIFIED BLOOD CELL POPULATIONS
# Unfiltered (SuppFig_2)
# df = read.csv('source_data/SuppFig_2_purified_cells_beta_U_unfiltered.csv')
# df2 = read.csv('source_data/SuppFig_2_purified_cells_beta_M_unfiltered.csv')
# Filtered (Fig_1a)
df = read.csv('source_data/Fig_1a_purified_cells_beta_U_filtered.csv')
df2 = read.csv('source_data/Fig_1a_purified_cells_beta_M_filtered.csv')

df = df[,2:(dim(df)[2])] %>%
  pivot_longer(cols = everything(), names_to = 'category', values_to = 'values')
df2 = df2[,2:(dim(df2)[2])] %>%
  pivot_longer(cols = everything(), names_to = 'category', values_to = 'values')

bplot = ggplot(df, aes(x = category, y = values)) + theme_classic() +
  geom_boxplot(fill = '#1f78b4') + theme(
    axis.text.y = element_text(size=12), #14
    axis.text.x = element_text(size=12, angle = 45, hjust = 1), # 14 30
    legend.position = 'none',
    plot.title = element_text(size=15, face='bold'),
    axis.title.x = element_text(size=15, face='bold'),
    axis.title.y = element_text(size=15, face='bold')
  ) + 
  ylab('β-value') + xlab('Cell Type') +
  scale_x_discrete(labels = c('B', 'CD4 T', 'CD8 T', 'Granulocyte', 'Monocyte')) +
  ggtitle(glue("Unmethylated (n={dim(df[df$category == 'B',])[1]})"))

bplot2 = ggplot(df2, aes(x = category, y = values)) + theme_classic() +
  geom_boxplot(fill = '#e31a1c') + theme(
    axis.text.y = element_text(size=12), #14
    axis.text.x = element_text(size=12, angle = 45, hjust = 1), # 14 30
    legend.position = 'none',
    plot.title = element_text(size=15, face='bold'),
    axis.title.x = element_text(size=15, face='bold'),
    axis.title.y = element_text(size=15, face='bold')
  ) + 
  ylab(NULL) + xlab('Cell Type') +
  scale_x_discrete(labels = c('B', 'CD4 T', 'CD8 T', 'Granulocyte', 'Monocyte')) +
  ggtitle(glue("Methylated (n={dim(df2[df2$category == 'B',])[1]})"))
          
fig = ggarrange(bplot, bplot2, nrow=1, ncol=2)
fig

# ggsave('figures/SuppFig_2_purified_cells_beta_unfiltered.png', fig, width = 7, height = 5)
ggsave('figures/Fig_1a_purified_cells_beta_filtered.png', fig, width = 8, height = 4.5)



################################# Fig_1b_GSE105018_mean_ESL_beta.png ########################################
### ESL DENSITY PLOT IN DISCOVERY COHORT
df = read.csv('source_data/Fig_1b_GSE105018_mean_ESL_beta.csv')
names(df) = c('cpg','value')

# Compute density
density_data <- density(df$value)
# Convert to data frame for ggplot
density_df <- data.frame(x = density_data$x, y = density_data$y)
# Split the data into two halves: left and right of the mean (or any threshold)
threshold <- 0.5
density_df$color <- ifelse(density_df$x < threshold, "#1f78b4", "#e31a1c")
# Plot with different colors for each half
p = ggplot(density_df, aes(x = x, y = y)) + theme_classic() +
  geom_area(aes(fill = color), alpha = 0.5) +  # Use fill based on the condition
  geom_line(color = "black", size = 0.5) +
  scale_fill_identity() +  # Use the specified colors (red, blue)
  theme(axis.title.x = element_text(size=15, face='bold'),
        axis.title.y = element_text(size=15, face='bold'),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  labs(x = "β-value", y = "Density")
p

ggsave('figures/Fig_1b_GSE105018_mean_ESL_beta.png', p, width = 7, height = 4)



################################# Fig_1c_methylation_atlas_beta.png ################################
### METHYL-ATLAS HUMAN CELL TYPES VIOLIN PLOT
df = read.csv('source_data/Fig_1c_methylation_atlas_beta_U.csv')
df2 = read.csv('source_data/Fig_1c_methylation_atlas_beta_M.csv')

df = df[,2:(dim(df)[2])] %>%
  pivot_longer(cols = everything(), names_to = 'category', values_to = 'values')
df2 = df2[,2:(dim(df2)[2])] %>%
  pivot_longer(cols = everything(), names_to = 'category', values_to = 'values')

violin_plot = ggplot(df2, aes(x = category, y = values)) + theme_classic() +
  geom_violin(fill = "#e31a1c", color='#e31a1c', width = 0.8, adjust=0.2, bw=0.2) +
  theme(axis.title.y = element_text(size=14, face='bold'),
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "β-value")

violin_plot2 = ggplot(df, aes(x = category, y = values)) + theme_classic() +
  geom_violin(fill = "#1f78b4", color='#1f78b4', width = 0.8, adjust=0.2, bw=0.2) +
  theme(axis.title.x = element_text(size=15, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        axis.text.x = element_text(size=11, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12)) +
  labs(x = "Cell Type", y = "β-value")
fig = violin_plot / violin_plot2

ggsave('figures/Fig_1c_methylation_atlas_beta.png', fig, width = 15, height = 6)



################################# Fig_2_ESLs_in_cancer_U.png ####################################
################################# ExtFig_1_ESLs_in_cancer_U.png #################################
### MEAN BETA VALUES IN DIFFERENT BLOOD CANCER COHORTS
# df = fread('source_data/Fig_2_ESLs_in_cancer_U.csv')
df = fread('source_data/ExtFig_1_ESLs_in_cancer_M.csv')
df = as.data.frame(df)

levels = unique(df$group)
df$group = factor(df$group, levels = levels)
levels = unique(df$label)
df$label = factor(df$label, levels = levels)

set.seed(678)
p <- ggplot(df, aes(x = factor(group), y = beta, color = label, shape = label)) +  # 'x' mapped to factor for categorical axis
  geom_jitter(width = 0.35, height = 0, size = 0.6) +  # Add jitter with small width
  theme_classic() + ylim(0, 1) +
  guides(color = guide_legend(override.aes = list(size = 2)), shape = guide_legend(override.aes = list(size = 2))) +
  labs(x = NULL, y = "Mean β-value for each ESL", color = NULL, shape = NULL) +
  theme(axis.title.y = element_text(size=12, face='bold'),
        legend.position = 'inside',
        # legend.position.inside = c(0.15, 0.66), # for unmethylated
        legend.position.inside = c(0.15, 0.34), # for methylated
        legend.key.height = unit(0.35, 'cm'),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = NA, color = NA),  
        legend.box.background = element_rect(fill = NA, color = NA) ) +
  
  # geom_hline(yintercept = 0.0586, linetype = "dashed", color = "red", size = 0.7) + # for unmethylated
  geom_hline(yintercept = 0.8364, linetype = "dashed", color = "red", size = 0.7) + # for methylated
  
  scale_color_manual(
    values = c('black', 'black', 'black', 'grey', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c',
               '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#8dd3c7')) +
    scale_shape_manual(values = c(16, 15, 15, 15, 15, 16, 16, 15, 15, 16, 17, 16, 16, 16, 16, 16))
p

# ggsave('figures/Fig_2_ESLs_in_cancer_U.png', p, width = 10, height = 4)
ggsave('figures/ExtFig_1_ESLs_in_cancer_M.png', p, width = 10, height = 4)



################################# SuppFig_3_ESLs_in_cancer_median_destabilized_U.png #################################
### MEDIAN NUMBER OF DESTABILIZED ESLs IN EACH CANCER COHORT
df = read.csv('source_data/SuppFig_3_ESLs_in_cancer_median_destabilized_U.csv')
names(df) = c('cohort', 'median')
levels = unique(df$cohort)
df$cohort = factor(df$cohort, levels = levels)

p <- ggplot(df, aes(x = cohort, y = median, fill = cohort)) +
  geom_bar(stat = "identity") + guides(fill = 'none') +
  geom_text(aes(label = median),  # Add text labels
            vjust = -0.5,        # Position the text above the bars
            size = 3.2) +
  labs(x = NULL, y = "Median Number of Destabilized ESLs") +
  theme_classic() + ylim(0, max(df$median) + 50) +
  theme(axis.title.y = element_text(size=11, face='bold'),
        axis.text.y = element_text(size=10)) + 
  scale_fill_manual(
    values = c('black', 'black', 'black', 'grey', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c',
               '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#8dd3c7'))
p

ggsave('figures/SuppFig_3_ESLs_in_cancer_median_destabilized_U.png', p, width = 10, height = 3)



################################# Fig_6a_ESL_genomic_distribution.png #################################
esl.u = read.csv('source_data/ESL_U.csv', header=F)[,1]
esl.m = read.csv('source_data/ESL_M.csv', header=F)[,1]
valid_probes = read.csv('../data/general/valid_probes.csv', header=F)[,1]
manifest450k = read.csv('../data/general/humanmethylation450_15017482_v1-2.csv',
                        skip = 7)
manifest450k[manifest450k$Relation_to_UCSC_CpG_Island == '', 'Relation_to_UCSC_CpG_Island'] = 'Open Sea'
manifest450k$Relation_to_UCSC_CpG_Island = gsub('N_', '', manifest450k$Relation_to_UCSC_CpG_Island)
manifest450k$Relation_to_UCSC_CpG_Island = gsub('S_', '', manifest450k$Relation_to_UCSC_CpG_Island)
rownames(manifest450k) = manifest450k$IlmnID
manifest450k = manifest450k[valid_probes, ]

df.u = manifest450k[esl.u,]
df.m = manifest450k[esl.m,]
non.esl = setdiff(setdiff(manifest450k$IlmnID, esl.u), esl.m)
df.non = manifest450k[non.esl,]

counts.u = df.u %>% count(Relation_to_UCSC_CpG_Island, name = 'value')
counts.u <- counts.u %>%
  mutate(percentage = round(value / sum(value) * 100, 1))
counts.non = df.non %>% count(Relation_to_UCSC_CpG_Island, name = 'value')
counts.non <- counts.non %>%
  mutate(percentage = round(value / sum(value) * 100, 1))

p1 <- ggplot(counts.u, aes(x = "", y = value, fill = Relation_to_UCSC_CpG_Island)) +
  geom_bar(stat = "identity", width = 1, color='black',) +
  coord_polar("y") +
  theme_void() +
  labs(fill = "", title='ESL (n=31,569)') + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  geom_text(aes(x = 1.7,
                label = ifelse(percentage > 3, paste0(percentage, "%"), "")),
                position = position_stack(vjust = 0.5), size=3) +
  scale_fill_manual(
    values = c('#e31a1c', "#6a3d9a", "#33a02c", "#1f78b4"))

p2 <- ggplot(counts.non, aes(x = "", y = value, fill = Relation_to_UCSC_CpG_Island)) +
  geom_bar(stat = "identity", width = 1, color='black',) +
  coord_polar("y") +
  theme_void() +
  labs(fill = "", title='Non-ESL (n=357,915)') + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  geom_text(aes(x = 1.7,
                label = ifelse(percentage > 3, paste0(percentage, "%"), "")),
            position = position_stack(vjust = 0.5), size=3) +
  scale_fill_manual(
    values = c('#e31a1c', "#6a3d9a", "#33a02c", "#1f78b4")) 

fig = ggarrange(p1, p2, nrow=1, ncol=2,common.legend = T, legend = 'bottom')
fig
ggsave('figures/Fig_6a_ESL_genomic_distribution.png', fig, width = 7, height = 3)



################################# Fig_6b_ESL_proximity_to_TSS.png #################################
esl.u = read.csv('source_data/Fig_6b_ESL-U_TSS_distances.csv')
non.esl = read.csv('source_data/Fig_6b_NonESL_TSS_distances.csv')

# proportion within 1500bp of TSS:
dim(esl.u[esl.u$distTSS < 1500 & esl.u$distTSS > -1500,])
dim(esl.u[esl.u$distTSS < 0 & esl.u$distTSS > -1500,])

df = data.frame(type = c(rep('ESL', dim(esl.u)[1]), rep('Non-ESL', dim(non.esl)[1])),
                value = c(esl.u$distTSS, non.esl$distTSS))

dplot = ggplot(df, aes(x = value, color = type)) +  # Map 'fill' for the legend
  theme_classic() +
  geom_density(lwd = 1) +
  xlim(-1500, 1500) +
  labs(x = "Distance to TSS", y = "Density", color = NULL) +  # Set the legend title

  theme(axis.title.x = element_text(size=17, face='bold'),
        axis.title.y = element_text(size=17, face='bold'),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        legend.position = c(0.9, 0.85),
        legend.text = element_text(size = 13)) +
  scale_color_manual(values = c('#1b9e77', '#fc8d62'))
dplot

ggsave('figures/Fig_6b_ESL_proximity_to_TSS.png', dplot, width = 7, height = 4.5)



#################################### Fig_3cd_GSE139369_chromVAR_scores.png ####################################
df = read.csv('source_data/Fig_3cd_GSE139369_chromVAR_scores.csv', row.names = 'X')
df$Comparison[df$Comparison == 'HSC/CMP/GMP'] = 'HSC/\nCMP/GMP'
features = c('LymphoidEnriched.1072.hg19.100bp.frags', 'MyeloidEnriched.81.hg19.100bp.frags')
comparisons <- list(c("HSC/\nCMP/GMP", "Lymphoid"))

for (i in 1:length(features)) {
  ### Boxplot for cells of interest (mature myeloid lymphoid)
  figfile = paste0('Fig_3cd_GSE139369_chromVAR_', strsplit(features[[i]], split='\\.')[[1]][[1]], '.png')
  bplot = ggplot(df, aes(x = Comparison, y = get(features[[i]]), fill = Comparison)) +
    geom_boxplot() + 
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",
                       comparisons = comparisons) +
    theme_classic() +
    ylab('chromVAR score') + xlab(NULL) +
    theme(
      axis.text.y = element_text(size=14),
      axis.text.x = element_text(size=14),
      legend.position = 'none',
      axis.title.y = element_text(size=15, face='bold')
    )
  ggsave( paste0('figures/', figfile), bplot, width = 3, height = 5)
}



#################################### ExtFig_2cd_GSE203251_chromVAR_scores.png ####################################
# df = read.csv('source_data/ExtFig_2c_GSE203251_chromVAR_scores_LymphoidEnriched.csv')
# scorecol = 'LymphoidEnriched'
df = read.csv('source_data/ExtFig_2d_GSE203251_chromVAR_scores_MyeloidEnriched.csv')
scorecol = 'MyeloidEnriched'

df$group[df$group == 'HSC/CMP/GMP'] = 'HSC/\nCMP/GMP'
comparisons = list(c('HSC/\nCMP/GMP', 'Lymphoid'))
bplot = ggplot(df, aes(x = group, y = score, fill = group)) + geom_boxplot() +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     comparisons = comparisons) +
  theme_classic() + ylab('chromVAR score') + xlab(NULL) +
  theme(
    axis.text.y = element_text(size=14),
    axis.text.x = element_text(size=14),
    legend.position = 'none',
    axis.title.y = element_text(size=15, face='bold')
  )
ggsave(paste0('figures/ExtFig_2cd_GSE203251_chromVAR_', scorecol, '.png'), bplot, width = 3, height = 5)



#################################### Fig_6d_AML_promoter_hypermethylation.png ####################################
df = read.csv('../data/gene_expression/SupplementaryTable3_FULL.csv')
df = df %>% arrange(desc(mean_recurrence))
df <- df %>% distinct(gene, .keep_all = TRUE)
genes234 = unique(read.csv('../data/gene_expression/SupplementaryTable3.csv')$gene)
df = df[df$gene %in% genes234,]

# replace inf with max
max_value <- max(df$WGBS_.log10.Pvalue.[is.finite(df$WGBS_.log10.Pvalue.)], na.rm = TRUE)
df$WGBS_.log10.Pvalue.[!is.finite(df$WGBS_.log10.Pvalue.)] <- max_value

p_thresh = -log10(0.05 / nrow(df)) # Bonferroni corrected
# Subset data based on conditions
passed_red <- df[df$WGBS_log2FC > 0 & df$WGBS_.log10.Pvalue. > p_thresh, ]
passed_blue <- df[df$WGBS_log2FC < 0 & df$WGBS_.log10.Pvalue. > p_thresh, ]
# Get indices of rows not in passed_red or passed_blue
badIdx <- setdiff(rownames(df), c(rownames(passed_red), rownames(passed_blue)))
# Subset the failed rows
failed <- df[badIdx, ]

# label PRDM5
x = df[df$gene == 'PRDM5', 'WGBS_log2FC']
y = df[df$gene == 'PRDM5', 'WGBS_.log10.Pvalue.']

# Create scatter plots
vplot = ggplot() +
  geom_point(data = failed, aes(x = WGBS_log2FC, y = WGBS_.log10.Pvalue.), color = "grey", size = 2) +
  geom_point(data = passed_red, aes(x = WGBS_log2FC, y = WGBS_.log10.Pvalue.), color = "red", size = 2) +
  geom_point(data = passed_blue, aes(x = WGBS_log2FC, y = WGBS_.log10.Pvalue.), color = "#4292c6", size = 2) +
  
  geom_point(data = subset(df, gene == 'PRDM5'), 
             aes(x = WGBS_log2FC, y = WGBS_.log10.Pvalue.), 
             color = "black", size = 2) + 
  
  labs(x = "log2(Fold Change)", y = "-log10(P value)") +
  theme_classic() + 
  theme(axis.title.x = element_text(size=17, face='bold'),
        axis.title.y = element_text(size=17, face='bold'),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15)) + xlim(-3.79, 3.79) +
  annotate('text', x = df[df$gene == 'PRDM5', 'WGBS_log2FC'],
            y = df[df$gene == 'PRDM5', 'WGBS_.log10.Pvalue.'], label = "PRDM5",
  vjust = 0, hjust = 0, size = 6, color = "black")
vplot

ggsave('figures/Fig_6d_AML_promoter_hypermethylation.png', vplot, width = 7, height = 4.5)
