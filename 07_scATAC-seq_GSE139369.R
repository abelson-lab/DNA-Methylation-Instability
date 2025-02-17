.libPaths('/../R_4.2.2_package_library_ubuntu20')
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scater)
library(Signac)
library(Seurat)
library(Matrix)
library(BPCells)
library(ArchR)
library(ggplot2)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RColorBrewer)

signac_obj = readRDS('/data/scATAC/scATAC-All-Hematopoiesis-MPAL-191120_signac.rds')
myeloid_cells = c('01_HSC', '05_CMP.LMPP', '07_GMP', '08_GMP.Neut')
lymphoid_cells = c('15_CLP.2', '16_Pre.B', '17_B', '19_CD8.N', '20_CD4.N1', '21_CD4.N2', '22_CD4.M', '23_CD8.EM', '24_CD8.CM', '25_NK')
compare_cells = c(myeloid_cells, lymphoid_cells)
signac_healthy = subset(signac_obj, subset = BioClassification %in% compare_cells)

### PROVIDE FOLDER; CREATE LIST OF FEATURES (DIFFERENT CPG SETS TO USE)
folder_name = '/data/scATAC/frags/'
create_feature_list <- function(folder_name) {
  items <- list.files(folder_name, pattern = "\\.RDS$", full.names = TRUE, ignore.case = TRUE)
  named_list <- list()
  for (file_path in items) {
    file_name <- basename(file_path)
    name_without_extension <- sub("\\.RDS$", "", file_name, ignore.case = TRUE)
    named_list[[name_without_extension]] <- gsub('_', '-', unique(readRDS(file_path)))
  }
  return(named_list)
}
features <- create_feature_list(folder_name)
signac_healthy <- AddChromatinModule(signac_healthy, features = features, genome = BSgenome.Hsapiens.UCSC.hg19, verbose = T)

## AVERAGE SCORES FROM AGGREGATED GROUPS
cell_groups = list('HSC' = c('01_HSC'),
                   'CMP' = c('05_CMP.LMPP'),
                   'GMP' = c('07_GMP', '08_GMP.Neut'),
                   'CLP' = c('15_CLP.2'),
                   'Pre.B' = c('16_Pre.B'),
                   'B' = c('17_B'),
                   'CD4' = c('20_CD4.N1', '21_CD4.N2', '22_CD4.M'),
                   'CD8' = c('19_CD8.N', '23_CD8.EM', '24_CD8.CM'),
                   'NK' = c('25_NK'))

for (n in names(features)) {
  newcol = paste0(n, '-mean')
  signac_healthy@meta.data[, newcol] = NA
  
  for (celltypes in cell_groups) {
    scores = signac_healthy@meta.data[signac_healthy@meta.data$BioClassification %in% celltypes, n]
    scores = scores[!is.na(scores)]
    mean_score = mean(scores)
    signac_healthy@meta.data[signac_healthy@meta.data$BioClassification %in% celltypes, newcol] = mean_score
  }
}

### PLOT UMAPS
figures_dir = '/plots/figures/'

h_meta = signac_healthy@meta.data
h_meta = h_meta[h_meta$BioClassification %in% compare_cells,]
h_meta['Comparison'] = ''
h_meta[h_meta$BioClassification %in% myeloid_cells, 'Comparison'] = 'HSC/CMP/GMP'
h_meta[h_meta$BioClassification %in% lymphoid_cells, 'Comparison'] = 'Lymphoid'
h_meta = h_meta[,c('Comparison', 'LymphoidEnriched-1072-hg19-100bp-frags', 'MyeloidEnriched-81-hg19-100bp-frags')]
### save chromVAR scores for plotting
write.csv(h_meta, 'plots/inputs/GSE139369_chromVAR_scores.csv')

### Add embeddings for UMAP based on a subset of the full dataset
myumap = readRDS('/data/scATAC/my_umap_coords_2.RDS')
mycells = myumap[!is.na(myumap$MY_UMAP1), 'Barcode']

og_meta = signac_healthy@meta.data
all_cells = signac_healthy@meta.data
all_cells$og_row = rownames(all_cells)
og_row_order = rownames(all_cells)

merge.meta = merge(all_cells, myumap, by=c('Group','Barcode'))
rownames(merge.meta) = merge.meta$og_row
merge.meta = merge.meta[og_row_order,]
signac_healthy@meta.data = as.data.frame(merge.meta)

coords = as.matrix(signac_healthy@meta.data[, c('MY_UMAP1','MY_UMAP2')])
signac_healthy[['my_umap']] <- CreateDimReducObject(
  embeddings = coords,
  key = 'myumap_',
  assay = DefaultAssay(signac_healthy)
)

### UMAP colored by chromVAR score
xrng = range(signac_healthy@meta.data$MY_UMAP1[!is.na(signac_healthy@meta.data$MY_UMAP1)])
yrng = range(signac_healthy@meta.data$MY_UMAP2[!is.na(signac_healthy@meta.data$MY_UMAP2)])

panels = list()
for (i in 1:length(features)) {
  figfile = gsub('frags', 'umap.png', names(features)[[i]])
  this_fig = FeaturePlot(signac_healthy, reduction = "my_umap", pt.size = 0.5, raster = F, features = c( gsub('-', '.', paste0(names(features)[[i]], '-mean')) )) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
  ) + ggtitle(NULL) + guides(color = "none") +
    coord_cartesian(
      xlim = xrng + c(-0, 0),  # Expand x-axis by 1 unit on both sides
      ylim = yrng + c(-0, 0)  # Expand y-axis by 0.5 units on both sides
    )
  panels[[i]] = this_fig
}
panelfig = panels[[1]] | panels[[2]]
ggsave( paste0(figures_dir, 'HCpanels-12-14.png'), panelfig, width = 20, height = 10 )
