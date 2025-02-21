library(Seurat)
library(Signac)
library(Matrix)
library(data.table)
library(BPCells)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(patchwork)

test = readRDS('data/scATAC/GSE203251_combined_mat_array_peaks_hg38_100bp.RDS')

# for all cells
metadata = fread('data/scATAC/input_data_gotcha//GSE203251_MetadataPatients.csv.gz')
counts = test

# create chrom assay for signac to use
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = 'data/scATAC/input_data_gotcha/combined.fragments.tsv.gz',
  min.cells = 1, # 5,
  min.features = 1) # 200)

# load metadata related to cells and standardize cell names to match
frags = open_fragments_dir('data/scATAC/input_data_gotcha/combined_atac_filtered')
gotcha_metadata = metadata
gotcha_metadata = as.data.frame(gotcha_metadata)
row.names(gotcha_metadata) = gotcha_metadata$Whitelist
gotcha_metadata$Patient = gsub(' ', '_', gotcha_metadata$Patient)
gotcha_metadata$Patient = gsub('\\(|\\)', '', gotcha_metadata$Patient)
gotcha_metadata = gotcha_metadata[frags@cell_names,]

# create signac object
signac_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = gotcha_metadata
)

signac_obj <- FindTopFeatures(signac_obj, min.cutoff = 1)
signac_obj <- RunTFIDF(signac_obj)
# signac_obj <- RunSVD(signac_obj)
signac_obj

signac_obj2 = subset(signac_obj, subset = JAK2_Genotype == 'WT')

#########
### PROVIDE FOLDER; CREATE LIST OF FEATURES (DIFFERENT CPG SETS TO USE)
folder_name = 'data/scATAC/hg38/'
create_feature_list <- function(folder_name) {
  # Create a list of all items within the folder
  items <- list.files(folder_name, pattern = "\\.RDS$", full.names = TRUE, ignore.case = TRUE)
  # Initialize an empty named list
  named_list <- list()
  # Iterate through each file
  for (file_path in items) {
    # Extract the file name without the extension
    file_name <- basename(file_path)
    name_without_extension <- sub("\\.RDS$", "", file_name, ignore.case = TRUE)
    # Add the file path to the named list with the extracted name
    # named_list[[name_without_extension]] <- unique(readRDS(file_path))
    
    # weird error where a few rows not found in signac_obj
    named_list[[name_without_extension]] <- intersect(unique(readRDS(file_path)), rownames(signac_obj))
  }
  return(named_list)
}
features <- create_feature_list(folder_name)

signac_obj2 <- AddChromatinModule(signac_obj2, features = features, genome = BSgenome.Hsapiens.UCSC.hg38, verbose = T)
# All metadata including the scores for each module is found in signac_obj@meta.data (a data frame)

## compare scores between myeloid and lymphoid cells 
lymphoid_cells = c('CLP', 'B_cells', 'CD4_effector', 'CD4_Tcells', 'CD8_Tcells', 'NK')
myeloid_cells = c('HSC', 'HSC_GM', 'CMP', 'GMP')

## AVERAGE SCORES FOR CELL TYPES
for (n in names(features)) {
  newcol = paste0(n, '-mean')
  signac_obj2@meta.data[, newcol] = signac_obj2@meta.data[, n]
  
  for (celltype in unique(signac_obj@meta.data$ClusterAnnotated)) {
    scores = signac_obj2@meta.data[signac_obj2@meta.data$ClusterAnnotated == celltype, n]
    scores = scores[!is.na(scores)]
    mean_score = mean(scores)
    signac_obj2@meta.data[signac_obj2@meta.data$ClusterAnnotated == celltype, newcol] = mean_score
  }
}

jak2_wt = rownames(signac_obj@meta.data[signac_obj@meta.data$JAK2_Genotype == 'WT',])

mydf = signac_obj2@meta.data

for (scorecol in c("LymphoidEnriched-1072-hg38-100bp-frags", "MyeloidEnriched-81-hg38-100bp-frags")) {
  lymph_wt = mydf[(mydf$ClusterAnnotated %in% lymphoid_cells) & mydf$JAK2_Genotype == 'WT', scorecol]
  lymph_mut = mydf[(mydf$ClusterAnnotated %in% lymphoid_cells) & mydf$JAK2_Genotype == 'MUT', scorecol]
  myelo_wt = mydf[(mydf$ClusterAnnotated %in% myeloid_cells) & mydf$JAK2_Genotype == 'WT', scorecol]
  myelo_mut = mydf[(mydf$ClusterAnnotated %in% myeloid_cells) & mydf$JAK2_Genotype == 'MUT', scorecol]
  
  wtdf <- data.frame(
    score = c(lymph_wt, myelo_wt),
    group = factor(rep(c("Lymphoid", "HSC/CMP/GMP"),
                       c(length(lymph_wt), length(myelo_wt)))))
  
  cpgset = strsplit(scorecol, '-')[[1]][1]
  write.csv(wtdf, paste0('plots/source_data/ExtFig_2_GSE203251_chromVAR_scores_', cpgset, '.csv'))
}

# Plot via umap coordinates from pub provided metadata
umap_pub = as.matrix(signac_obj2@meta.data[,c('UMAP1', 'UMAP2')])
signac_obj2[["umap_pub"]] <- CreateDimReducObject(embeddings = umap_pub, key = "UMAP_", assay = DefaultAssay(signac_obj2))

### PLOT UMAPS
panels = list()
# UMAPs colored by module score (chromVAR deviation score)
for (i in 1:length(features)) {
  ### UMAP colored by chromVAR score
  figfile = gsub('frags', 'umap.png', names(features)[[i]])
  
    this_fig = FeaturePlot(signac_obj2, reduction = "umap_pub",  cells = jak2_wt, pt.size = 1.5,  raster = F, features = c( paste0(names(features)[[i]], '-mean') )) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
    ) + ggtitle(NULL)  + coord_cartesian(expand = FALSE) + guides(color = "none")
  panels[[i]] = this_fig
}

panelfig = panels[[1]] | panels[[2]]
ggsave('plots/figures/ExtFig_2ab_scATAC_additional.png', panelfig, width = 20, height = 10 )
