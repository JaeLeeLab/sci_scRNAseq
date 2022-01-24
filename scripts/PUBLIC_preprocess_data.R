
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)

tmp_dir <- 'D:/MiamiProject/sci_scRNAseq/'

# Extract sci expression and dimreduc data -----------------------------------

sci <- readRDS(file = paste0(tmp_dir, 'data/sci.rds'))
DefaultAssay(sci) <- 'RNA'
sci <- NormalizeData(sci)
vars_sci <- sci[['RNA']]@meta.features
x_sci <- sci[['RNA']]@counts
log_x_sci <- sci[['RNA']]@data
obs_sci <- sci@meta.data
umap_sci <- sci[['umap']]@cell.embeddings
colnames(umap_sci) <- c('sci_UMAP_1','sci_UMAP_2')


# Extract myeloid metadata and dimreduc -----------------------------------

myeloid <- readRDS(file = paste0(tmp_dir, 'data/myeloid.rds'))
DefaultAssay(myeloid) <- 'RNA'
myeloid[['RNAcorrected']] <- NULL
myeloid[['integrated']] <- NULL
myeloid$Layer1_compartment <- 'Myeloid'
myeloid$Layer2_celltype <- myeloid$celltype
myeloid$Layer3_subtype <- myeloid$myeloid_subcluster
myeloid$preprint_subtype <- myeloid$old_subcluster
myeloid$myeloid_seurat_res.0.35 <- myeloid$integrated_snn_res.0.35
obs_myeloid <- c('Layer1_compartment','Layer2_celltype','Layer3_subtype','celltype','myeloid_subcluster','preprint_subtype','myeloid_seurat_res.0.35')
umap_myeloid <- myeloid[['umap']]@cell.embeddings
colnames(umap_myeloid) <- c('myeloid_UMAP_1','myeloid_UMAP_2')
obs_myeloid <- cbind(myeloid@meta.data[, obs_myeloid], umap_myeloid)
for (i in 1:ncol(obs_myeloid)) {
  if (class(obs_myeloid[[i]]) == 'factor') {
    obs_myeloid[[i]] <- as.character(obs_myeloid[[i]])
  }
}
rm(myeloid)

# Extract vascular metadata and dimreduc -----------------------------------

vascular <- readRDS(file = paste0(tmp_dir, 'data/vascular.rds'))
DefaultAssay(vascular) <- 'RNA'
vascular[['RNAcorrected']] <- NULL
vascular[['integrated']] <- NULL
vascular$Layer1_compartment <- 'Vascular'
vascular$Layer2_celltype <- vascular$celltype
vascular$Layer3_subtype <- vascular$vascular_subcluster
vascular$preprint_subtype <- plyr::mapvalues(
  x = vascular$integrated_snn_res.0.4,
  from = 0:8,
  to = c('C1-Endothelial', 'C2-Endothelial','Tip Cell','A-Endothelial','U-Vascular','V-Endothelial','Fibroblast','Pericyte','VSMC')
)
vascular$vascular_seurat_res.0.4 <- vascular$integrated_snn_res.0.4
obs_vascular <- c('Layer1_compartment','Layer2_celltype','Layer3_subtype','celltype','vascular_subcluster','preprint_subtype','vascular_seurat_res.0.4')
umap_vascular <- vascular[['umap']]@cell.embeddings
colnames(umap_vascular) <- c('vascular_UMAP_1','vascular_UMAP_2')
obs_vascular <- cbind(vascular@meta.data[, obs_vascular], umap_vascular)
for (i in 1:ncol(obs_vascular)) {
  if (class(obs_vascular[[i]]) == 'factor') {
    obs_vascular[[i]] <- as.character(obs_vascular[[i]])
  }
}
rm(vascular)

# Extract macroglia metadata and dimreduc -----------------------------------

macroglia <- readRDS(file = paste0(tmp_dir, 'data/macroglia.rds'))
DefaultAssay(macroglia) <- 'RNA'
macroglia[['RNAcorrected']] <- NULL
macroglia[['integrated']] <- NULL
macroglia$Layer1_compartment <- 'Macroglia'
macroglia$Layer2_celltype <- macroglia$celltype
macroglia$Layer3_subtype <- macroglia$macroglia_subcluster
macroglia$preprint_subtype <- macroglia$macroglia_subcluster
macroglia$macroglia_seurat_res.0.4 <- macroglia$integrated_snn_res.0.4
obs_macroglia <- c('Layer1_compartment','Layer2_celltype','Layer3_subtype','celltype','macroglia_subcluster','preprint_subtype','macroglia_seurat_res.0.4')
umap_macroglia <- macroglia[['umap']]@cell.embeddings
colnames(umap_macroglia) <- c('macroglia_UMAP_1','macroglia_UMAP_2')
obs_macroglia <- cbind(macroglia@meta.data[, obs_macroglia], umap_macroglia)
for (i in 1:ncol(obs_macroglia)) {
  if (class(obs_macroglia[[i]]) == 'factor') {
    obs_macroglia[[i]] <- as.character(obs_macroglia[[i]])
  }
}
rm(macroglia)

# Combine metadata ---------------------------------------------------------

# convenience function to pull myeloid/vascular/macroglia metadata in order of 
# barcodes in sci@meta.data rows
match_metadata <- function(order.char, obs.ls, col.char) {
  meta_dat <- rep(x = NA, times = length(order.char))
  for (i in 1:length(obs.ls)) {
    if (col.char %in% colnames(obs.ls[[i]])) {
      meta_match <- match(x = rownames(obs.ls[[i]]),
                          table = order.char)
      meta_dat[meta_match] <- obs.ls[[i]][, col.char]
    }
  }
  return(meta_dat)
}
cols_transfer <- c('Layer1_compartment', 'Layer2_celltype', 'Layer3_subtype',
                   'preprint_subtype', 'myeloid_UMAP_1', 'myeloid_UMAP_2',
                   'vascular_UMAP_1', 'vascular_UMAP_2', 'macroglia_UMAP_1',
                   'macroglia_UMAP_2', 'myeloid_subcluster', 
                   'vascular_subcluster', 'macroglia_subcluster')
obs_transfer <- data.frame(
  x = lapply(
    X = cols_transfer,
    FUN = match_metadata,
    order.char = rownames(obs_sci),
    obs.ls = list(obs_myeloid, obs_vascular, obs_macroglia)
  )
)
colnames(obs_transfer) <- cols_transfer

sci_meta_retain <- c('sample_id','time','orig.ident','nCount_RNA','nFeature_RNA','S.Score','G2M.Score','Phase','CC.Difference','percent_mt','percent_rp','doublet_scores','dissociationMethod','chemistry', 'celltype','percent_hbb')
obs_sci <- do.call(
  what = cbind, 
  args = list(obs_transfer, obs_sci[, sci_meta_retain], umap_sci)
)


# Tidying metadata --------------------------------------------------------

# Rename neurons and lymphocyte compartments
obs_sci$Layer1_compartment[obs_sci$celltype == 'Neuron'] <- 'Neural'
obs_sci$Layer1_compartment[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'
obs_sci$Layer2_celltype[obs_sci$celltype == 'Neuron'] <- 'Neuron'
obs_sci$Layer2_celltype[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'
obs_sci$Layer3_subtype[obs_sci$celltype == 'Neuron'] <- 'Neuron'
obs_sci$Layer3_subtype[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'
obs_sci$preprint_subtype[obs_sci$celltype == 'Neuron'] <- 'Neuron'
obs_sci$preprint_subtype[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'

meta_order <- c(
  'celltype',
  'myeloid_subcluster',
  'vascular_subcluster',
  'macroglia_subcluster',
  'time',
  'sample_id',
  'orig.ident',
  'preprint_subtype',
  'nCount_RNA', 'nFeature_RNA', 'S.Score','G2M.Score','Phase','CC.Difference',
  'percent_mt', 'percent_rp', 'percent_hbb','doublet_scores',
  'dissociationMethod','chemistry',
  'Layer1_compartment', 'Layer2_celltype', 'Layer3_subtype',
  'sci_UMAP_1', 
  'sci_UMAP_2', 
  'myeloid_UMAP_1', 
  'myeloid_UMAP_2', 
  'vascular_UMAP_1', 
  'vascular_UMAP_2',
  'macroglia_UMAP_1', 
  'macroglia_UMAP_2'
)
obs_sci <- obs_sci[meta_order]



# SingleR prediction result import ----------------------------------------

Rosenberg_SingleR <- read.table(file = './results/sci_annotation_crossReference/rosenberg_singleR_results.tsv', sep = '\t', header = TRUE, row.names = 1)
colnames(Rosenberg_SingleR) <- paste('Rosenberg_SingleR', colnames(Rosenberg_SingleR), sep = '.')
Sathyamurthy_SingleR <- read.table(file = './results/sci_annotation_crossReference/Sathyamurthy_singleR_results.tsv', sep = '\t', header = TRUE, row.names = 1)
colnames(Sathyamurthy_SingleR) <- paste('Sathyamurthy_SingleR', colnames(Sathyamurthy_SingleR), sep = '.')
obs_sci <- Reduce(f = cbind, x = list(obs_sci, Rosenberg_SingleR, Sathyamurthy_SingleR))


# Save --------------------------------------------------------------------

write.csv(x = vars_sci, file = './data/PUBLIC_sci_scRNAseq/vars_sci.csv')
write.csv(x = obs_sci, file = './data/PUBLIC_sci_scRNAseq/obs_sci.csv')
saveRDS(object = x_sci, file = './data/PUBLIC_sci_scRNAseq/x_sci.rds')
saveRDS(object = sci, file = './data/PUBLIC_sci_scRNAseq/sci.rds')
