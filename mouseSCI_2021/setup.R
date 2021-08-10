


# Data extraction --------------------------------------------------

require('Seurat')

## Myeloid cell-level metadata ----
myeloid <- readRDS(file = './data/myeloid.rds')
DefaultAssay(myeloid) <- 'RNA'
myeloid$L1_taxon <- 'Myeloid'
myeloid$L2_taxon <- myeloid$celltype
myeloid$L3_taxon <- myeloid$myeloid_subcluster
myeloid$preprint_subtype <- myeloid$old_subcluster
myeloid$myeloid_seurat_res.0.35 <- myeloid$integrated_snn_res.0.35
obs_myeloid <- c('L1_taxon','L2_taxon','L3_taxon','celltype','myeloid_subcluster','preprint_subtype','myeloid_seurat_res.0.35')
umap_myeloid <- myeloid[['umap']]@cell.embeddings
colnames(umap_myeloid) <- c('myeloid_UMAP_1','myeloid_UMAP_2')
obs_myeloid <- cbind(myeloid@meta.data[, obs_myeloid], umap_myeloid)
for (i in 1:ncol(obs_myeloid)) {
  if (class(obs_myeloid[[i]]) == 'factor') {
    obs_myeloid[[i]] <- as.character(obs_myeloid[[i]])
  }
}
rm(myeloid, umap_myeloid)


## Vascular cell-level metadata ----
vascular <- readRDS(file = './data/vascular.rds')
DefaultAssay(vascular) <- 'RNA'
vascular$L1_taxon <- 'Vascular'
vascular$L2_taxon <- vascular$celltype
vascular$L3_taxon <- vascular$vascular_subcluster
vascular$preprint_subtype <- plyr::mapvalues(
  x = vascular$integrated_snn_res.0.4,
  from = 0:8,
  to = c('C1-Endothelial', 'C2-Endothelial','Tip Cell','A-Endothelial','U-Vascular','V-Endothelial','Fibroblast','Pericyte','VSMC')
)
vascular$vascular_seurat_res.0.4 <- vascular$integrated_snn_res.0.4
obs_vascular <- c('L1_taxon','L2_taxon','L3_taxon','celltype','vascular_subcluster','preprint_subtype','vascular_seurat_res.0.4')
umap_vascular <- vascular[['umap']]@cell.embeddings
colnames(umap_vascular) <- c('vascular_UMAP_1','vascular_UMAP_2')
obs_vascular <- cbind(vascular@meta.data[, obs_vascular], umap_vascular)
for (i in 1:ncol(obs_vascular)) {
  if (class(obs_vascular[[i]]) == 'factor') {
    obs_vascular[[i]] <- as.character(obs_vascular[[i]])
  }
}
rm(vascular, umap_vascular)


## Macroglia cell-level metadata ----
macroglia <- readRDS(file = './data/macroglia.rds')
DefaultAssay(macroglia) <- 'RNA'
macroglia$L1_taxon <- 'Macroglia'
macroglia$L2_taxon <- macroglia$celltype
macroglia$L3_taxon <- macroglia$macroglia_subcluster
macroglia$preprint_subtype <- macroglia$macroglia_subcluster
macroglia$macroglia_seurat_res.0.4 <- macroglia$integrated_snn_res.0.4
obs_macroglia <- c('L1_taxon','L2_taxon','L3_taxon','celltype','macroglia_subcluster','preprint_subtype','macroglia_seurat_res.0.4')
umap_macroglia <- macroglia[['umap']]@cell.embeddings
colnames(umap_macroglia) <- c('macroglia_UMAP_1','macroglia_UMAP_2')
obs_macroglia <- cbind(macroglia@meta.data[, obs_macroglia], umap_macroglia)
for (i in 1:ncol(obs_macroglia)) {
  if (class(obs_macroglia[[i]]) == 'factor') {
    obs_macroglia[[i]] <- as.character(obs_macroglia[[i]])
  }
}
rm(macroglia, umap_macroglia)



## SCI cell-level metadata ----
sci <- readRDS(file = './data/sci.rds')
DefaultAssay(sci) <- 'RNA'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL

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

cols_transfer <- c('celltype','myeloid_subcluster','vascular_subcluster','macroglia_subcluster','L1_taxon','L2_taxon','L3_taxon','preprint_subtype','myeloid_seurat_res.0.35','vascular_seurat_res.0.4','macroglia_seurat_res.0.4','myeloid_UMAP_1','myeloid_UMAP_2','vascular_UMAP_1','vascular_UMAP_2','macroglia_UMAP_1','macroglia_UMAP_2')
obs_transfer <- lapply(
  X = cols_transfer,
  FUN = match_metadata,
  order.char = rownames(sci@meta.data),
  obs.ls = list(obs_myeloid, obs_vascular, obs_macroglia)
)
obs_transfer <- data.frame(
  x = lapply(
    X = cols_transfer,
    FUN = match_metadata,
    order.char = rownames(sci@meta.data),
    obs.ls = list(obs_myeloid, obs_vascular, obs_macroglia)
  )
)
colnames(obs_transfer) <- cols_transfer
umap_sci <- sci[['umap']]@cell.embeddings
colnames(umap_sci) <- c('sci_UMAP_1','sci_UMAP_2')
cols_sci <- c('sample_id','time','orig.ident','nCount_RNA','nFeature_RNA','S.Score','G2M.Score','Phase','CC.Difference','percent_mt','percent_rp','percent_hbb','doublet_scores','dissociationMethod','chemistry')
obs_sci <- do.call(what = cbind, args = list(obs_transfer, sci@meta.data[, cols_sci], umap_sci))
obs_sci <- obs_sci[match(x = colnames(x_sci), table = rownames(obs_sci)),]
for (i in 1:ncol(obs_sci)) {
  if (class(obs_sci[[i]]) == 'factor') {
    obs_sci[[i]] <- as.character(obs_sci[[i]])
  }
}
obs_sci$L1_taxon[obs_sci$celltype == 'Neuron'] <- 'Neural'
obs_sci$L2_taxon[obs_sci$celltype == 'Neuron'] <- 'Neuron'
obs_sci$L2_taxon[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'

write.csv(x = obs_sci, file = './data/shinyAppData/obs_sci.csv', quote = FALSE)
obs_sci <- read.csv(file = './data/shinyAppData/obs_sci.csv', row.names = 1)



# Gene-level metadata -----------------------------------------------------

# Including which variables genes were used for clustering in each compoartment
# might be helpful information but it's pretty low priority. Will push back for
# a later date.

vars_sci <- sci@assays$RNA@meta.features

write.csv(x = vars_sci, file = './data/shinyAppData/vars_sci.csv', quote = FALSE, row.names = TRUE)
vars_sci <- read.csv(file = './data/shinyAppData/vars_sci.csv', header = FALSE, row.names = 1)



# Count matrix ------------------------------------------------------------

log_x_sci <- sci[['RNA']]@data
log_x_sci <- round(log_x_sci, digits = 4)

saveRDS(object = log_x_sci, file = './data/shinyAppData/log_x_sci.rds')
log_x_sci <- readRDS(file = './data/shinyAppData/log_x_sci.rds')



# Altogether now ----------------------------------------------------------

# Last minute finishes and data
genes <- rownames(log_x_sci)
numerical_data <- c(rownames(log_x_sci), colnames(obs_sci)[grepl('numeric', sapply(obs_sci, class))])
numerical_data <- numerical_data[!numerical_data %in% genes]
numerical_data <- numerical_data[!grepl('UMAP', numerical_data)]
categorical_data <- c(colnames(obs_sci)[grepl('character', sapply(obs_sci, class))])

s1 <- split(x = rownames(obs_sci), f = obs_sci$time)
s1 <- lapply(s1, function(x) {x[sample(1:length(x), size = length(x)/2)]})
s1 <- unlist(s1, use.names = FALSE)
s1 <- obs_sci[s1,]
s1 <- split(x = rownames(s1), f = s1$celltype)
s1 <- lapply(s1, function(x) {x[sample(1:length(x), size = length(x)/2)]})
s1 <- unlist(s1, use.names = FALSE)
sample_subset <- s1
saveRDS(sample_subset, file = './data/shinyAppData/sample_subset.rds')

log_x_sci <- log_x_sci[,sample_subset]
obs_sci <- obs_sci[sample_subset,]

save(log_x_sci, obs_sci, vars_sci, numerical_data, categorical_data, genes,
     file = './mouseSCI_2021/SCI_portal_data.RData')
     
# obs_sci <- obs_sci[,c(19,4,5,6,1,2,3,7:18,20:34)]
# genes <- genes[c(19443, 19515, 21138, 21213)]