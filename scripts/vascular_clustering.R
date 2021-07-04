
######### Vascular cluster analysis #########



# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('dplyr')
require('ggplot2')
require('clustree')
require('Seurat')
require('SingleCellExperiment')
# data_in <- './data/'
data_out <- './data/data_integration/'
results_out <- './results/vascular_clustering/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = data_out)
dir.create(path = results_out)

# sci <- readRDS(file = './data/sci.rds')
vascular <- readRDS(file = './data/vascular.rds')

# utils
firstup <- function(x) {x <- tolower(x); substr(x,1,1) <- toupper(substr(x,1,1)); return(x)}




# Batch Correction between replicates -------------------------------------


# Subset cells of the vascular compartment. Delete SCT and integrated assay slots
# for memory saving.
DefaultAssay(sci) <- 'RNA'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL
vascular_cells <- c('Fibroblast','Endothelial','Pericyte')
vascular <- subset(sci, idents = c('Fibroblast','Endothelial','Pericyte'))


# We can correct for batch effects (primarily due to sequencing depth) using
# the "batchelor" package. The function below is a wrapper for 
# batchelor::rescaleBatches(), which performs a linear regression-based 
# correction (ie downsampling counts). STRONG assumption: composition of cells
# in datasets are the same. Correction can fail if composition varies.
# Inputs:
#   obj_list: list of Seurat objects. These should be technical replicates of a 
#     given condition/group.
correct_replicates <- function(obj_list) {
  
  # Identify shared genes across datasets and take common subset
  gene_universe <- lapply(X = obj_list, 
                          FUN = function(x) {
                            return(rownames(slot(x[['RNA']], 'data')))
                          }
  )
  gene_universe <- Reduce(intersect, gene_universe)
  obj_list_sce <- lapply(X = obj_list, 
                         FUN = function(x) {
                           x <- as.SingleCellExperiment(x)
                           x <- x[gene_universe,]
                           return(x)
                         }
  )
  
  # Extract raw count and cell-level meta data. Merge all into
  # single matrix. These values derived directly from Seurat objects.
  raw_counts <- do.call(cbind, lapply(X = obj_list_sce, FUN = counts))
  log_raw_counts <- do.call(cbind, lapply(X = obj_list_sce, FUN = logcounts))
  cell_metadata <- do.call(rbind, lapply(X = obj_list_sce, FUN = colData))
  cell_metadata <- data.frame(cell_metadata)
  
  # Perform linear regression-based batch correction (ie scale down counts).
  # Automatically generates single, merged matrix. Recalculate the corrected 
  # "raw counts".
  corrected_log <- batchelor::rescaleBatches(obj_list_sce, 
                                             log.base = exp(1),
                                             pseudo.count = 1)
  corrected_log <- Matrix::Matrix(data = assays(corrected_log)[['corrected']],
                                  sparse = TRUE)
  
  # Set corrected values as new Seurat Assay slot
  corrected_log <- CreateAssayObject(counts = corrected_log)
  out_seurat <- CreateSeuratObject(counts = raw_counts,
                                   assay = 'RNA')
  out_seurat[['RNAcorrected']] <- corrected_log
  
  # Import log counts and metadata. Return assembled Seurat object.
  slot(out_seurat[['RNA']], 'data') <- log_raw_counts
  out_seurat@meta.data <- cell_metadata
  return(out_seurat)
}


# Set Seurat object slots
Idents(vascular) <- 'time'
DefaultAssay(vascular) <- 'RNA'

# Setup structures by condition
inj_groups <- levels(vascular$time)
vascular_corrected <- vector(mode = 'list', length = length(inj_groups))
names(vascular_corrected) <- inj_groups

# Perform the correction
for(ii in 1:length(inj_groups)) {
  inj_time <- inj_groups[ii]
  inj_cells <- unlist(x = CellsByIdentities(vascular, idents = inj_time),
                      use.names = FALSE)
  inj_cells <- subset(vascular, cells = inj_cells)
  inj_cells <- SplitObject(inj_cells, split.by = 'sample_id')
  inj_cells <- correct_replicates(obj_list = inj_cells)
  vascular_corrected[[inj_time]] <- inj_cells
  message(paste('Done with:', inj_time))
}

# clean up
# rm(inj_time, inj_cells, vascular, vascular_cells, sci)



# Vascular data integration --------------------------------------------------

memory.limit(64000)

# To identify shared cell-types/states across conditions, we perform data 
# integration according to Seurat's pipeline. We use RNA assay to identify 
# variable genes, but identify anchors with the RNAcorrected assay. 
vascular_corrected <- lapply(
  X = vascular_corrected,
  FUN = FindVariableFeatures,
  assay = 'RNA',
  nfeatures = 2000
)
vascular_corrected <- lapply(
  X = vascular_corrected,
  FUN = function(x) {
    slot(x[['RNAcorrected']], 'var.features') <- slot(x[['RNA']], 'var.features')
    return(x)
  }
)
vascular_anchors <- FindIntegrationAnchors(
  object.list = vascular_corrected, 
  assay = rep('RNAcorrected', length(vascular_corrected)),
  normalization.method = 'LogNormalize'
)

# Integrate into single dataset
vascular <- IntegrateData(vascular_anchors)
vascular <- ScaleData(vascular, vars.to.regress = 'CC.Difference')
vascular <- RunPCA(vascular, npcs = 20)
npcs <- 1:10
vascular <- RunUMAP(vascular, dims = npcs, min.dist = 0.3, n.neighbors = 30L, umap.method = 'uwot')
vascular <- FindNeighbors(vascular, dims = npcs)
vascular <- FindClusters(vascular, resolution = 0.8) # default resolution

tmp <- FindAllMarkers(vascular, assay = 'RNA', only.pos = TRUE, logfc.threshold = 0.75)
View(tmp %>% group_by(cluster) %>% top_n(n = 5, wt = -p_val_adj))


# Remove contaminating myeloid cluster (11) and macroglia cluster (14). 
# Cluster 11 enriched for Tyrobp, Csf1r, etc.
# Cluster 14 enriched for Mag, Mobp, Plp, etc.
myeloid_contamination <- unlist(CellsByIdentities(vascular, idents = 11), use.names = FALSE)
write.table(x = myeloid_contamination, 
            file = paste0(results_out, 'myeloid_contamination_in_vascular.tsv'),
            sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
macroglia_contamination <- unlist(CellsByIdentities(vascular, idents = 14), use.names = FALSE)
write.table(x = macroglia_contamination, 
            file = paste0(results_out, 'macroglia_contamination_in_vascular.tsv'),
            sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
remove_cells <- c(myeloid_contamination, macroglia_contamination)
retain <- colnames(vascular)[!colnames(vascular) %in% remove_cells]
vascular <- subset(vascular, cells = retain)


# Repeat integration:
vascular_corrected <- SplitObject(vascular, split.by = 'time')
vascular_corrected <- lapply(
  X = vascular_corrected,
  FUN = FindVariableFeatures,
  assay = 'RNA',
  nfeatures = 2000
)
vascular_corrected <- lapply(
  X = vascular_corrected,
  FUN = function(x) {
    slot(x[['RNAcorrected']], 'var.features') <- slot(x[['RNA']], 'var.features')
    return(x)
  }
)
vascular_anchors <- FindIntegrationAnchors(
  object.list = vascular_corrected, 
  assay = rep('RNAcorrected', length(vascular_corrected)),
  normalization.method = 'LogNormalize'
)

# Integrate into single dataset
vascular <- IntegrateData(vascular_anchors)
vascular <- ScaleData(vascular, vars.to.regress = 'CC.Difference')
vascular <- RunPCA(vascular, npcs = 20)
npcs <- 1:10
vascular <- RunUMAP(vascular, dims = npcs, min.dist = 0.3, n.neighbors = 30L, umap.method = 'uwot')
vascular <- FindNeighbors(vascular, dims = npcs)
vascular <- FindClusters(vascular, resolution = 0.8) # default resolution
vascular[['default_vascular_subcluster']] <- vascular[['seurat_clusters']]


# Set factors for metadata variables
vascular$time <- factor(x = vascular$time, levels = c('Uninjured','1dpi','3dpi','7dpi'))
vascular$sample_id <- factor(vascular$sample_id, levels = c('uninj_sample1', 'uninj_sample2', 'uninj_sample3', '1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2'))





# Determining stable clusters ---------------------------------------------


# Use IKAP package to determine stable clustering
source('./ref/IKAP_Seurat3.R')
DefaultAssay(vascular) <- 'integrated'
vascular <- IKAP(sobj = vascular, 
                 out.dir = paste0(results_out, 'IKAP'), 
                 confounders = c(),
                 pcs = 5:15,
                 scale.data = FALSE)
ikap_results <- vascular@meta.data[grepl('^PC.+[0-9]$', colnames(vascular@meta.data))]
write.table(x = ikap_results, file = paste0(results_out, 'IKAP/IKAP_results.tsv'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)
ikap_results <- read.table(file = paste0(results_out, 'IKAP/IKAP_results.tsv'),
                           sep = '\t', header = TRUE, row.names = 1)
vascular@meta.data[grepl('^PC.+[0-9]$', colnames(vascular@meta.data))] <- NULL


# Examine cluster results with clustree. We don't test FindClusters() "resolution"
# parameter because # of clusters is more meaningful (K in IKAP results).
check_pcs <- 6:10
cluster_results <- cbind(ikap_results, 
                         'time' = as.character(vascular@meta.data[['time']]),
                         'UMAP_1' = vascular[['umap']]@cell.embeddings[,1],
                         'UMAP_2' = vascular[['umap']]@cell.embeddings[,2])
clustree_plots <- vector(mode = 'list', length = length(check_pcs))
names(clustree_plots) <- paste0('PC', check_pcs)
for(ii in 1:length(check_pcs)) {
  clustree_plots[[ii]] <- clustree(cluster_results, prefix = paste0('PC', check_pcs[ii], 'K'))
}
clustree_out <- paste0(results_out, 'clustree/')
dir.create(path = clustree_out)
for(ii in 1:length(clustree_plots)) {
  ggsave(filename = paste0(clustree_out, names(clustree_plots)[ii], '.tiff'),
         plot = clustree_plots[[ii]], device = 'tiff', height = 9, width = 10)
}


# Elbow plot for PC selection. PCs selected based on IKAP results, clustree
# visualizations, and VizDimLoadings inspection.
DefaultAssay(vascular) <- 'integrated'
npcs <- 1:11
elbowplot <- ElbowPlot(vascular, ndims = 50) +
  geom_vline(mapping = aes(xintercept = max(npcs)), size = 1, color = 'red', linetype = 'dashed') +
  labs(title = paste0('PCA of vascular dataset (', ncol(vascular), ' cells)'),
       subtitle = 'PCs up to dashed line used for downstream analysis')
ggsave(plot = elbowplot, filename = paste0(results_out, 'pca_elbowplot.tiff'),
       device = 'tiff', height = 3.5, width = 5)


# Cluster vascular cells with stable parameters
npcs <- 1:11
DefaultAssay(vascular) <- 'integrated'
vascular <- RunUMAP(vascular, dims = npcs)
vascular <- FindNeighbors(vascular, dims = npcs)
vascular <- FindClusters(vascular, resolution = 0.2)
vascular[['vascular_subcluster']] <- vascular[['integrated_snn_res.0.2']]




# Preliminary UMAPs -------------------------------------------------------


# Global options
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

sample_cols <- RColorBrewer::brewer.pal(n = 10, name = 'Spectral')
sample_cols <- c('firebrick','dodgerblue','goldenrod','firebrick','dodgerblue','goldenrod','firebrick','dodgerblue','firebrick','dodgerblue')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')

dissociationMethod_cols <- c('firebrick','dodgerblue')
names(dissociationMethod_cols) <- c('Standard','Enriched')

chemistry_cols <- c('firebrick','dodgerblue')
names(chemistry_cols) <- c('v2','v3')


# set themes
umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 12, color = 'black'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 12, color = 'black'))


# Subclusters UMAP
tmp <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','vascular_subcluster')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = vascular_subcluster), size = 0.5, alpha = 0.5) +
  umap_theme +
  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
umap_subcluster <- LabelClusters(plot = tmp, id = 'vascular_subcluster', size = 5, repel = FALSE)
ggsave(filename = paste0(results_out, 'umap_subcluster.tiff'),
       plot = umap_subcluster, device = 'tiff', height = 4, width = 4.25)


# Injury time-point UMAP
umap_time <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = time), size = 0.5, alpha = 0.5) +
  scale_color_manual(values = time_cols) +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Time after\ninjury', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_time.tiff'),
       plot = umap_time, device = 'tiff', height = 4, width = 5.5)


# Original sample UMAP
umap_sample <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 0.5, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 1) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_sample.tiff'),
       plot = umap_sample, device = 'tiff', height = 3.5, width = 15)


# dissociationMethod UMAP
umap_dissociationMethod <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','dissociationMethod')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = dissociationMethod), size = 0.5, alpha = 0.5) +
  scale_color_manual(values = dissociationMethod_cols) +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Dissociation Method', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_dissociationMethod.tiff'),
       plot = umap_dissociationMethod, device = 'tiff', height = 4, width = 6.25)


# Chemistry UMAP
umap_chemistry <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','chemistry')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = chemistry), size = 0.5, alpha = 0.5) +
  scale_color_manual(values = chemistry_cols) +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = '10X Chromium\nchemistry', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_chemistry.tiff'),
       plot = umap_chemistry, device = 'tiff', height = 4, width = 5.75)



# Vascular UMAP by sample ------------------------------------------------------

vascular <- readRDS(file = './data/vascular.rds')

Idents(object = vascular) <- 'vascular_subcluster'

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 16, color = 'black'),
                    legend.title = element_text(size = 16, color = 'black'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 16, color = 'black'))

# cell-type annotation split by time UMAP
sample_cols <- c('#800000','#f58231','#4363d8','#800000','#f58231','#4363d8','#800000','#f58231','#800000','#f58231')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')
celltype_counts <- table(vascular$vascular_subcluster)
celltype_label <- paste0(names(celltype_counts), ' (', celltype_counts, ')')
names(celltype_label) <- names(celltype_counts)
umap_sample <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 1, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 2) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 16, color = 'black'),
        legend.title = element_text(size = 16, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))
umap_sample
ggsave(filename = paste0(results_out, 'vascular_bySample_umap.tiff'),
       plot = umap_sample, device = 'tiff', height = 6, width = 8.5)
rm(vascular); gc()



# Density UMAP ------------------------------------------------------------


# Take UMAP coordinates and for each cell calculate median distance to 50 
# closest cells. Take sample of cells from each time-point so that each group
# has same number of cells.
minCells_perTime <- min(table(vascular@meta.data[['time']]))
umap_sample <- split(x = vascular@meta.data, f = vascular@meta.data$time)
umap_sample <- lapply(
  X = umap_sample,
  FUN = function(x) {
    rownames(x)[sample(x = 1:nrow(x), size = minCells_perTime)]
  }
)
umap_coord <- vascular[['umap']]@cell.embeddings[unlist(umap_sample, use.names = FALSE),]
vascular_k_dist <- BiocNeighbors::findKNN(X = umap_coord, k = 50)[['distance']]
median_k_dist <- median(vascular_k_dist[, ncol(vascular_k_dist)])

# For each time point, calculate the neighborhood density of each cell by
# using a tri-cube kernel on UMAP coordinates. Median distance of all cells to 
# 50 nearest neighbors used to define maximum distance for the kernel. Method
# inspired by: https://www.nature.com/articles/s41586-019-0933-9#Sec7
times <- levels(vascular@meta.data[['time']])
dens <- c()
for (t in times) {
  tmp_coord <- vascular[['umap']]@cell.embeddings[umap_sample[[t]],]
  dists <- BiocNeighbors::findNeighbors(X = tmp_coord,
                                        threshold = median_k_dist)[['distance']]
  tmp_dens <- sapply(
    X = dists,
    FUN = function(x) {
      sum((1 - (x / median_k_dist)^3)^3)
    }
  )
  names(tmp_dens) <- rownames(tmp_coord)
  dens <- c(dens, tmp_dens)
}
umap_data <- expand.grid(unlist(umap_sample, use.names = FALSE),
                         levels(vascular@meta.data[['time']]),
                         stringsAsFactors = FALSE)
colnames(umap_data) <- c('barcode','time')
umap_data <- cbind(umap_data,vascular[['umap']]@cell.embeddings[umap_data$barcode,])
umap_data[['dens']] <- NA
umap_data[['dens']] <- ifelse(
  test = umap_data$time == vascular@meta.data[umap_data$barcode,'time'],
  yes = dens[umap_data$barcode[umap_data$time == vascular@meta.data[umap_data$barcode,'time']]],
  no = NA
)
umap_data[['dens']] <- log2(umap_data[['dens']])
umap_data[['time']] <- factor(
  x = umap_data[['time']],
  levels = levels(vascular@meta.data[['time']])
)

tmp_cols <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')[3:9])(50)
umap_density <- umap_data[order(umap_data$dens, na.last = FALSE, decreasing = FALSE),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = dens), size = 1.5, alpha = 0.7) +
  scale_color_gradientn(colors = tmp_cols,
                        na.value = 'grey70',
                        breaks = c(0, max(umap_data$dens, na.rm = TRUE)),
                        labels = c('Low','High')) +
  facet_wrap(. ~ time) +
  xlab(label = 'UMAP 1') + 
  ylab(label = 'UMAP 2') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.line = element_line(size = 0),
        axis.title = element_text(size = 20, color = 'black'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.text = element_text(size = 20, color = 'black'),
        legend.title = element_text(size = 20, angle = 90),
        legend.text = element_text(size = 20, color = 'black')) +
  guides(color = guide_colorbar(title = 'log2(density)',
                                title.position = 'left',
                                title.hjust = 0.5,
                                frame.colour = 'black',
                                ticks.colour = 'black'))
umap_density
ggsave(filename = paste0(results_out, 'vascular_density_umap.tiff'),
       plot = umap_density, device = 'tiff', height = 7, width = 8.5)



# save data
saveRDS(vascular, file = './data/vascular.rds')


# rm(list = ls()); gc()