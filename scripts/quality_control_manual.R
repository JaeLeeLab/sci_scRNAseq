
####### Manual Quality Control ########


# Through intial cluster analysis of each sample individually, we identified
# one large cluster of cells in 7dpi_sample1 and one large of cells in 
# 3dpi_sample1 which are distinct from all other cells and do not express
# known marker genes at high levels. These clusters are similar in that they
# express some markers associated with myeloids and form a circular topology 
# in the resulting UMAPs or tSNEs. Furthermore, differential expression tests 
# show that these clusters have very few DE genes compared to other clusters in
# their respective datasets. We suspect these are low quality cells and manually
# discard them here before downstream data integration. 
# 
# i.e. cluster-based approach to remove putative low quality cells.


# stochastic methods are used
set.seed(123)


# Data Import -------------------------------------------------------------

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
data_in <- './data/QC_filtered_feature_bc_matrix/'
data_out <- './data/QC_filtered_feature_bc_matrix/'
results_out <- './results/quality_control/'
ref_out <- './ref/'


needs_filter <- c('3dpi_sample1','7dpi_sample1', '1dpi_sample1') # Include 1dpi_sample1 for control

# Import data
counts_in <- paste0(data_in, list.files(path = data_in))
counts_in <- counts_in[grepl(pattern = paste(needs_filter, collapse = '|'), x = counts_in)]
counts <- vector(mode = 'list', length = length(counts_in))
for(ii in 1:length(counts_in)) {
  counts[[ii]] <- readRDS(file = counts_in[ii])
}
names(counts) <- gsub(pattern = paste(data_in, '.rds', sep = "|"), replacement = '', counts_in)



# Manual filtering via low-quality cluster identification -----------------


# Remove duplicated and blank gene names in count matrix (if names not available
# in Ensembl database, assume not informative for biological investigation).
remove_duplicated_genes <- function(x) {
  bad_gene <- duplicated(rownames(x)) | (nchar(rownames(x)) == 0)
  counts <- x[!bad_gene,]
  return(counts)
}
counts <- lapply(X = counts, FUN = remove_duplicated_genes)

seurat_workflow <- function(x) {
  x <- x %>%
    CreateSeuratObject() %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = 15) %>%
    RunUMAP(dims = 1:15) %>%
    FindNeighbors(dims = 1:15) %>%
    FindClusters(resolution = 0.4)
  return(x)
}
counts <- lapply(counts, FUN = seurat_workflow)
umaps <- lapply(counts, FUN = DimPlot, label = TRUE, label.size = 5)
umaps <- cowplot::plot_grid(plotlist = umaps, ncol = 1)

markers <- lapply(
  X = counts,
  FUN = FindAllMarkers,
  only.pos = TRUE,
  logfc.threshold = 0.5,
  assay = 'RNA',
  slot = 'data',
  test.use = 'wilcox'
)
markers <- lapply(
  X = markers,
  FUN = function(x) {
    marker_gene <- x %>%
      group_by(cluster) %>%
      top_n(n = 1, wt = avg_logFC) %>%
      .[['gene']]
    return(marker_gene)
  }
)
markers_heatmap <- vector(mode = 'list', length = length(needs_filter))
names(markers_heatmap) <- names(counts)
for(ii in 1:length(markers_heatmap)) {
  sample_id <- names(markers_heatmap)[ii]
  counts[[sample_id]] <- ScaleData(counts[[sample_id]], 
                                   features = markers[[sample_id]])
  markers_heatmap[[sample_id]] <- DoHeatmap(
    object = counts[[sample_id]],
    features = markers[[sample_id]],
    angle = 0
  ) +
    labs(title = paste(sample_id, 'marker heatmap')) +
    theme(plot.title = element_text(size = 10, color = 'black', face = 'bold'),
          legend.position = 'none')
}
markers_heatmap <- cowplot::plot_grid(plotlist = markers_heatmap, ncol = 1)
needs_filter_results <- cowplot::plot_grid(markers_heatmap, umaps, ncol = 2, rel_widths = c(1, 0.6))

# save result
ggsave(filename = paste0(results_out, 'manual_filter_summary.tiff'),
       plot = needs_filter_results, device = 'tiff',
       height = length(counts) * 3.5, width = 12)



# Filter and save data ----------------------------------------------------

manual_filter_out <- paste0(data_out, 'pre_manual_filter_data/')
dir.create(path = manual_filter_out)


# Filter putative low quality cells
# 3dpi_sample1:
tmp <- counts[['qc_filtered_feature_bc_matrix_3dpi_sample1']]
tmp <- tmp[, tmp$seurat_clusters != 0]
counts[['new_qc_filtered_feature_bc_matrix_3dpi_sample1']] <- tmp
# 7dpi_sample1:
tmp <- counts[['qc_filtered_feature_bc_matrix_7dpi_sample1']]
tmp <- tmp[, tmp$seurat_clusters != 1]
counts[['new_qc_filtered_feature_bc_matrix_7dpi_sample1']] <- tmp


saveRDS(object = counts[['qc_filtered_feature_bc_matrix_3dpi_sample1']]@assays$RNA@counts,
        file = paste0(manual_filter_out, 'premanual_qc_filtered_feature_bc_matrix_3dpi_sample1.rds'))
saveRDS(object = counts[['qc_filtered_feature_bc_matrix_7dpi_sample1']]@assays$RNA@counts,
        file = paste0(manual_filter_out, 'premanual_qc_filtered_feature_bc_matrix_7dpi_sample1.rds'))
saveRDS(object = counts[['new_qc_filtered_feature_bc_matrix_3dpi_sample1']]@assays$RNA@counts,
        file = paste0(data_out, 'qc_filtered_feature_bc_matrix_3dpi_sample1.rds'))
saveRDS(object = counts[['new_qc_filtered_feature_bc_matrix_7dpi_sample1']]@assays$RNA@counts,
        file = paste0(data_out, 'qc_filtered_feature_bc_matrix_7dpi_sample1.rds'))
saveRDS(object = counts, file = paste0(manual_filter_out, 'premanual_filter_processed_data.rds'))

# Note: file sizes slightly smaller because duplicate/unnamed genes were removed.


rm(list = ls())
gc()