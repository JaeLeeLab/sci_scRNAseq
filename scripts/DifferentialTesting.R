
######## Differential Expression testing  ########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('scran')
require('foreach')
require('parallel')
require('doParallel')
require('dplyr')
require('ggplot2')
require('Seurat')
results_out <- './results/differential_expression_testing/'
dir.create(path = results_out)

sci <- readRDS(file = './data/sci.rds')





# Pairwise wilcox DE test within cell-types b/w samples -------------------

DefaultAssay(sci) <- 'RNA'
Idents(sci) <- 'celltype'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL


# Based on ependymal uninj_sample1 vs uninj_sample2, good threshold for p_val_adj
# is ~1e-10 (72 hits). When max.cells.per.ident is set to same number of cells 
# per group (62), then the number of results below 1e-10 is reduced to 8.
celltypes <- levels(sci@meta.data[['celltype']])
sample_names <- levels(sci@meta.data[['sample_id']])
times <- levels(sci@meta.data[['time']])
cell_count <- table(sci@meta.data[['celltype']], 
                    sci@meta.data[['sample_id']])
group_pairs <- t(combn(sample_names, m = 2))

# Setup results list
results_byCell <- vector(mode = 'list', length = length(celltypes))
names(results_byCell) <- celltypes

# For each cell-type
for (c in 1:length(celltypes)) {
  
  # Extract expression data
  c_dat <- sci[,sci$celltype == celltypes[c]]
  Idents(c_dat) <- 'sample_id'
  
  # Wrapper for Seurat::FindMarkers
  test_pairs <- function(object, subset.ident, ident.1, ident.2) {
    
    # If less than 3 cells per group, return empty data.frame
    if (cell_count[subset.ident, ident.1] < 3 | 
        cell_count[subset.ident, ident.2] < 3) {
      de_result <- data.frame(row.names = row.names(object))
      de_result[['gene']] <- rownames(de_result)
      de_result[['avg_logFC']] <- NA
      de_result[['p_val_adj']] <- NA
      return(de_result)
    }
    de_result <- Seurat::FindMarkers(
      object = object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      assay = 'RNA',
      slot = 'data', 
      test.use = 'wilcox'
    )
    
    # Take only gene, p-vals, and fold-changes values.
    de_result <- de_result[c('avg_logFC', 'p_val_adj')]
    colnames(de_result) <- paste(colnames(de_result),
                                 paste(ident.1, ident.2, sep = '.'),
                                 sep = '.')
    de_result[['gene']] <- rownames(de_result)
    return(de_result)
  }
  
  # For tracking progress
  message(paste('Testing:', celltypes[c]))
  
  # For parallel computing, use all but 1 core.
  cores <- detectCores()
  cl <- makeCluster(spec = cores[1] - 1)
  registerDoParallel(cl = cl)
  
  # Use foreach package to loop through DE tests
  results_byCell[[celltypes[c]]] <- foreach(i = 1:nrow(group_pairs)) %dopar% {
    tmp <- test_pairs(
      object = c_dat,
      subset.ident = celltypes[c],
      ident.1 = group_pairs[i, 1],
      ident.2 = group_pairs[i, 2]
    )
  }
  
  # Remove setup for parallel
  stopCluster(cl = cl)
}
rm(cores, cl, test_pairs, c_dat); gc()

# merge results into single data.frame per cell-type
results_byCell <- lapply(
  X = results_byCell,
  FUN = function(x) {
    Reduce(
      f = function(x, y) {
        merge(x, y, all.x = TRUE)
      },
      x = x
    )
  }
)

# Determine number of overlapping DE genes and their correlation between
# all combinations of comparisons.
pval_thresh <- 1e-10
sect_byCell <- vector(mode = 'list', length = length(results_byCell))
cor_byCell <- vector(mode = 'list', length = length(results_byCell))
set_pairs <- expand.grid(1:nrow(group_pairs), 1:nrow(group_pairs))
# set_pairs <- t(combn(x = 1:nrow(group_pairs), m = 2))
for (c in 1:length(results_byCell)) {
  de_sect <- rep(NA, times = nrow(set_pairs))
  de_cor <- rep(NA, times = nrow(set_pairs))
  for (sp in 1:nrow(set_pairs)) {
    set_col <- c(paste(group_pairs[set_pairs[sp, 1],], collapse = '.'),
                 paste(group_pairs[set_pairs[sp, 2],], collapse = '.'))
    set_col <- paste('p_val_adj', set_col, sep = '.')
    set_col <- match(x = set_col, table = colnames(results_byCell[[c]]))
    if (any(is.na(set_col))) {
      de_sect[sp] <- NA
      de_cor[sp] <- NA
      next
    }
    if (set_col[1] == set_col[2]) {
      de_sect[sp] <- sum(results_byCell[[c]][set_col[1]] < pval_thresh)
      de_cor[sp] <- 1
      next
    } else {
      pv1 <- which(results_byCell[[c]][set_col[1]] < pval_thresh)
      pv2 <- which(results_byCell[[c]][set_col[2]] < pval_thresh)
      if (length(pv1) == 0 | length(pv2) == 0) {
        de_sect[sp] <- NA
        de_cor[sp] <- NA
        next
      } else {
        g1 <- results_byCell[[c]][['gene']][pv1]
        g2 <- results_byCell[[c]][['gene']][pv2]
        de_sect[sp] <- length(x = intersect(x = g1, y = g2))
        de_cor[sp] <- c(cor(
          x = results_byCell[[c]][set_col[1]],
          y = results_byCell[[c]][set_col[2]],
          use = 'complete.obs',
          method = 'spearman'
        ))
      }
    }
  }
  sect_byCell[[c]] <- matrix(
    data = de_sect,
    nrow = nrow(group_pairs),
    byrow = FALSE
  )
  cor_byCell[[c]] <- matrix(
    data = de_cor,
    nrow = nrow(group_pairs),
    byrow = FALSE
  )
  tmp_names <- apply(
    X = group_pairs, 
    MARGIN = 1,
    FUN = function(x) paste(x[1], x[2], sep = '.')
  )
  dimnames(sect_byCell[[c]]) <- list(tmp_names, tmp_names)
  dimnames(cor_byCell[[c]]) <- list(tmp_names, tmp_names)
}

# # Example heatmap
# ComplexHeatmap::Heatmap(
#   matrix = cor_byCell[[c]],
#   cluster_rows = FALSE,
#   cluster_columns = FALSE
# )

saveRDS(results_byCell, file = paste0(results_out, 'results_byCell_wilcox.rds'))
results_byCell_wilcox <- readRDS(file = paste0(results_out, 'results_byCell_wilcox.rds'))





# Pairwise MAST DE test within cell-types b/w samples -------------------

# setup
DefaultAssay(sci) <- 'RNA'
Idents(sci) <- 'celltype'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL

# metadata prep
celltypes <- levels(sci@meta.data[['celltype']])
sample_names <- levels(sci@meta.data[['sample_id']])
times <- levels(sci@meta.data[['time']])
cell_count <- table(sci@meta.data[['celltype']], 
                    sci@meta.data[['sample_id']])
group_pairs <- t(combn(sample_names, m = 2))

# Setup results list
results_byCell <- vector(mode = 'list', length = length(celltypes))
names(results_byCell) <- celltypes

# For each cell-type
for (c in 1:length(celltypes)) {
  
  # Extract expression data
  c_dat <- sci[,sci$celltype == celltypes[c]]
  Idents(c_dat) <- 'sample_id'
  
  # Wrapper for Seurat::FindMarkers
  test_pairs <- function(object, subset.ident, ident.1, ident.2) {
    
    # If less than 3 cells per group, return empty data.frame
    if (cell_count[subset.ident, ident.1] < 3 | 
        cell_count[subset.ident, ident.2] < 3) {
      de_result <- data.frame(row.names = row.names(object))
      de_result[['gene']] <- rownames(de_result)
      de_result[['avg_logFC']] <- NA
      de_result[['p_val_adj']] <- NA
      return(de_result)
    }
    de_result <- Seurat::FindMarkers(
      object = object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      assay = 'RNA',
      slot = 'data', 
      test.use = 'MAST',
      latent.vars = 'nFeature_RNA'
    )
    
    # Take only gene, p-vals, and fold-changes values.
    de_result <- de_result[c('avg_logFC', 'p_val_adj')]
    colnames(de_result) <- paste(colnames(de_result),
                                 paste(ident.1, ident.2, sep = '.'),
                                 sep = '.')
    de_result[['gene']] <- rownames(de_result)
    return(de_result)
  }
  
  # For tracking progress
  message(paste('Testing:', celltypes[c]))
  
  # For parallel computing, use all but 1 core.
  cores <- detectCores()
  cl <- makeCluster(spec = cores[1] - 1)
  registerDoParallel(cl = cl)
  
  # Use foreach package to loop through DE tests
  results_byCell[[celltypes[c]]] <- foreach(i = 1:nrow(group_pairs)) %dopar% {
    tmp <- test_pairs(
      object = c_dat,
      subset.ident = celltypes[c],
      ident.1 = group_pairs[i, 1],
      ident.2 = group_pairs[i, 2]
    )
  }
  
  # Remove setup for parallel
  stopCluster(cl = cl)
}
rm(cores, cl, test_pairs, c_dat); gc()

# merge results into single data.frame per cell-type
results_byCell <- lapply(
  X = results_byCell,
  FUN = function(x) {
    Reduce(
      f = function(x, y) {
        merge(x, y, all.x = TRUE)
      },
      x = x
    )
  }
)

# Determine number of overlapping DE genes and their correlation between
# all combinations of comparisons.
pval_thresh <- 1e-10
sect_byCell <- vector(mode = 'list', length = length(results_byCell))
cor_byCell <- vector(mode = 'list', length = length(results_byCell))
set_pairs <- expand.grid(1:nrow(group_pairs), 1:nrow(group_pairs))
# set_pairs <- t(combn(x = 1:nrow(group_pairs), m = 2))
for (c in 1:length(results_byCell)) {
  de_sect <- rep(NA, times = nrow(set_pairs))
  de_cor <- rep(NA, times = nrow(set_pairs))
  for (sp in 1:nrow(set_pairs)) {
    set_col <- c(paste(group_pairs[set_pairs[sp, 1],], collapse = '.'),
                 paste(group_pairs[set_pairs[sp, 2],], collapse = '.'))
    set_col <- paste('p_val_adj', set_col, sep = '.')
    set_col <- match(x = set_col, table = colnames(results_byCell[[c]]))
    if (any(is.na(set_col))) {
      de_sect[sp] <- NA
      de_cor[sp] <- NA
      next
    }
    if (set_col[1] == set_col[2]) {
      de_sect[sp] <- sum(results_byCell[[c]][set_col[1]] < pval_thresh)
      de_cor[sp] <- 1
      next
    } else {
      pv1 <- which(results_byCell[[c]][set_col[1]] < pval_thresh)
      pv2 <- which(results_byCell[[c]][set_col[2]] < pval_thresh)
      if (length(pv1) == 0 | length(pv2) == 0) {
        de_sect[sp] <- NA
        de_cor[sp] <- NA
        next
      } else {
        g1 <- results_byCell[[c]][['gene']][pv1]
        g2 <- results_byCell[[c]][['gene']][pv2]
        de_sect[sp] <- length(x = intersect(x = g1, y = g2))
        de_cor[sp] <- c(cor(
          x = results_byCell[[c]][set_col[1]],
          y = results_byCell[[c]][set_col[2]],
          use = 'complete.obs',
          method = 'spearman'
        ))
      }
    }
  }
  sect_byCell[[c]] <- matrix(
    data = de_sect,
    nrow = nrow(group_pairs),
    byrow = FALSE
  )
  cor_byCell[[c]] <- matrix(
    data = de_cor,
    nrow = nrow(group_pairs),
    byrow = FALSE
  )
  tmp_names <- apply(
    X = group_pairs, 
    MARGIN = 1,
    FUN = function(x) paste(x[1], x[2], sep = '.')
  )
  dimnames(sect_byCell[[c]]) <- list(tmp_names, tmp_names)
  dimnames(cor_byCell[[c]]) <- list(tmp_names, tmp_names)
}

# # Example heatmap
# ComplexHeatmap::Heatmap(
#   matrix = cor_byCell[[c]],
#   cluster_rows = FALSE,
#   cluster_columns = FALSE
# )

saveRDS(results_byCell, file = paste0(results_out, 'results_byCell_MAST.rds'))
results_byCell_mast <- readRDS(file = paste0(results_out, 'results_byCell_MAST.rds'))





# Pairwise MAST DE test within cell-types b/w samples w/ downsampling -------------------

# setup
DefaultAssay(sci) <- 'RNA'
Idents(sci) <- 'celltype'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL

# metadata prep
celltypes <- levels(sci@meta.data[['celltype']])
sample_names <- levels(sci@meta.data[['sample_id']])
times <- levels(sci@meta.data[['time']])
cell_count <- table(sci@meta.data[['celltype']], 
                    sci@meta.data[['sample_id']])
group_pairs <- t(combn(sample_names, m = 2))

# Setup results list
results_byCell <- vector(mode = 'list', length = length(celltypes))
names(results_byCell) <- celltypes

# For each cell-type
for (c in 1:length(celltypes)) {
  
  # Extract expression data
  c_dat <- sci[,sci$celltype == celltypes[c]]
  Idents(c_dat) <- 'sample_id'
  
  # Wrapper for Seurat::FindMarkers
  test_pairs <- function(object, subset.ident, ident.1, ident.2) {
    
    # If less than 3 cells per group, return empty data.frame
    if (cell_count[subset.ident, ident.1] < 3 | 
        cell_count[subset.ident, ident.2] < 3) {
      de_result <- data.frame(row.names = row.names(object))
      de_result[['gene']] <- rownames(de_result)
      de_result[['avg_logFC']] <- NA
      de_result[['p_val_adj']] <- NA
      return(de_result)
    }
    de_result <- Seurat::FindMarkers(
      object = object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      assay = 'RNA',
      slot = 'data', 
      test.use = 'MAST',
      latent.vars = 'nFeature_RNA',
      max.cells.per.ident = min(cell_count[subset.ident, ident.1],
                                cell_count[subset.ident, ident.2]),
      random.seed = 123
    )
    
    # Take only gene, p-vals, and fold-changes values.
    de_result <- de_result[c('avg_logFC', 'p_val_adj')]
    colnames(de_result) <- paste(colnames(de_result),
                                 paste(ident.1, ident.2, sep = '.'),
                                 sep = '.')
    de_result[['gene']] <- rownames(de_result)
    return(de_result)
  }
  
  # For tracking progress
  message(paste('Testing:', celltypes[c]))
  
  # For parallel computing, use all but 1 core.
  cores <- detectCores()
  cl <- makeCluster(spec = cores[1] - 1)
  registerDoParallel(cl = cl)
  
  # Use foreach package to loop through DE tests
  results_byCell[[celltypes[c]]] <- foreach(i = 1:nrow(group_pairs)) %dopar% {
    tmp <- test_pairs(
      object = c_dat,
      subset.ident = celltypes[c],
      ident.1 = group_pairs[i, 1],
      ident.2 = group_pairs[i, 2]
    )
  }
  
  # Remove setup for parallel
  stopCluster(cl = cl)
}
rm(cores, cl, test_pairs, c_dat); gc()

# merge results into single data.frame per cell-type
results_byCell <- lapply(
  X = results_byCell,
  FUN = function(x) {
    Reduce(
      f = function(x, y) {
        merge(x, y, all.x = TRUE)
      },
      x = x
    )
  }
)

# Determine number of overlapping DE genes and their correlation between
# all combinations of comparisons.
pval_thresh <- 1e-10
sect_byCell <- vector(mode = 'list', length = length(results_byCell))
cor_byCell <- vector(mode = 'list', length = length(results_byCell))
set_pairs <- expand.grid(1:nrow(group_pairs), 1:nrow(group_pairs))
# set_pairs <- t(combn(x = 1:nrow(group_pairs), m = 2))
for (c in 1:length(results_byCell)) {
  de_sect <- rep(NA, times = nrow(set_pairs))
  de_cor <- rep(NA, times = nrow(set_pairs))
  for (sp in 1:nrow(set_pairs)) {
    set_col <- c(paste(group_pairs[set_pairs[sp, 1],], collapse = '.'),
                 paste(group_pairs[set_pairs[sp, 2],], collapse = '.'))
    set_col <- paste('p_val_adj', set_col, sep = '.')
    set_col <- match(x = set_col, table = colnames(results_byCell[[c]]))
    if (any(is.na(set_col))) {
      de_sect[sp] <- NA
      de_cor[sp] <- NA
      next
    }
    if (set_col[1] == set_col[2]) {
      de_sect[sp] <- sum(results_byCell[[c]][set_col[1]] < pval_thresh)
      de_cor[sp] <- 1
      next
    } else {
      pv1 <- which(results_byCell[[c]][set_col[1]] < pval_thresh)
      pv2 <- which(results_byCell[[c]][set_col[2]] < pval_thresh)
      if (length(pv1) == 0 | length(pv2) == 0) {
        de_sect[sp] <- NA
        de_cor[sp] <- NA
        next
      } else {
        g1 <- results_byCell[[c]][['gene']][pv1]
        g2 <- results_byCell[[c]][['gene']][pv2]
        de_sect[sp] <- length(x = intersect(x = g1, y = g2))
        de_cor[sp] <- c(cor(
          x = results_byCell[[c]][set_col[1]],
          y = results_byCell[[c]][set_col[2]],
          use = 'complete.obs',
          method = 'spearman'
        ))
      }
    }
  }
  sect_byCell[[c]] <- matrix(
    data = de_sect,
    nrow = nrow(group_pairs),
    byrow = FALSE
  )
  cor_byCell[[c]] <- matrix(
    data = de_cor,
    nrow = nrow(group_pairs),
    byrow = FALSE
  )
  tmp_names <- apply(
    X = group_pairs, 
    MARGIN = 1,
    FUN = function(x) paste(x[1], x[2], sep = '.')
  )
  dimnames(sect_byCell[[c]]) <- list(tmp_names, tmp_names)
  dimnames(cor_byCell[[c]]) <- list(tmp_names, tmp_names)
}

# # Example heatmap
# ComplexHeatmap::Heatmap(
#   matrix = cor_byCell[[c]],
#   cluster_rows = FALSE,
#   cluster_columns = FALSE
# )

saveRDS(results_byCell, file = paste0(results_out, 'results_byCell_MAST_downsample.rds'))
results_byCell_downsample <- readRDS(file = paste0(results_out, 'results_byCell_MAST_downsample.rds'))




# MDS of cell-types and samples  ---------------------------------------------

sci_cols <- c('Neutrophil' = '#800000',
              'Monocyte' = '#9a6324',
              'Macrophage' = '#e6194b',
              'Dendritic' = '#f58231',
              'Microglia' = 'gold',
              'Div-Myeloid' = '#808000',
              'Fibroblast' = '#333300',
              'Endothelial' = '#3cb44b',
              'Pericyte' = '#008080',
              'OPC' = 'cyan3',
              'Oligodendrocyte' = '#000075',
              'Astrocyte' = '#4363d8',
              'Ependymal' = '#911eb4',
              'Lymphocyte' = '#f032e6',
              'Neuron' = '#e6beff')

DefaultAssay(sci) <- 'RNA'
meta_vars <- c('sample_id', 'time', 'dissociationMethod', 'chemistry', 'celltype')
expr_data <- FetchData(
  object = sci,
  vars = c(sci[['integrated']]@var.features, meta_vars)
)
expr_data <- expr_data %>%
  group_by_at(.vars = meta_vars) %>%
  summarise(across(.cols = where(is.numeric), .fns = mean))
expr_data_meta <- expr_data[c(meta_vars)]
expr_data_dist <- dist(
  x = expr_data[!colnames(expr_data) %in% meta_vars],
  method = 'euclidean'
)
expr_data_mds <- as.data.frame(cmdscale(
  d = expr_data_dist,
  k = 2
))
colnames(expr_data_mds) <- c('MDS_1','MDS_2')


# Generate the MDS plots. Change shape of points by different metadata.
sci_mds_1 <- cbind(expr_data_mds, expr_data_meta) %>%
  ggplot(mapping = aes(x = MDS_1, y = MDS_2)) + 
  geom_point(mapping = aes(fill = celltype, shape = time), size = 5) +
  scale_fill_manual(values = sci_cols) +
  scale_size_area(max_size = 10) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  labs(title = 'MDS of cell-type average expression', 
       subtitle = 'Points represent a cell-type from each sample') +
  theme(panel.background = element_rect(color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.title = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12),
        legend.box = 'horizontal',
        legend.key = element_rect(fill = NA)) +
  guides(fill = guide_legend(title = 'Cell-type',
                             override.aes = list(size = 5, color = sci_cols),
                             ncol = 2),
         shape = guide_legend(title = 'Time',
                              override.aes = list(fill = 'black')))

sci_mds_2 <- cbind(expr_data_mds, expr_data_meta) %>%
  ggplot(mapping = aes(x = MDS_1, y = MDS_2)) + 
  geom_point(mapping = aes(fill = dissociationMethod), size = 5, pch = 21) +
  scale_fill_manual(values = c('Standard' = 'dodgerblue',
                               'Enriched' = 'indianred')) +
  scale_size_area(max_size = 10) +
  theme(panel.background = element_rect(color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.title = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12),
        legend.box = 'horizontal',
        legend.key = element_rect(fill = NA)) +
  guides(fill = guide_legend(
    title = 'Dissociation\nMethod',
    override.aes = list(size = 5)))

sci_mds_3 <- cbind(expr_data_mds, expr_data_meta) %>%
  ggplot(mapping = aes(x = MDS_1, y = MDS_2)) + 
  geom_point(mapping = aes(fill = chemistry), size = 5, pch = 21) +
  scale_fill_manual(values = c('v2' = 'dodgerblue',
                               'v3' = 'indianred')) +
  scale_size_area(max_size = 10) +
  theme(panel.background = element_rect(color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.title = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12),
        legend.box = 'horizontal',
        legend.key = element_rect(fill = NA)) +
  guides(fill = guide_legend(
    title = '10X Chemistry',
    override.aes = list(size = 5)))

sci_mds_row1 <- cowplot::plot_grid(NULL, sci_mds_1, NULL, ncol = 3, rel_widths = c(0.15, 1, 0.15))
sci_mds_row2 <- cowplot::plot_grid(sci_mds_2, sci_mds_3, ncol = 2)
sci_mds <- cowplot::plot_grid(sci_mds_row1, sci_mds_row2, ncol = 1, rel_heights = c(1,0.95))
sci_mds

ggsave(filename = paste0(results_out, 'sci_mds.tiff'),
       plot = sci_mds, device = 'tiff', height = 8.5, width = 11.5)







# MDS of subclusters and samples  ---------------------------------------------


subclusters <- c('Neutrophil',
                 'Monocyte',
                 'Macrophage-A',
                 'Macrophage-B',
                 'Div-Myeloid',
                 'Microglia-A',
                 'Microglia-B',
                 'Microglia-C',
                 'Div-Microglia',
                 'BA-Macrophage',
                 'A-Endothelial',
                 'C1-Endothelial',
                 'C2-Endothelial',
                 'V-Endothelial',
                 'Tip Cell',
                 'Fibroblast',
                 'Pericyte',
                 'VSMC',
                 'Ependymal-A',
                 'Ependymal-B',
                 'Astroependymal',
                 'Astrocyte',
                 'OPC-A',
                 'OPC-B',
                 'Div-OPC',
                 'Oligodendrocyte')

subcluster_cols <- colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 11, name = 'Spectral')
)(length(subclusters))

names(subcluster_cols) <- subclusters


DefaultAssay(sci) <- 'RNA'
meta_vars <- c('sample_id', 'time', 'dissociationMethod', 'chemistry', 'subcluster')
sci_lite <- DietSeurat(
  object = sci,
  data = TRUE,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL
)
expr_data <- FetchData(
  object = sci_lite,
  vars = c(sci[['integrated']]@var.features, meta_vars)
)
expr_data[['subcluster']] <- factor(
  x = expr_data[['subcluster']],
  levels = subclusters
)
expr_data <- expr_data[expr_data$subcluster %in% subclusters,]
expr_data <- expr_data %>%
  group_by_at(.vars = meta_vars) %>%
  summarise(across(.cols = where(is.numeric), .fns = mean))
expr_data_meta <- expr_data[c(meta_vars)]
expr_data_dist <- dist(
  x = expr_data[!colnames(expr_data) %in% meta_vars],
  method = 'euclidean'
)
expr_data_mds <- as.data.frame(cmdscale(
  d = expr_data_dist,
  k = 2
))
colnames(expr_data_mds) <- c('MDS_1','MDS_2')


# Generate the MDS plots. Change shape of points by different metadata.
sci_mds_1 <- cbind(expr_data_mds, expr_data_meta) %>%
  ggplot(mapping = aes(x = MDS_1, y = MDS_2)) + 
  geom_point(mapping = aes(fill = subcluster, shape = time), size = 5) +
  scale_fill_manual(values = subcluster_cols) +
  scale_size_area(max_size = 10) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  labs(title = 'MDS of subcluster average expression', 
       subtitle = 'Points represent a subcluster from each sample') +
  theme(panel.background = element_rect(color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.title = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12),
        legend.box = 'horizontal',
        legend.key = element_rect(fill = NA)) +
  guides(fill = guide_legend(title = 'Subcluster', 
                             override.aes = list(fill = subcluster_cols,
                                                 color = 'black',
                                                 pch = 21)),
         shape = guide_legend(title = 'Time',
                              override.aes = list(fill = 'black')))

sci_mds_2 <- cbind(expr_data_mds, expr_data_meta) %>%
  ggplot(mapping = aes(x = MDS_1, y = MDS_2)) + 
  geom_point(mapping = aes(fill = dissociationMethod), size = 5, pch = 21) +
  scale_fill_manual(values = c('Standard' = 'dodgerblue',
                               'Enriched' = 'indianred')) +
  scale_size_area(max_size = 10) +
  theme(panel.background = element_rect(color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.title = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12),
        legend.box = 'horizontal',
        legend.key = element_rect(fill = NA)) +
  guides(fill = guide_legend(
    title = 'Dissociation\nMethod',
    override.aes = list(size = 5)))

sci_mds_3 <- cbind(expr_data_mds, expr_data_meta) %>%
  ggplot(mapping = aes(x = MDS_1, y = MDS_2)) + 
  geom_point(mapping = aes(fill = chemistry), size = 5, pch = 21) +
  scale_fill_manual(values = c('v2' = 'dodgerblue',
                               'v3' = 'indianred')) +
  scale_size_area(max_size = 10) +
  theme(panel.background = element_rect(color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.title = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12),
        legend.box = 'horizontal',
        legend.key = element_rect(fill = NA)) +
  guides(fill = guide_legend(
    title = '10X Chemistry',
    override.aes = list(size = 5)))

sci_mds_row1 <- cowplot::plot_grid(NULL, sci_mds_1, NULL, ncol = 3, rel_widths = c(0.15, 1, 0.15))
sci_mds_row2 <- cowplot::plot_grid(sci_mds_2, sci_mds_3, ncol = 2)
sci_mds <- cowplot::plot_grid(sci_mds_row1, sci_mds_row2, ncol = 1, rel_heights = c(1,0.95))
sci_mds

ggsave(filename = paste0(results_out, 'sci_mds_subcluster.tiff'),
       plot = sci_mds, device = 'tiff', height = 8.5, width = 11.5)





# Uninjured vs 1dpi DE pre-print data vs true biological replciates  -----------

# Setup data. First, subset seurat object with RNA assay for DE testing.
Idents(sci) <- 'celltype'
DefaultAssay(sci) <- 'RNA'
sci_lite <- DietSeurat(
  object = sci,
  assays = 'RNA',
  graphs = NULL,
  dimreducs = NULL
)

# We will test for DE between times
celltypes <- levels(sci_lite@meta.data[['celltype']])
study <- list(
  original = c('uninj_sample1','uninj_sample2','1dpi_sample1','1dpi_sample2'),
  trueRep = c('uninj_sample2','uninj_sample3','1dpi_sample2','1dpi_sample3')
)

# Iterate through cell-types, study
de_bycell <- vector(mode = 'list', length = length(celltypes))
names(de_bycell) <- celltypes
for (ct in celltypes) {
  dat <- sci_lite[,sci_lite$time %in% c('Uninjured','1dpi') & sci_lite$celltype == ct]
  count_byStudy <- list(
    'original' = table(dat$sample_id, dat$time)[study$original,1:2],
    'trueRep' = table(dat$sample_id, dat$time)[study$trueRep,1:2]
  )
  count_byStudy <- sapply(
    X = count_byStudy,
    FUN = function(x) {
      return(any(c(x[1:2, 1], x[3:4, 2]) == 0))
    }
  )
  if (any(count_byStudy)) {
    next
  }
  de_bystudy <- vector(mode = 'list', length = length(study))
  names(de_bystudy) <- names(study)
  for (s in 1:length(study)) {
    Idents(dat) <- 'sample_id'
    tmp_de <- FindMarkers(
      object = dat,
      assay = 'RNA',
      slot = 'data',
      ident.1 = study[[s]][1:2],
      ident.2 = study[[s]][3:4]
    )
    de_bystudy[[names(study)[s]]] <- tmp_de
    de_bystudy[[names(study)[s]]][['gene']] <- rownames(de_bystudy[[names(study)[s]]])
    de_bystudy[[names(study)[s]]][['study']] <- names(study)[s]
  }
  de_bystudy <- Reduce(f = rbind, x = de_bystudy)
  de_bycell[[ct]] <- de_bystudy
}
saveRDS(object = de_bycell,
        file = paste0(results_out, 'results_DE_wilcox_uninj_1dpi_byCell_byStudy.rds'))

# Retrieve summary of results: 1) table of up or down-regulated genes b/w uninj
# and 1dpi for each study, per cell-type, 2) names of genes that
# were up or down-reg for each comparison.
de_bycell_summary <- lapply(
  X = de_bycell,
  FUN = function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    tmp <- x %>%
      filter(p_val_adj < 1e-03) %>%
      mutate('dir' = ifelse(test = avg_logFC > 0,
                            yes = 'up',
                            no = 'down'))
    nde <- tmp %>%
      group_by(study, dir) %>%
      summarise('nde' = n_distinct(gene))
    de_genes <- split(
      x = tmp,
      f = list(tmp[['study']], tmp[['dir']])
    )
    de_genes <- lapply(
      X = de_genes,
      FUN = function(x) {
        x[['gene']]
      }
    )
    return(list('nde' = nde,
                'd_genes' = de_genes))
  }
)

# Generate Venn diagrams to show overlap 
de_venns <- vector(mode = 'list', length = length(de_bycell_summary))
names(de_venns) <- names(de_bycell_summary)
for (i in 1:length(de_bycell_summary)) {
  x <- de_bycell_summary[[i]]
  original <- grepl(pattern = 'original',
                    x = names(x[['d_genes']]))
  original <- Reduce(f = union, x = x[['d_genes']][original])
  trueRep <- grepl(pattern = 'trueRep',
                    x = names(x[['d_genes']]))
  trueRep <- Reduce(f = union, x = x[['d_genes']][trueRep])
  venn <- eulerr::euler(
    combinations = list(
      'Original (sample1 + sample2)' = original,
      'Biological replicates (sample2 + sample3)' = trueRep
    )
  )
  venn <- plot(
    x = venn,
    quantities = list(fontsize = 16),
    fills = c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)),
    main = paste('Uninjured vs 1dpi', names(de_bycell_summary)[i]),
    legend = list('side' = 'bottom'),
    # labels = list(font = 1, fontsize = 14),
  )
  de_venns[[i]] <- venn
}

# save venn diagrams
tiff(filename = paste0(results_out, 'StudyComparison_DE.tiff'),
     height = 10, width = 20, res = 440, units = 'in')
gridExtra::grid.arrange(grobs = de_venns, ncol = 5)
dev.off()

# Save venn diagrams of cell-types present in uninjured cord (ie remove cells
# that were enrichment-specific and peripheral immune cells)
keep <- c('Microglia','Pericyte','OPC','Ependymal','Endothelial')
tiff(filename = paste0(results_out, 'StudyComparison_DE_subset.tiff'),
     height = 5, width = 12, res = 440, units = 'in')
gridExtra::grid.arrange(
  grobs = de_venns[keep], 
  ncol = ceiling(sqrt(length(keep)))
)
dev.off()





# Uninjured vs 1dpi DE by dissociation method -----------------------------------

# Setup data. First, subset seurat object with RNA assay for DE testing.
Idents(sci) <- 'celltype'
sci_lite <- DietSeurat(
  object = sci,
  assays = 'RNA',
  graphs = NULL,
  dimreducs = NULL
)

# We will test for DE between times, for each dissociation method and cell-type.
celltypes <- levels(sci@meta.data[['celltype']])
diss <- unique(sci@meta.data[['dissociationMethod']])

# Iterate through cell-types, dissociation method
de_bycell <- vector(mode = 'list', length = length(celltypes))
names(de_bycell) <- celltypes
for (ct in celltypes) {
  dat <- sci_lite[,sci_lite$time %in% c('Uninjured','1dpi') & sci_lite$celltype == ct]
  if (any(table(dat$dissociationMethod, dat$time)[,1:2] < 3)) {
    next
  }
  de_bydiss <- vector(mode = 'list', length = length(diss))
  names(de_bydiss) <- diss
  for (m in diss) {
    Idents(dat) <- 'dissociationMethod'
    tmp_de <- FindMarkers(
      object = dat,
      assay = 'RNA',
      slot = 'data',
      subset.ident = m,
      group.by = 'time',
      ident.1 = 'Uninjured',
      ident.2 = '1dpi'
    )
    de_bydiss[[m]] <- tmp_de
    de_bydiss[[m]][['gene']] <- rownames(de_bydiss[[m]])
    de_bydiss[[m]][['d_method']] <- m
  }
  de_bydiss <- Reduce(f = rbind, x = de_bydiss)
  de_bycell[[ct]] <- de_bydiss
}
saveRDS(object = de_bycell,
        file = paste0(results_out, 'results_DE_wilcox_uninj_1dpi_byCell_byDissociationMethod.rds'))

# Retrieve summary of results: 1) table of up or down-regulated genes b/w uninj
# and 1dpi for each dissociation method, per cell-type, 2) names of genes that
# were up or down-reg for each comparison.
de_bycell_summary <- lapply(
  X = de_bycell,
  FUN = function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    tmp <- x %>%
      filter(p_val_adj < 1e-03) %>%
      mutate('dir' = ifelse(test = avg_logFC > 0,
                            yes = 'up',
                            no = 'down'))
    nde <- tmp %>%
      group_by(d_method, dir) %>%
      summarise('nde' = n_distinct(gene))
    de_genes <- split(
      x = tmp,
      f = list(tmp[['d_method']], tmp[['dir']])
    )
    de_genes <- lapply(
      X = de_genes,
      FUN = function(x) {
        x[['gene']]
      }
    )
    return(list('nde' = nde,
                'd_genes' = de_genes))
  }
)

# Generate Venn diagrams to show overlap 
de_venns <- vector(mode = 'list', length = length(de_bycell_summary))
names(de_venns) <- names(de_bycell_summary)
for (i in 1:length(de_bycell_summary)) {
  x <- de_bycell_summary[[i]]
  standard <- grepl(pattern = 'Standard',
                    x = names(x[['d_genes']]))
  standard <- Reduce(f = union, x = x[['d_genes']][standard])
  enriched <- grepl(pattern = 'Enriched',
                    x = names(x[['d_genes']]))
  enriched <- Reduce(f = union, x = x[['d_genes']][enriched])
  venn <- eulerr::euler(
    combinations = list(
      'Standard' = standard,
      'Enriched' = enriched
    )
  )
  venn <- plot(
    x = venn,
    quantities = list(fontsize = 13),
    fills = c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)),
    main = paste('Uninjured vs 1dpi', names(de_bycell_summary)[i]), 
    labels = list(font = 1, fontsize = 14)
  )
  de_venns[[i]] <- venn
}

# save venn diagrams
tiff(filename = paste0(results_out, 'DissociationComparison_DE.tiff'),
     height = 8, width = 8*1.9, res = 440, units = 'in')
gridExtra::grid.arrange(grobs = de_venns)
dev.off()

# Save venn diagrams of cell-types present in uninjured cord (ie remove cells
# that were enrichment-specific and peripheral immune cells)
keep <- c('Microglia','Pericyte','OPC','Ependymal','Endothelial')
tiff(filename = paste0(results_out, 'DissociationComparison_DE_subset.tiff'),
     height = 5, width = 10, res = 440, units = 'in')
gridExtra::grid.arrange(
  grobs = de_venns[keep], 
  ncol = ceiling(sqrt(length(keep)))
)
dev.off()







# Uninjured vs 1dpi downsampled DE by dissociation method ----------------------

# Setup data. First, subset seurat object with RNA assay for DE testing.
Idents(sci) <- 'celltype'
sci_lite <- DietSeurat(
  object = sci,
  assays = 'RNA',
  graphs = NULL,
  dimreducs = NULL
)

# We will test for DE between times, for each dissociation method and cell-type.
celltypes <- levels(sci@meta.data[['celltype']])
diss <- unique(sci@meta.data[['dissociationMethod']])

# Iterate through cell-types, dissociation method.
# !! We downsample the number of cells to have equal group sizes.
de_bycell_downsample <- vector(mode = 'list', length = length(celltypes))
names(de_bycell_downsample) <- celltypes
for (ct in celltypes) {
  dat <- sci_lite[,sci_lite$time %in% c('Uninjured','1dpi') & sci_lite$celltype == ct]
  # If not enough cells per group, cant do DE
  if (any(table(dat$dissociationMethod, dat$time)[,1:2] < 3)) {
    next
  }
  de_bydiss <- vector(mode = 'list', length = length(diss))
  names(de_bydiss) <- diss
  for (m in diss) {
    max_cells <- min(table(dat$dissociationMethod, dat$time)[m,1:2])
    Idents(dat) <- 'dissociationMethod'
    tmp_de <- FindMarkers(
      object = dat,
      assay = 'RNA',
      slot = 'data',
      subset.ident = m,
      group.by = 'time',
      ident.1 = 'Uninjured',
      ident.2 = '1dpi',
      max.cells.per.ident = max_cells
    )
    de_bydiss[[m]] <- tmp_de
    de_bydiss[[m]][['gene']] <- rownames(de_bydiss[[m]])
    de_bydiss[[m]][['d_method']] <- m
  }
  de_bydiss <- Reduce(f = rbind, x = de_bydiss)
  de_bycell_downsample[[ct]] <- de_bydiss
}
saveRDS(object = de_bycell_downsample,
        file = paste0(results_out, 'results_DE_wilcox_downsampled_uninj_1dpi_byCell_byDissociationMethod.rds'))

# Retrieve summary of results: 1) table of up or down-regulated genes b/w uninj
# and 1dpi for each dissociation method, per cell-type, 2) names of genes that
# were up or down-reg for each comparison.
de_bycell_downsample_summary <- lapply(
  X = de_bycell_downsample,
  FUN = function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    tmp <- x %>%
      filter(p_val_adj < 1e-03) %>% # Filter low p-values
      mutate('dir' = ifelse(test = avg_logFC > 0,
                            yes = 'up',
                            no = 'down'))
    nde <- tmp %>%
      group_by(d_method, dir) %>%
      summarise('nde' = n_distinct(gene))
    de_genes <- split(
      x = tmp,
      f = list(tmp[['d_method']], tmp[['dir']])
    )
    de_genes <- lapply(
      X = de_genes,
      FUN = function(x) {
        x[['gene']]
      }
    )
    return(list('nde' = nde,
                'd_genes' = de_genes))
  }
)

# Generate Venn diagrams to show overlap 
de_venns <- vector(mode = 'list', length = length(de_bycell_downsample_summary))
names(de_venns) <- names(de_bycell_downsample_summary)
for (i in 1:length(de_bycell_downsample_summary)) {
  x <- de_bycell_downsample_summary[[i]]
  standard <- grepl(pattern = 'Standard',
                    x = names(x[['d_genes']]))
  standard <- Reduce(f = union, x = x[['d_genes']][standard])
  enriched <- grepl(pattern = 'Enriched',
                    x = names(x[['d_genes']]))
  enriched <- Reduce(f = union, x = x[['d_genes']][enriched])
  venn <- eulerr::euler(
    combinations = list(
      'Standard' = standard,
      'Enriched' = enriched
    )
  )
  venn <- plot(
    x = venn,
    quantities = list(fontsize = 13),
    fills = c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)),
    main = paste('Uninjured vs 1dpi', names(de_bycell_downsample_summary)[i]), 
    labels = list(font = 1, fontsize = 14)
  )
  de_venns[[i]] <- venn
}

# save venn diagrams
tiff(filename = paste0(results_out, 'DissociationComparison_downsampled_DE.tiff'),
     height = 8, width = 8*1.9, res = 440, units = 'in')
gridExtra::grid.arrange(grobs = de_venns)
dev.off()

# Save venn diagrams of cell-types present in uninjured cord (ie remove cells
# that were enrichment-specific and peripheral immune cells)
keep <- c('Microglia','Pericyte','OPC','Ependymal','Endothelial')
tiff(filename = paste0(results_out, 'DissociationComparison_downsampled_DE_subset.tiff'),
     height = 5, width = 10, res = 440, units = 'in')
gridExtra::grid.arrange(
  grobs = de_venns[keep], 
  ncol = ceiling(sqrt(length(keep)))
)
dev.off()




# Uninjured vs 1dpi MAST-detection rate DE by dissociation method --------------

# Setup data. First, subset seurat object with RNA assay for DE testing.
Idents(sci) <- 'celltype'
sci_lite <- DietSeurat(
  object = sci,
  assays = 'RNA',
  graphs = NULL,
  dimreducs = NULL
)

# We will test for DE between times, for each dissociation method and cell-type.
celltypes <- levels(sci@meta.data[['celltype']])
diss <- unique(sci@meta.data[['dissociationMethod']])

# Iterate through cell-types, dissociation method.
# !! We downsample the number of cells to have equal group sizes.
de_bycell_mast <- vector(mode = 'list', length = length(celltypes))
names(de_bycell_mast) <- celltypes
for (ct in celltypes) {
  dat <- sci_lite[,sci_lite$time %in% c('Uninjured','1dpi') & sci_lite$celltype == ct]
  # If not enough cells per group, cant do DE
  if (any(table(dat$dissociationMethod, dat$time)[,1:2] < 3)) {
    next
  }
  de_bydiss <- vector(mode = 'list', length = length(diss))
  names(de_bydiss) <- diss
  for (m in diss) {
    max_cells <- min(table(dat$dissociationMethod, dat$time)[m,1:2])
    Idents(dat) <- 'dissociationMethod'
    tmp_de <- FindMarkers(
      object = dat,
      assay = 'RNA',
      slot = 'data',
      subset.ident = m,
      group.by = 'time',
      ident.1 = 'Uninjured',
      ident.2 = '1dpi',
      test.use = 'MAST',
      latent.vars = 'nFeature_RNA'
    )
    de_bydiss[[m]] <- tmp_de
    de_bydiss[[m]][['gene']] <- rownames(de_bydiss[[m]])
    de_bydiss[[m]][['d_method']] <- m
  }
  de_bydiss <- Reduce(f = rbind, x = de_bydiss)
  de_bycell_mast[[ct]] <- de_bydiss
}
saveRDS(object = de_bycell_mast,
        file = paste0(results_out, 'results_DE_MAST_nFeature_uninj_1dpi_byCell_byDissociationMethod.rds'))

# Retrieve summary of results: 1) table of up or down-regulated genes b/w uninj
# and 1dpi for each dissociation method, per cell-type, 2) names of genes that
# were up or down-reg for each comparison.
de_bycell_mast_summary <- lapply(
  X = de_bycell_mast,
  FUN = function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    tmp <- x %>%
      filter(p_val_adj < 1e-03) %>% # Filter low p-values
      mutate('dir' = ifelse(test = avg_logFC > 0,
                            yes = 'up',
                            no = 'down'))
    nde <- tmp %>%
      group_by(d_method, dir) %>%
      summarise('nde' = n_distinct(gene))
    de_genes <- split(
      x = tmp,
      f = list(tmp[['d_method']], tmp[['dir']])
    )
    de_genes <- lapply(
      X = de_genes,
      FUN = function(x) {
        x[['gene']]
      }
    )
    return(list('nde' = nde,
                'd_genes' = de_genes))
  }
)

# Generate Venn diagrams to show overlap 
de_venns <- vector(mode = 'list', length = length(de_bycell_mast_summary))
names(de_venns) <- names(de_bycell_mast_summary)
for (i in 1:length(de_bycell_mast_summary)) {
  x <- de_bycell_mast_summary[[i]]
  standard <- grepl(pattern = 'Standard',
                    x = names(x[['d_genes']]))
  standard <- Reduce(f = union, x = x[['d_genes']][standard])
  enriched <- grepl(pattern = 'Enriched',
                    x = names(x[['d_genes']]))
  enriched <- Reduce(f = union, x = x[['d_genes']][enriched])
  venn <- eulerr::euler(
    combinations = list(
      'Standard' = standard,
      'Enriched' = enriched
    )
  )
  venn <- plot(
    x = venn,
    quantities = list(fontsize = 13),
    fills = c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)),
    main = paste('Uninjured vs 1dpi', names(de_bycell_mast_summary)[i]), 
    labels = list(font = 1, fontsize = 14)
  )
  de_venns[[i]] <- venn
}

# save venn diagrams
tiff(filename = paste0(results_out, 'DissociationComparison_MAST_DE.tiff'),
     height = 8, width = 8*1.9, res = 440, units = 'in')
gridExtra::grid.arrange(grobs = de_venns)
dev.off()

# Save venn diagrams of cell-types present in uninjured cord (ie remove cells
# that were enrichment-specific and peripheral immune cells)
keep <- c('Microglia','Pericyte','OPC','Ependymal','Endothelial')
tiff(filename = paste0(results_out, 'DissociationComparison_MAST_DE_subset.tiff'),
     height = 5, width = 10, res = 440, units = 'in')
gridExtra::grid.arrange(
  grobs = de_venns[keep], 
  ncol = ceiling(sqrt(length(keep)))
)
dev.off()







# Focus on Microglia ------------------------------------------------------


# Setup data. First, subset seurat object with RNA assay for DE testing.
Idents(sci) <- 'celltype'
sci_lite <- DietSeurat(
  object = sci,
  assays = 'RNA',
  graphs = NULL,
  dimreducs = NULL
)

standard_control_de <- FindMarkers(
  object = sci_lite,
  assay = 'RNA',
  slot = 'data',
  test.use = 'wilcox',
  subset.ident = 'Microglia',
  group.by = 'sample_id',
  ident.1 = 'uninj_sample1',
  ident.2 = 'uninj_sample2'
)
standard_control_de_genes <- standard_control_de %>%
  filter(p_val_adj < 1e-03) %>%
  rownames(.)

enriched_control_de <- FindMarkers(
  object = sci_lite,
  assay = 'RNA',
  slot = 'data',
  test.use = 'wilcox',
  subset.ident = 'Microglia',
  group.by = 'sample_id',
  ident.1 = '1dpi_sample2',
  ident.2 = '1dpi_sample3'
)
enriched_control_de_genes <- enriched_control_de %>%
  filter(p_val_adj < 1e-03) %>%
  rownames(.)


de_bycell <- readRDS(file = paste0(results_out, 'results_DE_wilcox_uninj_1dpi_byCell_byDissociationMethod.rds'))

# Retrieve summary of results: 1) table of up or down-regulated genes b/w uninj
# and 1dpi for each dissociation method, per cell-type, 2) names of genes that
# were up or down-reg for each comparison.
de_bycell_summary <- lapply(
  X = de_bycell,
  FUN = function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    tmp <- x %>%
      filter(p_val_adj < 1e-03) %>%
      mutate('dir' = ifelse(test = avg_logFC > 0,
                            yes = 'up',
                            no = 'down'))
    nde <- tmp %>%
      group_by(d_method, dir) %>%
      summarise('nde' = n_distinct(gene))
    de_genes <- split(
      x = tmp,
      f = list(tmp[['d_method']], tmp[['dir']])
    )
    de_genes <- lapply(
      X = de_genes,
      FUN = function(x) {
        x[['gene']]
      }
    )
    return(list('nde' = nde,
                'd_genes' = de_genes))
  }
)


mg_venn <- {
  x <- de_bycell_summary[['Microglia']]
  standard <- grepl(pattern = 'Standard',
                    x = names(x[['d_genes']]))
  standard <- Reduce(f = union, x = x[['d_genes']][standard])
  enriched <- grepl(pattern = 'Enriched',
                    x = names(x[['d_genes']]))
  enriched <- Reduce(f = union, x = x[['d_genes']][enriched])
  venn <- eulerr::euler(
    combinations = list(
      'Standard' = standard,
      'Enriched' = enriched,
      'Standard_ctrl' = standard_control_de_genes,
      'Enriched_ctrl' = enriched_control_de_genes
    )
  )
  venn <- plot(
    x = venn,
    quantities = list(fontsize = 13),
    fills = c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)),
    main = paste('Uninjured vs 1dpi', 'Microglia'), 
    labels = list(font = 1, fontsize = 14)
  )
}

mg_upset <- {
  x <- de_bycell_summary[['Microglia']]
  standard <- grepl(pattern = 'Standard',
                    x = names(x[['d_genes']]))
  standard <- Reduce(f = union, x = x[['d_genes']][standard])
  enriched <- grepl(pattern = 'Enriched',
                    x = names(x[['d_genes']]))
  enriched <- Reduce(f = union, x = x[['d_genes']][enriched])
  
  comb_mat <- ComplexHeatmap::make_comb_mat(
    list(
      'Standard' = standard,
      'Enriched' = enriched,
      'Standard_ctrl' = standard_control_de_genes,
      'Enriched_ctrl' = enriched_control_de_genes
    ),
    mode = 'intersect'
  )
  ComplexHeatmap::UpSet(m = comb_mat)
}

tiff(filename = paste0(results_out, 'microglia_DE_AllComparisons.tiff'),
     height = 4, width = 12, res = 440, units = 'in')
cowplot::plot_grid(
  cowplot::as_grob(mg_venn),
  grid::grid.grabExpr(ComplexHeatmap::draw(mg_upset)),
  ncol = 2,
  rel_widths = c(0.6,1)
)
dev.off()






# Uninjured vs 1dpi DE by original data vs appended ---------------------------

# For clarification, this is to specifically compare DE test results using 
# either the original data used for the manuscript or 
# whether the data used to
# generate the results in the original manuscript 
# presented in our original manuscr

# Setup data. First, subset seurat object with RNA assay for DE testing.
Idents(sci) <- 'celltype'
sci_lite <- DietSeurat(
  object = sci,
  assays = 'RNA',
  graphs = NULL,
  dimreducs = NULL
)

# We will test for DE between times, for each dissociation method and cell-type.
celltypes <- levels(sci@meta.data[['celltype']])
diss <- unique(sci@meta.data[['dissociationMethod']])

# Iterate through cell-types, dissociation method
de_bycell <- vector(mode = 'list', length = length(celltypes))
names(de_bycell) <- celltypes
for (ct in celltypes) {
  dat <- sci_lite[,sci_lite$time %in% c('Uninjured','1dpi') & sci_lite$celltype == ct]
  if (any(table(dat$dissociationMethod, dat$time)[,1:2] < 3)) {
    next
  }
  de_bydiss <- vector(mode = 'list', length = length(diss))
  names(de_bydiss) <- diss
  for (m in diss) {
    Idents(dat) <- 'dissociationMethod'
    tmp_de <- FindMarkers(
      object = dat,
      assay = 'RNA',
      slot = 'data',
      subset.ident = m,
      group.by = 'time',
      ident.1 = 'Uninjured',
      ident.2 = '1dpi'
    )
    de_bydiss[[m]] <- tmp_de
    de_bydiss[[m]][['gene']] <- rownames(de_bydiss[[m]])
    de_bydiss[[m]][['d_method']] <- m
  }
  de_bydiss <- Reduce(f = rbind, x = de_bydiss)
  de_bycell[[ct]] <- de_bydiss
}

# Retrieve summary of results: 1) table of up or down-regulated genes b/w uninj
# and 1dpi for each dissociation method, per cell-type, 2) names of genes that
# were up or down-reg for each comparison.
de_bycell_summary <- lapply(
  X = de_bycell,
  FUN = function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    tmp <- x %>%
      filter(p_val_adj < 1e-03) %>%
      mutate('dir' = ifelse(test = avg_logFC > 0,
                            yes = 'up',
                            no = 'down'))
    nde <- tmp %>%
      group_by(d_method, dir) %>%
      summarise('nde' = n_distinct(gene))
    de_genes <- split(
      x = tmp,
      f = list(tmp[['d_method']], tmp[['dir']])
    )
    de_genes <- lapply(
      X = de_genes,
      FUN = function(x) {
        x[['gene']]
      }
    )
    return(list('nde' = nde,
                'd_genes' = de_genes))
  }
)

# Generate Venn diagrams to show overlap 
de_venns <- vector(mode = 'list', length = length(de_bycell_summary))
names(de_venns) <- names(de_bycell_summary)
for (i in 1:length(de_bycell_summary)) {
  x <- de_bycell_summary[[i]]
  standard <- grepl(pattern = 'Standard',
                    x = names(x[['d_genes']]))
  standard <- Reduce(f = union, x = x[['d_genes']][standard])
  enriched <- grepl(pattern = 'Enriched',
                    x = names(x[['d_genes']]))
  enriched <- Reduce(f = union, x = x[['d_genes']][enriched])
  venn <- eulerr::euler(
    combinations = list(
      'Standard' = standard,
      'Enriched' = enriched
    )
  )
  venn <- plot(
    x = venn,
    quantities = list(fontsize = 13),
    fills = c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)),
    main = paste('Uninjured vs 1dpi', names(de_bycell_summary)[i]), 
    labels = list(font = 1, fontsize = 14)
  )
  de_venns[[i]] <- venn
}
tiff(filename = paste0(results_out, 'DissociationComparison_uninj_1dpi_DE.tiff'),
     height = 8, width = 8*1.9, res = 440, units = 'in')
gridExtra::grid.arrange(grobs = de_venns)
dev.off()














# scratch -----------------------------------------------------------------


c <- 5
c_dat <- sci[,sci$celltype == celltypes[c]]

Idents(c_dat) <- 'sample_id'
uninj1_1dpi1 <- FindMarkers(
  object = c_dat,
  ident.1 = c('uninj_sample1'),
  ident.2 = c('1dpi_sample1'),
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  max.cells.per.ident = min(sum(c_dat$sample_id == 'uninj_sample1'),
                            sum(c_dat$sample_id == '1dpi_sample1')),
  random.seed = 123
)
uninj2_1dpi1 <- FindMarkers(
  object = c_dat,
  ident.1 = c('uninj_sample2'),
  ident.2 = c('1dpi_sample1'),
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  max.cells.per.ident = min(sum(c_dat$sample_id == 'uninj_sample2'),
                            sum(c_dat$sample_id == '1dpi_sample1')),
  random.seed = 123
)
uninj3_1dpi2 <- FindMarkers(
  object = c_dat,
  ident.1 = c('uninj_sample3'),
  ident.2 = c('1dpi_sample2'),
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  max.cells.per.ident = min(sum(c_dat$sample_id == 'uninj_sample3'),
                            sum(c_dat$sample_id == '1dpi_sample2')),
  random.seed = 123
)
uninj3_1dpi3 <- FindMarkers(
  object = c_dat,
  ident.1 = c('uninj_sample3'),
  ident.2 = c('1dpi_sample3'),
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  max.cells.per.ident = min(sum(c_dat$sample_id == 'uninj_sample3'),
                            sum(c_dat$sample_id == '1dpi_sample3')),
  random.seed = 123
)

get_gene_order <- function(x, pval_thresh) {
  tmp <- x %>%
    filter(p_val_adj < pval_thresh) %>%
    tibble::rownames_to_column(var = 'gene') %>%
    select(gene, p_val_adj)
  return(tmp)
}
set1 <- get_gene_order(uninj1_1dpi1, 1e-10)
colnames(set1) <- c('gene','set1')
set2 <- get_gene_order(uninj2_1dpi1, 1e-10)
colnames(set2) <- c('gene','set2')
set3 <- get_gene_order(uninj3_1dpi2, 1e-10)
colnames(set3) <- c('gene','set3')
set4 <- get_gene_order(uninj3_1dpi3, 1e-10)
colnames(set4) <- c('gene','set4')

test1 <- eulerr::euler(list('uninj1_1dpi1' = set1$gene,
                            'uninj2_1dpi1' = set2$gene,
                            'uninj3_1dpi2' = set3$gene,
                            'uninj3_1dpi3' = set4$gene))
plot(test1)



uninjA_1dpiA <- FindMarkers(
  object = c_dat,
  ident.1 = c('uninj_sample1','uninj_sample2'),
  ident.2 = c('1dpi_sample1'),
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  random.seed = 123
)
uninjB_1dpiB <- FindMarkers(
  object = c_dat,
  ident.1 = c('uninj_sample3'),
  ident.2 = c('1dpi_sample2','1dpi_sample3'),
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  random.seed = 123
)
uninj_1dpi <- FindMarkers(
  object = c_dat,
  ident.1 = c('uninj_sample1','uninj_sample2','uninj_sample3'),
  ident.2 = c('1dpi_sample1','1dpi_sample2','1dpi_sample3'),
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  random.seed = 123
)
setA <- get_gene_order(uninjA_1dpiA, 1e-10)
colnames(setA) <- c('gene','setA')
setB <- get_gene_order(uninjB_1dpiB, 1e-10)
colnames(setB) <- c('gene','setB')
setC <- get_gene_order(uninj_1dpi, 1e-10)
colnames(setC) <- c('gene','setC')
test2 <- eulerr::euler(list('setA' = setA$gene,
                            'setB' = setB$gene,
                            'setC' = setC$gene))
plot(test2)





gplots::venn(data = list('uninj3_1dpi2' = set3$gene,
                         'uninj3_1dpi3' = set4$gene))
gplots::venn(data = list('uninj1_1dpi1' = set1$gene,
                         'uninj2_1dpi1' = set2$gene,
                         'uninj3_1dpi2' = set3$gene,
                         'uninj3_1dpi3' = set4$gene))

set3 <- get_gene_order(tmp3, 1e-50)
colnames(set3) <- c('gene','set3')
set <- Reduce(f = function(x, y) merge(x, y, all.x = TRUE), x = list(set1, set2, set3))
cor(x = set$set1, set$set2, method = 'spearman', use = 'complete.obs')
cor(x = set$set2, set$set3, method = 'spearman', use = 'complete.obs')
cor(x = set$set1, set$set3, method = 'spearman', use = 'complete.obs')
par(mfrow = c(1,3))
gplots::venn(data = list('set1' = set$gene[!is.na(set$set1)], 'set2' = set$gene[!is.na(set$set2)]))
gplots::venn(data = list('set2' = set$gene[!is.na(set$set2)], 'set3' = set$gene[!is.na(set$set3)]))
gplots::venn(data = list('set1' = set$gene[!is.na(set$set1)], 'set3' = set$gene[!is.na(set$set3)]))








sci <- readRDS(file = './data/sci.rds')
DefaultAssay(sci) <- 'RNA'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL


rep_12 <- c('uninj_sample1','uninj_sample2','1dpi_sample1','1dpi_sample2')
rep_23 <- c('uninj_sample2','uninj_sample3','1dpi_sample2','1dpi_sample3')
sci_12 <- sci[,sci$sample_id %in% rep_12]
sci_23 <- sci[,sci$sample_id %in% rep_23]


mg_de_12 <- FindMarkers(
  object = sci_12,
  subset.ident = 'Microglia',
  group.by = 'time',
  ident.1 = '1dpi',
  ident.2 = 'Uninjured',
  assay = 'RNA',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA'
)
mg_de_23 <- FindMarkers(
  object = sci_23,
  subset.ident = 'Microglia',
  group.by = 'time',
  ident.1 = '1dpi',
  ident.2 = 'Uninjured',
  assay = 'RNA',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA'
)


mg_12_nde <- mg_de_12 %>%
  filter(p_val_adj < 1e-05) %>%
  rownames(.)
mg_23_nde <- mg_de_23 %>%
  filter(p_val_adj < 1e-05) %>%
  rownames(.)

length(intersect(mg_12_nde, mg_23_nde))
length(setdiff(mg_12_nde, mg_23_nde))
gplots::venn(data = list(mg_12_nde, mg_23_nde))


mg_12_nde <- mg_de_12 %>%
  filter(p_val_adj < 1e-05) %>%
  top_n(n = 100, wt = avg_logFC)
x <- mg_12_nde[['p_val_adj']]
names(x) <- rownames(mg_12_nde)
x <- x[order(names(x))]

mg_12_de <- rownames(mg_12_nde)
mg_23_nde <- mg_de_23 %>%
  filter(p_val_adj < 1e-05) %>%
  top_n(n = 100, wt = avg_logFC) %>%
  .[['avg_logFC']]


cor(x = mg_12_nde, y = mg_23_nde, method = 'spearman')




