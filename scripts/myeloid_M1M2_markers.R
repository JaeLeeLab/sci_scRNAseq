
######## Myeloid expression of M1/M2 genes ########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
require('ComplexHeatmap')
# data_in <- './data/'
# data_out <- './data/data_integration/'
results_out <- './results/myeloid_M1M2_markers/'
ref_in <- './ref/'
ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)

myeloid <- readRDS(file = './data/myeloid.rds')





# M1 M2 marker heatmap  ---------------------------------------------------


DefaultAssay(myeloid) <- 'RNAcorrected'
Idents(myeloid) <- 'myeloid_subcluster'

# Preset values
myeloid_cols <- c('Neutrophil' = '#800000',    
                  'Monocyte' = '#9a6324',
                  'Macrophage-A' = '#e6194b',
                  'Macrophage-B' = '#f58231',
                  'BA-Macrophage' = '#CCCC00',
                  'Dendritic' = '#808000',
                  'Div-Myeloid' = '#3cb44b',
                  'Microglia-A' = '#008080',
                  'Microglia-B' = 'cyan3',
                  'Microglia-C' = '#000075',
                  'Div-Microglia' = '#4363d8',
                  'IFN-Myeloid' = '#911eb4')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')
m1m2_cols <- c('M1' = 'indianred', 'M2' = 'dodgerblue')


# https://www.annualreviews.org/doi/full/10.1146/annurev-physiol-022516-034339
# M1/M2 gene expression marker panel derived from Table 2 of above.
m1_markers <- c('Il1a','Il1b','Il6','Il12a','Il12b','Il23a','Il27','Tnf','Csf3','Csf2','Nfkbiz','Ccl1','Cxcl13','Ccl11','Cxcl2','Tnfaip3','Socs3','Peli1','Nos2','Marco')
m2_markers <- c('Retnla','Clec10a','Ccl17','Ccl24','Irf4','Chil3','Mrc1','Arg1','Rnase2a','Ear2','Ccl8','Pdcd1lg2','Socs2','Cdh1','Ppard','Pparg','Ccl22')


# Select which genes to annotate
label_genes <- c(m1_markers[1:5], m2_markers[1:5])

# Extract scaled expression values
expr_data <- t(ScaleData(
  object = myeloid[['RNAcorrected']]@data, 
  features = c(m1_markers, m2_markers),
  scale.max = 3
))

# metadata
myeloid_meta <- myeloid@meta.data[c('myeloid_subcluster')]

# Merge metadata and exprs values, then calculate avg across cell/time
m1m2_data <- cbind(expr_data, myeloid_meta) %>%
  group_by(myeloid_subcluster) %>%
  summarise(across(.fns = mean)) %>%
  tibble::column_to_rownames(var = 'myeloid_subcluster') %>%
  as.matrix()


# Cell-level annotations
m1m2_cell_anno <- rowAnnotation(
  'Subcluster' = rownames(m1m2_data),
  col = list('Subcluster' = myeloid_cols)
)

# Gene-level annotations
m1m2_gene_anno <- HeatmapAnnotation(
  'M1/M2' = c(rep('M1', length(m1_markers)), rep('M2', length(m2_markers))),
  col = list('M1/M2' = m1m2_cols),
  annotation_legend_param = list(
    'M1/M2' = list(
      direction = 'horizontal'
    )
  )
  # 'M1/M2' = anno_mark(
  #   at = match(label_genes, colnames(m1m2_data)),
  #   labels = label_genes,
  #   side = 'right',
  #   labels_gp = gpar(fontsize = 10)
  # )
)

# generate heatmap
m1m2_heatmp <- Heatmap(
  matrix = m1m2_data,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 1.75)(100),
  rect_gp = gpar(col = 'black', lwd = 0.5),
  heatmap_legend_param = list(
    title = 'Scaled Expression',
    title_gp = gpar(fontsize = 12),
    title_position = 'topcenter',
    labels = c('Low',0, 'High'),
    at = c(min(m1m2_data), 0, max(m1m2_data)),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(2.5, units = 'cm'),
    grid_width = unit(0.5, units = 'cm'),
    border = 'black',
    title_gap = unit(3, units = 'cm'),
    column_gap = unit(5, "mm"),
    row_gap = unit(5, 'mm'),
    direction = 'horizontal'),
  use_raster = TRUE,
  top_annotation = m1m2_gene_anno,
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  clustering_distance_columns = 'spearman',
  clustering_distance_rows = 'spearman',
  # right_annotation = m1m2_cell_anno,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10)
)

tiff(filename = paste0(results_out, 'myeloid_M1M2_heatmap_joint.tiff'),
     height = 4, width = 10, units = 'in', res = 440)
draw(m1m2_heatmp, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()






# M1 marker only panel ----------------------------------------------------


# M1 genes
m1_markers <- c('Il1a','Il1b','Il6','Il12a','Il12b','Il23a','Il27','Tnf','Csf3','Csf2','Nfkbiz','Ccl1','Cxcl13','Ccl11','Cxcl2','Tnfaip3','Socs3','Peli1','Nos2','Marco')

# Extract scaled expression values
expr_data <- t(ScaleData(
  object = myeloid[['RNAcorrected']]@data, 
  features = c(m1_markers),
  scale.max = 3
))

# metadata
myeloid_meta <- myeloid@meta.data[c('myeloid_subcluster')]

# Merge metadata and exprs values, then calculate avg across cell/time
m1_data <- cbind(expr_data, myeloid_meta) %>%
  group_by(myeloid_subcluster) %>%
  summarise(across(.fns = mean)) %>%
  tibble::column_to_rownames(var = 'myeloid_subcluster') %>%
  as.matrix()


# generate heatmap
m1_heatmap <- Heatmap(
  matrix = m1_data,
  column_title = 'M1 gene expression',
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 1.3)(100),
  rect_gp = gpar(col = 'black', lwd = 0.5),
  heatmap_legend_param = list(
    title = 'Scaled Expression',
    title_gp = gpar(fontsize = 12),
    title_position = 'topcenter',
    labels = c('Low',0, 'High'),
    at = c(min(m1_data), 0, max(m1_data)),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(2.5, units = 'cm'),
    grid_width = unit(0.5, units = 'cm'),
    border = 'black',
    title_gap = unit(3, units = 'cm'),
    column_gap = unit(5, "mm"),
    row_gap = unit(5, 'mm'),
    direction = 'horizontal'),
  use_raster = TRUE,
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  clustering_distance_columns = 'spearman',
  clustering_distance_rows = 'spearman',
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 12)
)

# save
tiff(filename = paste0(results_out, 'myeloid_M1_heatmap.tiff'),
     height = 4.5, width = 6, units = 'in', res = 440)
draw(m1_heatmap, heatmap_legend_side = 'bottom')
dev.off()




# M1 marker only panel (functional annotation) ---------------------------------


myeloid@meta.data[['myeloid_functional']] <- plyr::mapvalues(
  x = myeloid@meta.data[['integrated_snn_res.0.35']],
  from = levels(myeloid@meta.data[['integrated_snn_res.0.35']]),
  to = c('Homeostatic Microglia',
         'Inflammatory Microglia',
         'Chemotaxis-Inducing Mac',
         'Dividing Microglia',
         'Inflammatory Mac',
         'Monocyte',
         'Neutrophil',
         'Migrating Microglia',
         'Dendritic',
         'Interferon Myeloid',
         'Dividing Myeloid',
         'Border-Associated Mac')
)
myeloid@meta.data[['myeloid_functional']] <- factor(
  x = myeloid@meta.data[['myeloid_functional']],
  levels = c('Neutrophil',
             'Monocyte',
             'Chemotaxis-Inducing Mac',
             'Inflammatory Mac',
             'Border-Associated Mac',
             'Dendritic',
             'Dividing Myeloid',
             'Homeostatic Microglia',
             'Inflammatory Microglia',
             'Dividing Microglia',
             'Migrating Microglia',
             'Interferon Myeloid')
)

# M1 genes
m1_markers <- c('Il1a','Il1b','Il6','Il12a','Il12b','Il23a','Il27','Tnf','Csf3','Csf2','Nfkbiz','Ccl1','Cxcl13','Ccl11','Cxcl2','Tnfaip3','Socs3','Peli1','Nos2','Marco')

# Extract scaled expression values
expr_data <- t(ScaleData(
  object = myeloid[['RNAcorrected']]@data, 
  features = c(m1_markers),
  scale.max = 3
))

# metadata
myeloid_meta <- myeloid@meta.data[c('myeloid_functional')]

# Merge metadata and exprs values, then calculate avg across cell/time
m1_data <- cbind(expr_data, myeloid_meta) %>%
  group_by(myeloid_functional) %>%
  summarise(across(.fns = mean)) %>%
  tibble::column_to_rownames(var = 'myeloid_functional') %>%
  as.matrix()

Idents(myeloid) <- 'myeloid_functional'
DefaultAssay(myeloid) <- 'RNA'
tmp <- FindAllMarkers(
  object = myeloid,
  features = m1_markers,
  only.pos = TRUE,
  assay = 'RNA',
  slot = 'data',
  logfc.threshold = 0,
  min.pct = 0
)

# generate heatmap
m1_heatmap <- Heatmap(
  matrix = m1_data,
  column_title = 'M1 gene expression',
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
  rect_gp = gpar(col = 'black', lwd = 0.5),
  heatmap_legend_param = list(
    title = 'z-score',
    title_gp = gpar(fontsize = 12),
    title_position = 'topcenter',
    labels = c('Low',0, 'High'),
    at = c(-max(m1_data), 0, max(m1_data)),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(2.5, units = 'cm'),
    grid_width = unit(0.5, units = 'cm'),
    border = 'black',
    title_gap = unit(3, units = 'cm'),
    column_gap = unit(5, "mm"),
    row_gap = unit(5, 'mm'),
    direction = 'horizontal'),
  use_raster = TRUE,
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  clustering_distance_columns = 'spearman',
  clustering_distance_rows = 'spearman',
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 12)
)

# save
tiff(filename = './results/revision_figures/functional_names/myeloid_M1_heatmap.tiff',
     height = 4.5, width = 6.75, units = 'in', res = 440)
draw(m1_heatmap, heatmap_legend_side = 'bottom')
dev.off()







# M2 marker only panel ----------------------------------------------------

# M2 genes
m2_markers <- c('Retnla','Clec10a','Ccl17','Ccl24','Irf4','Chil3','Mrc1','Arg1','Rnase2a','Ear2','Ccl8','Pdcd1lg2','Socs2','Cdh1','Ppard','Pparg','Ccl22')

# Extract scaled expression values
expr_data <- t(ScaleData(
  object = myeloid[['RNAcorrected']]@data, 
  features = c(m2_markers),
  scale.max = 3
))

# metadata
myeloid_meta <- myeloid@meta.data[c('myeloid_subcluster')]

# Merge metadata and exprs values, then calculate avg across cell/time
m2_data <- cbind(expr_data, myeloid_meta) %>%
  group_by(myeloid_subcluster) %>%
  summarise(across(.fns = mean)) %>%
  tibble::column_to_rownames(var = 'myeloid_subcluster') %>%
  as.matrix()

# generate heatmap
m2_heatmap <- Heatmap(
  matrix = m2_data,
  column_title = 'M2 gene expression',
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 2.2)(100),
  rect_gp = gpar(col = 'black', lwd = 0.5),
  heatmap_legend_param = list(
    title = 'Scaled Expression',
    title_gp = gpar(fontsize = 12),
    title_position = 'topcenter',
    labels = c('Low',0, 'High'),
    at = c(min(m2_data), 0, max(m2_data)),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(2.5, units = 'cm'),
    grid_width = unit(0.5, units = 'cm'),
    border = 'black',
    title_gap = unit(3, units = 'cm'),
    column_gap = unit(5, "mm"),
    row_gap = unit(5, 'mm'),
    direction = 'horizontal'),
  use_raster = TRUE,
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  clustering_distance_columns = 'spearman',
  clustering_distance_rows = 'spearman',
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 12)
)

# save
tiff(filename = paste0(results_out, 'myeloid_M2_heatmap.tiff'),
     height = 4.5, width = 6, units = 'in', res = 440)
draw(m2_heatmap, heatmap_legend_side = 'bottom')
dev.off()



# M2 marker only panel (functional annotation) --------------------------------

myeloid@meta.data[['myeloid_functional']] <- plyr::mapvalues(
  x = myeloid@meta.data[['integrated_snn_res.0.35']],
  from = levels(myeloid@meta.data[['integrated_snn_res.0.35']]),
  to = c('Homeostatic Microglia',
         'Inflammatory Microglia',
         'Chemotaxis-Inducing Mac',
         'Dividing Microglia',
         'Inflammatory Mac',
         'Monocyte',
         'Neutrophil',
         'Migrating Microglia',
         'Dendritic',
         'Interferon Myeloid',
         'Dividing Myeloid',
         'Border-Associated Mac')
)
myeloid@meta.data[['myeloid_functional']] <- factor(
  x = myeloid@meta.data[['myeloid_functional']],
  levels = c('Neutrophil',
             'Monocyte',
             'Chemotaxis-Inducing Mac',
             'Inflammatory Mac',
             'Border-Associated Mac',
             'Dendritic',
             'Dividing Myeloid',
             'Homeostatic Microglia',
             'Inflammatory Microglia',
             'Dividing Microglia',
             'Migrating Microglia',
             'Interferon Myeloid')
)


# M2 genes
m2_markers <- c('Retnla','Clec10a','Ccl17','Ccl24','Irf4','Chil3','Mrc1','Arg1','Rnase2a','Ear2','Ccl8','Pdcd1lg2','Socs2','Cdh1','Ppard','Pparg','Ccl22')

# Extract scaled expression values
expr_data <- t(ScaleData(
  object = myeloid[['RNAcorrected']]@data, 
  features = c(m2_markers),
  scale.max = 3
))

# metadata
myeloid_meta <- myeloid@meta.data[c('myeloid_functional')]

# Merge metadata and exprs values, then calculate avg across cell/time
m2_data <- cbind(expr_data, myeloid_meta) %>%
  group_by(myeloid_functional) %>%
  summarise(across(.fns = mean)) %>%
  tibble::column_to_rownames(var = 'myeloid_functional') %>%
  as.matrix()

Idents(myeloid) <- 'myeloid_functional'
DefaultAssay(myeloid) <- 'RNA'
tmp <- FindAllMarkers(
  object = myeloid,
  features = m2_markers,
  only.pos = TRUE,
  assay = 'RNA',
  slot = 'data',
  logfc.threshold = 0,
  min.pct = 0
)

# generate heatmap
m2_heatmap <- Heatmap(
  matrix = m2_data,
  column_title = 'M2 gene expression',
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 1.7)(100),
  rect_gp = gpar(col = 'black', lwd = 0.5),
  heatmap_legend_param = list(
    title = 'z-score',
    title_gp = gpar(fontsize = 12),
    title_position = 'topcenter',
    labels = c('Low',0, 'High'),
    at = c(-max(m2_data), 0, max(m2_data)),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(2.5, units = 'cm'),
    grid_width = unit(0.5, units = 'cm'),
    border = 'black',
    title_gap = unit(3, units = 'cm'),
    column_gap = unit(5, "mm"),
    row_gap = unit(5, 'mm'),
    direction = 'horizontal'),
  use_raster = TRUE,
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  clustering_distance_columns = 'spearman',
  clustering_distance_rows = 'spearman',
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 12)
)

# save
tiff(filename = './results/revision_figures/functional_names/myeloid_M2_heatmap.tiff',
     height = 4.5, width = 6.75, units = 'in', res = 440)
draw(m2_heatmap, heatmap_legend_side = 'bottom')
dev.off()




# M1 M2 marker heatmap (split by time) ----------------------------------------


DefaultAssay(myeloid) <- 'RNAcorrected'
Idents(myeloid) <- 'myeloid_subcluster'

# Preset values
myeloid_cols <- c('Neutrophil' = '#800000',    
                  'Monocyte' = '#9a6324',
                  'Macrophage-A' = '#e6194b',
                  'Macrophage-B' = '#f58231',
                  'BA-Macrophage' = '#CCCC00',
                  'Dendritic' = '#808000',
                  'Div-Myeloid' = '#3cb44b',
                  'Microglia-A' = '#008080',
                  'Microglia-B' = 'cyan3',
                  'Microglia-C' = '#000075',
                  'Div-Microglia' = '#4363d8',
                  'IFN-Myeloid' = '#911eb4')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')
m1m2_cols <- c('M1' = 'indianred', 'M2' = 'dodgerblue')


# https://www.annualreviews.org/doi/full/10.1146/annurev-physiol-022516-034339
# M1/M2 gene expression marker panel derived from Table 2 of above.
m1_markers <- c('Il1a','Il1b','Il6','Il12a','Il12b','Il23a','Il27','Tnf','Csf3','Csf2','Nfkbiz','Ccl1','Cxcl13','Ccl11','Cxcl2','Tnfaip3','Socs3','Peli1','Nos2','Marco')
m2_markers <- c('Retnla','Clec10a','Ccl17','Ccl24','Irf4','Chil3','Mrc1','Arg1','Rnase2a','Ear2','Ccl8','Pdcd1lg2','Socs2','Cdh1','Ppard','Pparg','Ccl22')


# Select which genes to annotate
label_genes <- c(m1_markers[1:5], m2_markers[1:5])

# Extract scaled expression values
expr_data <- t(ScaleData(
  object = myeloid[['RNAcorrected']]@data, 
  features = c(m1_markers, m2_markers),
  scale.max = 2
))

# metadata
myeloid_meta <- myeloid@meta.data[c('myeloid_subcluster','time')]

# Merge metadata and exprs values, then calculate avg across cell/time
m1m2_data <- cbind(expr_data, myeloid_meta) %>%
  group_by(myeloid_subcluster, time) %>%
  summarise(across(.fns = mean))

# Remove groups with low counts (e.g. Uninjured Macrophage-B)
cell_count <- table(myeloid@meta.data[['myeloid_subcluster']],
                     myeloid@meta.data[['time']])
cell_prop <- prop.table(cell_count, margin = 1)
remove_cell <- cell_prop < 0.02
remove_rows <- c()
for (i in 1:nrow(m1m2_data)) {
  if (remove_cell[m1m2_data$myeloid_subcluster[i], m1m2_data$time[i]]) {
    remove_rows <- c(remove_rows, i)
  }
}
m1m2_data <- m1m2_data[-remove_rows,]

# Create new ID of celltype + time
m1m2_data[['id']] <- paste(m1m2_data[['myeloid_subcluster']],
                           m1m2_data[['time']],
                           sep = '_')

# Cell-level annotations
m1m2_cell_anno <- rowAnnotation(
  'Subcluster' = m1m2_data[['myeloid_subcluster']],
  'Time' = m1m2_data[['time']],
  col = list('Subcluster' = myeloid_cols,
             'Time' = time_cols)
)

# Gene-level annotations
m1m2_gene_anno <- HeatmapAnnotation(
  'M1/M2' = c(rep('M1', length(m1_markers)), rep('M2', length(m2_markers))),
  col = list('Genes' = m1m2_cols),
  'M1/M2' = anno_mark(
    at = match(label_genes, colnames(m1m2_data)),
    labels = label_genes,
    side = 'right',
    labels_gp = gpar(fontsize = 10)
  )
)

# convert to matrix
m1m2_mat <- m1m2_data[!names(m1m2_data) %in% c('myeloid_subcluster','time')] %>%
  tibble::column_to_rownames(var = 'id') %>%
  as.matrix()


# Generate heatmap
m1m2_split_heatmap <- Heatmap(
  matrix = m1m2_mat,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 1.6)(100),
  heatmap_legend_param = list(title = 'Scaled Expression',
                              title_gp = gpar(fontsize = 12),
                              title_position = 'leftcenter-rot',
                              labels = c('Low',0, 'High'),
                              at = c(min(expr_data), 0, max(expr_data)),
                              labels_gp = gpar(fontsize = 10),
                              legend_height = unit(2.5, units = 'cm'),
                              grid_width = unit(0.5, units = 'cm'),
                              border = 'black',
                              title_gap = unit(1, units = 'cm'),
                              direction = 'vertical'),
  use_raster = TRUE,
  # cluster_columns = FALSE,
  top_annotation = m1m2_gene_anno,
  right_annotation = m1m2_cell_anno
)

# save
tiff(filename = paste0(results_out, 'myeloid_M1M2_timeSplit_heatmap.tiff'),
     height = 8, width = 11, units = 'in', res = 440)
m1m2_split_heatmap
dev.off()


