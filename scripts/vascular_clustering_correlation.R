
######## vascular clustering correlation to assess experimental design #########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
require('ComplexHeatmap')
results_out <- './results/vascular_clustering_correlation/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)


vascular <- readRDS(file = './data/vascular.rds')


Idents(vascular) <- 'vascular_subcluster'
vascular_cols <- c('A-Endothelial' = '#800000',
                   'C1-Endothelial' = '#e6194b',
                   'C2-Endothelial' = '#f58231',
                   'V-Endothelial' = '#808000',
                   'Tip Cell' = '#3cb44b',
                   'Pericyte' = '#008080',
                   'VSMC' = '#000075',
                   'Fibroblast' = '#4363d8',
                   'U-Vascular' = '#f032e6')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')





# Correlation matrix by subcluster and replicate --------------------------------


# Calculate average expression of genes by subcluster + replicate. Calculate avg
# for the variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
vascular_avg_byReplicate <- AverageExpression(
  object = vascular,
  assays = 'RNA',
  features = vascular[['integrated']]@var.features,
  add.ident = 'sample_id',
  slot = 'data'
)
vascular_avg_byReplicate <- vascular_avg_byReplicate[[1]]

# Calculate spearman correlation between subcluster + replicate
vascular_cor_byReplicate <- cor(
  x = vascular_avg_byReplicate,
  method = 'spearman'
)

# Heatmap annotation (cell counts barplot and vascular_subcluster color)
tmp_count <- table('sample_id' = vascular$sample_id,
                   'vascular_subcluster' = vascular$vascular_subcluster) %>%
  reshape2::melt() %>%
  mutate('group' = paste(vascular_subcluster, sample_id, sep = '_'))
tmp_time <- unique(vascular@meta.data[c('sample_id','time')])
tmp_count[['time']] <- tmp_time$time[match(tmp_count$sample_id, tmp_time$sample_id)]
tmp_count <- tmp_count[tmp_count$value != 0,]
tmp_count_index <- match(x = colnames(vascular_cor_byReplicate), table = tmp_count$group)
tmp_count_val <- tmp_count$value[tmp_count_index]
tmp_vascular_subcluster <- tmp_count$vascular_subcluster[tmp_count_index]
tmp_time <- tmp_count$time[tmp_count_index]
vascular_cor_anno <- rowAnnotation(
  'Cell count' = anno_barplot(
    x = tmp_count_val,
    axis_param = list(
      side = 'bottom',
      labels_rot = 90,
      direction = 'reverse'
    )
  ),
  'Subcluster' = tmp_vascular_subcluster,
  'Time' = tmp_time,
  col = list('Subcluster' = vascular_cols,
             'Time' = time_cols),
  annotation_name_rot = 90,
  annotation_name_side = 'bottom'
)

# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'vascular_spearmanCor_byReplicate_RNA.tiff'),
     height = 10, width = 12, units = 'in', res = 480)
Heatmap(
  matrix = vascular_cor_byReplicate,
  column_title = 'Correlation by vascular subcluster and replicate',
  column_title_gp = gpar(fontsize = 28),
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 0.8)(100),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 18),
                              labels_gp = gpar(fontsize = 14),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(4, units = 'mm')),
  left_annotation = vascular_cor_anno
)
dev.off()




# Correlation matrix by subcluster and replicate (batchelor corrected) -----------------------


# Calculate average expression of genes by subcluster + replicate. Calculate avg
# for the variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
vascular_avg_byReplicate <- AverageExpression(
  object = vascular,
  assays = 'RNAcorrected',
  features = vascular[['integrated']]@var.features,
  add.ident = 'sample_id',
  slot = 'data'
)
vascular_avg_byReplicate <- vascular_avg_byReplicate[[1]]

# Calculate spearman correlation between subcluster + replicate
vascular_cor_byReplicate <- cor(
  x = vascular_avg_byReplicate,
  method = 'spearman'
)

# Heatmap annotation (cell counts barplot and vascular_subcluster color)
tmp_count <- table('sample_id' = vascular$sample_id,
                   'vascular_subcluster' = vascular$vascular_subcluster) %>%
  reshape2::melt() %>%
  mutate('group' = paste(vascular_subcluster, sample_id, sep = '_'))
tmp_time <- unique(vascular@meta.data[c('sample_id','time')])
tmp_count[['time']] <- tmp_time$time[match(tmp_count$sample_id, tmp_time$sample_id)]
tmp_count <- tmp_count[tmp_count$value != 0,]
tmp_count_index <- match(x = colnames(vascular_cor_byReplicate), table = tmp_count$group)
tmp_count_val <- tmp_count$value[tmp_count_index]
tmp_vascular_subcluster <- tmp_count$vascular_subcluster[tmp_count_index]
tmp_time <- tmp_count$time[tmp_count_index]
vascular_cor_anno <- rowAnnotation(
  'Cell count' = anno_barplot(
    x = tmp_count_val,
    axis_param = list(
      side = 'bottom',
      labels_rot = 90,
      direction = 'reverse'
    )
  ),
  'Subcluster' = tmp_vascular_subcluster,
  'Time' = tmp_time,
  col = list('Subcluster' = vascular_cols,
             'Time' = time_cols),
  annotation_name_rot = 90,
  annotation_name_side = 'bottom'
)

# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'vascular_spearmanCor_byReplicate_corrected.tiff'),
     height = 10, width = 12, units = 'in', res = 480)
Heatmap(
  matrix = vascular_cor_byReplicate,
  column_title = 'Correlation by vascular subcluster and replicate',
  column_title_gp = gpar(fontsize = 28),
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 0.8)(100),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 18),
                              labels_gp = gpar(fontsize = 14),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(4, units = 'mm')),
  left_annotation = vascular_cor_anno
)
dev.off()




# Correlation matrix by subcluster and dissociation method ---------------------

# Calculate average expression of genes by subcluster + dissociation method.
# Calculate avg for the 2000 variable genes used for clustering. Drop-out values
# inflate the spearman correlation between unrelated groups.
vascular_avg_byDissociationMethod <- AverageExpression(
  object = vascular,
  assays = 'RNAcorrected',
  features = vascular[['integrated']]@var.features,
  add.ident = 'dissociationMethod',
  slot = 'data'
)
vascular_avg_byDissociationMethod <- vascular_avg_byDissociationMethod[[1]]

# Calculate spearman correlation between subcluster + replicate
vascular_cor_byDissociationMethod <- cor(
  x = vascular_avg_byDissociationMethod,
  method = 'spearman'
)

# Heatmap annotation
tmp_count <- table('dissociationMethod' = vascular$dissociationMethod,
                   'vascular_subcluster' = vascular$vascular_subcluster) %>%
  reshape2::melt() %>%
  mutate('group' = paste(vascular_subcluster, dissociationMethod, sep = '_'))
tmp_count <- tmp_count[tmp_count$value != 0,]
tmp_count_index <- match(x = colnames(vascular_cor_byDissociationMethod), table = tmp_count$group)
tmp_count_val <- tmp_count$value[tmp_count_index]
tmp_vascular_subcluster <- tmp_count$vascular_subcluster[tmp_count_index]
vascular_cor_anno <- rowAnnotation(
  'Cell count' = anno_barplot(
    x = tmp_count_val,
    axis_param = list(
      side = 'bottom',
      labels_rot = 90,
      direction = 'reverse'
    )
  ),
  'Subcluster' = tmp_vascular_subcluster,
  col = list('Subcluster' = vascular_cols),
  annotation_name_rot = 90,
  annotation_name_side = 'bottom'
)

# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'vascular_spearmanCor_byDissociationMethod_corrected.tiff'),
     height = 6.25, width = 7.75, units = 'in', res = 480)
Heatmap(
  matrix = vascular_cor_byDissociationMethod,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
  column_title = 'Correlation by subcluster and dissociation method',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 14),
                              labels_gp = gpar(fontsize = 12),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(1, units = 'cm')),
  left_annotation = vascular_cor_anno
)
dev.off()




# Correlation matrix by subcluster and 10X Chemsitry ---------------------

# Calculate average expression of genes by subcluster + 10X Chemsitry. Calculate avg
# for the 2000 variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
vascular_avg_byChemistry <- AverageExpression(
  object = vascular,
  assays = 'RNAcorrected',
  features = vascular[['integrated']]@var.features,
  add.ident = 'chemistry',
  slot = 'data'
)
vascular_avg_byChemistry <- vascular_avg_byChemistry[[1]]

# Calculate spearman correlation between subcluster + replicate
vascular_cor_byChemistry <- cor(
  x = vascular_avg_byChemistry,
  method = 'spearman'
)

# Heatmap annotation
tmp_count <- table('chemistry' = vascular$chemistry,
                   'vascular_subcluster' = vascular$vascular_subcluster) %>%
  reshape2::melt() %>%
  mutate('group' = paste(vascular_subcluster, chemistry, sep = '_'))
tmp_count <- tmp_count[tmp_count$value != 0,]
tmp_count_index <- match(x = colnames(vascular_cor_byChemistry), table = tmp_count$group)
tmp_count_val <- tmp_count$value[tmp_count_index]
tmp_vascular_subcluster <- tmp_count$vascular_subcluster[tmp_count_index]
vascular_cor_anno <- rowAnnotation(
  'Cell count' = anno_barplot(
    x = tmp_count_val,
    axis_param = list(
      side = 'bottom',
      labels_rot = 90,
      direction = 'reverse'
    )
  ),
  'Subcluster' = tmp_vascular_subcluster,
  col = list('Subcluster' = vascular_cols),
  annotation_name_rot = 90,
  annotation_name_side = 'bottom'
)


# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'vascular_spearmanCor_byChemistry_corrected.tiff'),
     height = 7.5, width = 8.5, units = 'in', res = 480)
Heatmap(
  matrix = vascular_cor_byChemistry,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
  column_title = 'Correlation by subcluster and 10X Chemistry',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 14),
                              labels_gp = gpar(fontsize = 12),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(1, units = 'cm')),
  left_annotation = vascular_cor_anno
)
dev.off()




# Correlation matrix by subcluster and Time-point ---------------------

# Calculate average expression of genes by subcluster + Time-point. Calculate avg
# for the 2000 variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
vascular_avg_byTime <- AverageExpression(
  object = vascular,
  assays = 'RNAcorrected',
  features = vascular[['integrated']]@var.features,
  add.ident = 'time',
  slot = 'data'
)
vascular_avg_byTime <- vascular_avg_byTime[[1]]

# Calculate spearman correlation between subcluster + replicate
vascular_cor_byTime <- cor(
  x = vascular_avg_byTime,
  method = 'spearman'
)

# Heatmap annotation
tmp_count <- table('time' = vascular$time,
                   'vascular_subcluster' = vascular$vascular_subcluster) %>%
  reshape2::melt() %>%
  mutate('group' = paste(vascular_subcluster, time, sep = '_'))
tmp_count <- tmp_count[tmp_count$value != 0,]
tmp_count_index <- match(x = colnames(vascular_cor_byTime), table = tmp_count$group)
tmp_count_val <- tmp_count$value[tmp_count_index]
tmp_vascular_subcluster <- tmp_count$vascular_subcluster[tmp_count_index]
tmp_time <- tmp_count$time[tmp_count_index]
vascular_cor_anno <- rowAnnotation(
  'Cell count' = anno_barplot(
    x = tmp_count_val,
    axis_param = list(
      side = 'bottom',
      labels_rot = 90,
      direction = 'reverse'
    )
  ),
  'Subcluster' = tmp_vascular_subcluster,
  'Time' = tmp_time,
  col = list('Subcluster' = vascular_cols,
             'Time' = time_cols),
  annotation_name_rot = 90,
  annotation_name_side = 'bottom'
)


# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'vascular_spearmanCor_byTime_corrected.tiff'),
     height = 8.5, width = 10.5, units = 'in', res = 480)
Heatmap(
  matrix = vascular_cor_byTime,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 0.8)(100),
  column_title = 'Correlation by subcluster and Injury time-point',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 14),
                              labels_gp = gpar(fontsize = 12),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(1, units = 'cm')),
  left_annotation = vascular_cor_anno
)
dev.off()
