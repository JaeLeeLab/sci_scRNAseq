
######## SCI clustering correlation to assess experimental design #########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
require('ComplexHeatmap')
results_out <- './results/sci_clustering_correlation/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)

sci <- readRDS(file = './data/sci.rds')

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



# Correlation matrix by cell-type and replicate ---------------------------


# Calculate average expression of genes by cell-type + replicate. Calculate avg
# for the 2000 variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
sci_avg_byReplicate <- AverageExpression(
  object = sci,
  assays = 'RNA',
  features = sci[['integrated']]@var.features,
  add.ident = 'sample_id',
  slot = 'data'
)
sci_avg_byReplicate <- sci_avg_byReplicate[['RNA']]

# Calculate spearman correlation between cell-type + replicate
sci_cor_byReplicate <- cor(
  x = sci_avg_byReplicate,
  method = 'spearman'
)

# Heatmap annotation (cell counts barplot and celltype color)
tmp_count <- table('sample_id' = sci$sample_id,
                  'celltype' = sci$celltype) %>%
  reshape2::melt() %>%
  mutate('group' = paste(celltype, sample_id, sep = '_'))
tmp_count <- tmp_count[tmp_count$value != 0,]
tmp_count_index <- match(x = colnames(sci_cor_byReplicate), table = tmp_count$group)
tmp_count_val <- tmp_count$value[tmp_count_index]
tmp_celltype <- tmp_count$celltype[tmp_count_index]
tmp_batch <- tmp_count$sample_id[tmp_count_index]
sci_cor_count_anno <- rowAnnotation(
  'Cell count' = anno_barplot(
    x = tmp_count_val,
    axis_param = list(
      side = 'bottom',
      labels_rot = 90,
      direction = 'reverse'
    )
  ),
  'Cell-type' = tmp_celltype,
  'Batch' = tmp_batch,
  col = list('Cell-type' = sci_cols),
  annotation_name_rot = 90,
  annotation_name_side = 'bottom'
)


# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'sci_spearmanCor_byReplicate_RNA.tiff'),
         height = 15, width = 17, units = 'in', res = 480)
Heatmap(
  matrix = sci_cor_byReplicate,
  column_title = 'Correlation by cell-type and replicate',
  column_title_gp = gpar(fontsize = 32),
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
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
  left_annotation = sci_cor_count_anno
)
dev.off()



# Using SCTransform normalized values.
# Calculate average expression of genes by cell-type + replicate. Calculate avg
# for the 2000 variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
Idents(sci) <- 'celltype'
sci_avg_byReplicate <- AverageExpression(
  object = sci,
  assays = 'SCT',
  features = sci[['integrated']]@var.features,
  add.ident = 'sample_id',
  slot = 'data'
)
sci_avg_byReplicate <- sci_avg_byReplicate[['SCT']]

# Calculate spearman correlation between cell-type + replicate
sci_cor_byReplicate <- cor(
  x = sci_avg_byReplicate,
  method = 'spearman'
)

# Heatmap annotation (cell counts barplot and celltype color)
tmp_count <- table('sample_id' = sci$sample_id,
                   'celltype' = sci$celltype) %>%
  reshape2::melt() %>%
  mutate('group' = paste(celltype, sample_id, sep = '_'))
tmp_count <- tmp_count[tmp_count$value != 0,]
tmp_count_index <- match(x = colnames(sci_cor_byReplicate), table = tmp_count$group)
tmp_count_val <- tmp_count$value[tmp_count_index]
tmp_celltype <- tmp_count$celltype[tmp_count_index]
tmp_batch <- tmp_count$sample_id[tmp_count_index]
sci_cor_count_anno <- rowAnnotation(
  'Cell count' = anno_barplot(
    x = tmp_count_val,
    axis_param = list(
      side = 'bottom',
      labels_rot = 90,
      direction = 'reverse'
    )
  ),
  'Cell-type' = tmp_celltype,
  'Batch' = tmp_batch,
  col = list('Cell-type' = sci_cols),
  annotation_name_rot = 90,
  annotation_name_side = 'bottom'
)


# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'sci_spearmanCor_byReplicate_SCT.tiff'),
     height = 15, width = 17, units = 'in', res = 480)
Heatmap(
  matrix = sci_cor_byReplicate,
  column_title = 'Correlation by cell-type and replicate',
  column_title_gp = gpar(fontsize = 32),
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
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
  left_annotation = sci_cor_count_anno
)
dev.off()




# Correlation matrix by cell-type and dissociation method ---------------------

# Calculate average expression of genes by cell-type + dissociation method.
# Calculate avg for the 2000 variable genes used for clustering. Drop-out values
# inflate the spearman correlation between unrelated groups.
sci_avg_byDissociationMethod <- AverageExpression(
  object = sci,
  assays = 'RNA',
  features = sci[['integrated']]@var.features,
  add.ident = 'dissociationMethod',
  slot = 'data'
)
sci_avg_byDissociationMethod <- sci_avg_byDissociationMethod[['RNA']]

# Calculate spearman correlation between cell-type + replicate
sci_cor_byDissociationMethod <- cor(
  x = sci_avg_byDissociationMethod,
  method = 'spearman'
)

# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'sci_spearmanCor_byDissociationMethod_RNA.tiff'),
     height = 7.5, width = 8.5, units = 'in', res = 480)
Heatmap(
  matrix = sci_cor_byDissociationMethod,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
  column_title = 'Correlation by cell-type and dissociation method',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 18),
                              labels_gp = gpar(fontsize = 14),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(4, units = 'mm'))
)
dev.off()




# Correlation matrix by cell-type and 10X Chemsitry ---------------------

# Calculate average expression of genes by cell-type + 10X Chemsitry. Calculate avg
# for the 2000 variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
sci_avg_byChemistry <- AverageExpression(
  object = sci,
  assays = 'RNA',
  features = sci[['integrated']]@var.features,
  add.ident = 'chemistry',
  slot = 'data'
)
sci_avg_byChemistry <- sci_avg_byChemistry[['RNA']]

# Calculate spearman correlation between cell-type + replicate
sci_cor_byChemistry <- cor(
  x = sci_avg_byChemistry,
  method = 'spearman'
)

# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'sci_spearmanCor_byChemistry_RNA.tiff'),
     height = 7.5, width = 8.5, units = 'in', res = 480)
Heatmap(
  matrix = sci_cor_byChemistry,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
  column_title = 'Correlation by cell-type and 10X Chemistry',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 18),
                              labels_gp = gpar(fontsize = 14),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(4, units = 'mm'))
)
dev.off()




# Correlation matrix by cell-type and Time-point ---------------------

# Calculate average expression of genes by cell-type + Time-point. Calculate avg
# for the 2000 variable genes used for clustering. Drop-out values inflate the
# spearman correlation between unrelated groups.
sci_avg_byTime <- AverageExpression(
  object = sci,
  assays = 'RNA',
  features = sci[['integrated']]@var.features,
  add.ident = 'time',
  slot = 'data'
)
sci_avg_byTime <- sci_avg_byTime[['RNA']]

# Calculate spearman correlation between cell-type + replicate
sci_cor_byTime <- cor(
  x = sci_avg_byTime,
  method = 'spearman'
)

# Spearman correlation heatmap
tiff(filename = paste0(results_out, 'sci_spearmanCor_byTime_RNA.tiff'),
     height = 10, width = 11, units = 'in', res = 480)
Heatmap(
  matrix = sci_cor_byTime,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
  column_title = 'Correlation by cell-type and Injury time-point',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = 'Spearman\nCorrelation',
                              title_gp = gpar(fontsize = 18),
                              labels_gp = gpar(fontsize = 14),
                              legend_height = unit(4, units = 'cm'),
                              grid_width = unit(0.75, units = 'cm'),
                              border = 'black',
                              title_gap = unit(4, units = 'mm'))
)
dev.off()
