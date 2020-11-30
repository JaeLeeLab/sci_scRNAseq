

######## Macroglia expression of axon regeneration modulators ########


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
results_out <- './results/macroglia_axonregen_molecules/'
ref_in <- './ref/'
ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)

macroglia <- readRDS(file = './data/macroglia.rds')





# Genes of interest ------------------------------------------------------------


# We look at the differential expression of various permissive and inhibitory axon growth modulators. These genes were derived from [Anderson et al. Nature 2016](https://www.nature.com/articles/nature17623).  


# inhibitory_genes <- c('Acan','Bcan','Ncan','Vcan','Ptprz1','Xylt1','Tnr','Epha4','Ephb2','Efnb3','Ntn1','Sema3a','Sema3f','Plxna1','Plxnb1','Nrp1','Unc5b','Dcc','Neo1','Rgma','Rgmb','Slit1','Slit2','Slitrk1','Robo1','Robo2','Robo3','Draxin')
inhibitory_ligands <- c('Acan','Bcan','Ncan','Vcan','Xylt1','Tnr','Epha4','Ephb2','Efnb3','Ntn1','Sema3a','Sema3f','Rgma','Rgmb','Slit1','Slit2','Slitrk1','Draxin')
# inhibitory_genes <- inhibitory_genes[inhibitory_genes %in% rownames(macroglia[['RNA']]@counts)]
inhibitory_ligands <- inhibitory_ligands[inhibitory_ligands %in% rownames(macroglia[['RNA']]@counts)]

# permissive_genes <- c('Cspg4','Cspg5','Tnc','Sdc1','Sdc2','Sdc3','Sdc4','Bdnf','Ntf3','Gdnf','Lif','Cntf','Igf1','Fgf2','Tgfa','Lama1','Lama2','Lama4','Lama5','Lamb1','Lamc1','Col4a1','Fn1','Hspg2','Gpc1','Gpc3','Gpc5','Dcn','Lgals1','Ncam1','Matn2')
permissive_ligands <- c('Cspg4','Cspg5','Tnc','Sdc1','Sdc2','Sdc3','Sdc4','Bdnf','Ntf3','Lif','Cntf','Igf1','Fgf2','Tgfa','Lama1','Lama2','Lama4','Lama5','Lamb1','Lamc1','Col4a1','Fn1','Hspg2','Gpc1','Gpc3','Gpc5','Dcn','Lgals1','Ncam1','Matn2')
# permissive_genes <- permissive_genes[permissive_genes %in% rownames(macroglia[['RNA']]@counts)]
permissive_ligands <- permissive_ligands[permissive_ligands %in% rownames(macroglia[['RNA']]@counts)]






# Axon regen molecules heatmap (split by time) ----------------------------------------


DefaultAssay(macroglia) <- 'RNAcorrected'
Idents(macroglia) <- 'macroglia_subcluster'

# Preset values
macroglia_cols <- c('Ependymal-A' = '#800000',
                    'Ependymal-B' = '#e6194b',
                    'Astroependymal' = '#f58231',
                    'Astrocyte' = 'goldenrod',
                    'OPC-A' = '#3cb44b',
                    'OPC-B' = '#008080',
                    'Div-OPC' = '#4363d8',
                    'Pre-Oligo' = '#911eb4',
                    'Oligodendrocyte' = '#f032e6')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')
axon_cols <- c('Inhibitory' = 'indianred', 'Permissive' = 'dodgerblue')


# Select which genes to annotate
label_genes <- c(inhibitory_ligands[1:5], permissive_ligands[1:5])

# Extract scaled expression values
expr_data <- t(ScaleData(
  object = macroglia[['RNAcorrected']]@data, 
  features = c(inhibitory_ligands, permissive_ligands),
  scale.max = 3
))

# metadata
macroglia_meta <- macroglia@meta.data[c('macroglia_subcluster','time')]

# Merge metadata and exprs values, then calculate avg across cell/time
axon_data <- cbind(expr_data, macroglia_meta) %>%
  group_by(macroglia_subcluster, time) %>%
  summarise(across(.fns = mean))

# Remove groups with low counts (e.g. Uninjured Macrophage-B)
cell_count <- table(macroglia@meta.data[['macroglia_subcluster']],
                    macroglia@meta.data[['time']])
cell_prop <- prop.table(cell_count, margin = 1)
remove_cell <- cell_prop < 0.02
remove_rows <- c()
for (i in 1:nrow(axon_data)) {
  if (remove_cell[axon_data$macroglia_subcluster[i], axon_data$time[i]]) {
    remove_rows <- c(remove_rows, i)
  }
}
axon_data <- axon_data[-remove_rows,]

# Create new ID of celltype + time
axon_data[['id']] <- paste(axon_data[['macroglia_subcluster']],
                           axon_data[['time']],
                           sep = '_')

# Cell-level annotations
axon_cell_anno <- rowAnnotation(
  'Subcluster' = axon_data[['macroglia_subcluster']],
  'Time' = axon_data[['time']],
  col = list('Subcluster' = macroglia_cols,
             'Time' = time_cols)
)

# Gene-level annotations
axon_gene_anno <- HeatmapAnnotation(
  'Axon-regen effect' = c(rep('Inhibitory', length(inhibitory_ligands)), rep('Permissive', length(permissive_ligands))),
  col = list('Axon-regen effect' = axon_cols),
  'Axon-regen effect' = anno_mark(
    at = match(label_genes, colnames(axon_data)),
    labels = label_genes,
    side = 'right',
    labels_gp = gpar(fontsize = 10)
  )
)

# convert to matrix
axon_mat <- axon_data[!names(axon_data) %in% c('macroglia_subcluster','time')] %>%
  tibble::column_to_rownames(var = 'id') %>%
  as.matrix()


# Generate heatmap
axon_split_heatmap <- Heatmap(
  matrix = axon_mat,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 1.6)(100),
  heatmap_legend_param = list(title = 'Scaled\nExpression',
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
  top_annotation = axon_gene_anno,
  right_annotation = axon_cell_anno
)

# save
tiff(filename = paste0(results_out, 'macroglia_axonregen_timeSplit_heatmap.tiff'),
     height = 8, width = 14, units = 'in', res = 440)
axon_split_heatmap
dev.off()






# Inhibitory molecules heatmap (split by time) ----------------------------------------


DefaultAssay(macroglia) <- 'RNAcorrected'
Idents(macroglia) <- 'macroglia_subcluster'

# Preset values
macroglia_cols <- c('Ependymal-A' = '#800000',
                    'Ependymal-B' = '#e6194b',
                    'Astroependymal' = '#f58231',
                    'Astrocyte' = 'goldenrod',
                    'OPC-A' = '#3cb44b',
                    'OPC-B' = '#008080',
                    'Div-OPC' = '#4363d8',
                    'Pre-Oligo' = '#911eb4',
                    'Oligodendrocyte' = '#f032e6')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')
axon_cols <- c('Inhibitory' = 'indianred', 'Permissive' = 'dodgerblue')


# Extract scaled expression values
expr_data <- t(ScaleData(
  object = macroglia[['RNAcorrected']]@data, 
  features = c(inhibitory_ligands),
  scale.max = 3
))

# metadata
macroglia_meta <- macroglia@meta.data[c('macroglia_subcluster','time')]

# Merge metadata and exprs values, then calculate avg across cell/time
axon_data <- cbind(expr_data, macroglia_meta) %>%
  group_by(macroglia_subcluster, time) %>%
  summarise(across(.fns = mean))

# Remove groups with low counts (e.g. Uninjured astroependymal)
cell_count <- table(macroglia@meta.data[['macroglia_subcluster']],
                    macroglia@meta.data[['time']])
cell_prop <- prop.table(cell_count, margin = 1)
remove_cell <- cell_prop < 0.02
remove_rows <- c()
for (i in 1:nrow(axon_data)) {
  if (remove_cell[axon_data$macroglia_subcluster[i], axon_data$time[i]]) {
    remove_rows <- c(remove_rows, i)
  }
}
axon_data <- axon_data[-remove_rows,]

# Create new ID of celltype + time
axon_data[['id']] <- paste(axon_data[['macroglia_subcluster']],
                           axon_data[['time']],
                           sep = '_')

# Cell-level annotations
axon_cell_anno <- rowAnnotation(
  'Subcluster' = axon_data[['macroglia_subcluster']],
  'Time' = axon_data[['time']],
  col = list('Subcluster' = macroglia_cols,
             'Time' = time_cols)
)


# convert to matrix
axon_mat <- axon_data[!names(axon_data) %in% c('macroglia_subcluster','time')] %>%
  tibble::column_to_rownames(var = 'id') %>%
  as.matrix()


# Generate heatmap
axon_split_heatmap <- Heatmap(
  matrix = axon_mat,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 1.6)(100),
  heatmap_legend_param = list(title = 'Scaled\nExpression',
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
  right_annotation = axon_cell_anno
)

# save
tiff(filename = paste0(results_out, 'macroglia_axonregeninhibitory_timeSplit_heatmap.tiff'),
     height = 8, width = 8, units = 'in', res = 440)
axon_split_heatmap
dev.off()





# Inhibitory molecules violin plot (subclusters split by time) ----------------------------


inhibitory_ligands <- c('Acan','Bcan','Ncan','Vcan','Xylt1','Tnr','Epha4','Ephb2','Efnb3','Ntn1','Sema3a','Sema3f','Rgma','Rgmb','Slit1','Slit2','Slitrk1','Draxin')
inhibitory_ligands <- inhibitory_ligands[inhibitory_ligands %in% rownames(macroglia[['RNA']]@counts)]

# Preset values
macroglia_cols <- c('Ependymal-A' = '#800000',
                    'Ependymal-B' = '#e6194b',
                    'Astroependymal' = '#f58231',
                    'Astrocyte' = 'goldenrod',
                    'OPC-A' = '#3cb44b',
                    'OPC-B' = '#008080',
                    'Div-OPC' = '#4363d8',
                    'Pre-Oligo' = '#911eb4',
                    'Oligodendrocyte' = '#f032e6')


# Gather expression data and metadata
expr_data <- Matrix::t(macroglia[['RNAcorrected']]@data[inhibitory_ligands,])
macroglia_meta <- macroglia@meta.data[c('macroglia_subcluster','time')]
axon_data <- cbind(expr_data, macroglia_meta)

# Remove genes not expressed at least one group at 10%
remove_low_expr <- axon_data %>%
  group_by(macroglia_subcluster, time) %>%
  summarise(across(all_of(inhibitory_ligands), .fns = function(x) mean(x > 0) < 0.1)) %>%
  ungroup() %>%
  summarise(across(all_of(inhibitory_ligands), .fns = all)) %>%
  unlist(c(.)) %>%
  which(.) %>%
  names(.)

# Long form
axon_data <- axon_data[,!colnames(axon_data) %in% remove_low_expr] %>%
  reshape2::melt(id.vars = c('macroglia_subcluster','time'))

# generate violin plot
axon_inhibitory_vln <- axon_data %>%
  ggplot(mapping = aes(x = time, y = value)) +
  geom_violin(mapping = aes(fill = macroglia_subcluster), scale = 'width') +
  facet_grid(variable ~ macroglia_subcluster, scales = 'free_y', switch = 'y') +
  scale_fill_manual(values = macroglia_cols) +
  scale_y_continuous(position = 'right') +
  ylab(label = 'Log-normalized expression') +
  xlab(label = 'Time') +
  theme(axis.text.x = element_text(size = 12, color = 'black', angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        # strip.background = element_rect(fill = NA, color = NA),
        panel.border = element_rect(fill = NA, color = 'black'),
        strip.text.y.left = element_text(size = 12, face = 'bold', color = 'black', angle = 0, hjust = 1),
        strip.text.x = element_text(size = 12, face = 'bold', color = 'black'),
        legend.position = 'none',
        axis.ticks = element_blank(),
        axis.text.y = element_blank())
# save
ggsave(filename =  paste0(results_out, 'macroglia_axonregeninhibitory_vln.tiff'),
       plot = axon_inhibitory_vln, device = 'tiff', height = 9, width = 14)


# Repeat with only Astrocytes, Astroependymal, OPC-A, and OPC-B
axon_inhibitory_vln <- axon_data %>%
  filter(macroglia_subcluster %in% c('Astroependymal','Astrocyte','OPC-A','OPC-B')) %>%
  ggplot(mapping = aes(x = time, y = value)) +
  geom_violin(mapping = aes(fill = macroglia_subcluster), scale = 'width') +
  facet_grid(variable ~ macroglia_subcluster, scales = 'free_y', switch = 'y') +
  scale_fill_manual(values = macroglia_cols) +
  scale_y_continuous(position = 'right') +
  ylab(label = 'Log-normalized expression') +
  xlab(label = 'Time') +
  theme(axis.text.x = element_text(size = 12, color = 'black', angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        # strip.background = element_rect(fill = NA, color = NA),
        panel.border = element_rect(fill = NA, color = 'black'),
        strip.text.y.left = element_text(size = 12, face = 'bold', color = 'black', angle = 0, hjust = 1),
        strip.text.x = element_text(size = 12, face = 'bold', color = 'black'),
        legend.position = 'none',
        axis.ticks = element_blank(),
        axis.text.y = element_blank())
ggsave(filename =  paste0(results_out, 'macroglia_axonregeninhibitory_selectSubclusters_vln.tiff'),
       plot = axon_inhibitory_vln, device = 'tiff', height = 8.5, width = 6.5)



# Inhibitory molecules violin plot (OPCs vs Astrocytes) ----------------------------


inhibitory_ligands <- c('Acan','Bcan','Ncan','Vcan','Xylt1','Tnr','Epha4','Ephb2','Efnb3','Ntn1','Sema3a','Sema3f','Rgma','Rgmb','Slit1','Slit2','Slitrk1','Draxin')
inhibitory_ligands <- inhibitory_ligands[inhibitory_ligands %in% rownames(macroglia[['RNA']]@counts)]

# Preset values
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

# Gather expression data and metadata
expr_data <- Matrix::t(macroglia[['RNAcorrected']]@data[inhibitory_ligands,])
macroglia_meta <- macroglia@meta.data[c('celltype','time')]
axon_data <- cbind(expr_data, macroglia_meta)

# Remove genes not expressed at least one group at 10%
remove_low_expr <- axon_data %>%
  group_by(celltype, time) %>%
  summarise(across(all_of(inhibitory_ligands), .fns = function(x) mean(x > 0) < 0.1)) %>%
  ungroup() %>%
  summarise(across(all_of(inhibitory_ligands), .fns = all)) %>%
  unlist(c(.)) %>%
  which(.) %>%
  names(.)

# Long form
axon_data <- axon_data[,!colnames(axon_data) %in% remove_low_expr] %>%
  reshape2::melt(id.vars = c('celltype','time'))

# generate violin plot
axon_inhibitory_vln <- axon_data %>%
  filter(celltype %in% c('Astrocyte','OPC')) %>%
  ggplot(mapping = aes(x = time, y = value)) +
  geom_violin(mapping = aes(fill = celltype), scale = 'width') +
  facet_grid(variable ~ celltype, scales = 'free_y', switch = 'y') +
  scale_fill_manual(values = sci_cols) +
  scale_y_continuous(position = 'right') +
  ylab(label = 'Log-normalized expression') +
  xlab(label = 'Time') +
  theme(axis.text.x = element_text(size = 12, color = 'black', angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'),
        axis.title.x = element_blank(),
        # panel.background = element_rect(fill = NA, color = NA),
        # strip.background = element_rect(fill = NA, color = NA),
        # panel.border = element_rect(fill = NA, color = 'grey60'),
        strip.text.y.left = element_text(size = 12, face = 'bold', color = 'black', angle = 0, hjust = 1),
        strip.text.x = element_text(size = 12, face = 'bold', color = 'black'),
        legend.position = 'none',
        axis.ticks = element_blank(),
        axis.text.y = element_blank())

# save
ggsave(filename =  paste0(results_out, 'macroglia_axonregeninhibitory_OPC_Astrocyte_vln.tiff'),
       plot = axon_inhibitory_vln, device = 'tiff', height = 9, width = 3.5)


