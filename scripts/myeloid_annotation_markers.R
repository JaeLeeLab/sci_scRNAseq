
######## Myeloid cell-type annotation via marker identification ########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('scran')
require('Seurat')
require('dplyr')
require('ggplot2')
require('ComplexHeatmap')
results_out <- './results/myeloid_annotation_markers/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)

myeloid <- readRDS(file = './data/myeloid.rds')





# Additional markers UMAP ----------------------------------------------------


DefaultAssay(myeloid) <- 'RNA'

# These are marker genes identified in previous studies. Multiple sources.
addtl_markers <- c('Adgre1','Itgam','P2ry12','Tmem119','Siglech','Cx3cr1','Csf1r','Irf8','Msr1','Igf1','Mki67','Egr1','Ramp1','Lyz2','Ly6c2','Ccr2','Plac8','Cxcl3','Cebpb','Arg1','Cd63','Cd68','Spp1','Apoe','Trem2','Gpnmb','Ms4a7','Ly6g','H2-Aa','Mrc1')

umap_coordinates <- FetchData(myeloid, vars = c('UMAP_1','UMAP_2'), slot = 'data')
label_x <- min(umap_coordinates[['UMAP_1']]) * 0.96
label_y <- max(umap_coordinates[['UMAP_2']]) * 0.85
legend_x <- 0.85
legend_y <- 0.85
expr_colors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = 'Spectral'))(100))

addtl_markers_umap <- vector(mode = 'list', length = length(addtl_markers))
names(addtl_markers_umap) <- addtl_markers

for (ii in 1:length(addtl_markers)) {
  gene <- addtl_markers[ii]
  # gather data
  expr <- FetchData(myeloid, vars = gene, slot = 'data')
  expr <- cbind(expr, umap_coordinates)
  expr <- expr[order(expr[[gene]], decreasing = FALSE),]
  max_expr <- ceiling(max(expr[[gene]]) * 10) / 10
  
  if (grepl('-', x = gene)) {
    names(addtl_markers_umap)[ii] <- gsub(pattern = '-', replacement = '.', x = names(addtl_markers_umap)[ii])
    colnames(expr) <- gsub(pattern = '-', replacement = '.', x = colnames(expr))
    gene <- gsub(pattern = '-', replacement = '.', x = gene)
  }
  
  # assemble gg
  addtl_markers_umap[[gene]] <- expr %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes_string(color = gene), size = 0.5, alpha = 0.3) +
    geom_text(data = data.frame('gene' = gsub(pattern = '\\.', replacement = '-', x = gene),
                                'x_pos' = label_x,
                                'y_pos' = label_y),
              mapping = aes(x = x_pos, y = y_pos, label = gene),
              fontface = 'italic',
              size = 7,
              hjust = 0) +
    scale_color_gradientn(colors = expr_colors,
                          breaks = c(0, max_expr),
                          limits = c(0, max_expr)) +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(legend_x, legend_y), 
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(0.25, units = 'cm'), 
          legend.spacing.x = unit(0.1, units = 'cm'),
          legend.text = element_text(size = 14, color = 'black')) +
    guides(color = guide_colorbar(frame.colour = 'black', 
                                  ticks = FALSE,
                                  barwidth = 0.8,
                                  frame.linewidth = 1))
}
# save
addtl_markers_umap_grid <- cowplot::plot_grid(plotlist = addtl_markers_umap, ncol = 6)
ggsave(filename = paste0(results_out, 'addtl_markers_umap.tiff'),
       plot = addtl_markers_umap_grid, device = 'tiff',
       height = 14, width = 18)






# DE marker genes calculation -----------------------------------------------

# Use wilcox test while blocking against Time (group) to identify marker genes
# per cluster. scran::findMarkers() performs pair-wise for all cluster 
# combinations. 
rna_sce <- as.SingleCellExperiment(myeloid, assay = 'RNA')
markers_time <- findMarkers(x = rna_sce,
                            groups = colData(rna_sce)[['integrated_snn_res.0.35']],
                            test.type = 'wilcox',
                            direction = 'up',
                            pval.type = 'all',
                            block = colData(rna_sce)[['time']])

# Retain genes with FDR < 0.001
min_fdr <- 1e-03
sig_markers_time <- vector(mode = 'list', length = length(markers_time))
names(sig_markers_time) <- paste0('C', names(markers_time))
for (ii in 1:length(markers_time)) {
  cluster_id <- paste0('C', names(markers_time)[ii])
  tmp <- markers_time[[ii]]
  tmp <- tmp[tmp[['FDR']] < min_fdr, c('p.value', 'FDR')]
  tmp[['gene']] <- rownames(tmp)
  tmp[['cluster']] <- cluster_id
  sig_markers_time[[cluster_id]] <- tmp
}
sig_markers_time <- data.frame(do.call(rbind, sig_markers_time))
write.table(x = sig_markers_time, 
            file = paste0(results_out, 'myeloid_defaultclusters_DE_wilcox_blockTime.tsv'),
            sep = '\t', row.names = TRUE, col.names = NA)


# Alternatively, perform naive DE wilcox with Seurat's wrapper. 
Idents(myeloid) <- 'integrated_snn_res.0.35'
seurat_markers <- FindAllMarkers(object = myeloid,
                                 assay = 'RNA', 
                                 logfc.threshold = 0.75,
                                 only.pos = TRUE,
                                 min.pct = 0.4)
write.table(seurat_markers, file = paste0(results_out, 'myeloid_defaultclusters_topDE_wilcox.tsv'),
            sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)





# Cell-type annotation ----------------------------------------------------
# 
# myeloid@meta.data[['old_subcluster']] <- plyr::mapvalues(
#   x = myeloid@meta.data[['integrated_snn_res.0.35']],
#   from = levels(myeloid@meta.data[['integrated_snn_res.0.35']]),
#   to = c('H-Microglia',
#          'DAM-A',
#          'Macrophage-A',
#          'DAM-B',
#          'Macrophage-B',
#          'Monocyte',
#          'Neutrophil',
#          'DAM-C',
#          'Dendritic',
#          'IFN-Myeloid',
#          'Div-Myeloid',
#          'BA-Macrophage')
# )
# myeloid@meta.data[['old_subcluster']] <- factor(
#   x = myeloid@meta.data[['old_subcluster']],
#   levels = c('Neutrophil',
#              'Monocyte',
#              'Macrophage-A',
#              'Macrophage-B',
#              'BA-Macrophage',
#              'Dendritic',
#              'Div-Myeloid',
#              'H-Microglia',
#              'DAM-A',
#              'DAM-B',
#              'DAM-C',
#              'IFN-Myeloid')
# )
# Idents(myeloid) <- 'myeloid_subcluster'
# 
# # Save barcode-identities data.frame
# write.table(x = myeloid@meta.data['myeloid_subcluster'],
#             file = paste0(results_out, 'myeloid_subcluster.tsv'),
#             sep = '\t', row.names = TRUE, col.names = NA)
# 
# 
#
# Functional subcluster annotation --------------------------------------


# dir.create(path = './results/revision_figures/functional_names/')

myeloid@meta.data[['myeloid_subcluster']] <- plyr::mapvalues(
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
myeloid@meta.data[['myeloid_subcluster']] <- factor(
  x = myeloid@meta.data[['myeloid_subcluster']],
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

# Save barcode-identities data.frame
write.table(x = myeloid@meta.data['myeloid_subcluster'],
            file = paste0(results_out, 'myeloid_subcluster.tsv'),
            sep = '\t', row.names = TRUE, col.names = NA)


# Cell-type UMAP across time ----------------------------------------------

myeloid_cols <- c('#800000', 
                  '#9a6324',
                  '#e6194b', 
                  '#f58231', 
                  '#CCCC00',
                  '#808000', 
                  '#3cb44b',
                  '#008080', 
                  'cyan3', 
                  '#4363d8',
                  '#000075', 
                  '#911eb4')
names(myeloid_cols) <- c('Neutrophil',
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
Idents(myeloid) <- 'myeloid_subcluster'

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 16, color = 'black'),
                    legend.title = element_text(size = 16, color = 'black'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 16, color = 'black'))

# cell-type annotation UMAP 
myeloid_subcluster_counts <- table(myeloid$myeloid_subcluster)
myeloid_subcluster_label <- paste0(names(myeloid_subcluster_counts), ' (', myeloid_subcluster_counts, ')')
names(myeloid_subcluster_label) <- names(myeloid_subcluster_counts)
myeloid_subcluster_umap <- FetchData(object = myeloid, vars = c('UMAP_1','UMAP_2','myeloid_subcluster')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = myeloid_subcluster), size = 0.2, alpha = 0.5) +
  scale_color_manual(values = myeloid_cols, 
                     breaks = names(myeloid_subcluster_label),
                     label = myeloid_subcluster_label) +
  xlab(label = 'UMAP 1') +
  ylab(label = 'UMAP 2') +
  umap_theme +
  theme(legend.text = element_text(size = 16, color = 'black')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
ggsave(filename = paste0(results_out, 'myeloid_subcluster_annotation_umap.tiff'),
       plot = myeloid_subcluster_umap, device = 'tiff', height = 6, width = 10.25)


# cell-type annotation split by time UMAP (figure 1b)
myeloid_subcluster_split_umap <- FetchData(myeloid, vars = c('UMAP_1','UMAP_2','myeloid_subcluster','time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = myeloid_subcluster), size = 0.2, alpha = 0.5) +
  facet_wrap(. ~ time, ncol = 2) +
  scale_color_manual(values = myeloid_cols, breaks = names(myeloid_subcluster_label), label = myeloid_subcluster_label) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
ggsave(filename = paste0(results_out, 'myeloid_subcluster_annotation_split_umap.tiff'),
       plot = myeloid_subcluster_split_umap, device = 'tiff', height = 6, width = 9)




# Microglia DE markers dot plot (functional annotation) -----------------------


### Using same markers as flow cytometry ###

# Subset cells
microglia <- c('Homeostatic Microglia',
               'Inflammatory Microglia',
               'Dividing Microglia',
               'Migrating Microglia')
Idents(myeloid) <- 'myeloid_subcluster'
microglia <- subset(myeloid, idents = microglia)

microglia_markers <- c('P2ry12','Igf1','Msr1','Cdk1')
tmp_de <- FindAllMarkers(
  object = microglia,
  assay = 'RNA',
  slot = 'data',
  features = microglia_markers,
  only.pos = TRUE
)
tmp_de$p_val_adj <- signif(tmp_de$p_val_adj, digits = 2)
tmp_de$variable <- tmp_de$gene
tmp_de$signif <- ifelse(test = tmp_de$p_val < 1e-10,
                        yes = '*',
                        no = '')

# Calculate average and percent expression (scaled to within microglia)
DefaultAssay(microglia) <- 'RNA'
avg_exp <- data.frame(t(ScaleData(microglia[['RNA']]@data, features = microglia_markers)))
avg_exp <- cbind(avg_exp, 'myeloid_subcluster' = as.character(microglia@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(t(microglia[['RNA']]@counts[microglia_markers,]))
pct_exp <- cbind(pct_exp, 'myeloid_subcluster' = as.character(microglia@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
# min_expr <- floor(min(avg_exp$avg.exp)*10)/10
max_expr <- ceiling(max(avg_exp$avg.exp)*10)/10
min_expr <- -1*max_expr

# Plot
microglia_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(myeloid_subcluster = factor(myeloid_subcluster, levels = rev(levels(myeloid$myeloid_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  geom_text(data = tmp_de, mapping = aes(x = variable, y = cluster, label = signif),
            size = 8, fontface = 'bold', nudge_y = -0.05, nudge_x = 0.01) +
  scale_size(range = c(0,14), limits = c(0,100)) +
  scale_size(range = c(0,14), limits = c(0,100)) +
  scale_fill_gradientn(colors = rev(colorRampPalette(
    colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100)),
    breaks = c(min_expr, 0, max_expr),
    labels = c(min_expr, 0, max_expr),
    limits = c(min_expr, max_expr)) +
  scale_y_discrete(position = 'right') +
  theme(plot.margin = margin(0, 0, 0, 20, unit = 'mm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = 'black', face = 'italic'),
        axis.text.y = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 12, color = 'black', hjust = 0),
        legend.title = element_text(size = 12, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(5,0,10,0),
        legend.box = 'horizontal',
        legend.direction = 'vertical',
        legend.position = 'right',
        legend.spacing.x = unit(x = 2, units = 'mm'),
        panel.border = element_rect(fill = NA, size = 1),
        panel.background = element_rect(fill = NA)) +
  guides(fill = guide_colorbar(title = 'z-score', 
                               barwidth = 1.25,
                               frame.colour = 'black', 
                               frame.linewidth = 1.25,
                               ticks.colour = 'black', 
                               ticks.linewidth = 1.25,
                               title.position = 'left'), 
         size = guide_legend(title = '% expression', 
                             override.aes = list(fill = 'black'),
                             title.position = 'left'))
microglia_markers_dotplot
ggsave(filename = paste0(results_out, 'microglia_markers_dotplot.tiff'),
       plot = microglia_markers_dotplot, device = 'tiff', height = 2.75, width = 6.5)




# Microglia DE markers heatmap (functional annotation) -------------------------

# Preset values
Idents(myeloid) <- 'myeloid_subcluster'
DefaultAssay(myeloid) <- 'RNA'
myeloid_cols <- c('#800000', 
                  '#9a6324',
                  '#e6194b', 
                  '#f58231', 
                  '#CCCC00',
                  '#808000', 
                  '#3cb44b',
                  '#008080', 
                  'cyan3', 
                  '#4363d8',
                  '#000075', 
                  '#911eb4')
names(myeloid_cols) <- c('Neutrophil',
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
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

# Subset cells
microglia_names <- c('Homeostatic Microglia',
                     'Inflammatory Microglia',
                     'Dividing Microglia',
                     'Migrating Microglia',
                     'Interferon Myeloid')
microglia <- subset(myeloid, idents = microglia_names)

# Use Seurat wilcox test
microglia_srat_genes <- FindAllMarkers(object = microglia,
                                       assay = 'RNA',
                                       logfc.threshold = 0.25,
                                       only.pos = TRUE)
write.table(file = paste0(results_out, 'microglia_DE_wilcox.tsv'),
            x = microglia_srat_genes, sep = '\t', row.names = TRUE, 
            col.names = NA)
microglia_srat_genes <- read.table(file = paste0(results_out, 'microglia_DE_wilcox.tsv'),
                                   sep = '\t', header = TRUE, row.names = 1)

# Build heatmap
microglia_de_genes <- microglia_srat_genes %>%
  group_by(cluster) %>%
  filter(avg_logFC > 0.5) %>%
  filter(p_val_adj < 1e-03) %>%
  top_n(n = 10, wt = avg_logFC) %>%
  .[['gene']] %>%
  unique()
expr_data <- t(ScaleData(object = microglia@assays$RNA@data[microglia_de_genes,])) %>%
  cbind(microglia@meta.data[c('myeloid_subcluster')]) %>%
  group_by(myeloid_subcluster) %>%
  summarise(across(.cols = all_of(microglia_de_genes), .fns = mean)) %>%
  reshape2::melt(id.vars = 'myeloid_subcluster') %>%
  mutate(variable = factor(variable, levels = unique(variable)),
         myeloid_subcluster = factor(myeloid_subcluster, rev(levels(myeloid_subcluster))))

expr_data$pval <- 1
for (i in 1:nrow(microglia_srat_genes)) {
  tmp <- which(as.character(expr_data$myeloid_subcluster) == as.character(microglia_srat_genes$cluster[i]) &
                 as.character(expr_data$variable) == as.character(microglia_srat_genes$gene[i]))
  expr_data$pval[tmp] <- microglia_srat_genes$p_val_adj[i]
}
expr_data$significant <- ifelse(test = expr_data$pval < 1e-10,
                                yes = '*',
                                no = ' ')
max_expr <- 1.5
min_expr <- -1*max_expr
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))

microglia_de_genes_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_tile(mapping = aes(fill = value), color = 'black', size = 0.3) +
  geom_text(mapping = aes(label = significant)) +
  scale_fill_gradientn(colors = myColors,
                       limits = c(min_expr, max_expr), 
                       na.value = myColors[100],
                       breaks = c(min_expr,0,max_expr)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = 'black', face = 'italic'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.title = element_text(color = 'black', size = 14, hjust = 0.5),
        legend.text = element_text(size = 12, color = 'black'),
        legend.margin=margin(-10, -10, 0, -10),
        legend.box.margin=margin(-5,0,0,0,unit='mm')) +
  guides(fill = guide_colorbar(title = 'z-score',
                               title.position = 'top',
                               frame.colour = 'black',
                               ticks.colour = 'black',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               title.hjust = 0.5))
microglia_de_genes_heatmap
ggsave(filename = paste0(results_out, 'microglia_de_genes_heatmap.tiff'),
       plot = microglia_de_genes_heatmap, device = 'tiff',
       height = 3, width = 12)





# Complex Microglia DE heatmap --------------------------------------------

# Preset values
Idents(myeloid) <- 'myeloid_subcluster'
DefaultAssay(myeloid) <- 'RNA'
myeloid_cols <- c('Neutrophil' = '#800000',    
                  'Monocyte' = '#9a6324',
                  'Macrophage-A' = '#e6194b',
                  'Macrophage-B' = '#f58231',
                  'BA-Macrophage' = '#CCCC00',
                  'Dendritic' = '#808000',
                  'Div-Myeloid' = '#3cb44b',
                  'DAM-A' = '#008080',
                  'DAM-B' = 'cyan3',
                  'DAM-C' = '#000075',
                  'Div-Microglia' = '#4363d8',
                  'IFN-Myeloid' = '#911eb4')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

# Subset cells
microglia_names <- c('DAM-A',
                     'DAM-B',
                     'DAM-C',
                     'H-Microglia')
microglia <- subset(myeloid, idents = microglia_names)

# Use Seurat wilcox test
microglia_srat_genes <- FindAllMarkers(object = microglia,
                                       assay = 'RNA',
                                       logfc.threshold = 0.25,
                                       only.pos = TRUE)
microglia_de_genes <- microglia_srat_genes %>%
  group_by(cluster) %>%
  filter(avg_logFC > 0.5) %>%
  .[['gene']] %>%
  unique()

# Select which genes to annotate
label_genes <- c('P2ry12','Siglech','Tmem119','Il1b','Msr1','Igf1','Apoe','Lgals3','Irf7','Isg15','Cdk1','Mki67','Lpl','Csf1r','Plin2','Mif','Cxcl2','Fabp5','Lyz2','Ccl6','Ccl9','Spp1','Cd68','Cd5l','Ldha')


# Extract scaled expression
expr_data <- ScaleData(object = microglia[['RNAcorrected']]@data[microglia_de_genes,], 
                       scale.max = 4)
microglia_meta <- microglia@meta.data[c('myeloid_subcluster','time','sample_id')]

# Cell-level annotations
microglia_cell_anno <- HeatmapAnnotation(
  'Subcluster' = microglia_meta[['myeloid_subcluster']],
  'Time' = microglia_meta[['time']],
  col = list('Subcluster' = myeloid_cols[names(myeloid_cols) %in% microglia_names],
             'Time' = time_cols)
)

# Gene-level annotations
microglia_gene_anno <- rowAnnotation(
  'Genes' = anno_mark(
    at = match(label_genes, rownames(expr_data)),
    labels = label_genes,
    side = 'right',
    labels_gp = gpar(fontsize = 10)
  )
)

# Build heatmap
tmp <- Heatmap(
  matrix = expr_data, 
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 0.6)(100),
  column_title = 'Differentially expressed genes among microglia subclusters',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(title = 'Scaled Expression',
                              title_gp = gpar(fontsize = 12),
                              title_position = 'leftcenter-rot',
                              labels = c('Low','High'),
                              at = c(min(expr_data), max(expr_data)),
                              labels_gp = gpar(fontsize = 10),
                              legend_height = unit(2.5, units = 'cm'),
                              grid_width = unit(0.5, units = 'cm'),
                              border = 'black',
                              title_gap = unit(1, units = 'cm'),
                              direction = 'vertical'),
  use_raster = TRUE,
  top_annotation = microglia_cell_anno,
  right_annotation = microglia_gene_anno
)

# Save heatmap
pdf(file = paste0(results_out, 'microglia_heatmap_complex_small.pdf'),
    height = 5, width = 9)
tmp
dev.off()
svg(file = paste0(results_out, 'microglia_heatmap_complex_small.svg'),
    height = 5, width = 9)
tmp
dev.off()





# Macrophage DE markers dot plot (functional annotation) ----------------------


### Using same markers as flow cytometry + preprint %%%%%%%%%

# Subset cells
macrophage <- c('Neutrophil',
                'Monocyte',
                'Chemotaxis-Inducing Mac',
                'Inflammatory Mac',
                'Border-Associated Mac',
                'Dendritic',
                'Dividing Myeloid')
Idents(myeloid) <- 'myeloid_subcluster'
macrophage <- subset(myeloid, idents = macrophage)

# Previously identified genes
macrophage_markers <- c('Cxcl3','Plac8','Arg1','Hmox1','Apoe','Cd63','Mrc1','Cd74','Cdk1')
tmp_de <- FindAllMarkers(
  object = macrophage,
  assay = 'RNA',
  slot = 'data',
  features = macrophage_markers,
  only.pos = TRUE
)
tmp_de$p_val_adj <- signif(tmp_de$p_val_adj, digits = 2)
tmp_de$variable <- tmp_de$gene
tmp_de$signif <- ifelse(test = tmp_de$p_val < 1e-10,
                        yes = '*',
                        no = '')

# Calculate average and percent expression (scaled to within macrophage)
DefaultAssay(macrophage) <- 'RNA'
avg_exp <- data.frame(t(ScaleData(macrophage[['RNA']]@data, features = macrophage_markers)))
avg_exp <- cbind(avg_exp, 'myeloid_subcluster' = as.character(macrophage@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(t(macrophage[['RNA']]@counts[macrophage_markers,]))
pct_exp <- cbind(pct_exp, 'myeloid_subcluster' = as.character(macrophage@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
my_cols <- rev(colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))

# Plot
max_expr <- 1.4
min_expr <- -1.4
macrophage_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(myeloid_subcluster = factor(myeloid_subcluster, levels = rev(levels(myeloid$myeloid_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  geom_text(data = tmp_de, mapping = aes(x = variable, y = cluster, label = signif),
            size = 8, fontface = 'bold', nudge_y = -0.09, nudge_x = 0.01) +
  scale_size(range = c(0,11), limits = c(0,100)) +
  scale_fill_gradientn(
    colors = my_cols,
    limits = c(min_expr, max_expr),
    labels = c(min_expr, 0, max_expr),
    breaks = c(min_expr, 0, max_expr),
    na.value = my_cols[100]) +
  scale_y_discrete(position = 'right') +
  theme(plot.margin = margin(0, 0, 0, 20, unit = 'mm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black', face = 'italic'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(5,0,10,0),
        legend.box = 'vertical',
        legend.box.margin = margin(0,0,-10,0),
        legend.position = 'right',
        legend.spacing.x = unit(x = 2, units = 'mm'),
        panel.border = element_rect(fill = NA, size = 1),
        panel.background = element_rect(fill = NA)) +
  guides(fill = guide_colorbar(title = 'z-score', 
                               barwidth = 1.25,
                               barheight = 4,
                               frame.colour = 'black', 
                               frame.linewidth = 1.25,
                               ticks.colour = 'black', 
                               ticks.linewidth = 1.25,
                               title.position = 'left'), 
         size = guide_legend(title = '% expression', 
                             override.aes = list(fill = 'black'),
                             title.position = 'left'))
macrophage_markers_dotplot
ggsave(filename = paste0(results_out, 'macrophage_markers_dotplot.tiff'),
       plot = macrophage_markers_dotplot, device = 'tiff', height = 3.5, width = 8)



### Generating new marker genes %%%%%%%%%

# Use Seurat wilcox test
macrophage_srat_genes <- FindAllMarkers(object = macrophage,
                                        assay = 'RNA',
                                        logfc.threshold = 0.25,
                                        only.pos = TRUE)
macrophage_markers <-  macrophage_srat_genes %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC) %>%
  .[['gene']]

# Calculate average and percent expression (scaled to within macrophage)
DefaultAssay(macrophage) <- 'RNA'
avg_exp <- data.frame(t(ScaleData(macrophage[['RNA']]@data, features = macrophage_markers)))
avg_exp <- cbind(avg_exp, 'myeloid_subcluster' = as.character(macrophage@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(t(macrophage[['RNA']]@counts[macrophage_markers,]))
pct_exp <- cbind(pct_exp, 'myeloid_subcluster' = as.character(macrophage@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
expr_colors <- colorRampPalette(colors = c('dodgerblue', 'white',' indianred'))(100)

# Plot
max_exp <- 2
macrophage_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(myeloid_subcluster = factor(myeloid_subcluster, levels = rev(levels(myeloid$myeloid_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  scale_size(range = c(0,10), limits = c(0,100)) +
  scale_fill_gradient2(low = "#2166AC", 
                       mid = "#F7F7F7",
                       high = "#B2182B", 
                       limits = c(NA, max_exp), 
                       na.value ="#B2182B",
                       breaks = seq(-5, 10, 1)) +
  scale_y_discrete(position = 'right') +
  theme(plot.margin = margin(0, 0, 0, 20, unit = 'mm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(5,0,10,0),
        legend.box = 'horizontal',
        legend.position = 'right',
        legend.spacing.x = unit(x = 2, units = 'mm'),
        panel.border = element_rect(fill = NA, size = 1),
        panel.background = element_rect(fill = NA)) +
  guides(fill = guide_colorbar(title = 'Scaled\nexpression', 
                               barwidth = 1.25,
                               frame.colour = 'black', 
                               frame.linewidth = 1.25,
                               ticks.colour = 'black', 
                               ticks.linewidth = 1.25,
                               title.position = 'left'), 
         size = guide_legend(title = '% expression', 
                             override.aes = list(fill = 'black'),
                             title.position = 'left'))
ggsave(filename = paste0(results_out, 'macrophage_markers_dotplot_new.tiff'),
       plot = macrophage_markers_dotplot, device = 'tiff', height = 3, width = 11)



# Macrophage DE markers heatmap (functional annotation) ------------------------

DefaultAssay(myeloid) <- 'RNA'
macrophage <- c('Neutrophil',
                'Monocyte',
                'Chemotaxis-Inducing Mac',
                'Inflammatory Mac',
                'Border-Associated Mac',
                'Dendritic',
                'Dividing Myeloid')
Idents(myeloid) <- 'myeloid_subcluster'
macrophage <- subset(myeloid, idents = macrophage)


# Use Seurat wilcox test
macrophage_srat_genes <- FindAllMarkers(object = macrophage,
                                        assay = 'RNA',
                                        logfc.threshold = 0.25,
                                        only.pos = TRUE)
write.table(file = paste0(results_out, 'macrophage_DE_wilcox.tsv'),
            x = macrophage_srat_genes, sep = '\t', row.names = TRUE, 
            col.names = NA)
macrophage_srat_genes <- read.table(file = paste0(results_out, 'macrophage_DE_wilcox.tsv'),
                                    sep = '\t', header = TRUE, row.names = 1)
macrophage_de_genes <- macrophage_srat_genes %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_logFC) %>%
  .[['gene']]

# Build heatmap
expr_data <- t(ScaleData(object = macrophage@assays$RNA@data[macrophage_de_genes,])) %>%
  cbind(macrophage@meta.data[c('myeloid_subcluster')]) %>%
  group_by(myeloid_subcluster) %>%
  summarise(across(.cols = all_of(macrophage_de_genes), .fns = mean)) %>%
  reshape2::melt(id.vars = 'myeloid_subcluster') %>%
  mutate(variable = factor(variable, levels = unique(variable)),
         myeloid_subcluster = factor(myeloid_subcluster, levels = rev(levels(myeloid_subcluster))))

expr_data$pval <- 1
for (i in 1:nrow(macrophage_srat_genes)) {
  tmp <- which(as.character(expr_data$myeloid_subcluster) == as.character(macrophage_srat_genes$cluster[i]) &
                 as.character(expr_data$variable) == as.character(macrophage_srat_genes$gene[i]))
  expr_data$pval[tmp] <- macrophage_srat_genes$p_val_adj[i]
}
expr_data$significant <- ifelse(test = expr_data$pval < 1e-10,
                                yes = '*',
                                no = ' ')

max_expr <- 2
# min_expr <- round(floor(min(expr_data$value)*10)/10,1)
min_expr <- -2
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))

macrophage_de_genes_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_tile(mapping = aes(fill = value), color = 'black', size = 0.3) +
  geom_text(mapping = aes(label = significant)) +
  scale_fill_gradientn(colors = myColors,
                       limits = c(min_expr, max_expr), 
                       na.value = myColors[100],
                       breaks = c(min_expr, 0, max_expr)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = 'black', face = 'italic'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.title = element_text(color = 'black', angle = 90, size = 14, hjust = 0),
        legend.text = element_text(size = 12, color = 'black'),
        legend.box.margin=margin(0,0,0,0),
        legend.margin = margin(0,0,0,-3, unit = 'mm')) +
  guides(fill = guide_colorbar(title = 'z-score',
                               title.position = 'left',
                               frame.colour = 'black',
                               ticks.colour = 'black',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               title.hjust = 0.5))
macrophage_de_genes_heatmap
ggsave(filename = paste0(results_out, 'macrophage_de_genes_heatmap.tiff'),
       plot = macrophage_de_genes_heatmap, device = 'tiff',
       height = 3.5, width = 14)





 # Complex Macrophage DE heatmap --------------------------------------------

# Subset cells
macrophage_names <- c('Neutrophil',
                      'Monocyte',
                      'Macrophage-A',
                      'Macrophage-B',
                      'BA-Macrophage',
                      'Dendritic',
                      'Div-Myeloid')
Idents(myeloid) <- 'myeloid_subcluster'
macrophage <- subset(myeloid, idents = macrophage_names)

# Use Seurat wilcox test
macrophage_srat_genes <- FindAllMarkers(object = macrophage,
                                        assay = 'RNA',
                                        logfc.threshold = 0.5,
                                        only.pos = TRUE)
macrophage_de_genes <- macrophage_srat_genes %>%
  group_by(cluster) %>%
  filter(avg_logFC > 0.5) %>%
  top_n(n = 40, wt = avg_logFC) %>%
  .[['gene']] %>%
  unique()

# Select which genes to annotate
label_genes <- c('G0s2','Mmp9','Chil3','Fn1','Ccl9','Tgfbi','Ccr2','F10','Ly6c2','Arg1','Spp1','Fabp5','Lgals3','Ccl2','Ms4a7','Apoe','Trem2','Mrc1','Cbr2','Lyve1','Folr2','Cd74','H2-Aa','Top2a','Birc5')


# Extract scaled expression
expr_data <- ScaleData(object = macrophage[['RNAcorrected']]@data[macrophage_de_genes,], 
                       scale.max = 4)
macrophage_meta <- macrophage@meta.data[c('myeloid_subcluster','time','sample_id')]

# Cell-level annotations
macrophage_cell_anno <- HeatmapAnnotation(
  'Subcluster' = macrophage_meta[['myeloid_subcluster']],
  'Time' = macrophage_meta[['time']],
  col = list('Subcluster' = myeloid_cols[names(myeloid_cols) %in% macrophage_names],
             'Time' = time_cols)
)

# Gene-level annotations
macrophage_gene_anno <- rowAnnotation(
  'Genes' = anno_mark(
    at = match(label_genes, rownames(expr_data)),
    labels = label_genes,
    side = 'right',
    labels_gp = gpar(fontsize = 10)
  )
)

# Build heatmap
tmp <- Heatmap(
  matrix = expr_data, 
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 0.6)(100),
  column_title = 'Differentially expressed genes among macrophage subclusters',
  column_title_gp = gpar(fontsize = 20),
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  column_dend_height = unit(20, units = 'mm'),
  row_dend_width = unit(20, units = 'mm'),
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(title = 'Scaled Expression',
                              title_gp = gpar(fontsize = 12),
                              title_position = 'leftcenter-rot',
                              labels = c('Low','High'),
                              at = c(min(expr_data), max(expr_data)),
                              labels_gp = gpar(fontsize = 10),
                              legend_height = unit(2.5, units = 'cm'),
                              grid_width = unit(0.5, units = 'cm'),
                              border = 'black',
                              title_gap = unit(1, units = 'cm'),
                              direction = 'vertical'),
  use_raster = TRUE,
  top_annotation = macrophage_cell_anno,
  right_annotation = macrophage_gene_anno
)

# Save heatmap
pdf(file = paste0(results_out, 'macrophage_heatmap_complex.pdf'),
    height = 8, width = 13)
tmp
dev.off()
svg(file = paste0(results_out, 'macrophage_heatmap_complex_small.svg'),
    height = 5, width = 9)
tmp
dev.off()




# Microglia population dynamics (functional annotation) -----------------------


counts <- data.frame(table(myeloid$myeloid_subcluster, myeloid$time))
names(counts) <- c('subcluster','time','count')
count_these <- c('Homeostatic Microglia',
                 'Inflammatory Microglia',
                 'Dividing Microglia',
                 'Migrating Microglia')
counts[['time']] <- plyr::mapvalues(
  x = counts[['time']],
  from = c('Uninjured'),
  to = c('Uninj')
)

microglia_prop <- counts %>%
  filter(subcluster %in% count_these) %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) + 
  geom_bar(aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black', 
           size = 0.75,
           position = 'fill') + 
  scale_y_continuous() + 
  scale_x_discrete() +
  scale_fill_manual(values = myeloid_cols) +
  ylab('Prop. of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = 'black'), 
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'),
        legend.spacing.y = unit(x = 2.5, units = 'mm')) + 
  guides(fill = guide_legend(ncol = 1))
microglia_prop
ggsave(filename = paste0(results_out, 'microglia_population_dynamics.tiff'),
       plot = microglia_prop, device = 'tiff', height = 2.75, width = 5.25)





# Macrophage population dynamics (functional annotation) ----------------------


counts <- data.frame(table(myeloid$myeloid_subcluster, myeloid$time))
names(counts) <- c('subcluster','time','count')
count_these <- c('Neutrophil',
                 'Monocyte',
                 'Chemotaxis-Inducing Mac',
                 'Inflammatory Mac',
                 'Border-Associated Mac',
                 'Dendritic',
                 'Dividing Myeloid')

macrophage_prop <- counts %>%
  filter(subcluster %in% count_these) %>%
  filter(time != 'Uninjured') %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) + 
  geom_bar(aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black', 
           size = 0.75,
           position = 'fill') + 
  scale_y_continuous() + 
  scale_x_discrete() +
  scale_fill_manual(values = myeloid_cols) +
  ylab('Prop. of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = 'black'), 
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'),
        legend.spacing.y = unit(x = 2.5, units = 'mm')) + 
  guides(fill = guide_legend(ncol = 1))
macrophage_prop
ggsave(filename = paste0(results_out, 'macrophage_population_dynamics.tiff'),
       plot = macrophage_prop, device = 'tiff', height = 2.75, width = 5.25)





# Abundance analysis ------------------------------------------------------

time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

myeloid_subcluster_counts <- data.frame(table(myeloid$myeloid_subcluster, myeloid$sample_id))
colnames(myeloid_subcluster_counts) <- c('myeloid_subcluster','sample_id','count')
myeloid_subcluster_counts[['time']] <- substr(x = myeloid_subcluster_counts[['sample_id']],
                                                start = 1,
                                                stop = regexpr(pattern = '_', text = myeloid_subcluster_counts[['sample_id']]) -1)
myeloid_subcluster_counts[['time']] <- plyr::mapvalues(x = myeloid_subcluster_counts[['time']], from ='uninj', to = 'Uninjured')
myeloid_subcluster_counts[['time']] <- factor(myeloid_subcluster_counts[['time']], levels = levels(myeloid$time))
myeloid_subcluster_counts_plot <- myeloid_subcluster_counts %>%
  ggplot(mapping = aes(x = sample_id, y = count)) +
  geom_bar(mapping = aes(fill = time), color = 'black', stat = 'identity') +
  labs(title = 'Cell abundance across time and sample') +
  ylab(label = 'Number of cells') +
  facet_wrap(. ~ myeloid_subcluster, scales = 'free_y', ncol = 5) +
  scale_fill_manual(values = time_cols) +
  theme(plot.title = element_text(size = 14, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(fill = guide_legend(title = 'Time after\ninjury'))
ggsave(filename = paste0(results_out, 'myeloid_subcluster_counts_bySample.tiff'),
       plot = myeloid_subcluster_counts_plot, device = 'tiff', height = 6.5, width = 15)



