
######## SCI cell-type annotation via marker identification ########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
require('scran')
# data_in <- './data/'
# data_out <- './data/data_integration/'
results_out <- './results/sci_annotation_markers/'
ref_in <- './ref/'
ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)

sci <- readRDS(file = './data/sci.rds')



# Canonical markers UMAP ----------------------------------------------------

DefaultAssay(sci) <- 'SCT'

# These are marker genes identified in previous studies. Multiple sources.
canonical_markers <- list('Microglia' = c('Cx3cr1','Tmem119'),
                          'Leukocyte' = c('Itgam','Ptprc'),
                          'Monocyte' = c('Ly6c2','Ccr2'),
                          'Neutrophil' = c('Ly6g','Mpo'),
                          'Dendritic' = c('Cd74'),
                          'Lymphocyte' = c('Cd3e'),
                          'Div-Myeloid' = c('Mki67','Cdk1'),
                          'Endothelial' = c('Pecam1','Cldn5'),
                          'Fibroblast' = c('Col1a1'),
                          'Pericyte' = c('Cspg4','Pdgfrb'),
                          'Ependymal' = c('Foxj1'),
                          'Astrocyte' = c('Aldoc','Slc1a2','Gfap'),
                          'Oligodendrocyte' = c('Plp1','Opalin'),
                          'OPC' = c('Pdgfra','Olig2'),
                          'Neuron' = c('Snap25'))
canonical_markers <- unlist(canonical_markers, use.names = FALSE)

umap_coordinates <- FetchData(sci, vars = c('UMAP_1','UMAP_2'), slot = 'data')
label_x <- min(umap_coordinates[['UMAP_1']]) * 1
label_y <- max(umap_coordinates[['UMAP_2']]) * 0.85
expr_colors <- colorRampPalette(colors = c('grey85', 'red3'))(100)

canonical_markers_umap <- vector(mode = 'list', length = length(canonical_markers))
names(canonical_markers_umap) <- canonical_markers

for(ii in 1:length(canonical_markers)) {
  gene <- canonical_markers[ii]
  # gather data
  expr <- FetchData(sci, vars = gene, slot = 'data')
  expr <- cbind(expr, umap_coordinates)
  expr <- expr[order(expr[[gene]], decreasing = FALSE),]
  max_expr <- ceiling(max(expr[[gene]])*10)/10
  colnames(expr) <- gsub(pattern = '-', replacement = '.', x = colnames(expr))
  gene <- gsub(pattern = '-', replacement = '.', x = gene)
  # assemble gg
  canonical_markers_umap[[gene]] <- expr %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes_string(color = gene), size = 0.1, alpha = 0.4) +
    geom_text(data = data.frame('gene' = gene,
                                'x_pos' = label_x,
                                'y_pos' = label_y),
              mapping = aes(x = x_pos, y = y_pos, label = gene),
              fontface = 'italic',
              size = 6,
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
          legend.position = c(0.875, 0.825), 
          legend.key.size = unit(0.25, units = 'cm'), 
          legend.spacing.x = unit(0.1, units = 'cm'),
          legend.text = element_text(size = 12, color = 'black')) +
    guides(color = guide_colorbar(frame.colour = 'black', 
                                  ticks = FALSE,
                                  barwidth = 0.8,
                                  frame.linewidth = 1))
}
# save
canonical_markers_umap_grid <- cowplot::plot_grid(
  plotlist = canonical_markers_umap, 
  ncol = 3, 
  byrow = TRUE
)
ggsave(filename = paste0(results_out, 'canonical_markers_umap.tiff'),
       plot = canonical_markers_umap_grid, device = 'tiff',
       height = 24, width = 9)




# DE marker genes calculation -----------------------------------------------


# Use wilcox test while blocking against Time (group) to identify marker genes
# per cluster. scran::findMarkers() performs pair-wise for all cluster 
# combinations. 
rna_sce <- as.SingleCellExperiment(sci, assay = 'RNA')
markers_time <- findMarkers(x = rna_sce,
                            groups = colData(rna_sce)[['default_cluster']],
                            test.type = 'wilcox',
                            direction = 'up',
                            pval.type = 'all',
                            block = colData(rna_sce)[['time']])
min_fdr <- 1e-03
sig_markers_time <- lapply(X = markers_time,
                           FUN = function(x) {
                             return(x[x[['FDR']] < min_fdr,])
                           })

Idents(sci) <- 'default_cluster'
seurat_markers <- FindAllMarkers(object = sci,
                                 assay = 'RNA', 
                                 logfc.threshold = 1,
                                 only.pos = TRUE,
                                 min.pct = 0.4)
write.table(seurat_markers, file = paste0(results_out, 'seurat_markers_default_clusters.tsv'),
            sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
seurat_markers_top <- seurat_markers %>%
  mutate('pct.diff' = pct.1 - pct.2) %>%
  filter(p_val_adj < 1e-10) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = pct.diff)



# Cell-type annotation ----------------------------------------------------

# Based seurat findmarkers results
sci$celltype <- plyr::mapvalues(
  x = sci$default_cluster,
  from = c(0:29),
  to = c('Macrophage', 
         'Microglia', 
         'Microglia', 
         'Endothelial',
         'Microglia',
         'Microglia',
         'Endothelial',
         'Ependymal',
         'Macrophage',
         'Ependymal',
         'Endothelial',
         'Monocyte',
         'OPC',
         'Dendritic',
         'Neutrophil', # ?
         'Div-Myeloid',
         'Oligodendrocyte',
         'Ependymal',
         'Pericyte',
         'OPC',
         'Astrocyte',
         'Fibroblast',
         'Microglia',
         'Pericyte',
         'Endothelial',
         'Neuron',
         'Div-Myeloid',
         'Lymphocyte',
         'Endothelial',
         'Div-Myeloid')
)
sci$celltype <- factor(x = sci$celltype, levels = c('Neutrophil',
                                                    'Monocyte',
                                                    'Macrophage',
                                                    'Dendritic',
                                                    'Microglia',
                                                    'Div-Myeloid',
                                                    'Fibroblast',
                                                    'Endothelial',
                                                    'Pericyte',
                                                    'OPC',
                                                    'Oligodendrocyte',
                                                    'Astrocyte',
                                                    'Ependymal',
                                                    'Lymphocyte',
                                                    'Neuron'))
Idents(sci) <- 'celltype'

# Cell-type UMAP across time ----------------------------------------------

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
Idents(object = sci) <- 'celltype'

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
                    legend.title = element_text(size = 14, color = 'black', face = 'bold'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 14, color = 'black'))

# cell-type annotation UMAP
celltype_counts <- table(sci$celltype)
celltype_label <- paste0(names(celltype_counts), ' (', celltype_counts, ')')
names(celltype_label) <- names(celltype_counts)
celltype_umap <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','celltype')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = celltype), size = 0.2, alpha = 0.5) +
  scale_color_manual(values = sci_cols, 
                     breaks = names(celltype_label),
                     label = celltype_label) +
  xlab(label = 'UMAP 1') +
  ylab(label = 'UMAP 2') +
  umap_theme +
  theme(legend.text = element_text(size = 16, color = 'black')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
ggsave(filename = paste0(results_out, 'celltype_annotation_umap.tiff'),
       plot = celltype_umap, device = 'tiff', height = 6, width = 9.25)

# cell-type annotation split by time UMAP
celltype_split_umap <- FetchData(sci, vars = c('UMAP_1','UMAP_2','celltype','time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = celltype), size = 0.2, alpha = 0.5) +
  facet_wrap(. ~ time, ncol = 2) +
  scale_color_manual(values = sci_cols, breaks = names(celltype_label), label = celltype_label) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
ggsave(filename = paste0(results_out, 'celltype_annotation_split_umap.tiff'),
       plot = celltype_split_umap, device = 'tiff', height = 6, width = 9)





# DE markers dot plot -----------------------------------------------------

# Import
de_markers <- list('Neutrophil' = c('S100a9','Mmp9','Ly6g','Cd177','Ltf'),
                   'Monocyte' = c('Ly6c2','Ccr2','F10','Plac8'),
                   'Macrophage' = c('Thbs1','Ms4a7','Pf4','Fabp4','Gpnmb'),
                   'Microglia' = c('Gpr84','Ptgs1','P2ry12','Gpr34','Lag3','Siglech','Cst7'),
                   'Div-Myeloid' = c('Top2a', 'Mki67', 'Cdk1', 'Ccnb2'),
                   'Fibroblast' = c('Col1a1','Col6a1','Dcn','Lum','Postn'),
                   'Endothelial' = c('Pecam1','Kdr','Flt1'),
                   'Pericyte' = c('Kcnj8', 'Gm13861', 'Higd1b', 'Abcc9', 'Rgs5'),
                   'OPC' = c('Lhfpl3', 'Tnr', 'Igsf21', 'Neu4', 'Gpr17'),
                   'Oligodendrocyte' = c('Gjb1', 'Ermn', 'Mog', 'Tmem125', 'Hapln2'),
                   'Astrocyte' = c('Agt','Ntsr2','Acsbg1','Slc6a11','Slc7a10','Fgfr3'),
                   'Ependymal' = c('Cfap126','Fam183b','Tmem212','Pifo','Tekt1','Dnah12'))

DefaultAssay(sci) <- 'SCT'
avg_exp <- ScaleData(sci[['SCT']]@data, features = unlist(de_markers, use.names = FALSE))
avg_exp <- cbind(t(avg_exp), sci@meta.data[,c('celltype','time')]) %>%
  reshape2::melt(id.vars = c('celltype','time')) %>%
  group_by(celltype, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- sci[['RNA']]@counts[unlist(de_markers, use.names = FALSE),]
pct_exp <- cbind(t(pct_exp), sci@meta.data[,c('celltype','time')]) %>%
  reshape2::melt(id.vars = c('celltype','time')) %>%
  group_by(celltype, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
expr_colors <- colorRampPalette(colors = c('grey85', 'red3'))(100)

de_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  filter(celltype %in% names(de_markers)) %>%
  mutate(time = factor(time, levels = rev(levels(time)))) %>%
  ggplot(mapping = aes(x = variable, y = time)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  facet_grid(celltype ~ ., drop = TRUE, switch = 'y') +
  scale_size(range = c(0,7), limits = c(0,100)) +
  scale_fill_gradientn(colors = expr_colors,
                       limits = c(NA, 3.5),
                       breaks = seq(-5, 10, 1),
                       na.value = expr_colors[length(expr_colors)]) +
  scale_y_discrete(position = 'right') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'bold', size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(5,0,10,0),
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
ggsave(filename = paste0(results_out, 'celltype_markers_dotplot.tiff'),
       plot = de_markers_dotplot, device = 'tiff', height = 10, width = 16.8)






# DE markers dot plot (U-vascular removed) ------------------------------------


# We first have to map subcluster identities back to full SCI dataset.
subcluster_paths <- list(
  'myeloid_subclusters' = './results/myeloid_annotation_markers/myeloid_subcluster.tsv',
  'vascular_subclusters' = './results/vascular_annotation_markers/vascular_subcluster.tsv',
  'macroglia_subclusters' = './results/macroglia_annotation_markers/macroglia_subcluster.tsv'
)
# Load
subclusters_bySet <- lapply(
  X = subcluster_paths,
  FUN = read.table,
  sep = '\t',
  header = TRUE,
  row.names = 1
)
# Extract labels and barcodes
subclusters <- Reduce(
  f = c, 
  x = sapply(
    X = subclusters_bySet,
    FUN = `[[`,
    1
  )
)
barcodes <- Reduce(
  f = c,
  x = sapply(
    X = subclusters_bySet,
    FUN = rownames
  )
)
names(subclusters) <- barcodes
rm(barcodes, subcluster_paths)

# Add to sci metadata
sci@meta.data[['subcluster']] <- subclusters[match(x = rownames(sci@meta.data), 
                                                   table = names(subclusters))]
# label unassigned cells
sci@meta.data[['subcluster']][is.na(sci@meta.data[['subcluster']])] <- 'Unassigned'
Idents(sci) <- 'subcluster'


# Import de markers
de_markers <- list('Neutrophil' = c('S100a9','Mmp9','Ly6g','Cd177','Ltf'),
                   'Monocyte' = c('Ly6c2','Ccr2','F10','Plac8'),
                   'Macrophage' = c('Thbs1','Ms4a7','Pf4','Fabp4','Gpnmb'),
                   'Microglia' = c('Gpr84','Ptgs1','P2ry12','Gpr34','Lag3','Siglech','Cst7'),
                   'Div-Myeloid' = c('Top2a', 'Mki67', 'Cdk1', 'Ccnb2'),
                   'Fibroblast' = c('Col1a1','Col6a1','Dcn','Lum','Postn'),
                   'Endothelial' = c('Pecam1','Kdr','Flt1'),
                   'Pericyte' = c('Kcnj8', 'Gm13861', 'Higd1b', 'Abcc9', 'Rgs5'),
                   'OPC' = c('Lhfpl3', 'Tnr', 'Igsf21', 'Neu4', 'Gpr17'),
                   'Oligodendrocyte' = c('Gjb1', 'Ermn', 'Mog', 'Tmem125', 'Hapln2'),
                   'Astrocyte' = c('Agt','Ntsr2','Acsbg1','Slc6a11','Slc7a10','Fgfr3'),
                   'Ependymal' = c('Cfap126','Fam183b','Tmem212','Pifo','Tekt1','Dnah12'))

DefaultAssay(sci) <- 'SCT'
remove_these <- c('Unassigned','U-Vascular')
avg_exp <- ScaleData(sci[['SCT']]@data[,!sci@meta.data[['subcluster']] %in% remove_these], features = unlist(de_markers, use.names = FALSE))
avg_exp <- cbind(t(avg_exp), sci@meta.data[!sci@meta.data[['subcluster']] %in% remove_these,c('celltype','time')]) %>%
  reshape2::melt(id.vars = c('celltype','time')) %>%
  group_by(celltype, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- sci[['RNA']]@counts[unlist(de_markers, use.names = FALSE),!sci@meta.data[['subcluster']] %in% remove_these]
pct_exp <- cbind(t(pct_exp), sci@meta.data[!sci@meta.data[['subcluster']] %in% remove_these, c('celltype','time')]) %>%
  reshape2::melt(id.vars = c('celltype','time')) %>%
  group_by(celltype, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
expr_colors <- colorRampPalette(colors = c('grey85', 'red3'))(100)

de_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  filter(celltype %in% names(de_markers)) %>%
  mutate(time = factor(time, levels = rev(levels(time)))) %>%
  ggplot(mapping = aes(x = variable, y = time)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  facet_grid(celltype ~ ., drop = TRUE, switch = 'y') +
  scale_size(range = c(0,7), limits = c(0,100)) +
  scale_fill_gradientn(colors = expr_colors,
                       limits = c(NA, 3.5),
                       breaks = seq(-5, 10, 1),
                       na.value = expr_colors[length(expr_colors)]) +
  scale_y_discrete(position = 'right') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'bold', size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(5,0,10,0),
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
ggsave(filename = paste0(results_out, 'celltype_markers_dotplot_filtered.tiff'),
       plot = de_markers_dotplot, device = 'tiff', height = 10, width = 16.8)






# Individual marker UMAP --------------------------------------------------


DefaultAssay(sci) <- 'SCT'

# single markers
single_markers <- list('Neutrophil' = 'Ltf',
                       'Monocyte' = 'Ccr2',
                       'Macrophage' = 'Ms4a7',
                       'Microglia' = 'Gpr84',
                       'Div-Myeloid' = 'Mki67',
                       'Fibroblast' = 'Col1a1',
                       'Endothelial' = 'Pecam1',
                       'Pericyte' = 'Higd1b',
                       'OPC' = 'Tnr',
                       'Oligodendrocyte' = 'Ermn',
                       'Astrocyte' = 'Slc6a11',
                       'Ependymal' = 'Pifo')
umap_coordinates <- FetchData(sci, vars = c('UMAP_1','UMAP_2'), slot = 'data')
label_x <- min(umap_coordinates[['UMAP_1']]) * 1
label_y <- max(umap_coordinates[['UMAP_2']]) * 0.85
expr_colors <- colorRampPalette(colors = c('grey85', 'red3'))(100)

single_markers_umap <- vector(mode = 'list', length = length(single_markers))
names(single_markers_umap) <- single_markers

for(ii in 1:length(single_markers)) {
  gene <- single_markers[[ii]]
  celltype <- names(single_markers)[ii]
  # gather data
  expr <- FetchData(sci, vars = gene, slot = 'data')
  expr <- cbind(expr, umap_coordinates)
  expr <- expr[order(expr[[gene]], decreasing = FALSE),]
  max_expr <- ceiling(max(expr[[gene]])*10)/10
  # assemble gg
  single_markers_umap[[gene]] <- expr %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes_string(color = gene), size = 0.1, alpha = 0.2) +
    geom_text(data = data.frame('gene' = gene,
                                'x_pos' = label_x,
                                'y_pos' = label_y),
              mapping = aes(x = x_pos, y = y_pos, label = gene),
              fontface = 'italic',
              size = 4.5,
              hjust = 0) +
    labs(title = celltype) +
    scale_color_gradientn(colors = expr_colors,
                          breaks = c(0, max_expr),
                          limits = c(0, max_expr)) +
    theme(plot.title = element_text(size = 14, color = 'black', face = 'bold'),
          panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.875, 0.825), 
          legend.key.size = unit(0.25, units = 'cm'), 
          legend.background = element_rect(fill = NA),
          legend.spacing.x = unit(0.1, units = 'cm'),
          legend.text = element_text(size = 10, color = 'black')) +
    guides(color = guide_colorbar(frame.colour = 'black', 
                                  ticks = FALSE,
                                  barwidth = 0.8,
                                  frame.linewidth = 1))
}
# save
single_markers_umap_grid <- cowplot::plot_grid(plotlist = single_markers_umap, ncol = 6)
ggsave(filename = paste0(results_out, 'celltype_single_markers_umap.tiff'),
       plot = single_markers_umap_grid, device = 'tiff',
       height = 5, width = 15)




# Abundance analysis ------------------------------------------------------

time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

celltype_counts <- data.frame(table(sci$celltype, sci$sample_id))
colnames(celltype_counts) <- c('celltype','sample_id','count')
celltype_counts[['time']] <- substr(x = celltype_counts[['sample_id']],
                                    start = 1,
                                    stop = regexpr(pattern = '_', text = celltype_counts[['sample_id']]) -1)
celltype_counts[['time']] <- plyr::mapvalues(x = celltype_counts[['time']], from ='uninj', to = 'Uninjured')
celltype_counts[['time']] <- factor(celltype_counts[['time']], levels = levels(sci$time))
celltype_counts_plot <- celltype_counts %>%
  ggplot(mapping = aes(x = sample_id, y = count)) +
  geom_bar(mapping = aes(fill = time), color = 'black', stat = 'identity') +
  labs(title = 'Cell abundance across time and sample') +
  ylab(label = 'Number of cells') +
  facet_wrap(. ~ celltype, scales = 'free_y', ncol = 5) +
  scale_fill_manual(values = time_cols) +
  theme(plot.title = element_text(size = 14, color = 'black', face = 'bold'),
        panel.border = element_rect(fill = NA, color = 'black'),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(fill = guide_legend(title = 'Time after\ninjury'))
ggsave(filename = paste0(results_out, 'celltype_counts.tiff'),
       plot = celltype_counts_plot, device = 'tiff', height = 6.5, width = 15)








