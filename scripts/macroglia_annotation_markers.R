
######## Macroglia cell-type annotation via marker identification ########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
require('scran')
require('parallel')
require('doParallel')
require('foreach')
require('ComplexHeatmap')
# data_in <- './data/'
# data_out <- './data/data_integration/'
results_out <- './results/macroglia_annotation_markers/'
ref_in <- './ref/'
ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)

macroglia <- readRDS(file = './data/macroglia.rds')






# Additional markers UMAP ----------------------------------------------------


DefaultAssay(macroglia) <- 'RNA'

# These are marker genes identified in previous studies. Multiple sources.
addtl_markers <- list('Astrocyte' = c('Slc1a2', 'Slc1a3', 'Atp1a2','Aqp4','Gfap','Aldh1l1'),
                      'Ependymal' = c('Foxj1','Tmem212','Cfap126'),
                      'OPC' = c('Olig2','Cspg4','Pdgfra','Sox10'),
                      'Oligodendrocyte' = c('Opalin','Mog','Plp1','Mag'),
                      'Pre-Oligo' = c('Mycl','Bmp4'),
                      'Div-OPC' = c('Cdk1','Top2a','Mki67'),
                      'Astroependymal' = c('Crym','Vim'))
addtl_markers <- unlist(addtl_markers, use.names = FALSE)

umap_coordinates <- FetchData(macroglia, vars = c('UMAP_1','UMAP_2'), slot = 'data')
label_x <- min(umap_coordinates[['UMAP_1']]) * 0.975
label_y <- min(umap_coordinates[['UMAP_2']]) * 0.85
legend_x <- 0.15
legend_y <- 0.85
expr_colors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = 'Spectral'))(100))

addtl_markers_umap <- vector(mode = 'list', length = length(addtl_markers))
names(addtl_markers_umap) <- addtl_markers

for (ii in 1:length(unlist(addtl_markers, use.names = FALSE))) {
  gene <- addtl_markers[ii]
  # gather data
  expr <- FetchData(macroglia, vars = gene, slot = 'data')
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
    geom_point(mapping = aes_string(color = gene), size = 0.5, alpha = 0.4) +
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
          legend.text = element_text(size = 12, color = 'black')) +
    guides(color = guide_colorbar(frame.colour = 'black', 
                                  ticks = FALSE,
                                  barwidth = 0.8,
                                  frame.linewidth = 1))
}
# save
addtl_markers_umap_grid <- cowplot::plot_grid(plotlist = addtl_markers_umap, ncol = 3)
ggsave(filename = paste0(results_out, 'addtl_markers_umap.tiff'),
       plot = addtl_markers_umap_grid, device = 'tiff',
       height = 19, width = 7.5)






# DE marker genes calculation -----------------------------------------------


# Use wilcox test while blocking against Time (group) to identify marker genes
# per cluster. scran::findMarkers() performs pair-wise for all cluster 
# combinations. 
rna_sce <- as.SingleCellExperiment(macroglia, assay = 'RNA')
markers_time <- findMarkers(x = rna_sce,
                            groups = colData(rna_sce)[['integrated_snn_res.0.3']],
                            test.type = 'wilcox',
                            direction = 'up',
                            pval.type = 'all',
                            block = colData(rna_sce)[['time']])
min_fdr <- 1e-03
sig_markers_time <- vector(mode = 'list', length = length(markers_time))
names(sig_markers_time) <- paste0('C', names(markers_time))
for(ii in 1:length(markers_time)) {
  cluster_id <- paste0('C', names(markers_time)[ii])
  tmp <- markers_time[[ii]]
  tmp <- tmp[tmp[['FDR']] < min_fdr, c('p.value', 'FDR')]
  tmp[['gene']] <- rownames(tmp)
  tmp[['cluster']] <- cluster_id
  sig_markers_time[[cluster_id]] <- tmp
}
sig_markers_time <- data.frame(do.call(rbind, sig_markers_time))
# sig_markers_time <- lapply(X = markers_time,
#                            FUN = function(x) {
#                              return(x[x[['FDR']] < min_fdr,])
#                            })
# sig_markers_time <- do.call(rbind, sig_markers_time)

Idents(macroglia) <- 'default_cluster'
seurat_markers <- FindAllMarkers(object = macroglia,
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

macroglia@meta.data[['macroglia_subcluster']] <- plyr::mapvalues(
  x = macroglia@meta.data[['integrated_snn_res.0.4']],
  from = levels(macroglia@meta.data[['integrated_snn_res.0.4']]),
  to = c('Ependymal-B',
         'Ependymal-A',
         'OPC-A',
         'Astroependymal',
         'Oligodendrocyte',
         'Astrocyte',
         'OPC-B',
         'Pre-Oligo',
         'Oligodendrocyte',
         'Div-OPC')
)
macroglia@meta.data[['macroglia_subcluster']] <- factor(
  x = macroglia@meta.data[['macroglia_subcluster']],
  levels = c('Ependymal-A',
             'Ependymal-B',
             'Astroependymal',
             'Astrocyte',
             'OPC-A',
             'OPC-B',
             'Div-OPC',
             'Pre-Oligo',
             'Oligodendrocyte')
)
Idents(macroglia) <- 'macroglia_subcluster'

# Save barcode-identities data.frame
write.table(x = macroglia@meta.data['macroglia_subcluster'],
            file = paste0(results_out, 'macroglia_subcluster.tsv'),
            sep = '\t', row.names = TRUE, col.names = NA)





# Cell-type UMAP across time ----------------------------------------------

# macroglia colors
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


Idents(object = macroglia) <- 'macroglia_subcluster'

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 18, color = 'black', face = 'bold'),
                    legend.title = element_text(size = 14, color = 'black', face = 'bold'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 14, color = 'black'))

# cell-type annotation UMAP 
macroglia_subcluster_counts <- table(macroglia$macroglia_subcluster)
macroglia_subcluster_label <- paste0(names(macroglia_subcluster_counts), ' (', macroglia_subcluster_counts, ')')
names(macroglia_subcluster_label) <- names(macroglia_subcluster_counts)
macroglia_subcluster_umap <- FetchData(object = macroglia, vars = c('UMAP_1','UMAP_2','macroglia_subcluster')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = macroglia_subcluster), size = 1, alpha = 0.5) +
  scale_color_manual(values = macroglia_cols, 
                     breaks = names(macroglia_subcluster_label),
                     label = macroglia_subcluster_label) +
  xlab(label = 'UMAP 1') +
  ylab(label = 'UMAP 2') +
  umap_theme +
  theme(legend.text = element_text(size = 16, color = 'black'),
        legend.box.margin = margin(0,0,0,-5,unit='mm')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
macroglia_subcluster_umap
ggsave(filename = paste0(results_out, 'macroglia_subcluster_annotation_umap.tiff'),
       plot = macroglia_subcluster_umap, device = 'tiff', height = 6, width = 9.25)

# cell-type annotation split by time UMAP (figure 4b)
macroglia_subcluster_split_umap <- FetchData(macroglia, vars = c('UMAP_1','UMAP_2','macroglia_subcluster','time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = macroglia_subcluster), size = 1, alpha = 0.5) +
  facet_wrap(. ~ time, ncol = 2) +
  scale_color_manual(values = macroglia_cols, breaks = names(macroglia_subcluster_label), label = macroglia_subcluster_label) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
ggsave(filename = paste0(results_out, 'macroglia_subcluster_annotation_split_umap.tiff'),
       plot = macroglia_subcluster_split_umap, device = 'tiff', height = 6, width = 9)







# DE markers dot plot -----------------------------------------------------

# Calculate DE
de_genes <- FindAllMarkers(object = macroglia,
                           assay = 'RNA',
                           logfc.threshold = 0.5,
                           only.pos = TRUE)
# Differences between Ependymal-A and Ependymal-B
ependymalA_ependymalB <- FindMarkers(macroglia,
                                     ident.1 = 'Ependymal-A',
                                     ident.2 = 'Ependymal-B',
                                     assay = 'RNA',
                                     logfc.threshold = 0.5)
# Differences between OPC-A and OPC-B
opcA_opcB <- FindMarkers(macroglia,
                         ident.1 = 'OPC-A',
                         ident.2 = 'OPC-B',
                         assay = 'RNA',
                         logfc.threshold = 0.25)
union(de_genes$gene[de_genes$cluster == 'OPC-B'], rownames(opcA_opcB)[opcA_opcB$avg_logFC < 0])


# Select marker genes
de_markers <- c('Pan-Ependymal' = 'Foxj1',
                'Ependymal-A' = 'Tmem212',
                'Ependymal-B' = 'Rbp1',
                'Astroependymal' = 'Crym',
                'Astrocyte' = 'Agt',
                'OPC-A1' = 'Tnr',
                'OPC-B' = 'Tnc',
                'Div-OPC' = 'Cdk1',
                'Pre-Oligo' = 'Bmp4',
                'Oligodendrocyte' = 'Mag')

DefaultAssay(macroglia) <- 'RNA'
avg_exp <- data.frame(t(ScaleData(macroglia[['RNA']]@data, features = de_markers)))
avg_exp <- cbind(avg_exp, 'macroglia_subcluster' = as.character(macroglia@meta.data[,c('macroglia_subcluster')])) %>%
  reshape2::melt(id.vars = c('macroglia_subcluster')) %>%
  group_by(macroglia_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(t(macroglia[['RNA']]@counts[unlist(de_markers, use.names = FALSE),]))
pct_exp <- cbind(pct_exp, 'macroglia_subcluster' = as.character(macroglia@meta.data[,c('macroglia_subcluster')])) %>%
  reshape2::melt(id.vars = c('macroglia_subcluster')) %>%
  group_by(macroglia_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)

max_expr <- 3
min_expr <- floor(min(avg_exp$avg.exp)*10)/10
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min_expr, 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))

de_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  filter(macroglia_subcluster != 'U-macroglia') %>%
  mutate(macroglia_subcluster = factor(macroglia_subcluster, levels = rev(levels(macroglia$macroglia_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = macroglia_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  scale_size(range = c(0,10), limits = c(0,100)) +
  scale_fill_gradientn(
    colors = myColors,
    values = scales::rescale(x = myBreaks, to = c(0,1)),
    limits = c(min_expr, max_expr), 
    na.value = myColors[100],
    breaks = c(min_expr, 0, max_expr)) +
  scale_y_discrete(position = 'right') +
  theme(plot.margin = margin(0, 0, 0, 20, unit = 'mm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'bold', size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 12, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.box.margin = margin(10,0,0,0,unit = 'mm'),
        legend.margin = margin(0,0,0,0, unit = 'mm'),
        legend.box = 'vertical',
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
de_markers_dotplot
ggsave(filename = paste0(results_out, 'macroglia_subcluster_markers_dotplot.tiff'),
       plot = de_markers_dotplot, device = 'tiff', height = 3.75, width = 6.75)





# DE markers dot plot (significant p-values) ----------------------------------

# Calculate DE
de_genes <- FindAllMarkers(object = macroglia,
                           assay = 'RNA',
                           logfc.threshold = 0.5,
                           only.pos = TRUE)
# Differences between Ependymal-A and Ependymal-B
ependymalA_ependymalB <- FindMarkers(macroglia,
                                     ident.1 = 'Ependymal-A',
                                     ident.2 = 'Ependymal-B',
                                     assay = 'RNA',
                                     logfc.threshold = 0.5)
# Differences between OPC-A and OPC-B
opcA_opcB <- FindMarkers(macroglia,
                         ident.1 = 'OPC-A',
                         ident.2 = 'OPC-B',
                         assay = 'RNA',
                         logfc.threshold = 0.25)
union(de_genes$gene[de_genes$cluster == 'OPC-B'], rownames(opcA_opcB)[opcA_opcB$avg_logFC < 0])


# Select marker genes
de_markers <- c('Pan-Ependymal' = 'Foxj1',
                'Ependymal-A' = 'Tmem212',
                'Ependymal-B' = 'Rbp1',
                'Astroependymal' = 'Crym',
                'Astrocyte' = 'Agt',
                'OPC-A1' = 'Tnr',
                'OPC-B' = 'Tnc',
                'Div-OPC' = 'Cdk1',
                'Pre-Oligo' = 'Bmp4',
                'Oligodendrocyte' = 'Mag')
tmp_de <- FindAllMarkers(
  object = macroglia,
  assay = 'RNA',
  slot = 'data',
  features = de_markers,
  only.pos = TRUE
)
tmp_de$p_val_adj <- signif(tmp_de$p_val_adj, digits = 2)
tmp_de$variable <- tmp_de$gene
tmp_de$signif <- ifelse(test = tmp_de$p_val < 1e-10,
                        yes = '*',
                        no = '')

DefaultAssay(macroglia) <- 'RNA'
avg_exp <- data.frame(t(ScaleData(macroglia[['RNA']]@data, features = de_markers)))
avg_exp <- cbind(avg_exp, 'macroglia_subcluster' = as.character(macroglia@meta.data[,c('macroglia_subcluster')])) %>%
  reshape2::melt(id.vars = c('macroglia_subcluster')) %>%
  group_by(macroglia_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(t(macroglia[['RNA']]@counts[unlist(de_markers, use.names = FALSE),]))
pct_exp <- cbind(pct_exp, 'macroglia_subcluster' = as.character(macroglia@meta.data[,c('macroglia_subcluster')])) %>%
  reshape2::melt(id.vars = c('macroglia_subcluster')) %>%
  group_by(macroglia_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)

# max_expr <- 3
# min_expr <- floor(min(avg_exp$avg.exp)*10)/10
max_expr <- 3
min_expr <- -3
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
# myBreaks <- c(seq(min_expr, 0, length.out = 50),
#               seq(max_expr/100, max_expr, length.out = 50))

de_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  filter(macroglia_subcluster != 'U-macroglia') %>%
  mutate(macroglia_subcluster = factor(macroglia_subcluster, levels = rev(levels(macroglia$macroglia_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = macroglia_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  geom_text(data = tmp_de, mapping = aes(x = variable, y = cluster, label = signif),
            size = 8, fontface = 'bold', nudge_y = -0.09, nudge_x = 0.01) +
  scale_size(range = c(0,10), limits = c(0,100)) +
  scale_fill_gradientn(
    colors = myColors,
    limits = c(min_expr, max_expr), 
    na.value = myColors[100],
    breaks = c(min_expr, 0, max_expr)) +
  scale_y_discrete(position = 'right') +
  theme(plot.margin = margin(0, 0, 0, 20, unit = 'mm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'bold', size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 12, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.box.margin = margin(10,0,0,0,unit = 'mm'),
        legend.margin = margin(0,0,0,0, unit = 'mm'),
        legend.box = 'vertical',
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
de_markers_dotplot
ggsave(filename = './results/revision_figures/macroglia_subcluster_markers_dotplot.tiff',
       plot = de_markers_dotplot, device = 'tiff', height = 3.75, width = 6.75)







# DE markers heatmap ------------------------------------------------------

Idents(macroglia) <- 'macroglia_subcluster'
DefaultAssay(macroglia) <- 'RNA'

# Use wilcox test while blocking against Time (group) to identify marker genes
# per cluster. scran::findMarkers() performs pair-wise for all cluster
# combinations.
rna_sce <- as.SingleCellExperiment(macroglia, assay = 'RNA')
macroglia_markers_scran <- findMarkers(x = rna_sce,
                                      groups = colData(rna_sce)[['macroglia_subcluster']],
                                      test.type = 'wilcox',
                                      direction = 'up',
                                      pval.type = 'all')
for(ii in 1:length(macroglia_markers_scran)) {
  cluster_id <- names(macroglia_markers_scran)[ii]
  tmp <- macroglia_markers_scran[[ii]]
  tmp <- tmp[tmp[['FDR']] < min_fdr, c('p.value', 'FDR')]
  tmp[['gene']] <- rownames(tmp)
  tmp[['cluster']] <- cluster_id
  macroglia_markers_scran[[cluster_id]] <- tmp
}
macroglia_markers_scran <- do.call(rbind, macroglia_markers_scran)

# Same statistical test but no pairwise comparison. See Seurat::FindMarkers() for
# more detail.
macroglia_markers_seurat <- FindAllMarkers(object = macroglia,
                                          assay = 'RNA',
                                          slot = 'data',
                                          only.pos = TRUE,
                                          min.pct = 0.2)


# Heatmap with gene results from scran test
scran_genes <- macroglia_markers_scran %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = -p.value) %>%
  select(gene)
expr_data <- t(ScaleData(object = macroglia@assays$RNA@data[scran_genes$gene,])) %>%
  cbind(macroglia@meta.data[c('macroglia_subcluster')]) %>%
  group_by(macroglia_subcluster) %>%
  summarise(across(.cols = scran_genes[['gene']], .fns = mean)) %>%
  reshape2::melt(id.vars = 'macroglia_subcluster') %>%
  mutate(variable = factor(variable, levels = rev(unique(variable)))) 
max_expr <- 2.5
myColors <- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min(expr_data$value), 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))
de_markers_scran_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = macroglia_subcluster, y = variable)) +
  geom_tile(mapping = aes(fill = value), color = 'black') +
  scale_fill_gradientn(colors = myColors,
                       values = scales::rescale(x = myBreaks, to = c(0,1)),
                       limits = c(NA, max_expr), 
                       na.value = myColors[100],
                       breaks = seq(-5,5,0.5)); de_markers_scran_heatmap


# Heatmap with gene results from Seurat
seurat_genes <- macroglia_markers_seurat %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_logFC) %>%
  # top_n(n = 10, wt = -p.value) %>%
  .[['gene']]
expr_data <- t(ScaleData(object = macroglia@assays$RNA@data[seurat_genes,])) %>%
  cbind(macroglia@meta.data[c('macroglia_subcluster')]) %>%
  group_by(macroglia_subcluster) %>%
  summarise(across(.cols = all_of(seurat_genes), .fns = mean)) %>%
  reshape2::melt(id.vars = 'macroglia_subcluster') %>%
  mutate(variable = factor(variable, levels = rev(unique(variable)))) 
max_expr <- 2.5
myColors <- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min(expr_data$value), 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))
de_markers_seurat_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = macroglia_subcluster, y = variable)) +
  geom_tile(mapping = aes(fill = value), color = 'black', size = 0.3) +
  scale_fill_gradientn(colors = myColors,
                       values = scales::rescale(x = myBreaks, to = c(0,1)),
                       limits = c(NA, max_expr), 
                       na.value = myColors[100],
                       breaks = seq(-5,5,1)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.position = 'right',
        legend.title = element_text(angle = 90, color = 'black', size = 14, vjust = 1, hjust = 0),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(fill = guide_colorbar(title = 'Scaled\nExpression',
                               title.position = 'left',
                               frame.colour = 'black',
                               ticks.colour = 'black',
                               frame.linewidth = 1,
                               ticks.linewidth = 1)); de_markers_seurat_heatmap
ggsave(filename = paste0(results_out, 'macroglia_de_markers_heatmap.tiff'),
       plot = de_markers_seurat_heatmap, device = 'tiff',
       height = 11, width = 5.25)





# Oligo lineage population dynamics -------------------------------------------

counts <- data.frame(table(macroglia$macroglia_subcluster, macroglia$time))
names(counts) <- c('subcluster','time','count')
counts[['time']] <- plyr::mapvalues(
  x = counts[['time']],
  from = c('Uninjured'),
  to = c('Uninj')
)
count_these <- c('OPC-A',
                 'OPC-B',
                 'Div-OPC',
                 'Pre-Oligo',
                 'Oligodendrocyte')

oligo_prop <- counts %>%
  filter(subcluster %in% count_these) %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) + 
  geom_bar(aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black', 
           size = 0.75,
           position = 'fill') + 
  scale_y_continuous() + 
  scale_x_discrete() +
  scale_fill_manual(values = macroglia_cols) +
  ylab('Prop. of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'),
        legend.spacing.y = unit(x = 2.5, units = 'mm'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(ncol = 1))
oligo_prop
ggsave(filename = paste0(results_out, 'oligo_population_prop_dynamics.tiff'),
       plot = oligo_prop, device = 'tiff', height = 4, width = 2.75)


oligo_counts <- counts %>%
  filter(subcluster %in% count_these) %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) +
  geom_bar(mapping = aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black',
           size = 0.75) +
  scale_y_continuous(breaks = seq(0, 5000, 500)) + 
  scale_x_discrete() +
  scale_fill_manual(values = macroglia_cols) +
  ylab(label = '# of cells') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 12, color = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.position = 'none', 
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = NA))

oligo_legend <- cowplot::get_legend(plot = oligo_prop)
oligo_pop <- cowplot::plot_grid(oligo_counts, oligo_prop + theme(legend.position = 'none'), ncol = 2, rel_widths = c(1, 1))
oligo_pop <- cowplot::plot_grid(oligo_pop, oligo_legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave(filename = paste0(results_out, 'oligo_population_dynamics.tiff'),
       plot = oligo_pop, device = 'tiff', height = 2.75, width = 7)








# Ependymal population dynamics -------------------------------------------

counts <- data.frame(table(macroglia$macroglia_subcluster, macroglia$time))
names(counts) <- c('subcluster','time','count')
counts[['time']] <- plyr::mapvalues(
  x = counts[['time']],
  from = c('Uninjured'),
  to = c('Uninj')
)
count_these <- c('Ependymal-A',
                 'Ependymal-B',
                 'Astroependymal',
                 'Astrocyte')

ependymal_prop <- counts %>%
  filter(subcluster %in% count_these) %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) + 
  geom_bar(aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black', 
           size = 0.75,
           position = 'fill') + 
  scale_y_continuous() + 
  scale_x_discrete() +
  scale_fill_manual(values = macroglia_cols) +
  ylab('Prop. of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'),
        legend.spacing.y = unit(x = 2.5, units = 'mm'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(ncol = 1))
ependymal_prop
ggsave(filename = paste0(results_out, 'ependymal_population_prop_dynamics.tiff'),
       plot = ependymal_prop, device = 'tiff', height = 4, width = 2.75)


ependymal_counts <- counts %>%
  filter(subcluster %in% count_these) %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) +
  geom_bar(mapping = aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black',
           size = 0.75) +
  scale_y_continuous(breaks = seq(0, 5000, 500)) + 
  scale_x_discrete() +
  scale_fill_manual(values = macroglia_cols) +
  ylab(label = '# of cells') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 12, color = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.position = 'none', 
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = NA))

ependymal_legend <- cowplot::get_legend(plot = ependymal_prop)
ependymal_pop <- cowplot::plot_grid(ependymal_counts, ependymal_prop + theme(legend.position = 'none'), ncol = 2, rel_widths = c(1, 1))
ependymal_pop <- cowplot::plot_grid(ependymal_pop, ependymal_legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave(filename = paste0(results_out, 'ependymal_population_dynamics.tiff'),
       plot = ependymal_pop, device = 'tiff', height = 2.75, width = 7)




# Joint pop diagram (for paper) ------------------------------------------------

# Requires patchwork package
joint_prop <- ependymal_prop + oligo_prop
ggsave(filename = paste0(results_out, 'joint_population_prop_dynamics.tiff'),
       plot = joint_prop, device = 'tiff', height = 4, width = 5.25)


# Ependymal cilia heatmap  ------------------------------------------------

# Goal here is to better characterize similarties/difference between the newly
# identified clusters (i.e. Astroependymal, OPC-A, and OPC-B). Exclude 
# oligodendrocytes, pre-oligos, and div-OPCs bc their relationship is clearer.
DefaultAssay(macroglia) <- 'RNA'
Idents(macroglia) <- 'macroglia_subcluster'

# # Setup groups (merge ependymal-A and B bc they appear to be along continuous 
# # spectrum with no strongly distinguishing markers)
# groups <- list(c('Ependymal-A','Ependymal-B'),
#                'Astroependymal',
#                'Astrocyte',
#                'OPC-A',
#                'OPC-B')
# # Matrix of pairs of groups
# comparisons <- t(x = combn(x = 1:length(groups), m = 2))
# # Set names of DE results list to compared groups. Merged subclusteres (e.g.
# # Ependymal-A and Ependymal-B) separated by underscore, groups separated by 
# # period.
# comparisons_names <- paste(
#   sapply(X = groups[comparisons[,1]],
#          FUN = function(x) paste(x, collapse = '_')), 
#   sapply(X = groups[comparisons[,2]],
#          FUN = function(x) paste(x, collapse = '_')),
#   sep = '.'
# )
# # Setup for parallel processing. Use all available cores minus 1.
# cores <- detectCores()
# cl <- makeCluster(spec = cores[1] - 1)
# doParallel::registerDoParallel(cl = cl)
# # Subset seurat data for lighter memory load
# macroglia_lite <- DietSeurat(
#   object = macroglia,
#   assays = 'RNA',
#   dimreducs = NULL,
#   graphs = NULL
# )
# # Use foreach package to iterate in parallel. `.final` sets the names of the 
# # resulting list to the previously defined names.
# de_results <- foreach(
#   i = 1:nrow(comparisons),
#   .final = function(x) setNames(x, comparisons_names)
# ) %dopar% {
#   tmp <- Seurat::FindMarkers(
#     object = macroglia_lite,
#     ident.1 = groups[[comparisons[i, 1]]],
#     ident.2 = groups[[comparisons[i, 2]]],
#     test.use = 'wilcox',
#     assay = 'RNA',
#     slot = 'data'
#   )
# }
# stopCluster(cl = cl); gc()
# 
# # Extract top DE genes per comparison. Allow downstream clustering to reveal the
# # relationships.
# top_de <- unique(
#   unlist(
#     lapply(
#       X = de_results,
#       FUN = function(x) {
#         x %>%
#           filter(p_val_adj < 1e-10) %>% # p-value threshold
#           arrange(p_val_adj, abs(avg_logFC)) %>%
#           rownames(.) %>%
#           .[1:5] # top by p-value, then magnitude of fold-change
#         }
#       ),
#     use.names = FALSE
#   )
# )
# 
# # shared genes between astrocytes and astroependymal cells compared to ependymal
# astro_genes <- FindMarkers(
#   object = macroglia,
#   ident.1 = c('Astroependymal','Astrocyte'),
#   ident.2 = c('Ependymal-A','Ependymal-B'),
#   test.use = 'wilcox',
#   assay = 'RNA',
#   slot = 'data'
# )
# shared_astro_genes <- intersect(
#   x = rownames(de_results$`Ependymal-A_Ependymal-B.Astroependymal`)[de_results$`Ependymal-A_Ependymal-B.Astroependymal`$avg_logFC < 0 & de_results$`Ependymal-A_Ependymal-B.Astroependymal`$p_val_adj < 1e-10],
#   y = rownames(astro_genes)[astro_genes$avg_logFC > 0 & astro_genes$p_val_adj < 1e-10]
# )
# top_astro_genes <- astro_genes[shared_astro_genes,] %>%
#   arrange(desc(avg_logFC)) %>%
#   rownames(.) %>%
#   .[1:10]
# 
# # cilia filament associated proteins
# cfap_genes <- grep(
#   pattern = 'Cfap',
#   x = rownames(macroglia[['RNA']]@counts),
#   value = TRUE
# )



# Find shared genes between astrocytes/astroependymal, and OPC-B/astroependymal
# de_genes <- FindAllMarkers(
#   object = macroglia,
#   assay = 'RNA',
#   slot = 'data',
#   only.pos = TRUE
# )
# write.table(x = de_genes, file = paste0(results_out, 'macroglia_de_genes.tsv'),
#             sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
de_genes <- read.table(file = paste0(results_out, 'macroglia_de_genes.tsv'), header = TRUE, sep = '\t', row.names = 1)
astro_astroepen_shared <- intersect(
  x = de_genes$gene[de_genes$cluster == 'Astrocyte' & de_genes$p_val_adj < 1e-10],
  y = de_genes$gene[de_genes$cluster == 'Astroependymal' & de_genes$p_val_adj < 1e-10]
)
astro_astroepen_shared <- de_genes[order(de_genes$avg_logFC, decreasing = TRUE),][astro_astroepen_shared, 'gene'][1:15]
opcB_astroepen_shared <- intersect(
  x = de_genes$gene[de_genes$cluster == 'OPC-B' & de_genes$p_val_adj < 1e-10],
  y = de_genes$gene[de_genes$cluster == 'Astroependymal' & de_genes$p_val_adj < 1e-10]
)
opcB_astroepen_shared <- de_genes[order(de_genes$avg_logFC, decreasing = TRUE),][opcB_astroepen_shared, 'gene'][1:15]

# Find genes distinguishing OPC-A and OPC-B
# opcA_opcB_de <- FindMarkers(
#   object = macroglia,
#   ident.1 = 'OPC-A',
#   ident.2 = 'OPC-B',
#   assay = "RNA",
#   slot = 'data',
#   logfc.threshold = 0,
#   only.pos = FALSE
# )
# opcA_opcB_de$which_high <- ifelse(
#   test = opcA_opcB_de$avg_logFC > 0,
#   yes = 'OPC-A',
#   no = 'OPC-B'
# )
# write.table(x = opcA_opcB_de, file = paste0(results_out, 'opcA_opcB_de_genes.tsv'),
#             sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
opcA_opcB_de <- read.table(file = paste0(results_out, 'opcA_opcB_de_genes.tsv'),
                           header = TRUE, sep = '\t', row.names = 1)
opcA_opcB_genes <- opcA_opcB_de %>%
  mutate('gene' = rownames(.)) %>%
  mutate('cell' = ifelse(avg_logFC > 0, 'OPC-A','OPC-B')) %>%
  group_by(cell) %>%
  top_n(n = 10, wt = abs(avg_logFC)) %>%
  .[['gene']]

# Find genes distinguishing Astroependymal/ependymal
# astroependymal_ependymal_de <- FindMarkers(
#   object = macroglia, 
#   ident.1 = 'Astroependymal',
#   ident.2 = c('Ependymal-A','Ependymal-B'),
#   assay = "RNA",
#   slot = 'data', 
#   logfc.threshold = 0.5,
#   only.pos = FALSE
# )
# astroependymal_ependymal_de$which_high <- ifelse(
#   test = astroependymal_ependymal_de$avg_logFC > 0,
#   yes = 'Astroependymal',
#   no = 'Ependymal'
# )
# write.table(x = astroependymal_ependymal_de, file = paste0(results_out, 'astroependymal_ependymal_de_genes.tsv'),
#             sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
astroependymal_ependymal_de <- read.table(
  file = paste0(results_out, 'astroependymal_ependymal_de_genes.tsv'),
  header = TRUE, sep = '\t', row.names = 1)
astroependymal_ependymal_genes <- astroependymal_ependymal_de %>%
  mutate('gene' = rownames(.)) %>%
  top_n(n = 10, wt = abs(avg_logFC)) %>%
  .[['gene']]

# Globally unique genes per cluster
unique_genes <- de_genes %>%
  group_by(cluster) %>%
  filter(p_val_adj < 1e-03) %>%
  top_n(n = 7, wt = -p_val_adj) %>%
  top_n(n = 7, wt = avg_logFC) %>%
  ungroup() %>%
  .[['gene']] %>%
  unique()

# Cilia-associated genes (see: )
cilia_genes <- c(paste0('Cfap', c(20,36,43,44,45,52,53,54,69,70,77,100,126,161,206)),'Ccdc180','Odf3b','Fam183b','Tekt1','Dnah12','Sntn','Dnah6','Tmem107','Acta2','Stmn2','Tmsb10')

# Combine all genes to plot
all_genes <- unique(c(astro_astroepen_shared, opcB_astroepen_shared, opcA_opcB_genes, astroependymal_ependymal_genes, unique_genes, cilia_genes))
write.table(x = all_genes, file = paste0(results_out, 'macroglia_cilia_heatmap_genes.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# Calculate scaled expression values
expr_data <- t(ScaleData(object = macroglia[['RNA']]@data[all_genes,], 
                         scale.max = 4)) %>%
  cbind(macroglia@meta.data[c('macroglia_subcluster')]) %>%
  group_by(macroglia_subcluster) %>%
  summarise(across(.cols = all_of(all_genes), .fns = mean)) %>%
  tibble::column_to_rownames(var = 'macroglia_subcluster')
expr_data <- t(as.matrix(x = expr_data))

# Generate heatmap
cilia_heatmap <- Heatmap(
  matrix = expr_data,
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  col = rev(
    colorRampPalette(
      colors = RColorBrewer::brewer.pal(
        n = 11,
        name = 'RdYlBu'
      ),
      bias = 0.65
    )(100)
  ),
  border = TRUE,
  clustering_distance_rows = 'pearson',
  clustering_distance_columns = 'pearson',
  column_dend_height = unit(x = 15, units = 'mm'),
  row_dend_width = unit(x = 15, units = 'mm'),
  column_names_rot = 45,
  rect_gp = gpar(col = 'black'),
  heatmap_legend_param = list(
    title = 'z-score',
    title_gp = gpar(fontsize = 12),
    title_position = 'leftcenter-rot',
    labels = c('Low', 0, 'High'),
    at = c(min(expr_data), 0, max(expr_data)),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(2.5, units = 'cm'),
    grid_width = unit(0.5, units = 'cm'),
    border = 'black',
    title_gap = unit(1, units = 'cm'),
    direction = 'vertical'),
  # column_km = 4,
)
draw(cilia_heatmap, heatmap_legend_side = 'right')
tiff(file = paste0(results_out, 'macroglia_cilia_heatmap.tiff'),
     height = 21, width = 5.5, res = 420, units = 'in')
draw(cilia_heatmap, heatmap_legend_side = 'right')
dev.off()




# Ependymal cilia heatmap with old list of genes --------------------------------

cilia.genes <- c(paste0('Cfap', c(20,36,43,44,45,52,53,54,69,70,77,100,126,161,206)),'Ccdc180','Odf3b','Fam183b','Tekt1','Dnah12','Sntn','Dnah6','Tmem107','Acta2','Stmn2','Tmsb10')
tmp1 <- FindMarkers(macroglia, 'Astrocyte', ident.2 = paste0('Ependymal-', LETTERS[1:2]), assay = 'RNA', slot = 'data', only.pos = TRUE)
tmp2 <- FindMarkers(macroglia, 'Astroependymal', ident.2 = paste0('Ependymal-', LETTERS[1:2]), assay = 'RNA', slot = 'data', only.pos = TRUE)
astrocyte.genes <- c('Gfap','Agt','Ntsr2','Acsbg1','Aldh1l1','Aqp4','Slc6a11','Slc7a10','Fgfr3')
astroependymal.genes <- c('Ccnd1','S100a10','Hbegf','Stmn2','Timp1','Vim','Hspb1','Crym')
oligo.genes <- c('Sox10','Olig1','Olig2','Ascl1','Bmp4','Mycl','Plp1','Mog','Mbp','Sox9','Mki67','Top2a')
tmp.thres <- 0.4
shared.genes <- intersect(rownames(tmp1)[tmp1$p_val_adj < 1e-10 & tmp1$avg_logFC > tmp.thres], rownames(tmp2)[tmp1$p_val_adj < 1e-10 & tmp1$avg_logFC > tmp.thres])

marker.genes <- c(cilia.genes, astrocyte.genes, astroependymal.genes, oligo.genes, shared.genes)
marker.genes <- unique(marker.genes)
all_genes <- marker.genes

expr_data <- t(ScaleData(object = macroglia[['RNA']]@data[all_genes,], 
                         scale.max = 4)) %>%
  cbind(macroglia@meta.data[c('macroglia_subcluster')]) %>%
  group_by(macroglia_subcluster) %>%
  summarise(across(.cols = all_of(all_genes), .fns = mean)) %>%
  tibble::column_to_rownames(var = 'macroglia_subcluster')
expr_data <- as.matrix(x = expr_data)

cilia_heatmap <- Heatmap(
  matrix = expr_data,
  clustering_method_columns = 'ward.D2',
  clustering_method_rows = 'ward.D2',
  col = rev(
    colorRampPalette(
      colors = RColorBrewer::brewer.pal(
        n = 9,
        name = 'RdYlBu'
      ),
      bias = 0.65
    )(100)
  ),
  border = TRUE,
  clustering_distance_rows = 'pearson',
  clustering_distance_columns = 'pearson',
  column_dend_height = unit(x = 15, units = 'mm'),
  row_dend_width = unit(x = 15, units = 'mm'),
  # column_names_rot = 45,
  rect_gp = gpar(col = 'black'),
  heatmap_legend_param = list(
    title = 'Scaled Expression',
    title_gp = gpar(fontsize = 12),
    title_position = 'topcenter',
    labels = c('Low', 0, 'High'),
    at = c(min(expr_data), 0, max(expr_data)),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(2.5, units = 'cm'),
    grid_width = unit(0.5, units = 'cm'),
    border = 'black',
    title_gap = unit(1, units = 'cm'),
    direction = 'horizontal'),
  # column_km = 4,
)
draw(cilia_heatmap, heatmap_legend_side = 'bottom')
tiff(file = paste0(results_out, 'astroependymal_cilia_heatmap_oldMethod.tiff'),
    height = 5, width = 20, res = 420, units = 'in')
draw(cilia_heatmap, heatmap_legend_side = 'bottom')
dev.off()
png(filename = paste0(results_out, 'astroependymal_cilia_heatmap_oldMethod.png'),
    height = 5, width = 20, res = 420, units = 'in')
draw(cilia_heatmap, heatmap_legend_side = 'bottom')
dev.off()



# Abundance analysis ------------------------------------------------------

time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

macroglia_subcluster_counts <- data.frame(table(macroglia$macroglia_subcluster, macroglia$sample_id))
colnames(macroglia_subcluster_counts) <- c('macroglia_subcluster','sample_id','count')
macroglia_subcluster_counts[['time']] <- substr(x = macroglia_subcluster_counts[['sample_id']],
                                              start = 1,
                                              stop = regexpr(pattern = '_', text = macroglia_subcluster_counts[['sample_id']]) -1)
macroglia_subcluster_counts[['time']] <- plyr::mapvalues(x = macroglia_subcluster_counts[['time']], from ='uninj', to = 'Uninjured')
macroglia_subcluster_counts[['time']] <- factor(macroglia_subcluster_counts[['time']], levels = levels(macroglia$time))
macroglia_subcluster_counts_plot <- macroglia_subcluster_counts %>%
  ggplot(mapping = aes(x = sample_id, y = count)) +
  geom_bar(mapping = aes(fill = time), color = 'black', stat = 'identity') +
  labs(title = 'Cell abundance across time and sample') +
  ylab(label = 'Number of cells') +
  facet_wrap(. ~ macroglia_subcluster, scales = 'free_y', ncol = 5) +
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
ggsave(filename = paste0(results_out, 'macroglia_subcluster_counts_bySample.tiff'),
       plot = macroglia_subcluster_counts_plot, device = 'tiff', height = 6.5, width = 15)








# OPC-A vs OPC-B -------------------------------------------------------

# Differences between OPC-A and OPC-B
opcA_opcB <- FindMarkers(macroglia,
                         ident.1 = 'OPC-A',
                         ident.2 = 'OPC-B',
                         assay = 'RNA',
                         logfc.threshold = 0.25)
write.table(opcA_opcB, file = paste0(results_out, 'OPC-A_vs_OPC-B_DEgenes.tsv'),
            sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
DefaultAssay(macroglia) <- 'RNA'

# Select DE genes
opcA_opcB_genes <- list('OPC' = c('Neu4','Ptprz1','Cspg5','Bcan','Cntn1'),
                      'U-OPC' = c('Vgf','Lgals1','Ccl2','Tnc','Klk8'))
opcA_opcB_genes <- unlist(opcA_opcB_genes, use.names = FALSE)

umap_coordinates <- FetchData(macroglia, vars = c('UMAP_1','UMAP_2'), slot = 'data')
label_x <- min(umap_coordinates[['UMAP_1']]) * 0.975
label_y <- min(umap_coordinates[['UMAP_2']]) * 0.85
legend_x <- 0.15
legend_y <- 0.85
expr_colors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = 'Spectral'))(100))

opcA_opcB_genes_umap <- vector(mode = 'list', length = length(opcA_opcB_genes))
names(opcA_opcB_genes_umap) <- opcA_opcB_genes

for (ii in 1:length(unlist(opcA_opcB_genes, use.names = FALSE))) {
  gene <- opcA_opcB_genes[ii]
  # gather data
  expr <- FetchData(macroglia, vars = gene, slot = 'data')
  expr <- cbind(expr, umap_coordinates)
  expr <- expr[order(expr[[gene]], decreasing = FALSE),]
  max_expr <- ceiling(max(expr[[gene]]) * 10) / 10
  
  if (grepl('-', x = gene)) {
    names(opcA_opcB_genes_umap)[ii] <- gsub(pattern = '-', replacement = '.', x = names(opcA_opcB_genes_umap)[ii])
    colnames(expr) <- gsub(pattern = '-', replacement = '.', x = colnames(expr))
    gene <- gsub(pattern = '-', replacement = '.', x = gene)
  }
  
  # assemble gg
  opcA_opcB_genes_umap[[gene]] <- expr %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes_string(color = gene), size = 0.5, alpha = 0.4) +
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
          legend.text = element_text(size = 12, color = 'black')) +
    guides(color = guide_colorbar(frame.colour = 'black', 
                                  ticks = FALSE,
                                  barwidth = 0.8,
                                  frame.linewidth = 1))
}
# save
opcA_opcB_genes_umap_grid <- cowplot::plot_grid(plotlist = opcA_opcB_genes_umap, nrow = 2, byrow = TRUE)
ggsave(filename = paste0(results_out, 'OPC-A_vs_OPC-B_DEgenes_umap.tiff'),
       plot = opcA_opcB_genes_umap_grid, device = 'tiff',
       height = 5, width = 12)


# Violin plots of transcription factors to confirm OPC identity
transcription_factors <- c('Ascl1','Id2','Id4','Olig1','Olig2','Nkx2-2','Sox9','Sox10','Myt1','Foxj1','Nog')
DefaultAssay(macroglia) <- 'RNA'
transcription_factor_vlnplot <- 
  FetchData(object = macroglia,
            vars = c('macroglia_subcluster', transcription_factors), 
            slot = 'data') %>%
  reshape2::melt(id.vars = c('macroglia_subcluster')) %>%
  ggplot(mapping = aes(x = macroglia_subcluster, y = value)) +
  geom_violin(mapping = aes(fill = macroglia_subcluster), scale = 'width') +
  facet_wrap(. ~ variable, scales = 'free_y', strip.position = 'left', ncol = 1) +
  scale_fill_manual(values = macroglia_cols) +
  scale_y_continuous(breaks = seq(0, 10, 1), position = 'right') +
  ylab(label = 'Normalized expression') +
  theme(axis.text.x = element_text(angle = 45, size = 14, hjust = 1),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = 'black', fill = NA),
        axis.title.x = element_blank(),
        axis.title.y.right = element_text(size = 14, color = 'black', margin = margin(0, 0, 0, 3, unit = 'mm')),
        axis.text.y = element_text(size = 12, color = 'black'),
        strip.text.y.left = element_text(size = 14, color = 'black', angle = 0),
        legend.position = 'none',
        plot.margin = margin(0, 1, 0, 1, unit = 'cm'))
transcription_factor_vlnplot
ggsave(filename = paste0(results_out, 'transcription_factors_vln.tiff'), 
       plot = transcription_factor_vlnplot,
       height = 9, width = 5, device = 'tiff')






# # Measuring relationships between macroglia subclusters ----------------------
# 
# # Preset values
# Idents(macroglia) <- 'macroglia_subcluster'
# DefaultAssay(macroglia) <- 'RNA'
# macroglia_cols <- c('Neutrophil' = '#800000',    
#                     'Monocyte' = '#9a6324',
#                     'Macrophage-A' = '#e6194b',
#                     'Macrophage-B' = '#f58231',
#                     'BA-Macrophage' = '#CCCC00',
#                     'Dendritic' = '#808000',
#                     'Div-Myeloid' = '#3cb44b',
#                     'Microglia-A' = '#008080',
#                     'Microglia-B' = 'cyan3',
#                     'Microglia-C' = '#000075',
#                     'Div-Microglia' = '#4363d8',
#                     'IFN-Myeloid' = '#911eb4')
# time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
# names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')
# 
# 
# # Use Seurat wilcox test
# microglia_srat_genes <- FindAllMarkers(object = microglia,
#                                        assay = 'RNA',
#                                        logfc.threshold = 0.25,
#                                        only.pos = TRUE)
# microglia_de_genes <- microglia_srat_genes %>%
#   group_by(cluster) %>%
#   filter(avg_logFC > 0.5) %>%
#   .[['gene']] %>%
#   unique()
# 
# # Select which genes to annotate
# label_genes <- c('P2ry12','Siglech','Tmem119','Il1b','Msr1','Igf1','Apoe','Lgals3','Irf7','Isg15','Cdk1','Mki67','Lpl','Csf1r','Plin2','Mif','Cxcl2','Fabp5','Lyz2','Ccl6','Ccl9','Spp1','Cd68','Cd5l','Ldha')
# 
# 
# # Extract scaled expression
# expr_data <- ScaleData(object = microglia[['RNAcorrected']]@data[microglia_de_genes,], 
#                        scale.max = 4)
# microglia_meta <- microglia@meta.data[c('macroglia_subcluster','time','sample_id')]
# 
# # Cell-level annotations
# microglia_cell_anno <- HeatmapAnnotation(
#   'Subcluster' = microglia_meta[['macroglia_subcluster']],
#   'Time' = microglia_meta[['time']],
#   col = list('Subcluster' = macroglia_cols[names(macroglia_cols) %in% microglia_names],
#              'Time' = time_cols)
# )
# 
# # Gene-level annotations
# microglia_gene_anno <- rowAnnotation(
#   'Genes' = anno_mark(
#     at = match(label_genes, rownames(expr_data)),
#     labels = label_genes,
#     side = 'right',
#     labels_gp = gpar(fontsize = 10)
#   )
# )
# 
# # Build heatmap
# tmp <- Heatmap(
#   matrix = expr_data, 
#   col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")), bias = 0.6)(100),
#   column_title = 'Differentially expressed genes among microglia subclusters',
#   column_title_gp = gpar(fontsize = 20),
#   clustering_method_columns = 'ward.D2',
#   clustering_method_rows = 'ward.D2',
#   column_dend_height = unit(20, units = 'mm'),
#   row_dend_width = unit(20, units = 'mm'),
#   show_row_names = FALSE,
#   show_column_names = FALSE,
#   heatmap_legend_param = list(title = 'Scaled Expression',
#                               title_gp = gpar(fontsize = 12),
#                               title_position = 'leftcenter-rot',
#                               labels = c('Low','High'),
#                               at = c(min(expr_data), max(expr_data)),
#                               labels_gp = gpar(fontsize = 10),
#                               legend_height = unit(2.5, units = 'cm'),
#                               grid_width = unit(0.5, units = 'cm'),
#                               border = 'black',
#                               title_gap = unit(1, units = 'cm'),
#                               direction = 'vertical'),
#   use_raster = TRUE,
#   top_annotation = microglia_cell_anno,
#   right_annotation = microglia_gene_anno
# )
# 
# # Save heatmap
# pdf(file = paste0(results_out, 'microglia_heatmap_complex_small.pdf'),
#     height = 5, width = 9)
# tmp
# dev.off()
# svg(file = paste0(results_out, 'microglia_heatmap_complex_small.svg'),
#     height = 5, width = 9)
# tmp
# dev.off()





# # Heatmap characterization of ependymal cells w/ cilia --------------------
# 
# library('foreach')
# library('ComplexHeatmap')
# 
# # First calculate some DE genes per subcluster
# macroglia_lite <- DietSeurat(
#   object = macroglia,
#   assays = 'RNA', 
#   dimreducs = NULL,
#   graphs = NULL
# )
# 
# # For parallel computing, use all but 1 core.
# cores <- detectCores()
# cl <- makeCluster(spec = cores[1] - 1)
# doParallel::registerDoParallel(cl = cl)
# 
# # Wilcox test
# subs <- levels(macroglia@meta.data[['macroglia_subcluster']])
# de_wilcox <- vector(mode = 'list', length = length(subs))
# de_wilcox <- foreach::foreach(i = 1:length(subs)) %dopar% {
#   tmp <- Seurat::FindMarkers(
#     object = macroglia_lite,
#     ident.1 = subs[i],
#     assay = 'RNA',
#     slot = 'data',
#     test.use = 'wilcox'
#   )
# }
# names(de_wilcox) <- subs
# stopCluster(cl = cl)
# 
# 
# # Heatmap
# macroglia_cols <- c('Ependymal-A' = '#800000',
#                     'Ependymal-B' = '#e6194b',
#                     'Astroependymal' = '#f58231',
#                     'Astrocyte' = 'goldenrod',
#                     'OPC-A' = '#3cb44b',
#                     'OPC-B' = '#008080',
#                     'Div-OPC' = '#4363d8',
#                     'Pre-Oligo' = '#911eb4',
#                     'Oligodendrocyte' = '#f032e6')
# time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
# names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')
# top_genes <- sapply(
#   X = de_wilcox,
#   FUN = function(x) {
#     tmp <- x %>%
#       filter(p_val_adj < 1e-10) %>%
#       filter(avg_logFC > 0) %>%
#       arrange(desc(avg_logFC))
#     return(rownames(tmp)[1:20])
#   }
# )
# label_genes <- apply(X = top_genes, FUN = function(x) x[1:3], MARGIN = 2)
# top_genes <- unique(c(top_genes))
# label_genes <- unique(c(label_genes))
# expr_data <- ScaleData(
#   object = macroglia_lite[['RNA']]@data,
#   features = top_genes
# )
# row_gene_anno <- rowAnnotation(
#   'Genes' = anno_mark(
#     at = match(label_genes, rownames(expr_data)),
#     labels = label_genes,
#     side = 'right',
#     labels_gp = gpar(fontsize = 10)
#   )
# )
# col_cell_anno <- HeatmapAnnotation(
#   'Subcluster' = macroglia@meta.data[['macroglia_subcluster']],
#   'Time' = macroglia@meta.data[['time']],
#   col = list(
#     'Subcluster' = macroglia_cols,
#     'Time' = time_cols
#   )
# )
# 
# pdf(file = paste0(results_out, 'Cilia_heatmap.pdf'),
#     height = 6, width = 8)
# Heatmap(
#   matrix = expr_data,
#   show_row_names = FALSE,
#   show_column_names = FALSE,
#   use_raster = TRUE,
#   clustering_method_rows = 'ward.D2',
#   clustering_method_columns = 'ward.D2',
#   top_annotation = col_cell_anno,
#   right_annotation = row_gene_anno
# )
# dev.off()