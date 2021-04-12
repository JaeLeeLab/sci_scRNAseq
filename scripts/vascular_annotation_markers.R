
######## Vascular cell-type annotation via marker identification ########


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
results_out <- './results/vascular_annotation_markers/'
ref_in <- './ref/'
ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)

vascular <- readRDS(file = './data/vascular.rds')




# Additional markers UMAP ----------------------------------------------------


DefaultAssay(vascular) <- 'RNA'

# These are marker genes identified in previous studies. Multiple sources.
addtl_markers <- list('Endothelial' = c('Cldn5','Podxl','Pecam1'),
                      'Pericyte' = c('Pdgfrb','Cspg4','Kcnj8'),
                      'VSMC' = c('Acta2','Tagln','Myh11'),
                      'Fibroblast' = c('Col1a1','Postn','Dcn'),
                      # 'VLMC' = c('Lum','Pdgfra','Il33'),
                      'Tip cell' = c('Apln','Pgf','Angpt2'),
                      'Zonation' = c('Sema3g','Vcam1','Slc38a5'))
addtl_markers <- unlist(addtl_markers, use.names = FALSE)

umap_coordinates <- FetchData(vascular, vars = c('UMAP_1','UMAP_2'), slot = 'data')
label_x <- min(umap_coordinates[['UMAP_1']]) * 0.975
label_y <- max(umap_coordinates[['UMAP_2']]) * 0.85
legend_x <- 0.85
legend_y <- 0.2
expr_colors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = 'Spectral'))(100))

addtl_markers_umap <- vector(mode = 'list', length = length(addtl_markers))
names(addtl_markers_umap) <- addtl_markers

for (ii in 1:length(unlist(addtl_markers, use.names = FALSE))) {
  gene <- addtl_markers[ii]
  # gather data
  expr <- FetchData(vascular, vars = gene, slot = 'data')
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
    geom_point(mapping = aes_string(color = gene), size = 0.5, alpha = 0.5) +
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
addtl_markers_umap_grid <- cowplot::plot_grid(plotlist = addtl_markers_umap, nrow = 6)
ggsave(filename = paste0(results_out, 'addtl_markers_umap.tiff'),
       plot = addtl_markers_umap_grid, device = 'tiff',
       height = 17, width = 9)




# DE marker genes calculation -----------------------------------------------


# Use wilcox test while blocking against Time (group) to identify marker genes
# per cluster. scran::findMarkers() performs pair-wise for all cluster 
# combinations. 
rna_sce <- as.SingleCellExperiment(vascular, assay = 'RNA')
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

Idents(vascular) <- 'default_cluster'
seurat_markers <- FindAllMarkers(object = vascular,
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

vascular@meta.data[['vascular_subcluster']] <- plyr::mapvalues(
  x = vascular@meta.data[['vascular_subcluster']],
  from = levels(vascular@meta.data[['vascular_subcluster']]),
  to = c('C1-Endothelial',
         'C2-Endothelial',
         'Tip Cell',
         'V-Endothelial',
         'U-Vascular',
         'A-Endothelial',
         'Fibroblast',
         'Pericyte',
         'VSMC')
)
vascular@meta.data[['vascular_subcluster']] <- factor(
  x = vascular@meta.data[['vascular_subcluster']],
  levels = c('A-Endothelial',
             'C1-Endothelial',
             'C2-Endothelial',
             'V-Endothelial',
             'Tip Cell',
             'Pericyte',
             'VSMC',
             'Fibroblast',
             'U-Vascular')
)
Idents(vascular) <- 'vascular_subcluster'

# Save barcode-identities data.frame
write.table(x = vascular@meta.data['vascular_subcluster'],
            file = paste0(results_out, 'vascular_subcluster.tsv'),
            sep = '\t', row.names = TRUE, col.names = NA)





# Cell-type UMAP across time ----------------------------------------------

# join fibroblasts
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


Idents(object = vascular) <- 'vascular_subcluster'

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
vascular_subcluster_counts <- table(vascular$vascular_subcluster)
vascular_subcluster_label <- paste0(names(vascular_subcluster_counts), ' (', vascular_subcluster_counts, ')')
names(vascular_subcluster_label) <- names(vascular_subcluster_counts)
vascular_subcluster_umap <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','vascular_subcluster')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = vascular_subcluster), size = 1, alpha = 0.5) +
  scale_color_manual(values = vascular_cols, 
                     breaks = names(vascular_subcluster_label),
                     label = vascular_subcluster_label) +
  xlab(label = 'UMAP 1') +
  ylab(label = 'UMAP 2') +
  umap_theme +
  theme(legend.text = element_text(size = 16, color = 'black'),
        legend.box.margin = margin(0,0,0,-5,unit='mm')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
vascular_subcluster_umap
ggsave(filename = paste0(results_out, 'vascular_subcluster_annotation_umap.tiff'),
       plot = vascular_subcluster_umap, device = 'tiff', height = 6, width = 9)

# cell-type annotation split by time UMAP (figure 4b)
vascular_subcluster_split_umap <- FetchData(vascular, vars = c('UMAP_1','UMAP_2','vascular_subcluster','time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = vascular_subcluster), size = 1, alpha = 0.5) +
  facet_wrap(. ~ time, ncol = 2) +
  scale_color_manual(values = vascular_cols, breaks = names(vascular_subcluster_label), label = vascular_subcluster_label) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
ggsave(filename = paste0(results_out, 'vascular_subcluster_annotation_split_umap.tiff'),
       plot = vascular_subcluster_split_umap, device = 'tiff', height = 6, width = 9)


# DE markers dot plot -----------------------------------------------------

de_markers <- c('A-Endothelial1' = 'Gkn3',
                'A-Endothelial2' = 'Stmn2',
                'Endothelial1' = 'Cldn5',
                'Endothelial2' = 'Ly6a',
                'V-Endothelial1' = 'Slc38a5',
                'V-Endothelial2' = 'Icam1',
                'Tip Cell' = 'Apln',
                'Pericyte' = 'Kcnj8',
                'VSMC' = 'Tagln',
                'Fibroblast'= 'Col1a1')


DefaultAssay(vascular) <- 'RNA'
avg_exp <- data.frame(t(ScaleData(vascular[['RNA']]@data, features = de_markers)))
avg_exp <- cbind(avg_exp, 'vascular_subcluster' = as.character(vascular@meta.data[,c('vascular_subcluster')])) %>%
  reshape2::melt(id.vars = c('vascular_subcluster')) %>%
  group_by(vascular_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(t(vascular[['RNA']]@counts[unlist(de_markers, use.names = FALSE),]))
pct_exp <- cbind(pct_exp, 'vascular_subcluster' = as.character(vascular@meta.data[,c('vascular_subcluster')])) %>%
  reshape2::melt(id.vars = c('vascular_subcluster')) %>%
  group_by(vascular_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
my_cols <- rev(colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))


min_expr <- floor(min(avg_exp$avg.exp)*10)/10
# max_exp <- ceiling(max(avg_exp$avg.exp)*10)/10
max_expr <- 3
de_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  filter(vascular_subcluster != 'U-Vascular') %>%
  mutate(vascular_subcluster = factor(vascular_subcluster, levels = rev(levels(vascular$vascular_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = vascular_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  # facet_grid(vascular_subcluster ~ ., drop = TRUE, switch = 'y') +
  scale_size(range = c(0,12), limits = c(0,100)) +
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
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'bold', size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(0,0,0,0),
        legend.box = 'vertical',
        legend.position = 'right',
        legend.spacing.x = unit(x = 2, units = 'mm'),
        panel.border = element_rect(fill = NA, size = 1),
        panel.background = element_rect(fill = NA)) +
  guides(fill = guide_colorbar(title = 'Scaled\nexpression', 
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
de_markers_dotplot
ggsave(filename = paste0(results_out, 'vascular_subcluster_markers_dotplot.tiff'),
       plot = de_markers_dotplot, device = 'tiff', height = 3.75, width = 7.25)





# DE markers dot plot (with significant p-values) ------------------------------


de_markers <- c('A-Endothelial1' = 'Gkn3',
                'A-Endothelial2' = 'Stmn2',
                'Endothelial1' = 'Cldn5',
                'Endothelial2' = 'Ly6a',
                'V-Endothelial1' = 'Slc38a5',
                'V-Endothelial2' = 'Icam1',
                'Tip Cell' = 'Apln',
                'Pericyte' = 'Kcnj8',
                'VSMC' = 'Tagln',
                'Fibroblast'= 'Col1a1')
tmp_de <- FindAllMarkers(
  object = vascular,
  assay = 'RNA',
  slot = 'data',
  features = de_markers,
  only.pos = TRUE
)
tmp_de <- tmp_de[tmp_de$cluster != 'U-Vascular',]
tmp_de$p_val_adj <- signif(tmp_de$p_val_adj, digits = 2)
tmp_de$variable <- tmp_de$gene
tmp_de$signif <- ifelse(test = tmp_de$p_val < 1e-10,
                        yes = '*',
                        no = '')

DefaultAssay(vascular) <- 'RNA'
avg_exp <- data.frame(t(ScaleData(vascular[['RNA']]@data, features = de_markers)))
avg_exp <- cbind(avg_exp, 'vascular_subcluster' = as.character(vascular@meta.data[,c('vascular_subcluster')])) %>%
  reshape2::melt(id.vars = c('vascular_subcluster')) %>%
  group_by(vascular_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(t(vascular[['RNA']]@counts[unlist(de_markers, use.names = FALSE),]))
pct_exp <- cbind(pct_exp, 'vascular_subcluster' = as.character(vascular@meta.data[,c('vascular_subcluster')])) %>%
  reshape2::melt(id.vars = c('vascular_subcluster')) %>%
  group_by(vascular_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
my_cols <- rev(colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))


# min_expr <- floor(min(avg_exp$avg.exp)*10)/10
# max_exp <- ceiling(max(avg_exp$avg.exp)*10)/10
min_expr <- -3
max_expr <- 3
de_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  filter(vascular_subcluster != 'U-Vascular') %>%
  mutate(vascular_subcluster = factor(vascular_subcluster, levels = rev(levels(vascular$vascular_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = vascular_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  geom_text(data = tmp_de, mapping = aes(x = variable, y = cluster, label = signif),
            size = 8, fontface = 'bold', nudge_y = -0.09, nudge_x = 0.01) +
  # facet_grid(vascular_subcluster ~ ., drop = TRUE, switch = 'y') +
  scale_size(range = c(0,12), limits = c(0,100)) +
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
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'bold', size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(0,0,0,0),
        legend.box = 'vertical',
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
de_markers_dotplot
ggsave(filename = './results/revision_figures/vascular_subcluster_markers_dotplot.tiff',
       plot = de_markers_dotplot, device = 'tiff', height = 3.75, width = 7.25)



# DE markers heatmap ------------------------------------------------------

Idents(vascular) <- 'vascular_subcluster'
DefaultAssay(vascular) <- 'RNA'

# Use wilcox test while blocking against Time (group) to identify marker genes
# per cluster. scran::findMarkers() performs pair-wise for all cluster
# combinations.
rna_sce <- as.SingleCellExperiment(vascular, assay = 'RNA')
vascular_markers_scran <- findMarkers(x = rna_sce,
                                      groups = colData(rna_sce)[['vascular_subcluster']],
                                      test.type = 'wilcox',
                                      direction = 'up',
                                      pval.type = 'all')
for(ii in 1:length(vascular_markers_scran)) {
  cluster_id <- names(vascular_markers_scran)[ii]
  tmp <- vascular_markers_scran[[ii]]
  tmp <- tmp[tmp[['FDR']] < min_fdr, c('p.value', 'FDR')]
  tmp[['gene']] <- rownames(tmp)
  tmp[['cluster']] <- cluster_id
  vascular_markers_scran[[cluster_id]] <- tmp
}
vascular_markers_scran <- do.call(rbind, vascular_markers_scran)

# Same statistical test but no pairwise comparison. See Seurat::FindMarkers() for
# more detail.
vascular_markers_seurat <- FindAllMarkers(object = vascular,
                                          assay = 'RNA',
                                          slot = 'data',
                                          only.pos = TRUE,
                                          min.pct = 0.2)
write.table(x = vascular_markers_seurat,
            file = paste0(results_out, 'vascular_markers_seurat.tsv'),
            sep = '\t', row.names = FALSE)

# Heatmap with gene results from scran test
scran_genes <- vascular_markers_scran %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = -p.value) %>%
  select(gene)
expr_data <- t(ScaleData(object = vascular@assays$RNA@data[scran_genes$gene,])) %>%
  cbind(vascular@meta.data[c('vascular_subcluster')]) %>%
  group_by(vascular_subcluster) %>%
  summarise(across(.cols = scran_genes[['gene']], .fns = mean)) %>%
  reshape2::melt(id.vars = 'vascular_subcluster') %>%
  mutate(variable = factor(variable, levels = rev(unique(variable)))) 
max_expr <- 2.5
myColors <- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min(expr_data$value), 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))
de_markers_scran_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = vascular_subcluster, y = variable)) +
  geom_tile(mapping = aes(fill = value), color = 'black') +
  scale_fill_gradientn(colors = myColors,
                       values = scales::rescale(x = myBreaks, to = c(0,1)),
                       limits = c(NA, max_expr), 
                       na.value = myColors[100],
                       breaks = seq(-5,5,0.5)); de_markers_scran_heatmap


# Heatmap with gene results from Seurat
vascular_markers_seurat <- read.table(file = paste0(results_out, 'vascular_markers_seurat.tsv'),
                                      sep = '\t', header = TRUE)
seurat_genes <- vascular_markers_seurat %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_logFC) %>%
  # top_n(n = 10, wt = -p.value) %>%
  .[['gene']]
seurat_genes <- c(seurat_genes, c('Cldn5','Podxl','Pecam1','Cspg4','Pdgfrb','Lum','Pdgfra'))
expr_data <- t(ScaleData(object = vascular@assays$RNA@data[seurat_genes,])) %>%
  cbind(vascular@meta.data[c('vascular_subcluster')]) %>%
  group_by(vascular_subcluster) %>%
  summarise(across(.cols = all_of(seurat_genes), .fns = mean)) %>%
  reshape2::melt(id.vars = 'vascular_subcluster') %>%
  mutate(variable = factor(variable, levels = rev(unique(variable)))) 
max_expr <- 2.5
min_expr <- floor(min(expr_data$value)*10)/10
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'), bias = 1.1)(100))
myBreaks <- c(seq(min(expr_data$value), 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))
de_markers_seurat_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = vascular_subcluster, y = variable)) +
  geom_tile(mapping = aes(fill = value), color = 'black', size = 0.3) +
  scale_fill_gradientn(colors = myColors,
                       values = scales::rescale(x = myBreaks, to = c(0,1)),
                       limits = c(min_expr, max_expr), 
                       na.value = myColors[100],
                       breaks = c(min_expr,0,max_expr)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.position = 'right',
        legend.box.margin = margin(0,0,0,-2,unit='mm'),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(angle = 90, color = 'black', size = 14, vjust = 1, hjust = 0),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(fill = guide_colorbar(title = 'Scaled Expression',
                               title.position = 'left',
                               frame.colour = 'black',
                               ticks.colour = 'black',
                               frame.linewidth = 1,
                               ticks.linewidth = 1)); de_markers_seurat_heatmap
ggsave(filename = paste0(results_out, 'vascular_de_markers_heatmap.tiff'),
        plot = de_markers_seurat_heatmap, device = 'tiff',
        height = 11, width = 4.5)




# DE markers heatmap (with significant p-values) ------------------------------

Idents(vascular) <- 'vascular_subcluster'
DefaultAssay(vascular) <- 'RNA'

vascular_markers_seurat <- FindAllMarkers(object = vascular,
                                          assay = 'RNA',
                                          slot = 'data',
                                          only.pos = TRUE,
                                          min.pct = 0.2)
write.table(x = vascular_markers_seurat,
            file = paste0(results_out, 'vascular_markers_seurat.tsv'),
            sep = '\t', row.names = FALSE)

# Heatmap with gene results from Seurat
vascular_markers_seurat <- read.table(file = paste0(results_out, 'vascular_markers_seurat.tsv'),
                                      sep = '\t', header = TRUE)
seurat_genes <- vascular_markers_seurat %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_logFC) %>%
  # top_n(n = 10, wt = -p.value) %>%
  .[['gene']]
seurat_genes <- c(seurat_genes, c('Cldn5','Podxl','Pecam1','Cspg4','Pdgfrb','Lum','Pdgfra'))
expr_data <- t(ScaleData(object = vascular@assays$RNA@data[seurat_genes,])) %>%
  cbind(vascular@meta.data[c('vascular_subcluster')]) %>%
  group_by(vascular_subcluster) %>%
  summarise(across(.cols = all_of(seurat_genes), .fns = mean)) %>%
  reshape2::melt(id.vars = 'vascular_subcluster') %>%
  mutate(variable = factor(variable, levels = rev(unique(variable)))) 

expr_data$pval <- 1
for (i in 1:nrow(vascular_markers_seurat)) {
  tmp <- which(as.character(expr_data$vascular_subcluster) == as.character(vascular_markers_seurat$cluster[i]) &
                 as.character(expr_data$variable) == as.character(vascular_markers_seurat$gene[i]))
  expr_data$pval[tmp] <- vascular_markers_seurat$p_val_adj[i]
}
expr_data$significant <- ifelse(test = expr_data$pval < 1e-10,
                                yes = '*',
                                no = ' ')

max_expr <- 3
min_expr <- -3
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'), bias = 1.1)(100))

de_markers_seurat_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = vascular_subcluster, y = variable)) +
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
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.position = 'right',
        legend.box.margin = margin(0,0,0,-2,unit='mm'),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(angle = 90, color = 'black', size = 14, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(fill = guide_colorbar(title = 'z-score',
                               title.position = 'left',
                               frame.colour = 'black',
                               ticks.colour = 'black',
                               frame.linewidth = 1,
                               ticks.linewidth = 1))
de_markers_seurat_heatmap
ggsave(filename = './results/revision_figures/vascular_de_markers_heatmap.tiff',
       plot = de_markers_seurat_heatmap, device = 'tiff',
       height = 11, width = 4.5)




# Population dynamics -----------------------------------------------------

counts <- data.frame(table(vascular$vascular_subcluster, vascular$time))
names(counts) <- c('subcluster','time','count')
counts[['time']] <- plyr::mapvalues(
  x = counts[['time']],
  from = c('Uninjured'),
  to = c('Uninj')
)
count_these <- c('A-Endothelial',
                 'C1-Endothelial',
                 'C2-Endothelial',
                 'V-Endothelial',
                 'Tip Cell',
                 'Pericyte',
                 'VSMC',
                 'Fibroblast')

vascular_prop <- counts %>%
  filter(subcluster %in% count_these) %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) + 
  geom_bar(aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black', 
           size = 0.75,
           position = 'fill') + 
  scale_y_continuous() + 
  scale_x_discrete() +
  scale_fill_manual(values = vascular_cols) +
  ylab('Prop. of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'),
        legend.spacing.y = unit(x = 2.5, units = 'mm')) + 
  guides(fill = guide_legend(ncol = 1))
vascular_prop
ggsave(filename = paste0(results_out, 'vascular_population_prop_dynamics.tiff'),
       plot = vascular_prop, device = 'tiff', height = 2.75, width = 4.75)

vascular_counts <- counts %>%
  filter(subcluster %in% count_these) %>%
  ggplot(mapping = aes(x = time, y = count, group = subcluster)) +
  geom_bar(mapping = aes(fill = subcluster), 
           stat = 'identity', 
           color = 'black',
           size = 0.75) +
  scale_y_continuous(breaks = seq(0, 5000, 500)) + 
  scale_x_discrete() +
  scale_fill_manual(values = vascular_cols) +
  ylab(label = '# of cells') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 12, color = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.position = 'none', 
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = NA))

vascular_legend <- cowplot::get_legend(plot = vascular_prop)
vascular_pop <- cowplot::plot_grid(vascular_counts, vascular_prop + theme(legend.position = 'none'), ncol = 2, rel_widths = c(1, 1))
vascular_pop <- cowplot::plot_grid(vascular_pop, vascular_legend, ncol = 1, rel_heights = c(1, 0.4))
ggsave(filename = paste0(results_out, 'vascular_population_dynamics.tiff'),
       plot = vascular_pop, device = 'tiff', height = 4, width = 5.25)





# Apelin percent expression -----------------------------------------------

DefaultAssay(vascular) <- 'RNA'
apln_expr <- FetchData(object = vascular,
                       vars = c('time','vascular_subcluster','sample_id','Apln'),
                       slot = 'data')
apln_pct_vascular <- apln_expr %>%
  group_by(sample_id) %>%
  summarise(pct = mean(Apln > 0),
            time = unique(time)) %>%
  ggplot(mapping = aes(x = time, y = pct)) +
  geom_boxplot(mapping = aes(x = time, y = pct, fill = time)) +
  geom_point(mapping = aes(fill = time), 
             stat = 'identity', 
             size = 3, 
             pch = 21) +
  ylab(label = 'Apln(+) cells / all vascular cells') +
  scale_fill_manual(values = time_cols) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(axis.text.x = element_text(size = 16, 
                                   color = 'black', 
                                   angle = 45, 
                                   hjust = 1),
        axis.text.y = element_text(size = 14,
                                   color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14,
                                    color = 'black',
                                    hjust = 0,
                                    margin = margin(0, 5, 0, 0, unit = 'pt')),
        legend.title = element_text(size = 14, 
                                    color = 'black'),
        legend.text = element_text(size = 12, 
                                   color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.position = 'none')
ggsave(filename = paste0(results_out, 'apln_pct_vascular.tiff'),
       plot = apln_pct_vascular, device = 'tiff', height = 4, width = 3)


subcluster_cols <- c('Tip Cell' = 'indianred',
                     'Non-Tip Cell' = 'dodgerblue')
apln_pct_tipcell <- apln_expr %>%
  mutate(subcluster = ifelse(vascular_subcluster == 'Tip Cell',
                             yes = 'Tip Cell',
                             no = 'Non-Tip Cell')) %>%
  group_by(sample_id, subcluster) %>%
  summarise(pct = mean(Apln > 0),
            time = unique(time)) %>%
  ggplot(mapping = aes(x = time, y = pct)) +
  geom_boxplot(mapping = aes(x = time, y = pct, fill = subcluster), 
               position = position_dodge(0.75)) +
  geom_point(mapping = aes(group = subcluster, fill = subcluster), 
             stat = 'identity', 
             size = 2,
             pch = 21,
             position = position_dodge(0.75)) +
  ylab(label = 'Apln(+) cells / all vascular cells') +
  scale_fill_manual(values = subcluster_cols) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(axis.text.x = element_text(size = 16, 
                                   color = 'black', 
                                   angle = 45, 
                                   hjust = 1),
        axis.text.y = element_text(size = 14,
                                   color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14,
                                    color = 'black',
                                    hjust = 0,
                                    margin = margin(0, 5, 0, 0, unit = 'pt')),
        legend.title = element_text(size = 14, 
                                    color = 'black'),
        legend.text = element_text(size = 12, 
                                   color = 'black'),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(1, units = 'cm'),
        panel.border = element_rect(color = 'black', fill = NA)) +
  guides(fill = guide_legend(title = 'Vascular\nsubcluster'),
         group = guide_legend(override.aes = list(size = 3))); apln_pct_tipcell
ggsave(filename = paste0(results_out, 'apln_pct_tipcell.tiff'),
       plot = apln_pct_tipcell, device = 'tiff', height = 4, width = 5)




# Abundance analysis ------------------------------------------------------

time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

vascular_subcluster_counts <- data.frame(table(vascular$vascular_subcluster, vascular$sample_id))
colnames(vascular_subcluster_counts) <- c('vascular_subcluster','sample_id','count')
vascular_subcluster_counts[['time']] <- substr(x = vascular_subcluster_counts[['sample_id']],
                                              start = 1,
                                              stop = regexpr(pattern = '_', text = vascular_subcluster_counts[['sample_id']]) -1)
vascular_subcluster_counts[['time']] <- plyr::mapvalues(x = vascular_subcluster_counts[['time']], from ='uninj', to = 'Uninjured')
vascular_subcluster_counts[['time']] <- factor(vascular_subcluster_counts[['time']], levels = levels(vascular$time))
vascular_subcluster_counts_plot <- vascular_subcluster_counts %>%
  ggplot(mapping = aes(x = sample_id, y = count)) +
  geom_bar(mapping = aes(fill = time), color = 'black', stat = 'identity') +
  labs(title = 'Cell abundance across time and sample') +
  ylab(label = 'Number of cells') +
  facet_wrap(. ~ vascular_subcluster, scales = 'free_y', ncol = 5) +
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
ggsave(filename = paste0(results_out, 'vascular_subcluster_counts.tiff'),
       plot = vascular_subcluster_counts_plot, device = 'tiff', height = 6.5, width = 15)



