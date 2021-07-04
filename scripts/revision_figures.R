
require('Seurat')
require('dplyr')
require('ggplot2')

results_out <- 'results/revision_figures/'
dir.create(path = results_out)

sci <- readRDS(file = 'data/sci.rds')
myeloid <- readRDS(file = 'data/myeloid.rds')

DefaultAssay(sci) <- 'RNA'
DefaultAssay(myeloid) <- 'RNA'



# Functional subcluster nomenclature --------------------------------------

# NOTE: Original analysis scripts were modified to reflect changes to myeloid
# subcluster annotations. This chunk here only shows how the annotations were
# relabelled. Code to generate figures shown in paper are in respective analysis
# files, but figure results are saved to the "revision_figures" directory.


require('Seurat')
require('ggplot2')
require('dplyr')

myeloid <- readRDS(file = './data/myeloid.rds')

dir.create(path = './results/revision_figures/functional_names/')

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
Idents(myeloid) <- 'myeloid_functional'

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
myeloid_subcluster_counts <- table(myeloid$myeloid_functional)
myeloid_subcluster_label <- paste0(names(myeloid_subcluster_counts), ' (', myeloid_subcluster_counts, ')')
names(myeloid_subcluster_label) <- names(myeloid_subcluster_counts)
myeloid_subcluster_umap <- FetchData(object = myeloid, vars = c('UMAP_1','UMAP_2','myeloid_functional')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = myeloid_functional), size = 0.2, alpha = 0.5) +
  scale_color_manual(values = myeloid_cols, 
                     breaks = names(myeloid_subcluster_label),
                     label = myeloid_subcluster_label) +
  xlab(label = 'UMAP 1') +
  ylab(label = 'UMAP 2') +
  umap_theme +
  theme(legend.text = element_text(size = 16, color = 'black')) +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))
myeloid_subcluster_umap
ggsave(filename = './results/revision_figures/functional_names/myeloid_annotation_umap.tiff',
       plot = myeloid_subcluster_umap, 
       device = 'tiff', height = 6, width = 10.25)



# DAM genes ---------------------------------------------------------------

dam_genes <- c('Apoe','Trem2','Igf1','Msr1','Cdk1','Il1b','P2ry12','Spp1')

dam_umap <- FeaturePlot(
  object = sci,
  features = dam_genes,
  combine = FALSE,
  order = TRUE,
  split.by = 'time'
)
dam_umap <- lapply(
  X = dam_umap,
  FUN = function(x) {
    x + theme_bw() + NoLegend() +
      theme(axis.title.x.bottom = element_blank(),
            axis.title.y.left = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank())
  }
)
dam_umap <- cowplot::plot_grid(
  plotlist = dam_umap, 
  ncol = length(dam_genes),
  rel_widths = c(rep(1, length(dam_genes)-1), 1.1),
  rel_heights = c(1.1, rep(1, length(levels(myeloid$time)) - 1))
)
ggsave(filename = paste0(results_out, 'DAM-genes-sci-time_umap.tiff'),
       plot = dam_umap, device = 'tiff', height = 10, width = 20)


dam_umap <- FeaturePlot(
  object = myeloid,
  features = dam_genes,
  combine = FALSE,
  # order = TRUE,
  split.by = 'time',
  pt.size = 0.5
)
dam_umap <- lapply(
  X = dam_umap,
  FUN = function(x) {
    x + theme_bw() + NoLegend() +
      theme(axis.title.x.bottom = element_blank(),
            axis.title.y.left = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank())
  }
)
dam_umap <- cowplot::plot_grid(
  plotlist = dam_umap, 
  ncol = length(dam_genes),
  rel_widths = c(rep(1, length(dam_genes)-1), 1.1),
  rel_heights = c(1.1, rep(1, length(levels(myeloid$time)) - 1))
)
ggsave(filename = paste0(results_out, 'DAM-genes-myeloid-time_umap.tiff'),
       plot = dam_umap, device = 'tiff', height = 10, width = 20)
ggsave(filename = paste0(results_out, 'DAM-genes-myeloid-time_umap.png'),
       plot = dam_umap, device = 'png', height = 10, width = 20)




# Msr1 expression over time -----------------------------------------------

myeloid <- readRDS(file = './data/myeloid.rds')
DefaultAssay(myeloid) <- 'RNAcorrected'

p1 <- FetchData(
  object = myeloid, 
  vars = c('Msr1','time','UMAP_1','UMAP_2'), 
  slot = 'data'
) %>%
  arrange(Msr1) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = Msr1),
             size = 1) +
  facet_grid(. ~ time) +
  scale_color_gradient(low = 'grey', high = 'red3') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14))
p1
ggsave(filename = './results/revision_figures/msr1-over-time_myeloid.tiff',
       plot = p1, device ='tiff', height = 3, width = 12)





# Microglia DE heatmap with z-score, p-values -------------------------------


mg_markers <- read.table(
  file = 'results/myeloid_annotation_markers/microglia_DE_wilcox.tsv',
  sep = '\t',
  header = TRUE,
  row.names = 1
)
top_mg_markers <- mg_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = -p_val_adj) %>%
  top_n(n = 10, wt = avg_logFC)
microglia_de_genes <- top_mg_markers$gene

mg_subclusters <- c('H-Microglia',
                    'DAM-A',
                    'DAM-B',
                    'DAM-C',
                    'IFN-Myeloid')
microglia <- myeloid[,myeloid$myeloid_subcluster %in% mg_subclusters]
de_set <- FindAllMarkers(
  object = microglia,
  assay = 'RNA',
  features = microglia_de_genes,
  only.pos = TRUE
)

expr_data <- t(ScaleData(object = microglia@assays$RNA@data[microglia_de_genes,])) %>%
  cbind(microglia@meta.data[c('myeloid_subcluster')]) %>%
  group_by(myeloid_subcluster) %>%
  summarise(across(.cols = all_of(microglia_de_genes), .fns = mean)) %>%
  reshape2::melt(id.vars = 'myeloid_subcluster') %>%
  mutate(variable = factor(variable, levels = unique(variable)),
         myeloid_subcluster = factor(myeloid_subcluster, rev(levels(myeloid_subcluster))))

expr_data$pval <- 1
for (i in 1:nrow(de_set)) {
  tmp <- which(as.character(expr_data$myeloid_subcluster) == as.character(de_set$cluster[i]) &
                 as.character(expr_data$variable) == as.character(de_set$gene[i]))
  expr_data$pval[tmp] <- de_set$p_val_adj[i]
}
expr_data$significant <- ifelse(test = expr_data$pval < 1e-10,
                                yes = '*',
                                no = ' ')

max_expr <- 1.5
min_expr <- round(floor(min(expr_data$value)*10)/10,1)
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min(expr_data$value), 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))
microglia_de_genes_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_tile(mapping = aes(fill = value), color = 'black', size = 0.3) +
  geom_text(mapping = aes(label = significant)) +
  scale_fill_gradientn(colors = myColors,
                       values = scales::rescale(x = myBreaks, to = c(0,1)),
                       limits = c(min_expr, max_expr), 
                       na.value = myColors[100],
                       breaks = c(min_expr,0,max_expr)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = 'black'),
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
ggsave(filename = './results/revision_figures/microglia_de_heatmap_pvals.tiff',
       plot = microglia_de_genes_heatmap, device = 'tiff', height = 3, width = 10)






# Microglia markers Vln plot w/ p-values -----------------------------------------


# Subset cells
microglia <- c('H-Microglia',
               'DAM-A',
               'DAM-B',
               'DAM-C')
Idents(myeloid) <- 'myeloid_subcluster'
microglia <- subset(myeloid, idents = microglia)

microglia_markers <- c('P2ry12','Il1b','Igf1','Msr1','Cdk1')

DefaultAssay(microglia) <- 'RNAcorrected'
expr_dat <- FetchData(
  object = microglia,
  vars = c('myeloid_subcluster', microglia_markers),
  slot = 'data'
)

tmp_de <- FindAllMarkers(
  object = microglia, 
  assay = 'RNA',
  slot = 'data',
  features = microglia_markers,
  logfc.threshold = 0,
  only.pos = TRUE
)
tmp_de$p_val_adj <- signif(tmp_de$p_val_adj, digits = 2)
tmp_de$variable <- tmp_de$gene

tmp_ypos <- expr_dat %>%
  reshape2::melt(id.vars = 'myeloid_subcluster') %>%
  group_by(variable) %>%
  summarise('y.position' = max(value))
tmp_de$y.position <- tmp_ypos$y.position[match(tmp_de$gene, tmp_ypos$variable)]+0.25

expr_dat %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  ggplot(mapping = aes(x = myeloid_subcluster, y = value)) +
  geom_violin(mapping = aes(fill = myeloid_subcluster), scale = 'width') +
  geom_label(data = tmp_de, mapping = aes(x = cluster, y = y.position, label = p_val_adj), size = 3) +
  # stat_pvalue_manual(data = tmp_de, 
  #                    x = 'group1',
  #                    y.position = 'y.position',
  #                    label = 'p_val_adj',
  #                    size = 3,
  #                    hide.ns = TRUE,
  #                    group.by = 'facet') +
  facet_wrap(. ~ variable, ncol = 1, scales = 'free_y', strip.position = 'right') +
  # scale_y_continuous(limits = c(NA, max(test$y.position) + 4)) +
  xlab(label = 'Microglia subtype') + 
  ylab(label = 'log-normalized expression') + 
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(strip.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))





# Microglia DE dot plot ---------------------------------------------------


# Subset cells
microglia <- c('H-Microglia',
               'DAM-A',
               'DAM-B',
               'DAM-C')
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
avg_exp <- data.frame(Matrix::t(ScaleData(microglia[['RNA']]@data, features = microglia_markers)))
avg_exp <- cbind(avg_exp, 'myeloid_subcluster' = as.character(microglia@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- data.frame(Matrix::t(microglia[['RNA']]@counts[microglia_markers,]))
pct_exp <- cbind(pct_exp, 'myeloid_subcluster' = as.character(microglia@meta.data[,c('myeloid_subcluster')])) %>%
  reshape2::melt(id.vars = c('myeloid_subcluster')) %>%
  group_by(myeloid_subcluster, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
min_expr <- floor(min(avg_exp$avg.exp)*10)/10
max_expr <- ceiling(max(avg_exp$avg.exp)*10)/10

# Plot
microglia_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(myeloid_subcluster = factor(myeloid_subcluster, levels = rev(levels(myeloid$myeloid_subcluster)))) %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  geom_text(data = tmp_de, mapping = aes(x = variable, y = cluster, label = signif),
            size = 8, fontface = 'bold', nudge_y = -0.09, nudge_x = 0.01) +
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
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'bold', size = 20, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(5,0,10,0),
        legend.box = 'horizontal',
        legend.direction = 'vertical',
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
microglia_markers_dotplot
ggsave(filename = './results/revision_figures/microglia_markers_dotplot_igf1.tiff',
       plot = microglia_markers_dotplot, device = 'tiff', height = 2.75, width = 6.25)


# Macrophage DE markers heatmap  with z-score, p-values ------------------------

Idents(myeloid) <- 'myeloid_subcluster'
DefaultAssay(myeloid) <- 'RNA'

# Subset cells
macrophage <- c('Neutrophil',
                'Monocyte',
                'Macrophage-A',
                'Macrophage-B',
                'BA-Macrophage',
                'Dendritic',
                'Div-Myeloid')
Idents(myeloid) <- 'myeloid_subcluster'
macrophage <- subset(myeloid, idents = macrophage)

# Use Seurat wilcox test
mac_markers <- read.table(
  file = 'results/myeloid_annotation_markers/macrophage_DE_wilcox.tsv',
  sep = '\t',
  header = TRUE,
  row.names = 1
)
top_mac_markers <- mac_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = -p_val_adj) %>%
  top_n(n = 8, wt = avg_logFC)
macrophage_de_genes <- top_mac_markers$gene
de_set <- FindAllMarkers(
  object = macrophage,
  assay = 'RNA',
  features = macrophage_de_genes,
  only.pos = TRUE
)

# Build heatmap
expr_data <- t(ScaleData(object = macrophage@assays$RNA@data[macrophage_de_genes,])) %>%
  cbind(macrophage@meta.data[c('myeloid_subcluster')]) %>%
  group_by(myeloid_subcluster) %>%
  summarise(across(.cols = all_of(macrophage_de_genes), .fns = mean)) %>%
  reshape2::melt(id.vars = 'myeloid_subcluster') %>%
  mutate(variable = factor(variable, levels = unique(variable)),
         myeloid_subcluster = factor(myeloid_subcluster, levels = rev(levels(myeloid_subcluster))))

expr_data$pval <- 1
for (i in 1:nrow(de_set)) {
  tmp <- which(as.character(expr_data$myeloid_subcluster) == as.character(de_set$cluster[i]) &
                 as.character(expr_data$variable) == as.character(de_set$gene[i]))
  expr_data$pval[tmp] <- de_set$p_val_adj[i]
}
expr_data$significant <- ifelse(test = expr_data$pval < 1e-10,
                                yes = '*',
                                no = ' ')

max_expr <- 2
min_expr <- round(floor(min(expr_data$value)*10)/10,1)
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min(expr_data$value), 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))
macrophage_de_genes_heatmap <- expr_data %>%
  ggplot(mapping = aes(x = variable, y = myeloid_subcluster)) +
  geom_tile(mapping = aes(fill = value), color = 'black', size = 0.3) +
  geom_text(mapping = aes(label = significant)) +
  scale_fill_gradientn(colors = myColors,
                       values = scales::rescale(x = myBreaks, to = c(0,1)),
                       limits = c(min_expr, max_expr), 
                       na.value = myColors[100],
                       breaks = c(min_expr,0,max_expr)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.title = element_text(color = 'black', angle = 90, size = 14, hjust = 0),
        legend.text = element_text(size = 12, color = 'black'),
        legend.box.margin=margin(0,0,0,0),
        legend.margin = margin(0,0,0,-3, unit = 'mm')) +
  guides(fill = guide_colorbar(title = 'Scaled Expression',
                               title.position = 'left',
                               frame.colour = 'black',
                               ticks.colour = 'black',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               title.hjust = 1))
macrophage_de_genes_heatmap
ggsave(filename = 'results/revision_figures/macrophage_de_genes_heatmap.tiff',
       plot = macrophage_de_genes_heatmap, device = 'tiff',
       height = 3.5, width = 14)





# Microglia GO ------------------------------------------------------------


mg_markers <- read.table(
  file = 'results/myeloid_annotation_markers/microglia_DE_wilcox.tsv', 
  sep = '\t',
  header = TRUE,
  row.names = 1
)
go_genes <- mg_markers %>%
  mutate(pct.diff = pct.1 - pct.2) %>%
  filter(pct.diff > 0) %>%
  group_by(cluster)
ensembl_convert <- read.table('ref/gene_name_conversion.tsv', sep = '\t', header = TRUE)
go_genes$ensembl <- plyr::mapvalues(
  x = go_genes$gene,
  from = ensembl_convert$mgi_symbol,
  to = ensembl_convert$ensembl_gene_id,
  warn_missing = FALSE
)
write.csv(
  x = go_genes,
  file = 'results/revision_figures/microglia_go_genes.csv'
)

# Use gProfiler, use output
mg_go_results <- readxl::read_xlsx(path = 'results/revision_figures/microglia_GO_terms.xlsx')
h_mg <- mg_go_results[c(1,2)]
a_mg <- mg_go_results[c(4,5)]
b_mg <- mg_go_results[c(7,8)]
c_mg <- mg_go_results[c(10,11)]
mg_go_results <- list('Homeostatic Microglia' = h_mg, 
                      'Inflammatory Microglia' = a_mg, 
                      'Dividing Microglia' = b_mg,
                      'Migrating Microglia' = c_mg)
for (i in 1:length(mg_go_results)) {
  colnames(mg_go_results[[i]]) <- c('GO term', 'adj. p-value')
  mg_go_results[[i]]$Cluster = names(mg_go_results)[i]
  mg_go_results[[i]] <- mg_go_results[[i]][!apply(is.na(mg_go_results[[i]]), any, MARGIN = 1),]
}
mg_go_results <- Reduce(rbind, mg_go_results)
mg_go_results$log_pval <- -log10(mg_go_results$`adj. p-value`)
mg_go_results$tmp_id <- paste(
  mg_go_results$Cluster,
  mg_go_results$`GO term`,
  sep = '_'
)
mg_go_results$Cluster <- factor(
  x = mg_go_results$Cluster,
  levels = c('Homeostatic Microglia','Inflammatory Microglia','Dividing Microglia','Migrating Microglia')
)

mg_go_plot <- mg_go_results %>%
  group_by(Cluster) %>%
  top_n(n = 15, wt = log_pval) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ Cluster, 
             scales = 'free_y', 
             drop = TRUE, 
             ncol = 2, 
             labeller = label_wrap_gen()) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(adj. p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 16, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
mg_go_plot
ggsave(filename = './results/revision_figures/microglia_GO_plot.tiff',
       plot = mg_go_plot, device = 'tiff', height = 8, width = 14)



# SCI Cellcycle regression ----------------------------------------------------

# Dont ask me to do this again. Takes 3 hours...

sci_tmp <- readRDS(file = './data/sci.rds')
sci_celltype <- sci_tmp@meta.data['celltype']
rm(sci_tmp);gc()

# Set unique gene names (i.e. remove duplicates)
remove_duplicated_genes <- function(x) {
  bad_gene <- duplicated(rownames(x)) | (nchar(rownames(x)) == 0)
  counts <- x[!bad_gene,]
  return(counts)
}
counts <- lapply(X = counts, FUN = remove_duplicated_genes)


# Initialize Seurat objects
for(ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  sample_name <- gsub(pattern = 'qc_filtered_feature_bc_matrix_',
                      replacement = '',
                      x = sample_id)
  counts[[sample_id]] <- counts[[sample_id]] %>%
    CreateSeuratObject(project = sample_name) %>%
    AddMetaData(metadata = sample_metadata[[sample_id]]) %>%
    NormalizeData() %>%
    CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
}

# Cell-cycle module scoring
counts <- lapply(X = counts, 
                 FUN = function(x) {
                   x$CC.Difference <- x$S.Score - x$G2M.Score
                   return(x)
                 })

counts <- lapply(X = counts,
                 FUN = FindVariableFeatures,
                 nfeatures = 2500)

sci_anchors <- FindIntegrationAnchors(
  object.list = counts, 
  normalization.method = 'LogNormalize'
)

# Integrate into single dataset
sci <- IntegrateData(sci_anchors)
DefaultAssay(sci) <- 'integrated'
sci_cells <- colnames(sci)[colnames(sci) %in% rownames(sci_celltype)]
sci <- sci[, sci_cells]

sci <- ScaleData(sci)
sci <- RunPCA(sci, npcs = 40)
ElbowPlot(sci, ndims = 40)
sci <- sci %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

sci$celltype <- plyr::mapvalues(
  x = rownames(sci@meta.data),
  from = rownames(sci_celltype),
  to = as.character(sci_celltype$celltype)
)

p1 <- lapply(
  X = DimPlot(
    object = sci, 
    group.by = c('celltype', 'integrated_snn_res.0.8', 'Phase'),
    label = TRUE,
    repel = TRUE,
    combine = FALSE
  ),
  FUN = function(x) x + theme_bw()
)
p1[1:2] <- lapply(X = p1[1:2], FUN = function(x) x + NoLegend())
p1 <- Reduce(
  f = `+`,
  x = p1
)

sci <- ScaleData(sci, vars.to.regress = 'CC.Difference')
sci <- RunPCA(sci, npcs = 40)
ElbowPlot(sci, ndims = 40)
sci <- sci %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

p2 <- lapply(
  X = DimPlot(
    object = sci, 
    group.by = c('celltype', 'integrated_snn_res.0.8', 'Phase'),
    label = TRUE,
    repel = TRUE,
    combine = FALSE
  ),
  FUN = function(x) x + theme_bw()
)
p2[1:2] <- lapply(X = p2[1:2], FUN = function(x) x + NoLegend())
p2 <- Reduce(
  f = `+`,
  x = p2
)
p3 <- (p1 / p2) + patchwork::plot_annotation(
  title = 'Effect on cell-cycle regression on myeloid clustering',
  subtitle = 'Top row: cell-cycle regression applied. Bottom row: cell-cycle regression NOT applied.'
)
ggsave(filename = './results/revision_figures/cellCycleRegressionComparison_sci-log_cluster.tiff',
       plot = p3, device = 'tiff', height = 8, width = 13.33)




# Myeloid Cellcycle regression ----------------------------------------------------

myeloid <- readRDS('data/myeloid.rds')
DefaultAssay(myeloid) <- 'integrated'
myeloid <- FindClusters(myeloid, resolution = 0.35)
myeloid$myeloid_subcluster <- plyr::mapvalues(
  x = myeloid$myeloid_subcluster,
  from = c('Neutrophil',
                'Monocyte',
                'Macrophage-A',
                'Macrophage-B',
                'BA-Macrophage',
                'Dendritic',
                'Div-Myeloid',
                'H-Microglia',
                'DAM-A',
                'DAM-B',
                'DAM-C',
                'IFN-Myeloid'),
  to = c('Neutrophil',
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

myeloid_unscaled <- ScaleData(myeloid)
myeloid_unscaled <- RunPCA(myeloid_unscaled, npcs = 30)
ElbowPlot(myeloid_unscaled, ndims = 30)
myeloid_unscaled <- FindNeighbors(myeloid_unscaled, dims = 1:11)
myeloid_unscaled <- RunUMAP(myeloid_unscaled, dims = 1:11)
myeloid_unscaled <- FindClusters(myeloid_unscaled, resolution = 0.35)

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
p1 <- lapply(
  X = DimPlot(
    object = myeloid, 
    group.by = c('myeloid_subcluster','integrated_snn_res.0.35','Phase'),
    label = TRUE,
    combine = FALSE
  ),
  FUN = function(x) x + theme_bw()
)
p1[[1]] <- p1[[1]] + scale_color_manual(values = myeloid_cols)
p1[1:2] <- lapply(X = p1[1:2], FUN = function(x) x + NoLegend())
p1 <- Reduce(
  f = `+`,
  x = p1
)
p2 <- lapply(
  X = DimPlot(
    object = myeloid_unscaled, 
    group.by = c('myeloid_subcluster','integrated_snn_res.0.35','Phase'),
    label = TRUE,
    combine = FALSE
  ),
  FUN = function(x) x + theme_bw()
)
p2[[1]] <- p2[[1]] + scale_color_manual(values = myeloid_cols)
p2[1:2] <- lapply(X = p2[1:2], FUN = function(x) x + NoLegend())
p2 <- Reduce(
  f = `+`,
  x = p2
)

DefaultAssay(myeloid) <- 'RNA'
DefaultAssay(myeloid_unscaled) <- 'RNA'
genes <- c('P2ry12','Mki67','Ms4a7')
p3 <- Reduce(
  f = `+`,
  x = lapply(
    X = FeaturePlot(object = myeloid, features = genes, order = TRUE, combine = FALSE),
    FUN = function(x) x + theme_bw()
  )
)
p4 <-  Reduce(
  f = `+`,
  x = lapply(
    X = FeaturePlot(myeloid_unscaled, features = genes, order = TRUE, combine = FALSE),
    FUN = function(x) x + theme_bw()
  )
)
p5 <- (p1 / p2) + patchwork::plot_annotation(
  title = 'Effect on cell-cycle regression on myeloid clustering',
  subtitle = 'Top row: cell-cycle regression applied. Bottom row: cell-cycle regression NOT applied.'
)
p6 <- (p3 / p4) + patchwork::plot_annotation(
  title = 'Effect on cell-cycle regression on myeloid clustering',
  subtitle = 'Top row: cell-cycle regression applied. Bottom row: cell-cycle regression NOT applied.'
)

ggsave(filename = './results/revision_figures/cellCycleRegressionComparison_myeloid_cluster.tiff',
       plot = p5, device = 'tiff', height = 8, width = 13.33)
ggsave(filename = './results/revision_figures/cellCycleRegressionComparison_myeloid_genes.tiff',
       plot = p6, device = 'tiff', height = 6, width = 10)


rm(myeloid_unscaled, p1, p2, p3, p4, p5, p6)



# Macroglia Cellcycle regression ----------------------------------------------------

macroglia <- readRDS(file = './data/macroglia.rds')
DefaultAssay(macroglia)

macroglia_unscaled <- ScaleData(macroglia)
macroglia_unscaled <- RunPCA(macroglia_unscaled, npcs = 30)
ElbowPlot(macroglia_unscaled, ndims = 30)
npcs <- 1:8
macroglia_unscaled <- RunUMAP(macroglia_unscaled, dims = npcs)
macroglia_unscaled <- FindNeighbors(macroglia_unscaled, dims = npcs)
macroglia_unscaled <- FindClusters(macroglia_unscaled, resolution = 0.8)
DimPlot(macroglia_unscaled, label = TRUE, label.size = 6)

p1 <- lapply(
  X = DimPlot(
    object = macroglia, 
    group.by = c('macroglia_subcluster','integrated_snn_res.0.8','Phase'),
    label = TRUE,
    combine = FALSE
  ),
  FUN = function(x) x + theme_bw()
)
p1[1:2] <- lapply(X = p1[1:2], FUN = function(x) x + NoLegend())
p1 <- Reduce(
  f = `+`,
  x = p1
)
p2 <- lapply(
  X = DimPlot(
    object = macroglia_unscaled, 
    group.by = c('macroglia_subcluster','integrated_snn_res.0.8','Phase'),
    label = TRUE,
    combine = FALSE
  ),
  FUN = function(x) x + theme_bw()
)
p2[1:2] <- lapply(X = p2[1:2], FUN = function(x) x + NoLegend())
p2 <- Reduce(
  f = `+`,
  x = p2
)

DefaultAssay(macroglia) <- 'RNA'
DefaultAssay(macroglia_unscaled) <- 'RNA'
genes <- c('Cspg4','Mki67','Gfap')
p3 <- Reduce(
  f = `+`,
  x = lapply(
    X = FeaturePlot(object = macroglia, features = genes, order = TRUE, combine = FALSE),
    FUN = function(x) x + theme_bw()
  )
)
p4 <-  Reduce(
  f = `+`,
  x = lapply(
    X = FeaturePlot(macroglia_unscaled, features = genes, order = TRUE, combine = FALSE),
    FUN = function(x) x + theme_bw()
  )
)
p5 <- (p1 / p2) + patchwork::plot_annotation(
  title = 'Effect on cell-cycle regression on macroglia clustering',
  subtitle = 'Top row: cell-cycle regression applied. Bottom row: cell-cycle regression NOT applied.'
)
p6 <- (p3 / p4) + patchwork::plot_annotation(
  title = 'Effect on cell-cycle regression on macroglia clustering',
  subtitle = 'Top row: cell-cycle regression applied. Bottom row: cell-cycle regression NOT applied.'
)
ggsave(filename = './results/revision_figures/cellCycleRegressionComparison_macroglia_cluster.tiff',
       plot = p5, device = 'tiff', height = 8, width = 13.33)
ggsave(filename = './results/revision_figures/cellCycleRegressionComparison_macroglia_genes.tiff',
       plot = p6, device = 'tiff', height = 6, width = 10)

rm(macroglia_unscaled, p1, p2, p3, p4, p5, p6)



# Unbiased LR analysis ----------------------------------------------------


source('./scripts/ligand_receptor_analysis_functions.R')

sci <- readRDS(file = './data/sci.rds')


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

previous_genes <- c('Pdgfrb','Pdgfra','Pdgfb','Pdgfa')
lr_ref <- read.csv(file = './ref/fantom_PairsLigRec_mouse.csv')
keep_rows <- lr_ref$Ligand.ApprovedSymbol %in% previous_genes | 
  lr_ref$Receptor.ApprovedSymbol %in% previous_genes 
lr_ref <- lr_ref[keep_rows,]
lr_ref <- lr_ref[lr_ref$Pair.Evidence == 'literature supported' & lr_ref$Pair.Source == 'known',]


Idents(sci) <- 'celltype'

# Make lightweight data structures
DefaultAssay(sci) <- 'RNA'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL

# Set identities and calculate setup values (avg_exp, pct_exp, cell counts)
sci_lr_setup <- setupLR(seurat_object = sci,
                        lr_ref = lr_ref,
                        split_by = 'time',
                        assay = 'RNA',
                        slot = 'data')
gc()

# Permutation test
sci_result <- calculateLR(setup = sci_lr_setup,
                          resample = 1000,
                          adjust_pval = FALSE)

# Sample plot
tmp <- plotLR(results = sci_result)
ggsave(filename = 'results/revision_figures/pdgf_signaling_gold-standard.tiff',
       plot = tmp, device = 'tiff', height = 25, width = 8)





# Axon-tract Microglia DEG comparison -------------------------------------

mg_markers <- read.table(
  file = 'results/myeloid_annotation_markers/microglia_DE_wilcox.tsv', 
  sep = '\t',
  header = TRUE,
  row.names = 1
)

ks_dam <- read.table(file = '../Schachtrup_collab/ref/KerenShaul2017_TableS3_HomeostaticMG_vs_DAM.txt',
                     sep = '\t', header = TRUE)
li_p7c1 <- read.table(file = '../Schachtrup_collab/ref/Li2019_TableS5_MG_P7C1_DEgenes.txt',
                      sep = '\t', header = TRUE)
hammond_irm <- read.table(file = '../Schachtrup_collab/ref/Hammond2019_TableS1_IRM2_DEgenes.txt',
                          sep = '\t', header = TRUE)
hammond_c4 <- read.table(file = '../Schachtrup_collab/ref/Hammond2019_TableS1_MGcluster4_DEgenes.txt',
                         sep = '\t', header = TRUE)


# Keren-Shaul, 2017, DAMs vs Homeostatic Microglia
ks_dam_set <- with(
  data = ks_dam,
  expr = Fold.change..DAM.to.homeostatic.microglia. > 0 &
    DAM.FDR.p.value >= 3 &
    !is.na(DAM.FDR.p.value))
ks_dam_set <- unique(ks_dam[ks_dam_set,'Gene.name'])
length(ks_dam_set)

# Li, 2019, P7-C2 (phagocytic?) microglia vs other P7 clusters
li_p7c1_set <- with(
  data = li_p7c1,
  expr = avg_logFC > 0 &
    p_val_adj <= 0.001
)
li_p7c1_set <- unique(li_p7c1[li_p7c1_set,'gene'])
length(li_p7c1_set)

# Hammond, 2019, Injury-responsive microglia cluster 2 (IRM2) vs other microglia
# (which come from control and saline-injected)
hammond_irm_set <- with(
  data = hammond_irm,
  expr = Fold.Change > 0 & # all are positive bc theyre marker genes
    Padj..FDR. <= 0.001
)
hammond_irm_set <- unique(hammond_irm[hammond_irm_set,'Gene'])
length(hammond_irm_set)

# Hammond, 2019, P4/5 cluster 4 vs other microglia (Axon tract-associated MG)
hammond_c4_set <- with(
  data = hammond_c4,
  expr = Fold.Change > 0 & # all are positive bc theyre marker genes
    Padj..FDR. <= 0.001
)
hammond_c4_set <- unique(hammond_c4[hammond_c4_set,'Gene'])
length(hammond_c4_set)


par(mfrow = c(2,2))
for (id in unique(mg_markers$cluster)) {
  gene_set <- with(
    data = mg_markers,
    expr = cluster == id &
      avg_logFC > 0 &
      p_val_adj <= 1e-10
  )
  gene_set <- unique(mg_markers$gene[gene_set])
  message(length(gene_set))
  
  gene_lists <- list(
    'DAM (KerenShaul)' = ks_dam_set,
    # 'PAM (Li)' = li_p7c1_set,
    # 'IRM (Hammond)' = hammond_irm_set,
    'ATM (Hammond)' = hammond_c4_set,
    gene_set
  )
  names(gene_lists)[length(gene_lists)] <- id
  tiff(filename = paste0('results/revision_figures/DEG-overlap_', id, '.tiff'),
       height = 4, width = 5, units = 'in', res = 440)
  print(plot(euler(gene_lists), quantities = TRUE))
  dev.off()
}


# Taking top 100 DEG only
ks_dam_set <- ks_dam %>%
  filter(Fold.change..DAM.to.homeostatic.microglia. > 0) %>%
  filter(DAM.FDR.p.value >= 3) %>%
  filter(!is.na(DAM.FDR.p.value)) %>%
  top_n(n = 100, wt = DAM.FDR.p.value)
ks_dam_set <- unique(ks_dam_set$Gene.name)
length(ks_dam_set)

li_p7c1_set <- li_p7c1 %>%
  filter(avg_logFC > 0 & p_val_adj < 1e-03) %>%
  top_n(n = 100, wt = -p_val_adj)
li_p7c1_set <- unique(li_p7c1_set$gene)
length(li_p7c1_set)

hammond_irm_set <- hammond_irm %>%
  filter(Fold.Change > 0 & Padj..FDR. <= 0.001) %>%
  top_n(n = 100, wt = -Padj..FDR.) %>%
  top_n(n = 100, wt = Fold.Change)
hammond_irm_set <- unique(hammond_irm_set$Gene)
length(hammond_irm_set)

hammond_c4_set <- hammond_c4 %>%
  filter(Fold.Change > 0 & Padj..FDR. <= 0.001) %>%
  top_n(n = 100, wt = -Padj..FDR.) %>%
  top_n(n = 100, wt = Fold.Change)
hammond_c4_set <- unique(hammond_c4_set$Gene)
length(hammond_c4_set)


par(mfrow = c(2,2))
for (id in unique(mg_markers$cluster)) {
  gene_set <- mg_markers %>%
    filter(cluster == id & avg_logFC > 0 & p_val_adj < 1e-10) %>%
    top_n(n = 100, wt = -p_val_adj) %>%
    top_n(n = 100, wt = avg_logFC)
  gene_set <- unique(gene_set$gene)
  message(length(gene_set))
  
  gene_lists <- list(
    'DAM (KerenShaul)' = ks_dam_set,
    # 'PAM (Li)' = li_p7c1_set,
    # 'IRM (Hammond)' = hammond_irm_set,
    'ATM (Hammond)' = hammond_c4_set,
    gene_set
  )
  names(gene_lists)[length(gene_lists)] <- id
  tiff(filename = paste0('results/revision_figures/DEG-overlap_100genes_', id, '.tiff'),
       height = 4, width = 5, units = 'in', res = 440)
  print(plot(euler(gene_lists), quantities = TRUE))
  dev.off()
}




# IEG Stress-induced myeloid activation signature  ------------------------


# AddModuleScore
feats <- c('Rgs1', 'Hist2h2aa1', 'Hist1h4i', 'Nfkbiz', 'Klf2', 'Junb', 'Dusp1', 'Ccl3', 'Hspa1a', 'Hsp90aa1', 'Fos', 'Hspa1b', 'Jun', 'Jund', 'Nfkbid', 'Gem', 'Ccl4', 'Ier5', 'Txnip', 'Hist1h2bc', 'Zfp36', 'Hist1h1c', 'Egr1', 'Atf3', 'Rhob')
feat.use <- feats[feats %in% rownames(myeloid[['RNA']]@data)]
tmp_myeloid <- SplitObject(
  object = myeloid,
  split.by = 'orig.ident'
)
tmp_ieg_score <- lapply(
  X = tmp_myeloid,
  FUN = function(x) {
    x <- AddModuleScore(
      object = x,
      features = list('IEG' = feat.use),
      assay = 'RNA',
      name = 'IEG'
    )
    return(x$IEG1)
  }
)
# tmp_ieg_score <- tmp_ieg_score[sort(names(tmp_ieg_score))]
tmp_ieg_score <- unlist(tmp_ieg_score, use.names = FALSE)
myeloid$IEG <- tmp_ieg_score


umap <- function(obj, gene, split.var) {
  tmp <- FeaturePlot(obj, features = gene, split.by = split.var, order = TRUE,
                     combine = FALSE, pt.size = 0.5)
  tmp <- lapply(tmp, function(x) {
    x + 
      # scale_color_gradientn(
      #   colors = rev(RColorBrewer::brewer.pal(n=9,name='Spectral'))
      # ) +
      scale_color_viridis_c(option = 'A') +
      theme_bw() +
      theme(axis.text = element_blank()) +
      NoLegend()
  }
  )
  tmp[[5]] <- patchwork::plot_layout(ncol = 4)
  tmp <- Reduce(x = tmp, f = `+`)
  return(tmp)
}
p1 <- umap(myeloid, 'IEG', 'time')
p2 <- umap(myeloid, 'Fosb', 'time')
p3 <- VlnPlot(myeloid, features = c('IEG','Junb'), pt.size = 0, ncol = 2,
              group.by = 'myeloid_subcluster')
p4 <- p1 / p2 / p3
p4
ggsave(filename = 'results/revision_figures/IEG_microglia.tiff',
       plot = p4, device = 'tiff', height = 10, width = 12.5)




# SCI UMAP by sample ------------------------------------------------------

sci <- readRDS(file = './data/sci.rds')

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

# cell-type annotation split by time UMAP
sample_cols <- c('#800000','#f58231','#4363d8','#800000','#f58231','#4363d8','#800000','#f58231','#800000','#f58231')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')
celltype_counts <- table(sci$celltype)
celltype_label <- paste0(names(celltype_counts), ' (', celltype_counts, ')')
names(celltype_label) <- names(celltype_counts)
umap_sample <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 0.5, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 2) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))
umap_sample
ggsave(filename = './results/revision_figures/SCI_bySample_umap.tiff',
       plot = umap_sample, device = 'tiff', height = 7, width = 9.5)
rm(sci); gc()


# Myeloid UMAP by sample ------------------------------------------------------

myeloid <- readRDS(file = './data/myeloid.rds')

Idents(object = myeloid) <- 'myeloid_subcluster'

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
                    legend.title = element_text(size = 14, color = 'black', face = 'bold'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 14, color = 'black'))

# cell-type annotation split by time UMAP
sample_cols <- c('#800000','#f58231','#4363d8','#800000','#f58231','#4363d8','#800000','#f58231','#800000','#f58231')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')
celltype_counts <- table(myeloid$myeloid_subcluster)
celltype_label <- paste0(names(celltype_counts), ' (', celltype_counts, ')')
names(celltype_label) <- names(celltype_counts)
umap_sample <- FetchData(object = myeloid, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 0.5, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 2) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))
umap_sample
ggsave(filename = './results/revision_figures/myeloid_bySample_umap.tiff',
       plot = umap_sample, device = 'tiff', height = 7, width = 9.5)
rm(myeloid); gc()


# Vascular UMAP by sample ------------------------------------------------------

vascular <- readRDS(file = './data/vascular.rds')

Idents(object = vascular) <- 'vascular_subcluster'

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
                    legend.title = element_text(size = 14, color = 'black', face = 'bold'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 14, color = 'black'))

# cell-type annotation split by time UMAP
sample_cols <- c('#800000','#f58231','#4363d8','#800000','#f58231','#4363d8','#800000','#f58231','#800000','#f58231')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')
celltype_counts <- table(vascular$vascular_subcluster)
celltype_label <- paste0(names(celltype_counts), ' (', celltype_counts, ')')
names(celltype_label) <- names(celltype_counts)
umap_sample <- FetchData(object = vascular, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 1, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 2) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))
umap_sample
ggsave(filename = './results/revision_figures/vascular_bySample_umap.tiff',
       plot = umap_sample, device = 'tiff', height = 7, width = 9.5)
rm(vascular); gc()


# Macroglia UMAP by sample ------------------------------------------------------

macroglia <- readRDS(file = './data/macroglia.rds')

Idents(object = macroglia) <- 'macroglia_subcluster'

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
                    legend.title = element_text(size = 14, color = 'black', face = 'bold'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 14, color = 'black'))

# cell-type annotation split by time UMAP
sample_cols <- c('#800000','#f58231','#4363d8','#800000','#f58231','#4363d8','#800000','#f58231','#800000','#f58231')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')
celltype_counts <- table(macroglia$macroglia_subcluster)
celltype_label <- paste0(names(celltype_counts), ' (', celltype_counts, ')')
names(celltype_label) <- names(celltype_counts)
umap_sample <- FetchData(object = macroglia, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 1, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 2) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))
umap_sample
ggsave(filename = './results/revision_figures/macroglia_bySample_umap.tiff',
       plot = umap_sample, device = 'tiff', height = 7, width = 9.5)
rm(macroglia); gc()




# Macrophage WGCNA (by cell) --------------------------------------------------


# Summary: Total fail. Don't try WGCNA on scSeq data again.

library('Seurat')
library('dplyr')
library('ggplot2')
library('WGCNA')
require('doParallel')

myeloid <- readRDS('./data/myeloid.rds')

DefaultAssay(myeloid) <- 'integrated'
macrophage <- c('Monocyte','Macrophage-A','Macrophage-B')
macrophage <- myeloid[,myeloid$myeloid_subcluster %in% macrophage]

high_genes <- Matrix::rowSums(macrophage[['RNA']]@counts) > sqrt(ncol(macrophage))
high_genes <- Matrix::rowMeans(macrophage[['RNA']]@counts) > 0.1
high_genes <- VariableFeatures(myeloid)
# macrophage <- macrophage[high_genes,]

test <- macrophage[high_genes, sample(x = 1:ncol(macrophage), size = 1000)]
mac_mat <- data.frame(Matrix::t(test[['RNA']]@data))
# mac_mat <- Matrix::t(macrophage[['RNA']]@data)
# mac_cor <- cor(
#   x = mac_mat,
#   method = 'spearman'
# )

powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2))

# cl <- parallel::makeCluster(2)
# doParallel::registerDoParallel(cl = cl)
# disableWGCNAThreads()
sft <- pickSoftThreshold(
  data = mac_mat, 
  # dataIsExpr = FALSE,
  powerVector = powers,
  # corFnc = stats::cor,
  # corOptions = list(method = 'pearson'),
  networkType = "signed",
  verbose = 5
)


# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 1.3

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", 
     main = paste("Scale independence"))
text(x = sft$fitIndices[,1], 
     y = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex = cex1, 
     col="red")

# Red line corresponds to using an R^2 cut-off
abline(h = 0.80, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(x = sft$fitIndices[,1], 
     y = sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", 
     type = "n",
     main = paste("Mean connectivity"))
text(x = sft$fitIndices[,1], 
     y = sft$fitIndices[,5], 
     labels = powers, 
     cex = cex1,
     col = "red")


softPower = 12
adj <- adjacency(
  datExpr = mac_mat,
  type = 'signed',
  power = softPower
)

# Turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM <- TOMsimilarity(
  adjMat = adj, 
  TOMType = "signed"
)
SubGeneNames <- colnames(adj)
colnames(TOM) <- rownames(TOM) <- SubGeneNames
dissTOM <- 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average");

par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
tiff(filename = './results/revision_figures/myeloid_wgcna_dendrogram.tiff', height = 3.5, width = 8, units = 'in', res = 440)
# sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(mac_mat, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.8
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(mac_mat, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs


# Recalculate MEs with color labels

nGenes = ncol(mac_mat);
nSamples = nrow(mac_mat);





# Macrophage WGCNA (by sample) -------------------------------------------------

# From WGCNA, we identified modules of genes that were strongly weighted by individual samples. Gene Ontology analysis of a majority of these modules yielded no significant terms for some modules and broader terms, such as "regulation of biological quality", in others. 

# RNA pearson
# RNAcorrected spearman
# RNA spearman

require('Seurat')
require('dplyr')
require('ggplot2')
require('WGCNA')
require('doParallel')

myeloid <- readRDS('./data/myeloid.rds')

DefaultAssay(myeloid) <- 'integrated'
macrophage <- c('Monocyte','Macrophage-A','Macrophage-B')
macrophage <- myeloid[,myeloid$myeloid_subcluster %in% macrophage]
rm(myeloid); gc()

DefaultAssay(macrophage) <- 'RNA'
macrophage[['RNAcorrected']] <- NULL
macrophage[['integrated']] <- NULL

macrophage_bulk <- FetchData(
  object = macrophage,
  vars = c(rownames(macrophage[['RNA']]@counts), 'orig.ident'),
  slot = 'counts'
)
macrophage_bulk <- macrophage_bulk[!macrophage_bulk$orig.ident %in% c('uninj_sample1','uninj_sample2','uninj_sample3'),]
macrophage_bulk <- macrophage_bulk %>%
  group_by(orig.ident) %>%
  summarise(across(everything(), .fns = sum))
orig.ident <- macrophage_bulk$orig.ident
macrophage_bulk <- macrophage_bulk[,colnames(macrophage_bulk) != 'orig.ident']
macrophage_bulk <- t(as.matrix(macrophage_bulk))
colnames(macrophage_bulk) <- orig.ident
lib_size <- apply(macrophage_bulk, 2, sum)
macrophage_bulk <- sweep(x = macrophage_bulk, MARGIN = 2, FUN = '/', lib_size)

high_exprs_genes <- apply(
  X = macrophage_bulk,
  FUN = function(x) all(x > 0),
  MARGIN = 1
)
high_exprs_genes <- names(which(high_exprs_genes > 0))
macrophage_bulk <- macrophage_bulk[high_exprs_genes, ]
macrophage_bulk <- t(macrophage_bulk)

powers <- c(seq(4, 10, by = 1), seq(8, 40, by = 2))

sft <- pickSoftThreshold(
  data = macrophage_bulk, 
  # dataIsExpr = FALSE,
  powerVector = powers,
  # corFnc = stats::cor,
  corOptions = list(method = 'spearman'),
  networkType = "signed",
  verbose = 5
)

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 1.3

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", 
     main = paste("Scale independence"))
text(x = sft$fitIndices[,1], 
     y = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex = cex1, 
     col="red")

# Red line corresponds to using an R^2 cut-off
abline(h = 0.80, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(x = sft$fitIndices[,1], 
     y = sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", 
     type = "n",
     main = paste("Mean connectivity"))
text(x = sft$fitIndices[,1], 
     y = sft$fitIndices[,5], 
     labels = powers, 
     cex = cex1,
     col = "red")


softPower = 24
adj <- adjacency(
  datExpr = macrophage_bulk,
  corOptions = list(method = 'spearman'),
  type = 'signed',
  power = softPower
)

# Turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM <- TOMsimilarity(
  adjMat = adj, 
  TOMType = "signed"
)
SubGeneNames <- colnames(adj)
colnames(TOM) <- rownames(TOM) <- SubGeneNames
dissTOM <- 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
tiff(filename = './results/revision_figures/myeloid_wgcna_bySample_dendrogram.tiff', height = 3.5, width = 8, units = 'in', res = 440)
# sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(macrophage_bulk, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.75
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(macrophage_bulk, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs


# Export genes in each module
intModules <- unique(moduleColors)
dir.create(path = './results/revision_figures/WGCNA_results/')
for (module in intModules) {
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their gene names
  modGenes <- colnames(macrophage_bulk)[modGenes]
  # Write them into a file
  fileName = paste("./results/revision_figures/WGCNA_results/ModuleGenes_", module, ".txt", sep="")
  write.table(as.data.frame(modGenes), file = fileName,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}




# Monocle Trajectory Inference --------------------------------------------

dir.create('results/revision_figures/TrajectoryAnalysis/')

# Monocle installation
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
devtools::install_github('cole-trapnell-lab/leidenbase', force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3', force = TRUE)
remotes::install_github('satijalab/seurat-wrappers')

require('monocle3')
require('Seurat')
require('SeuratWrappers')
require('ggplot2')
require('dplyr')

myeloid <- readRDS('./data/myeloid.rds')
DefaultAssay(myeloid) <- 'integrated'

myeloid_cds <- as.cell_data_set(myeloid)
myeloid_cds <- cluster_cells(
  cds = myeloid_cds, 
  # reduction_method = 'PCA'
  # Uses UMAP coordinates
)
p1 <- plot_cells(myeloid_cds, 
                 show_trajectory_graph = FALSE) + 
  ggtitle(label = 'Monocle PCA clustering') +
  theme_bw()
p2 <- plot_cells(myeloid_cds, 
                 color_cells_by = "partition",
                 show_trajectory_graph = FALSE) + 
  ggtitle(label = 'Monocle partitions') +
  theme_bw()
p3 <- patchwork::wrap_plots(p1, p2)
ggsave(filename = 'results/revision_figures/TrajectoryAnalysis/Monocle_clustering.tiff',
       plot = p3, height = 3.5, width = 8.25, device = 'tiff')


tmp <- as.Seurat(myeloid_cds)
mac_cds <- tmp[,tmp$monocle3_partitions == 2]
mac_cds <- as.cell_data_set(mac_cds)
tmp_names <- names(mac_cds@clusters$UMAP$partitions)
mac_cds@clusters$UMAP$partitions <- factor(rep(1, length(mac_cds@clusters$UMAP$partitions)))
# mac_cds@clusters$UMAP$partitions <- as.character(mac_cds@clusters$UMAP$partitions)
names(mac_cds@clusters$UMAP$partitions) <- tmp_names
# colData(mac_cds)$monocle3_clusters <- as.character(colData(mac_cds)$monocle3_clusters)
mac_cds <- learn_graph(mac_cds)
plot_cells(mac_cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
mac_cds <- order_cells(mac_cds)
p1 <- plot_cells(mac_cds,
           color_cells_by = "pseudotime", 
           label_cell_groups = FALSE,
           label_leaves = FALSE, 
           label_branch_points = FALSE) +
  ggtitle(label = 'Macrophage pseudotime trajectory') + 
  theme_bw()

tmp <- rep(NA, times = nrow(myeloid@meta.data))
tmp[match(names(pseudotime(mac_cds)), rownames(myeloid@meta.data))] <- pseudotime(mac_cds)
# tmp <- pseudotime(mac_cds)[match(names(pseudotime(mac_cds)), rownames(myeloid@meta.data))]
myeloid$monocle_pseudotime <- tmp
p1 <- FeaturePlot(myeloid, features = 'monocle_pseudotime') + 
  scale_color_viridis_c(option = 'B') + 
  theme_bw() + 
  guides(color = guide_colorbar(frame.colour = 'black', 
                                ticks = FALSE,
                                title = 'Pseudotime'))
ggsave(filename = 'results/revision_figures/TrajectoryAnalysis/Monocle_pseudotime_umap.tiff',
       plot = p1, height = 3.5, width = 4.5, device = 'tiff')

# mac_cds_test <- graph_test(mac_cds, neighbor_graph = 'principal_graph')
# ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)



# Slingshot Trajectory Inference ------------------------------------------

BiocManager::install("slingshot")
install.packages('gam')

require('Seurat')
require('SeuratWrappers')
require('ggplot2')
require('dplyr')
require('slingshot')
require('monocle3')


myeloid <- readRDS('./data/myeloid.rds')
DefaultAssay(myeloid) <- 'integrated'

# monocle3 clustering of cells in PCA-space, PCA embeddings from Seurat
myeloid_cds <- as.cell_data_set(
  x = myeloid,
  assay = 'integrated'
)
myeloid_cds <- cluster_cells(
  cds = myeloid_cds, 
  # reduction_method = 'PCA'
  # Uses UMAP coordinates
)

# extract (partition & subcluster) macrophages, get fine clustering.
DefaultAssay(myeloid) <- 'integrated'
myeloid@meta.data$monocle3_clusters <- myeloid_cds@clusters$UMAP$clusters
myeloid@meta.data$monocle3_partitions <- myeloid_cds@clusters$UMAP$partitions
DimPlot(myeloid, group.by = 'monocle3_partitions')

extract_this <- 1
mac <- myeloid[, myeloid$monocle3_partitions == extract_this]
mac <- mac[,mac$myeloid_subcluster %in% c('Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid')]
mac <- FindNeighbors(mac, dims = 1:11, reduction = 'pca')
mac <- FindClusters(mac, resolution = 0.4)
Idents(mac) <- 'integrated_snn_res.0.4'
DimPlot(mac, label = TRUE)
mac_sce <- as.SingleCellExperiment(mac, assay = 'integrated')
mac_umap <- mac[['umap']]@cell.embeddings

# run slingshot on UMAP coordinates. MST drawn b/w fine clusters.
reducedDim(mac_sce, type = 'UMAP') <- mac_umap
mac_sce <- slingshot(data = mac_sce, 
                     clusterLabels = 'integrated_snn_res.0.4', 
                     reducedDim = 'UMAP',
                     start.clus = '4')
mac_sds <- SlingshotDataSet(mac_sce)


# extract curves/lineages coordinates and plot with subcluster
macrophage_cols <- c('Monocyte' = '#9a6324',
                     'Macrophage-A' = '#e6194b',
                     'Macrophage-B' = '#f58231',
                     'BA-Macrophage' = '#CCCC00',
                     'Dendritic' = '#808000',
                     'Div-Myeloid' = '#3cb44b')
curve1 <- data.frame(mac_sds@curves$curve1$s)[mac_sds@curves$curve1$ord,]
curve2 <- data.frame(mac_sds@curves$curve2$s)[mac_sds@curves$curve2$ord,]
curve3 <- data.frame(mac_sds@curves$curve3$s)[mac_sds@curves$curve3$ord,]

p1 <- DimPlot(object = mac,
        group.by = 'myeloid_subcluster',
        pt.size = 1) +
  scale_color_manual(values = macrophage_cols) +
  geom_path(data = curve1, mapping = aes(x = UMAP_1, y = UMAP_2), lwd = 1) +
  geom_text(data = data.frame(UMAP_1 = 4, UMAP_2 = 1.5),
            mapping = aes(label = 'Trajectory 1', x = UMAP_1, y = UMAP_2),
            size = 5) + 
  geom_path(data = curve2, mapping = aes(x = UMAP_1, y = UMAP_2), lwd = 1) + 
  geom_text(data = data.frame(UMAP_1 = 5.25, UMAP_2 = 5.5), 
            mapping = aes(label = 'Trajectory 2', x = UMAP_1, y = UMAP_2),
            size = 5) + 
  geom_path(data = curve3, mapping = aes(x = UMAP_1, y = UMAP_2), lwd = 1) +
  geom_text(data = data.frame(UMAP_1 = 3.5, UMAP_2 = -4.5),
            mapping = aes(label = 'Trajectory 3', x = UMAP_1, y = UMAP_2),
            size = 5) +
  ggtitle(label = 'Peripheral myeloid trajectory inference using Slingshot') +
  theme_bw() +
  theme(legend.text = element_text(size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6)))
p1
ggsave(filename = 'results/revision_figures/TrajectoryAnalysis/Slingshot-curves_umap.tiff',
       plot = p1, device = 'tiff', height = 5, width = 7)

# lineages + pseudotime
mac <- AddMetaData(
  object = mac, 
  metadata = slingPseudotime(x = mac_sds),
  col.name = c('sling_curve1', 'sling_curve2', 'sling_curve3')
)
p2 <- cbind(mac@meta.data[c('sling_curve1','sling_curve2','sling_curve3')], mac[['umap']]@cell.embeddings) %>%
  reshape2::melt(id.vars = c('UMAP_1','UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = value), size = 1) +
  scale_color_viridis_c(option = 'B',
                        limits = c(0,17),
                        breaks = c(0,17),
                        labels = c('  Early','Late  ')) + 
  facet_wrap(. ~ variable) +
  ggtitle(label = 'Slingshot pseudotime') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.text = element_text(size = 12),
        legend.title = element_text(angle = 90)) + 
  guides(color = guide_colorbar(title = 'Pseudotime',
                                frame.colour = 'black', 
                                ticks = FALSE,
                                title.position = 'left',
                                title.hjust = 0.5))
p2
ggsave(filename = 'results/revision_figures/TrajectoryAnalysis/slingshot_pseudotime_umap.tiff',
       plot = p2, device = 'tiff', height = 4.25, width = 12)

saveRDS(mac, file = 'results/revision_figures/TrajectoryAnalysis/mac.rds')
saveRDS(mac_sds, file = 'results/revision_figures/TrajectoryAnalysis/mac_sds.rds')


pt <- mac$sling_curve3 # mono to mac lineage
mac_expr <- mac[['RNA']]@data
fit_gam <- function(exp, t) {
  d <- data.frame(exp = exp, t = t)
  tmp <- gam::gam(exp ~ gam::lo(t), data = d)
  p <- summary(tmp)
  return(p)
}
test <- fit_gam(mac_expr['Apoe',], t = pt)

gam_pval <- apply(
  X = mac_expr,
  MARGIN = 1,
  FUN = nope
)


t <- sub_slingshot$slingPseudotime_1 #  lineage 1: KRT5-/KRT17

# Extract the gene expression matrix
Y <-assay(sub_slingshot)

# Fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 


install.packages('tradeSeq')
require('tradeSeq')
mac_sim <- fitGAM(counts = mac@assays$RNA@data)
mac_tmp <- myeloid[,colnames(mac)]





# Unbiased vascular LR analysis (Compute score)  -------------------------------

devtools::install_github(
  repo = "https://github.com/JamesChoi94/SingleCellTools.git", 
  INSTALL_opts = '--no-multiarch'
)

require('Seurat')
require('dplyr')
require('ggplot2')
require('SingleCellTools')

dir.create(path = './results/revision_figures/LR_results/')


sci <- readRDS(file = 'data/sci.rds')
DefaultAssay(sci) <- 'RNA'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL
gc()

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
# Reassign neurons
sci@meta.data[['subcluster']][sci$celltype == 'Neuron'] <- 'Neuron'

# label unassigned cells
sci@meta.data[['subcluster']][is.na(sci@meta.data[['subcluster']])] <- 'Unassigned'

# Clusters to merge
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('Ependymal-A','Ependymal-B'),
  to = c('Ependymal','Ependymal')
)
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('C1-Endothelial','C2-Endothelial','A-Endothelial','V-Endothelial'),
  to = c('Endothelial', 'Endothelial','Endothelial','Endothelial')
)
sci$subcluster <- plyr::mapvalues(
  x = sci$subcluster, from = 'Tip Cell', to = 'Tip-Cell'
)
sci <- sci[, !sci$subcluster %in% c('U-Vascular','Unassigned')]
Idents(sci) <- 'subcluster'
table(Idents(sci))
gc()


endo_lr <- StandardLR2(
  seurat.object = sci,
  lr.ref.path = './ref/fantom_PairsLigRec_mouse.csv',
  split.by = 'time',
  min.pct = 0.1,
  assay = 'RNA',
  slot = 'data',
  resample = 1000,
  adjust.pval = TRUE,
  receptor.cells = c('Endothelial','Tip-Cell'),
  split.vars = '1dpi',
  BPPARAM = BiocParallel::SnowParam(workers = 2, tasks = 21, progressbar = TRUE, type = 'SOCK')
)
write.csv(x = endo_lr, file = './results/revision_figures/LR_results/endothelial_1dpi_LR.rds', row.names = FALSE)


endo_lr <- StandardLR2(
  seurat.object = sci,
  lr.ref.path = './ref/fantom_PairsLigRec_mouse.csv',
  split.by = 'time',
  min.pct = 0.1,
  assay = 'RNA',
  slot = 'data',
  resample = 1000,
  adjust.pval = TRUE,
  receptor.cells = c('Endothelial','Tip-Cell'),
  split.vars = 'Uninjured',
  BPPARAM = BiocParallel::SnowParam(workers = 2, tasks = 21, progressbar = TRUE, type = 'SOCK')
)
write.csv(x = endo_lr, file = './results/revision_figures/LR_results/endothelial_uninj_LR.rds', row.names = FALSE)



# Unbiased vascular LR analysis (Plotting)  ----------------------------------

endo_uninj <- read.csv(file = './results/revision_figures/LR_results/endothelial_uninj_LR.rds')
endo_1dpi <- read.csv(file = './results/revision_figures/LR_results/endothelial_1dpi_LR.rds')
lr_ref <- read.csv(file = './ref/fantom_PairsLigRec_mouse.csv')

sci_subclusters <- c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid','H-Microglia','DAM-A','DAM-B','DAM-C','IFN-Myeloid','Endothelial','Tip-Cell','Pericyte','VSMC','Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC-A','OPC-B','Div-OPC','Pre-Oligo','Oligodendrocyte','Neuron')

endo_lr <- rbind(endo_uninj, endo_1dpi)

endo_lr$log_pval <- -log10(endo_lr$adj_pval)
endo_lr$log_pval[is.infinite(endo_lr$log_pva)] <- 3
endo_lr$Ligand_cell <- factor(endo_lr$Ligand_cell, levels = sci_subclusters)
endo_lr$Ligand_cell <- plyr::mapvalues(
  x = endo_lr$Ligand_cell,
  from = c('Neutrophil',
           'Monocyte',
           'Macrophage-A',
           'Macrophage-B',
           'BA-Macrophage',
           'Dendritic',
           'Div-Myeloid',
           'H-Microglia',
           'DAM-A',
           'DAM-B',
           'DAM-C',
           'IFN-Myeloid'),
  to = c('Neutrophil',
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
endo_lr$Receptor_cell <- factor(endo_lr$Receptor_cell, levels = sci_subclusters)
endo_lr$Receptor_cell <- plyr::mapvalues(
  x = endo_lr$Receptor_cell,
  from = c('Neutrophil',
           'Monocyte',
           'Macrophage-A',
           'Macrophage-B',
           'BA-Macrophage',
           'Dendritic',
           'Div-Myeloid',
           'H-Microglia',
           'DAM-A',
           'DAM-B',
           'DAM-C',
           'IFN-Myeloid'),
  to = c('Neutrophil',
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
endo_lr$compartment <- plyr::mapvalues(
  x = endo_lr$Ligand_cell,
  from = c('Neutrophil', 'Monocyte', 'Chemotaxis-Inducing Mac', 'Inflammatory Mac', 'Border-Associated Mac', 'Dendritic', 'Dividing Myeloid', 'Homeostatic Microglia', 'Inflammatory Microglia', 'Dividing Microglia', 'Migrating Microglia', 'Interferon Myeloid', 'Endothelial', 'Tip-Cell', 'Pericyte', 'VSMC', 'Fibroblast', 'Ependymal', 'Astroependymal', 'Astrocyte', 'OPC-A', 'OPC-B', 'Div-OPC', 'Pre-Oligo', 'Oligodendrocyte', 'Neuron'),
  to = c('Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Vascular','Vascular','Vascular','Vascular','Vascular','Macroglia','Macroglia','Macroglia','Macroglia','Macroglia','Macroglia','Macroglia','Macroglia','Neuron')
)
endo_lr$split.by <- factor(endo_lr$split.by, levels = c('Uninjured','1dpi'))

remove_pairs <- c('Calr_Scarf1','Gpi1_Amfr','Ltb_Cd40','Nampt_Insr','Npnt_Itgb1','Pf4_Ldlr','Pf4_Procr','Pf4_Thbd','Ptdss1_Scarb1','Vtn_Pvr','Vtn_Plaur','Rtn4_Tnfrsf19','Ptgs2_Cav1','Ngf_Sort1','Ncam1_Fgfr1','Kng2_Plaur','Kng2_Cd93','Hras_Cav1','Hbegf_Cd9','Hbegf_Cd82','Gpi1_Amfr','Gpc3_Igf1r','Gpc3_Cd81','Gnai2_Igf1r','Gnai2_Cav1','Efnb2_Pecam1','Efemp2_Plscr4','Calr_Itga2b','Calca_Ramp2','Calca_Calr1','Bgn_Tlr4','Arg1_Insr','App_Cav1','Anxa1_Dysf','Adm_Calcrl','Nampy_Insr','Lamb2_Rpsa','Lama2_Rpsa','Ccn2_Lrp6')
remove_cells <- c('Border-Associated Mac','Dendritic','Dividing Myeloid','Interferon Myeloid')

endo_lr$Ligand <- plyr::mapvalues(
  x = endo_lr$Ligand,
  from =  c('Tek','Flt1','Kdr','Flt4','Pgf'),
  to = c('Tie2','Vegfr1','Vegfr2','Vegfr3','Plgf')
)
endo_lr$Receptor <- plyr::mapvalues(
  x = endo_lr$Receptor,
  from =  c('Tek','Flt1','Kdr','Flt4','Pgf'),
  to = c('Tie2','Vegfr1','Vegfr2','Vegfr3','Plgf')
)
endo_lr$Pair_name <- paste(endo_lr$Ligand, endo_lr$Receptor, sep = '_')

lr_ref$Ligand.ApprovedSymbol <- plyr::mapvalues(
  x = lr_ref$Ligand.ApprovedSymbol,
  from =  c('Tek','Flt1','Kdr','Flt4','Pgf'),
  to = c('Tie2','Vegfr1','Vegfr2','Vegfr3','Plgf')
)
lr_ref$Receptor.ApprovedSymbol <- plyr::mapvalues(
  x = lr_ref$Receptor.ApprovedSymbol,
  from =  c('Tek','Flt1','Kdr','Flt4','Pgf'),
  to = c('Tie2','Vegfr1','Vegfr2','Vegfr3','Plgf')
)
lr_ref$Pair.Name <- paste(lr_ref$Ligand.ApprovedSymbol,
                          lr_ref$Receptor.ApprovedSymbol,
                          sep = '_')


filter_lr <- function(x, comp, time, pct = 1) {
  tmp <- x %>%
    filter(Pair_name %in% lr_ref$Pair.Name[lr_ref$Pair.Evidence == 'literature supported' & lr_ref$Pair.Source == 'known']) %>%
    filter(!Pair_name %in% remove_pairs) %>%
    filter(!Ligand_cell %in% remove_cells) %>%
    filter(Ligand_pct >= 0.1 & Receptor_pct >= 0.1) %>%
    filter(Ligand_cell_pct > pct & Receptor_cell_pct > pct) %>%
    filter(adj_pval < 0.05) %>%
    filter(compartment %in% comp) %>%
    filter(split.by %in% time)
  return(tmp)
}

plot_lr <- function(x, split = FALSE, max_score = NULL) {
  max_score <- max(x[['Score']], max_score)
  tmp <- ggplot(x) + 
    geom_point(mapping = aes(x = Ligand_cell, y = Pair_name, size = log_pval, fill = Score),
               color = 'black', pch = 21)
  if (split) {
    tmp <- tmp + facet_wrap(. ~ split.by + Receptor_cell,
                            strip.position = 'bottom', 
                            ncol = 4)
  }
  tmp <- tmp + 
    scale_x_discrete(position = 'top') +
    scale_size_continuous(limits = c(0,3)) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, 'Spectral')),
                         limits = c(0, max_score)) +
    ylab(label = 'Ligand_Receptor Pair') +
    xlab(label = 'Ligand Cell') +
    theme_bw() +
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0, size = 12),
          strip.text = element_text(size = 12),
          strip.background = element_rect(fill = NA, color = NA),
          axis.text.y.left = element_text(size = 11, face = 'italic')) +
    guides(fill = guide_colorbar(frame.colour = 'black', 
                                  ticks.colour = 'black'),
           size = guide_legend(title = '-Log10(p-value)',
                               override.aes = list(color = 'black',
                                                   fill = 'black')))
}

endo_mye_1dpi <- filter_lr(endo_lr, comp = 'Myeloid', time = '1dpi') %>% 
  plot_lr(split = TRUE, max_score = 2.7)
ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_myeloid.tiff',
       plot = endo_mye_1dpi, device = 'tiff', height = 9, width = 6)
# ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_myeloid.png',
#        plot = endo_mye_1dpi, device = 'png', height = 9, width = 6)
endo_vasc_1dpi <- filter_lr(endo_lr, comp = 'Vascular', time = '1dpi') %>% 
  plot_lr(split = TRUE, max_score = 2.7)
ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_vascular.tiff',
       plot = endo_vasc_1dpi, device = 'tiff', height = 16, width = 5)
# ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_vascular.png',
#        plot = endo_vasc_1dpi, device = 'png', height = 16, width = 5)
endo_macro_1dpi <- filter_lr(endo_lr, comp = 'Macroglia', time = '1dpi') %>% 
  plot_lr(split = TRUE, max_score = 2.7)
ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_macroglia.tiff',
       plot = endo_macro_1dpi, device = 'tiff', height = 14, width = 6)
# ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_macroglia.png',
#        plot = endo_macro_1dpi, device = 'png', height = 14, width = 6)


endo_vasc_all <- filter_lr(endo_lr, comp = 'Vascular', time = c('1dpi','Uninjured')) %>% 
  plot_lr(split = TRUE)
ggsave(filename = './results/revision_figures/LR_results/endothelial_uninj_1dpi_vascular.tiff',
       plot = endo_vasc_all, device = 'tiff', height = 16, width = 7)
ggsave(filename = './results/revision_figures/LR_results/endothelial_uninj_1dpi_vascular.png',
       plot = endo_vasc_all, device = 'png', height = 16, width = 7)


endo_vasc_all <- filter_lr(endo_lr, comp = 'Vascular', time = c('1dpi','Uninjured'), pct = 2) %>% 
  plot_lr(split = TRUE)
ggsave(filename = './results/revision_figures/LR_results/endothelial_uninj_1dpi_vascular.jpeg',
       plot = endo_vasc_all, device = 'jpeg', height = 16, width = 5.75)



# IL6 family LR analysis (Compute score)  ------------------------------------

devtools::install_github(
  repo = "https://github.com/JamesChoi94/SingleCellTools.git", 
  INSTALL_opts = '--no-multiarch'
)

require('Seurat')
require('dplyr')
require('ggplot2')
require('SingleCellTools')

dir.create(path = './results/revision_figures/LR_results/')


sci <- readRDS(file = 'data/sci.rds')
DefaultAssay(sci) <- 'RNA'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL
gc()

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
# Reassign neurons
sci@meta.data[['subcluster']][sci$celltype == 'Neuron'] <- 'Neuron'

# label unassigned cells
sci@meta.data[['subcluster']][is.na(sci@meta.data[['subcluster']])] <- 'Unassigned'

# Clusters to merge
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('Ependymal-A','Ependymal-B'),
  to = c('Ependymal','Ependymal')
)
sci$subcluster <- plyr::mapvalues(
  x = sci$subcluster, from = 'Tip Cell', to = 'Tip-Cell'
)
sci <- sci[, !sci$subcluster %in% c('U-Vascular','Unassigned')]
Idents(sci) <- 'subcluster'
table(Idents(sci))
gc()

il6.genes <- list('Ligands' = c('Il6','Cntf','Lif','Ctf1','Ctf2','Osm','Il11','Il31','Clcf1','Il27'),
                  'Receptors' = c('Il6st','Il6ra','Lifr','Cntfr','Osmr','Il11ra1'))
il6_ligs <- il6.genes$Ligands
il6_ligs <- c('Osm','Il6','Clcf1')
il6_recs <- il6.genes$Receptors
il6_genes <- c(il6_ligs, il6_recs)

lr_ref <- lr_ref[lr_ref$Ligand.ApprovedSymbol %in% il6_genes & lr_ref$Receptor.ApprovedSymbol %in% il6_genes,]

il6_lr <- StandardLR2(
  seurat.object = sci,
  # lr.ref.path = './ref/fantom_PairsLigRec_mouse.csv',
  lr.ref = lr_ref,
  split.by = 'time',
  min.pct = 0.1,
  assay = 'RNA',
  slot = 'data',
  resample = 1000,
  adjust.pval = TRUE,
  ligand.cells = c('Monocyte','Macrophage-A','Macrophage-B','H-Microglia','DAM-A','DAM-B','DAM-C','Fibroblast','Astroependymal','OPC-A','OPC-B'),
  receptor.cells = c('Pericyte','Fibroblast','Astroependymal','Astrocyte','OPC-A','OPC-B','Oligodendrocyte'),
  # split.vars = '1dpi',
  BPPARAM = BiocParallel::SnowParam(workers = 2, tasks = 21, progressbar = TRUE, type = 'SOCK')
)
write.csv(x = il6_lr, file = './results/revision_figures/LR_results/macroglia_il6_LR.csv', row.names = FALSE)



# Macroglia IL6 LR analysis (Plotting)  ----------------------------------

il6_lr <- read.csv(file = './results/revision_figures/LR_results/macroglia_il6_LR.rds')

sci_subclusters <- c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid','H-Microglia','DAM-A','DAM-B','DAM-C','IFN-Myeloid','Endothelial','Tip-Cell','Pericyte','VSMC','Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC-A','OPC-B','Div-OPC','Pre-Oligo','Oligodendrocyte','Neuron')

il6_lr$log_pval <- -log10(il6_lr$adj_pval)
il6_lr$log_pval[is.infinite(il6_lr$log_pva)] <- 3
il6_lr$Ligand_cell <- factor(il6_lr$Ligand_cell, levels = sci_subclusters)
il6_lr$Ligand_cell <- plyr::mapvalues(
  x = il6_lr$Ligand_cell,
  from = c('Neutrophil',
           'Monocyte',
           'Macrophage-A',
           'Macrophage-B',
           'BA-Macrophage',
           'Dendritic',
           'Div-Myeloid',
           'H-Microglia',
           'DAM-A',
           'DAM-B',
           'DAM-C',
           'IFN-Myeloid'),
  to = c('Neutrophil',
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
il6_lr$Receptor_cell <- factor(il6_lr$Receptor_cell, levels = sci_subclusters)
il6_lr$Receptor_cell <- plyr::mapvalues(
  x = il6_lr$Receptor_cell,
  from = c('Neutrophil',
           'Monocyte',
           'Macrophage-A',
           'Macrophage-B',
           'BA-Macrophage',
           'Dendritic',
           'Div-Myeloid',
           'H-Microglia',
           'DAM-A',
           'DAM-B',
           'DAM-C',
           'IFN-Myeloid'),
  to = c('Neutrophil',
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

il6_lr %>%
  mutate(split.by = factor(split.by, levels = c('Uninjured','1dpi','3dpi','7dpi'))) %>%
  filter(Ligand_cell_pct > 1 & Receptor_cell_pct > 1) %>%
  filter(pval < 0.05) %>%
  filter(!Receptor %in% 'Lifr') %>%
  ggplot() +
  geom_point(mapping = aes(x = Ligand_cell, y = Pair_name, size = log_pval, fill = Score),
             color = 'black',
             pch = 21) +
  facet_grid(split.by ~ Receptor_cell, drop = TRUE, switch = 'x') +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral'))) +
  scale_radius(limits = c(0,NA), range = c(1,6)) +
  scale_y_discrete(position = 'left') +
  scale_x_discrete(position = 'top') +
  xlab(label = 'Ligand cell') +
  ylab(label = 'Ligand_Receptor Pair')


filter_lr <- function(x, time, pct = 1) {
  tmp <- x %>%
    filter(Ligand_pct >= 0.1 & Receptor_pct >= 0.1) %>%
    filter(Ligand_cell_pct > pct & Receptor_cell_pct > pct) %>%
    filter(adj_pval < 0.05) %>%
    filter(split.by %in% time)
  return(tmp)
}

plot_lr <- function(x, split = FALSE, max_score = NULL) {
  max_score <- max(x[['Score']], max_score)
  tmp <- ggplot(x) + 
    geom_point(mapping = aes(x = Ligand_cell, y = Pair_name, size = log_pval, fill = Score),
               color = 'black', pch = 21)
  if (split) {
    tmp <- tmp + facet_wrap(. ~ split.by + Receptor_cell,
                            strip.position = 'bottom', 
                            ncol = 4)
  }
  tmp <- tmp + 
    scale_x_discrete(position = 'top') +
    scale_size_continuous(limits = c(0,3)) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, 'Spectral')),
                         limits = c(0, max_score)) +
    ylab(label = 'Ligand_Receptor Pair') +
    xlab(label = 'Ligand Cell') +
    theme_bw() +
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0),
          strip.text = element_text(size = 12),
          strip.background = element_rect(fill = NA, color = NA)) +
    guides(fill = guide_colorbar(frame.colour = 'black', 
                                 ticks.colour = 'black'),
           size = guide_legend(title = '-Log10(p-value)',
                               override.aes = list(color = 'black',
                                                   fill = 'black')))
}

