
require('Seurat')
require('dplyr')
require('ggplot2')

results_out <- 'results/revision_figures/'
dir.create(path = results_out)

sci <- readRDS(file = 'data/sci.rds')
myeloid <- readRDS(file = 'data/myeloid.rds')

DefaultAssay(sci) <- 'RNA'
DefaultAssay(myeloid) <- 'RNA'




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
mg_go_results <- list('H-Microglia' = h_mg, 
                      'Inflammatory-Microglia' = a_mg, 
                      'Dividing-Microglia' = b_mg,
                      'Migrating-Microglia' = c_mg)
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
  levels = c('H-Microglia','Inflammatory-Microglia','Dividing-Microglia','Migrating-Microglia')
)

mg_go_plot <- mg_go_results %>%
  group_by(Cluster) %>%
  top_n(n = 15, wt = log_pval) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ Cluster, scales = 'free_y', drop = TRUE, ncol = 2) +
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



# cellcycle regression ----------------------------------------------------

myeloid <- readRDS('data/myeloid.rds')




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
       plot = umap_sample, device = 'tiff', height = 9, width = 11)


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
       plot = umap_sample, device = 'tiff', height = 7, width = 8.5)


# Vascular UMAP by sample ------------------------------------------------------

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
       plot = umap_sample, device = 'tiff', height = 7, width = 8.5)





# Macrophage WGCNA --------------------------------------------------------


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
  corFnc = stats::cor,
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
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

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

MEDissThres = 0.25
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

MEs0 = moduleEigengenes(mac_mat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, macrophage@meta.data, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)




# -------------------------------------------------------------------------


