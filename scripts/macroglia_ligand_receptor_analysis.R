

###### Ligand-receptor interaction analysis of IL6-mediate gliosis #######


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
results_in <- './results/ligand_receptor_analysis/'
results_out <- './results/macroglia_ligand_receptor_analysis/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)


# data
sci <- readRDS(file = './data/sci.rds')
lr_results <- read.table(file = paste0(results_in, 'sci_cluster_LigandReceptorResults_manuscript_interactions.tsv'))
lr_ref <- read.csv(file = paste0(ref_in, 'fantom_PairsLigRec_mouse.csv'))




# Data setup --------------------------------------------------------------


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


# Clusters to merge
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('OPC-A','OPC-B'),
  to = c('OPC','OPC')
)
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('Ependymal-A','Ependymal-B'),
  to = c('Ependymal','Ependymal')
)
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('C1-Endothelial','C2-Endothelial'),
  to = c('C-Endothelial', 'C-Endothelial')
)

# Convert p-values into numbers
lr_results[['pval']] <- as.numeric(lr_results[['pval']])
lr_results[['pval']][lr_results[['Ligand_pct']] < 0.1] <- NA
lr_results[['pval']][lr_results[['Receptor_pct']] < 0.1] <- NA
lr_results[['log_pval']] <- -log10(lr_results[['pval']])
lr_results[['log_pval']][is.infinite(lr_results[['log_pval']])] <- -log10(1/1000) # 1000 permutations
lr_results[['time']] <- factor(lr_results[['split_by']], levels = levels(sci@meta.data[['time']]))
lr_results[['Ligand_cell']] <- factor(
  x = lr_results[['Ligand_cell']],
  levels = c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid','Microglia-A','Microglia-B','Microglia-C','Div-Microglia','A-Endothelial', 'C-Endothelial','V-Endothelial','Tip Cell','Pericyte','VSMC', 'Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC','Div-OPC','Pre-Oligo', 'Oligodendrocyte')
)
lr_results[['Receptor_cell']] <- factor(
  x = lr_results[['Receptor_cell']],
  levels = c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid','Microglia-A','Microglia-B','Microglia-C','Div-Microglia','A-Endothelial', 'C-Endothelial','V-Endothelial','Tip Cell','Pericyte','VSMC', 'Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC','Div-OPC','Pre-Oligo', 'Oligodendrocyte')
)





# IL6 signaling (Ligand-receptor plot) ------------------------------

il6.genes <- list('Ligands' = c('Il6','Cntf','Lif','Ctf1','Ctf2','Osm','Il11','Il31','Clcf1','Il27'),
                  'Receptors' = c('Il6st','Il6ra','Lifr','Cntfr','Osmr','Il11ra1'))
il6_ligs <- il6.genes$Ligands
il6_ligs <- c('Osm','Il6','Clcf1')
il6_recs <- il6.genes$Receptors
il6_genes <- c(il6_ligs, il6_recs)

# Take "known" interactions only
known_lr <- lr_ref$Pair.Name[lr_ref$Pair.Source == 'known']

ligand_celltypes <- c('Monocyte','Macrophage-A','Macrophage-B','H-Microglia','DAM-A','DAM-B','DAM-C','Fibroblast','Astroependymal','OPC')
receptor_celltypes <- c('Astroependymal','Astrocyte','Oligodendrocyte','OPC','Pericyte','Fibroblast')


# Reproduce manuscript figure
il6_data <- lr_results %>%
  # rename microglia to new cluster names
  mutate('Receptor_cell' = plyr::mapvalues(
    x = Receptor_cell, 
    from = c('Microglia-A','Microglia-B','Microglia-C','Div-Microglia'),
    to = c('H-Microglia','DAM-A','DAM-C','DAM-B')
  )) %>%
  mutate('Ligand_cell' = plyr::mapvalues(
    x = Ligand_cell, 
    from = c('Microglia-A','Microglia-B','Microglia-C','Div-Microglia'),
    to = c('H-Microglia','DAM-A','DAM-C','DAM-B')
  )) %>%
  filter(Ligand_cell %in% ligand_celltypes) %>%
  filter(Receptor_cell %in% receptor_celltypes) %>%
  filter(Ligand %in% il6_genes & Receptor %in% il6_genes) %>%
  filter(Ligand_cell_pct >= 2.5 & Receptor_cell_pct >= 2.5) %>%
  filter(Pair_name %in% known_lr) %>%
  filter(!Pair_name %in% c('Osm_Lifr')) %>%
  mutate(Receptor_cell = factor(Receptor_cell, 
                                levels = c('Pericyte','Fibroblast','Astroependymal','Astrocyte','OPC','Oligodendrocyte')),
         Ligand_cell = factor(Ligand_cell, 
                              levels = c('Monocyte','Macrophage-A','Macrophage-B','H-Microglia','DAM-A','DAM-B','DAM-C','Fibroblast','Astroependymal','Astrocyte','OPC'))) %>%
  mutate(Receptor_cell = plyr::mapvalues(x = Receptor_cell, from = 'Astroependymal', to = 'Astroepen.')) %>%
  mutate(Ligand_cell = plyr::mapvalues(x = Ligand_cell, from = 'Astroependymal', to = 'Astroepen.'))
max_score <- ceiling(max(il6_data$Score)*10)/10
min_score <- floor(min(il6_data$Score)*10)/10
il6_plot <- il6_data %>%
  ggplot() +
  geom_point(mapping = aes(x = Receptor_cell, y = Pair_name, size = log_pval, fill = Score),
             color = 'black',
             pch = 21) +
  facet_grid(time ~ Ligand_cell, drop = TRUE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral')),
                       breaks = c(min_score, 0.5, max_score),
                       limits = c(min_score, max_score),
                       labels = c(min_score,0.5, max_score)) +
  scale_radius(limits = c(0,NA), range = c(1,6)) +
  scale_y_discrete(position = 'left') +
  scale_x_discrete(position = 'bottom') +
  xlab(label = 'Ligand cell') +
  ylab(label = 'Ligand_Receptor Pair') +
  theme(strip.text = element_text(size = 12, color = 'black', face = 'bold'),
        strip.background = element_rect(fill = NA, color = NA),
        axis.title.x.bottom = element_text(size = 12, color = 'black', face = 'bold'),
        axis.title.y.left = element_text(size = 12, color = 'black', face = 'bold'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black', face = 'bold'),
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70'),
        legend.position = 'top',
        legend.direction = 'horizontal') +
  guides(fill = guide_colorbar(title = 'Score',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               frame.colour = 'black',
                               ticks.colour = 'black'),
         size = guide_legend(title = '-log10(p-value)', 
                             override.aes = list(fill = 'black'))); il6_plot
ggsave(filename = paste0(results_out, 'IL6_LigandReceptorplot.tiff'),
       plot = il6_plot, device = 'tiff', height = 7.75, width = 16)


# Flipped facet orientation
il6_plot <- il6_data %>%
  ggplot() +
  geom_point(mapping = aes(x = Ligand_cell, y = Pair_name, size = log_pval, fill = Score),
             color = 'black',
             pch = 21) +
  facet_grid(time ~ Receptor_cell, switch = 'x', drop = TRUE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral')),
                       breaks = c(min_score, 0.5, max_score),
                       limits = c(min_score, max_score),
                       labels = c(min_score,0.5, max_score)) +
  scale_radius(limits = c(0,NA), range = c(1,6)) +
  scale_y_discrete(position = 'left') +
  scale_x_discrete(position = 'top') +
  xlab(label = 'Ligand cell') +
  ylab(label = 'Ligand_Receptor Pair') +
  theme(strip.text = element_text(size = 12, color = 'black', face = 'bold'),
        strip.background = element_rect(fill = NA, color = NA),
        axis.title.x.top = element_text(size = 12, color = 'black', face = 'bold'),
        axis.title.y.left = element_text(size = 12, color = 'black', face = 'bold'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black', angle = 90, face = 'bold', hjust = 0.5),
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70'),
        legend.position = 'right',
        legend.direction = 'vertical') +
  guides(fill = guide_colorbar(title = 'Score',
                               title.position = 'left',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               frame.colour = 'black',
                               ticks.colour = 'black'),
         size = guide_legend(title = '-log10(p-value)', 
                             title.position = 'left',
                             override.aes = list(fill = 'black'))); il6_plot
ggsave(filename = paste0(results_out, 'IL6_LigandReceptorplot_flipped.tiff'),
       plot = il6_plot, device = 'tiff', height = 7.25, width = 16.5)






# Validates scores with expression dot plot ------------------------------------


# IL6 genes
il6_genes <- list('Ligands' = c('Il6','Cntf','Lif','Ctf1','Ctf2','Osm','Il11','Il31','Clcf1','Il27'),
                  'Receptors' = c('Il6st','Il6ra','Lifr','Cntfr','Osmr','Il11ra1'))
# il6_ligands <- il6_genes$Ligands
il6_ligands <- c('Osm','Clcf1','Il6')
il6_receptors <- il6_genes$Receptors


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

# Subset subcluster identities of interest
select_sci_subclusters <- c('Neutrophil',
                            'Monocyte',
                            'Macrophage-A',
                            'Macrophage-B',
                            'Div-Myeloid',
                            'H-Microglia',
                            'DAM-A',
                            'DAM-B',
                            'DAM-C',
                            'BA-Macrophage',
                            'A-Endothelial',
                            'C1-Endothelial',
                            'C2-Endothelial',
                            'V-Endothelial',
                            'Tip Cell',
                            'Fibroblast',
                            'Pericyte',
                            'VSMC',
                            'Ependymal-A',
                            'Ependymal-B',
                            'Astroependymal',
                            'Astrocyte',
                            'OPC-A',
                            'OPC-B',
                            'Div-OPC',
                            'Oligodendrocyte')
Idents(sci) <- 'subcluster'
select_sci <- subset(sci, idents = select_sci_subclusters)
Idents(select_sci) <- 'subcluster'


# Calculate values for ligand dot plot
DefaultAssay(select_sci) <- 'RNA'
avg_exp <- ScaleData(select_sci[['RNA']]@data, features = c(il6_ligands))
avg_exp <- cbind(t(avg_exp), select_sci@meta.data[,c('subcluster','time')]) %>%
  reshape2::melt(id.vars = c('subcluster','time')) %>%
  group_by(subcluster, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- select_sci[['RNA']]@counts[c(il6_ligands),]
pct_exp <- cbind(Matrix::t(pct_exp), select_sci@meta.data[,c('subcluster','time')]) %>%
  reshape2::melt(id.vars = c('subcluster','time')) %>%
  group_by(subcluster, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)

# Filter subclusters with low presence in a time point (prevent skew of scaled
# expression)
low_prop <- prop.table(table(select_sci$subcluster, select_sci$time), margin = 1) < 0.015
for (i in 1:nrow(avg_exp)) {
  if (low_prop[avg_exp$subcluster[i], avg_exp$time[i]]) {
    avg_exp$avg.exp[i] <- NA
  }
}
for (i in 1:nrow(pct_exp)) {
  if (low_prop[pct_exp$subcluster[i], pct_exp$time[i]]) {
    pct_exp$pct.exp[i] <- NA
  }
}

max_expr <- ceiling(max(avg_exp$avg.exp, na.rm = TRUE)*10)/10
min_expr <- floor(min(avg_exp$avg.exp, na.rm = TRUE)*10)/10
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min_expr, 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))

# il6 ligand dot plot
il6_ligand_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(subcluster = factor(subcluster, levels = select_sci_subclusters)) %>%
  ggplot(mapping = aes(x = subcluster, y = variable)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), pch = 21, color = 'black') +
  geom_vline(xintercept = 10.5, linetype = 'dashed', color = 'black', size = 1) +
  geom_vline(xintercept = 18.5, linetype = 'dashed', color = 'black', size = 1) +
  facet_grid(time ~ .) +
  scale_size(range = c(0,8), limits = c(0,100), breaks = seq(25,100,25)) +
  scale_fill_gradientn(
    colors = myColors,
    values = scales::rescale(x = myBreaks, to = c(0,1)),
    limits = c(min_expr, max_expr), 
    na.value = myColors[100],
    breaks = c(0, max_expr)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        strip.background = element_rect(fill = NA, colour = 'black'),
        strip.text.y = element_text(angle = 270, hjust = 0.5, face = 'bold', size = 12, color = 'black'),
        strip.text.x = element_text(size = 12, face = 'bold', color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 12, color = 'black', hjust = 0.5),
        legend.title = element_text(size = 12, color = 'black'),
        legend.box.margin = margin(0, 0, 10, 0),
        legend.box = 'horizontal',
        legend.position = 'bottom',
        panel.border = element_rect(fill = NA, size = 1, color = 'grey60'),
        panel.background = element_rect(fill = NA)) +
  guides(fill = guide_colorbar(title = 'Scaled Expression', 
                               barwidth = 4,
                               barheight = 1,
                               frame.colour = 'black', 
                               frame.linewidth = 1.25,
                               ticks.colour = 'black',
                               ticks.linewidth = 1.25,
                               title.position = 'top'), 
         size = guide_legend(title = '% Expression', 
                             override.aes = list(fill = 'black'), 
                             title.position = 'top')); il6_ligand_dotplot
ggsave(filename = paste0(results_out, 'IL6_ligand_dotplot.tiff'),
       plot = il6_ligand_dotplot, device = 'tiff', height = 6, width = 6.75)




# Calculate values for receptor dot plot
DefaultAssay(select_sci) <- 'RNA'
avg_exp <- ScaleData(select_sci[['RNA']]@data, features = c(il6_receptors), scale.max = 3.5)
avg_exp <- cbind(t(avg_exp), select_sci@meta.data[,c('subcluster','time')]) %>%
  reshape2::melt(id.vars = c('subcluster','time')) %>%
  group_by(subcluster, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- select_sci[['RNA']]@counts[c(il6_receptors),]
pct_exp <- cbind(Matrix::t(pct_exp), select_sci@meta.data[,c('subcluster','time')]) %>%
  reshape2::melt(id.vars = c('subcluster','time')) %>%
  group_by(subcluster, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)

# Filter subclusters with low presence in a time point (prevent skew of scaled
# expression)
low_prop <- prop.table(table(select_sci$subcluster, select_sci$time), margin = 1) < 0.015
for (i in 1:nrow(avg_exp)) {
  if (low_prop[avg_exp$subcluster[i], avg_exp$time[i]]) {
    avg_exp$avg.exp[i] <- NA
  }
}
for (i in 1:nrow(pct_exp)) {
  if (low_prop[pct_exp$subcluster[i], pct_exp$time[i]]) {
    pct_exp$pct.exp[i] <- NA
  }
}

max_expr <- ceiling(max(avg_exp$avg.exp, na.rm = TRUE)*10)/10
min_expr <- floor(min(avg_exp$avg.exp, na.rm = TRUE)*10)/10
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min_expr, 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))

# il6 receptor dot plot
il6_receptor_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(subcluster = factor(subcluster, levels = select_sci_subclusters)) %>%
  ggplot(mapping = aes(x = subcluster, y = variable)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), pch = 21, color = 'black') +
  geom_vline(xintercept = 10.5, linetype = 'dashed', color = 'black', size = 1) +
  geom_vline(xintercept = 18.5, linetype = 'dashed', color = 'black', size = 1) +
  facet_grid(time ~ .) +
  scale_size(range = c(0,8), limits = c(0,100), breaks = seq(25,100,25)) +
  scale_fill_gradientn(
    colors = myColors,
    values = scales::rescale(x = myBreaks, to = c(0,1)),
    limits = c(min_expr, max_expr), 
    na.value = myColors[100],
    breaks = c(0, max_expr)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        strip.background = element_rect(fill = NA, colour = 'black'),
        strip.text.y = element_text(angle = 270, hjust = 0.5, face = 'bold', size = 12, color = 'black'),
        strip.text.x = element_text(size = 12, face = 'bold', color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 12, color = 'black', hjust = 0.5),
        legend.title = element_text(size = 12, color = 'black', angle = 90, hjust = 0.5),
        legend.box.margin = margin(0, 0, 10, 0),
        legend.box = 'vertical',
        legend.position = 'right',
        panel.border = element_rect(fill = NA, size = 1, color = 'grey60'),
        panel.background = element_rect(fill = NA)) +
  guides(fill = guide_colorbar(title = 'Scaled\nExpression', 
                               barwidth = 1,
                               barheight = 4,
                               frame.colour = 'black', 
                               frame.linewidth = 1.25,
                               ticks.colour = 'black',
                               ticks.linewidth = 1.25,
                               title.position = 'left'), 
         size = guide_legend(title = '% Expression', 
                             override.aes = list(fill = 'black'), 
                             title.position = 'left')); il6_receptor_dotplot
ggsave(filename = paste0(results_out, 'IL6_receptor_dotplot.tiff'),
       plot = il6_receptor_dotplot, device = 'tiff', height = 6, width = 8.75)
