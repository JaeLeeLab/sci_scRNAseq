
###### Astrocyte vs Fibroblast ligand-receptor interactions comparison ######


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
results_out <- './results/AstrocyteFibroblast_ligand_receptor_analysis/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)

source('./scripts/ligand_receptor_analysis_functions.R')

sci <- readRDS(file = './data/sci.rds')





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


# Interaction analysis at subcluster level --------------------------------


# Subset subcluster identities of interest
select_sci_subclusters <- c('Neutrophil',
                            'Monocyte',
                            'Macrophage-A',
                            'Macrophage-B',
                            'Div-Myeloid',
                            'Microglia-A',
                            'Microglia-B',
                            'Microglia-C',
                            'Div-Microglia',
                            'BA-Macrophage',
                            'A-Endothelial',
                            'C-Endothelial',
                            'V-Endothelial',
                            'Tip Cell',
                            'Fibroblast',
                            'Pericyte',
                            'VSMC',
                            'Ependymal',
                            'Astroependymal',
                            'Astrocyte',
                            'OPC',
                            'Div-OPC',
                            'Oligodendrocyte')
Idents(sci) <- 'subcluster'
select_sci <- subset(sci, idents = select_sci_subclusters)
rm(sci); gc()

# Make lightweight data structures
DefaultAssay(select_sci) <- 'RNA'
select_sci[['SCT']] <- NULL
select_sci[['integrated']] <- NULL

# Set identities and calculate setup values (avg_exp, pct_exp, cell counts)
Idents(select_sci) <- 'subcluster'
astrofibro_setup <- setupLR(seurat_object = select_sci,
                            ref_path = './ref/fantom_PairsLigRec_mouse.csv',
                            split_by = 'time',
                            assay = 'RNA',
                            slot = 'data',
                            source_filter = 'known')
gc()

# Permutation test
astrofibro_result <- calculateLR(setup = astrofibro_setup,
                                 resample = 1000,
                                 adjust_pval = TRUE,
                                 receptor_subsets = c('Astrocyte','Fibroblast'),
                                 ligand_subsets = c('Macrophage-A','Macrophage-B'))

# Save results
write.table(x = astrofibro_result,
            file = paste0(results_out, 'sci_cluster_LigandReceptorResults_AstrocyteFibroblast_macrophage.tsv'),
            sep = '\t',
            col.names = TRUE)








# Plotting results --------------------------------------------------------

astrofibro_result <- read.table(file = paste0(results_out, 'sci_cluster_LigandReceptorResults_AstrocyteFibroblast_macrophage.tsv'))


# Interactions previously vetted by Pantelis
lr_vetted <- c('Sema4d_Plxnb1','Trf_Tfrc','Apoe_Ldlr','Apoe_Lrp8','F7_F3','Fn1_Itga5','Fn1_Itgb1','Fn1_Itgb8','Fn1_Sdc2','Gas6_Axl','Igf1_Igf1r','Osm_Il6st','Osm_Osmr','Pros1_Axl','Sema4d_Plxnb2','Spp1_Cd44','Spp1_Itgav','Thbs1_Cd47','Tnfsf12_Tnfrsf12a','Vegfa_Nrp1','Vegfa_Nrp2','Vegfb_Nrp1','Cxcl12_Ackr3','Il1a_Il1r1','Il1a_Il1r2','Il1b_Il1r1','Il1b_Il1r2','Il1rn_Il1r1','Il1rn_Il1r2','Pdgfa_Pdgfra','Pdgfb_Pdgfra','Pdgfb_Pdgfrb','Thbs1_Sdc1','Tgfb1_Tgfbr1','Vcam1_Itgb1','Vegfb_Flt1', 'Apoe_Lrp5','Tgfb1_Tgfbr2','Tgfb1_Tgfbr3')


# Convert p-values into numbers
astrofibro_result[['pval']] <- as.numeric(astrofibro_result[['pval']])
astrofibro_result[['pval']][astrofibro_result[['Ligand_pct']] < 0.1] <- NA
astrofibro_result[['pval']][astrofibro_result[['Receptor_pct']] < 0.1] <- NA
astrofibro_result[['log_pval']] <- -log10(astrofibro_result[['pval']])
astrofibro_result[['log_pval']][is.infinite(astrofibro_result[['log_pval']])] <- -log10(1/1000) # 1000 permutations

astrofibro_result[['adj_pval']][astrofibro_result[['Ligand_pct']] < 0.1] <- NA
astrofibro_result[['adj_pval']][astrofibro_result[['Receptor_pct']] < 0.1] <- NA
astrofibro_result[['log_adj_pval']] <- -log10(as.numeric(astrofibro_result[['adj_pval']]))
astrofibro_result[['log_adj_pval']][is.infinite(astrofibro_result[['log_adj_pval']])] <- -log10(1/1000) # 1000 permutations
astrofibro_result[['time']] <- factor(astrofibro_result[['split_by']], 
                                      levels = c('Uninjured','1dpi','3dpi','7dpi'))


# Take significant scores from 3dpi/7dpi
astrofibro_result_subset <- astrofibro_result %>%
  filter(time %in% c('3dpi','7dpi')) %>%
  .[!is.na(.[['pval']]),] %>%
  filter(pval < 0.01)

# Which LR_pairs are exclusive to each receptor_cell? 
tmp <- apply(X = table(astrofibro_result_subset[['LR_pair']], 
                       astrofibro_result_subset[['Receptor_cell']]), 
             FUN = function(x) which(x != 0), 
             MARGIN = 1)
tmp <- sapply(X = tmp, FUN = function(x) paste(names(x), collapse = '_'))
tmp_order <- names(sort(tmp))


# Generature LR plot
astrofibro_plot <- astrofibro_result_subset %>%
  filter(LR_pair %in% lr_vetted) %>%
  # mutate(LR_pair = factor(x = LR_pair, levels = sort(unique(LR_pair), decreasing = TRUE))) %>%
  mutate(LR_pair = factor(x = LR_pair, levels = rev(tmp_order))) %>%
  mutate(Ligand_cell = plyr::mapvalues(x = Ligand_cell,
                                       from = c('Macrophage-A','Macrophage-B'),
                                       to = c('Mac-A','Mac-B'))) %>%
  ggplot(mapping = aes(x = Ligand_cell, y = LR_pair)) +
  geom_point(mapping = aes(size = log_pval, fill = Score), color = 'black', pch = 21) +
  ggh4x::facet_nested(. ~ time + Receptor_cell, switch = 'x') +
  xlab(label = 'Ligand-expressing cell') + 
  ylab(label = 'Ligand-Receptor pair') +
  scale_y_discrete(drop = TRUE, position = 'right') +
  scale_x_discrete(position = 'top') +
  scale_radius(limits = c(0,3)) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral')),
                       limits = c(0,NA)) +
  theme(strip.text = element_text(size = 12, color = 'black', face = 'bold'),
        strip.background = element_rect(fill = NA, color = 'black'),
        axis.title.x.top = element_text(size = 12, color = 'black', face = 'bold'),
        axis.title.y.right = element_text(size = 12, color = 'black', face = 'bold', hjust = 0.25),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0, face = 'bold'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black', face = 'bold', angle = 90),
        legend.key = element_rect(fill = NA),
        legend.position = 'right',
        legend.justification = c(1,0),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,-10,unit = 'mm'),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70')) +
  guides(fill = guide_colorbar(title = 'Score',
                               title.position = 'left',
                               title.hjust = 0.5,
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               barheight = 5,
                               frame.colour = 'black',
                               ticks.colour = 'black'),
         size = guide_legend(title = '-log10(p-value)', 
                             title.position = 'left',
                             title.hjust = 0.5,
                             override.aes = list(fill = 'black')))
astrofibro_plot <- cowplot::plot_grid(
  astrofibro_plot,
  cowplot::ggdraw() + 
    cowplot::draw_text(text = 'Receptor-expressing cell',
                       fontface = 'bold',
                       size = 12,
                       hjust = 1.05),
  ncol = 1,
  rel_heights = c(1,0.05)
)
astrofibro_plot
ggsave(filename = paste0(results_out, 'AstrocyteFibroblast_MacrophageLigand_LRplot.tiff'),
       plot = astrofibro_plot, device = 'tiff', height = 7.5, width = 5.75)






# working space -----------------------------------------------------------


# Take pearson correlation between LRpairs across cell_pairs, then cluster by
# 1-correlation distance.
astrofibro_result <- read.table(file = paste0(results_out, 'sci_cluster_LigandReceptorResults_AstrocyteFibroblast_macrophage.tsv'))
astrofibro_result[['time']] <- factor(astrofibro_result[['split_by']], 
                                      levels = c('Uninjured','1dpi','3dpi','7dpi'))


tmp <- astrofibro_result[c('Cell_pair','LR_pair','Score','time')]
tmp <- unique(subset(tmp, subset = time %in% c('3dpi','7dpi') & LR_pair %in% lr_vetted))
tmp[['id']] <- paste(tmp[['Cell_pair']], tmp[['time']], sep = '_')
tmp <- tmp[c('id', 'LR_pair', 'Score')]
tmp <- reshape(tmp, idvar = 'id', timevar = 'LR_pair', direction = 'wide')
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
tmp_cor <- cor(x = tmp, method = 'pearson')
tmp_dist <- as.dist(m = 1 - tmp_cor)
tmp_hc <- hclust(d = tmp_dist, method = 'ward.D2')
tmp <- data.matrix(frame = tmp)
tmp_hc <- hclust(d = tmp_cor, method = 'ward.D2')
tmp_order <- gsub(pattern = 'Score.',
                  replacement = '',
                  x = tmp_hc$labels[tmp_hc$order])


astrofibro_plot <- astrofibro_result %>%
  mutate(LR_pair = factor(x = LR_pair, levels = sort(unique(LR_pair), decreasing = TRUE))) %>%
  mutate(Ligand_cell = plyr::mapvalues(x = Ligand_cell,
                                       from = c('Macrophage-A','Macrophage-B'),
                                       to = c('Mac-A','Mac-B'))) %>%
  ggplot(mapping = aes(x = Ligand_cell, y = LR_pair)) +
  geom_point(mapping = aes(size = log_pval, color = Score)) +
  scale_y_discrete(drop = TRUE, position = 'right') +
  scale_x_discrete(position = 'top') +
  ggh4x::facet_nested(. ~ time + Receptor_cell, switch = 'x') +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral'))) +
  theme(strip.text = element_text(size = 12, color = 'black', face = 'bold'),
        strip.background = element_rect(fill = NA, color = 'black'),
        axis.title.x.top = element_text(size = 12, color = 'black', face = 'bold'),
        axis.title.y.right = element_text(size = 12, color = 'black', face = 'bold'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black', face = 'bold'),
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70')) +
  guides(color = guide_colorbar(title = 'Score',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               frame.colour = 'black',
                               ticks.colour = 'black'),
         size = guide_legend(title = '-log10(p-value)', 
                             override.aes = list(fill = 'black')))



Idents(sci) <- 'subcluster'
macroA_macroB_de <- FindMarkers(
  object = sci,
  ident.1 = 'Macrophage-A',
  ident.2 = 'Macrophage-B',
  assay = 'RNA'
)


# Complex heatmap attempt
tmp <- subset(astrofibro_result, subset = split_by == '3dpi' & adj_pval < 0.05)
tmp <- tmp[c('Cell_pair','LR_pair','Score')]
tmp <- reshape(tmp, idvar = 'Cell_pair', timevar = 'LR_pair', direction = 'wide')
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
tmp <- data.matrix(frame = tmp)
tmp[is.na(tmp)] <- 0

ComplexHeatmap::Heatmap(
  matrix = data.matrix(frame = tmp),
  cluster_rows = FALSE,
  col = rev(RColorBrewer::brewer.pal(n = 9, name = 'Spectral')),
  clustering_method_columns = 'ward.D2'
)


# attempt to cluster first, then plot via ggplot
tmp <- subset(astrofibro_result, subset = adj_pval < 0.01 & time == '3dpi')
tmp <- tmp[c('Cell_pair','LR_pair','Score')]
tmp <- reshape(tmp, idvar = 'Cell_pair', timevar = 'LR_pair', direction = 'wide')
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
tmp <- t(data.matrix(frame = tmp))
tmp_cor <- cor(x = tmp, use = 'pairwise.complete.obs')
tmp_dist <- dist(x = t(data.matrix(frame = tmp)), method = '')


tmp_dist <- dist(x = tmp)
tmp_hc <- hclust(d = )

tmp <- data.matrix(frame = tmp, rownames.force = )
heatmap(x = as.matrix(tmp))

test <- tmp %>% reshape2::acast(Cell_pair ~ LR_pair)
tmp %>%
  reshape2::acast(formula = 'Cell_pair' ~ 'LR_pair', value.var = 'Score')
head(tmp)
