

########## Ligand-receptor interaction analysis of angiogenesis ###########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
results_in <- './results/ligand_receptor_analysis/'
results_out <- './results/vascular_ligand_receptor_analysis/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)


# data
sci <- readRDS(file = './data/sci.rds')
lr_results <- read.table(file = paste0(results_in, 'sci_cluster_LigandReceptorResults_manuscript_interactions.tsv'))




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





# Angiopoeitin signaling (Ligand-receptor plot) ------------------------------


angpt_genes <- c('Angpt1','Angpt2','Tie1','Tek')
ligand_celltypes <- c('A-Endothelial','C-Endothelial','V-Endothelial','Tip Cell', 'VSMC','Fibroblast','Astrocyte')
receptor_celltypes <- c('A-Endothelial','C-Endothelial','V-Endothelial','Tip Cell', 'Pericyte','Fibroblast')
angpt_plot <- lr_results %>%
  filter(Ligand_cell %in% ligand_celltypes) %>%
  filter(Receptor_cell %in% receptor_celltypes) %>%
  filter(Ligand %in% angpt_genes & Receptor %in% angpt_genes) %>%
  filter(Pair_name != c('Angpt2_Tie1')) %>%
  filter(Ligand_cell_pct >= 2.5 & Receptor_cell_pct >= 2.5) %>%
  mutate('Ligand_cell_adj' = ifelse(
    test = Ligand_cell %in% c('A-Endothelial','C-Endothelial','V-Endothelial'),
    yes = 'Endothelial',
    no = as.character(Ligand_cell)
  )) %>%
  mutate('Receptor_cell_adj' = ifelse(
    test = Receptor_cell %in% c('A-Endothelial','C-Endothelial','V-Endothelial'),
    yes = 'Endothelial',
    no = as.character(Receptor_cell)
  )) %>%
  group_by(Ligand_cell_adj, Receptor_cell_adj, split_by, Pair_name, time, Receptor, Ligand) %>%
  summarise('Score_adj' = mean(Score),
            'log_pval' = min(log_pval),
            'time' = unique(time)) %>%
  mutate(Receptor_cell_adj = factor(Receptor_cell_adj, 
                                    levels = c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid','Microglia-A','Microglia-B','Microglia-C','Div-Microglia','Endothelial','Tip Cell','Pericyte','VSMC', 'Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC','Div-OPC','Pre-Oligo', 'Oligodendrocyte')),
         Ligand_cell_adj = factor(Ligand_cell_adj, 
                                  levels = c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid','Microglia-A','Microglia-B','Microglia-C','Div-Microglia','Endothelial','Tip Cell','Pericyte','VSMC', 'Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC','Div-OPC','Pre-Oligo', 'Oligodendrocyte')),
         Pair_name = plyr::mapvalues(x = Pair_name, 
                                     from = c('Angpt1_Tek','Angpt2_Tek'),
                                     to = c('Angpt1_Tie2','Angpt2_Tie2'))) %>%
  ggplot() +
  geom_point(mapping = aes(x = Ligand_cell_adj, y = Pair_name, size = log_pval, fill = Score_adj),
             color = 'black',
             pch = 21) +
  facet_grid(time ~ Receptor_cell_adj, switch = 'x', drop = TRUE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral'))) +
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
        legend.title = element_text(size = 12, color = 'black', face = 'bold'),
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70')) +
  guides(fill = guide_colorbar(title = 'Score',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               frame.colour = 'black',
                               ticks.colour = 'black'),
         size = guide_legend(title = '-log10(p-value)', 
                             override.aes = list(fill = 'black'))); angpt_plot
ggsave(filename = paste0(results_out, 'angiopoeitin_LigandReceptorplot.tiff'),
       plot = angpt_plot, device = 'tiff', height = 4.75, width = 8)





# VEGF signaling (Ligand-receptor plot) --------------------------------------


vegf_genes <- c('Flt1','Kdr','Nrp1','Nrp2','Vegfa','Vegfb','Pgf')
ligand_celltypes <- c('Monocyte','Macrophage-A','Macrophage-B','Microglia-A','Microglia-B','Microglia-C','Div-Microglia', 'Fibroblast','Endothelial','Tip Cell','VSMC','Astrocyte','OPC')
receptor_celltypes <- c('A-Endothelial','C-Endothelial','V-Endothelial')
vegf_plot <- lr_results %>%
  filter(Ligand_cell %in% ligand_celltypes) %>%
  filter(Receptor_cell %in% receptor_celltypes) %>%
  filter(Ligand %in% vegf_genes & Receptor %in% vegf_genes) %>%
  filter(Ligand_cell_pct >= 2.5 & Receptor_cell_pct >= 2.5) %>%
  filter(!Pair_name %in% c('Pgf_Nrp2')) %>%
  # merge endothelial cells into one
  mutate('Ligand_cell_adj' = ifelse(
    test = Ligand_cell %in% c('A-Endothelial','C-Endothelial','V-Endothelial'),
    yes = 'Endothelial',
    no = as.character(Ligand_cell)
  )) %>%
  mutate('Receptor_cell_adj' = ifelse(
    test = Receptor_cell %in% c('A-Endothelial','C-Endothelial','V-Endothelial'),
    yes = 'Endothelial',
    no = as.character(Receptor_cell)
  )) %>%
  # rename microglia to new cluster names
  mutate('Receptor_cell_adj' = plyr::mapvalues(
    x = Receptor_cell_adj, 
    from = c('Microglia-A','Microglia-B','Microglia-C','Div-Microglia'),
    to = c('H-Microglia','DAM-A','DAM-C','DAM-B')
  )) %>%
  mutate('Ligand_cell_adj' = plyr::mapvalues(
    x = Ligand_cell_adj, 
    from = c('Microglia-A','Microglia-B','Microglia-C','Div-Microglia'),
    to = c('H-Microglia','DAM-A','DAM-C','DAM-B')
  )) %>%
  group_by(Ligand_cell_adj, Receptor_cell_adj, split_by, Pair_name, time, Receptor, Ligand) %>%
  summarise('Score_adj' = mean(Score),
            'log_pval' = min(log_pval),
            'time' = unique(time)) %>%
  mutate(Receptor_cell_adj = factor(Receptor_cell_adj, 
                                    levels = c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','H-Microglia','DAM-A','DAM-B','DAM-C','Endothelial','Tip Cell','Pericyte','VSMC', 'Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC','Div-OPC','Pre-Oligo', 'Oligodendrocyte')),
         Ligand_cell_adj = factor(Ligand_cell_adj, 
                                  levels = c('Neutrophil','Monocyte','Macrophage-A','Macrophage-B','BA-Macrophage','Dendritic','Div-Myeloid','H-Microglia','DAM-A','DAM-B','DAM-C','Endothelial','Tip Cell','Pericyte','VSMC', 'Fibroblast','Ependymal','Astroependymal','Astrocyte','OPC','Div-OPC','Pre-Oligo', 'Oligodendrocyte')),
         Receptor = plyr::mapvalues(x = Receptor, from = c('Flt1','Kdr'), to = c('Vegfr1','Vegfr2')),
         Ligand = plyr::mapvalues(x = Ligand, from = c('Pgf'), to = c('Plgf'))) %>%
  mutate(Pair_name = paste(Ligand, Receptor, sep = '_')) %>%
  ggplot() +
  geom_point(mapping = aes(x = Ligand_cell_adj, y = Pair_name, size = log_pval, fill = Score_adj),
             color = 'black',
             pch = 21) +
  facet_grid(time ~ Receptor_cell_adj, switch = 'x', drop = TRUE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = 'Spectral'))) +
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
        legend.title = element_text(size = 12, color = 'black', face = 'bold'),
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70')) +
  guides(fill = guide_colorbar(title = 'Score',
                               frame.linewidth = 1,
                               ticks.linewidth = 1,
                               frame.colour = 'black',
                               ticks.colour = 'black'),
         size = guide_legend(title = '-log10(p-value)', 
                             override.aes = list(fill = 'black'))); vegf_plot
ggsave(filename = paste0(results_out, 'vegf_LigandReceptorplot.tiff'),
       plot = vegf_plot, device = 'tiff', height = 8.5, width = 6.25)






# Validates scores with expression dot plot ------------------------------------


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
Idents(select_sci) <- 'celltype'

# Extract interesting subclusters
select_sci@meta.data[['tmp_lr']] <- ifelse(
  test = select_sci@meta.data[['subcluster']] == 'Tip Cell',
  yes = 'Tip Cell',
  no = as.character(select_sci@meta.data[['celltype']])
)
select_sci@meta.data[['tmp_lr']] <- ifelse(
  test = select_sci@meta.data[['subcluster']] == 'VSMC',
  yes = 'VSMC',
  no = as.character(select_sci@meta.data[['tmp_lr']])
)


select_sci@meta.data[['tmp_lr']] <- factor(
  x = select_sci@meta.data[['tmp_lr']],
  levels = c("Neutrophil", "Monocyte","Macrophage","Dendritic","Microglia","Div-Myeloid", "Fibroblast","Endothelial","Tip Cell","VSMC", "Pericyte","OPC", "Oligodendrocyte", "Astrocyte","Ependymal","Lymphocyte", "Neuron")
)
Idents(select_sci) <- 'tmp_lr'


# angiopoeitin genes
angpt_ligands <- c('Angpt1','Angpt2')
angpt_receptors <- c('Tie1','Tek')

DefaultAssay(select_sci) <- 'RNA'
avg_exp <- ScaleData(select_sci[['RNA']]@data, features = c(angpt_ligands, angpt_receptors))
avg_exp <- cbind(t(avg_exp), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- select_sci[['RNA']]@counts[c(angpt_ligands, angpt_receptors),]
pct_exp <- cbind(Matrix::t(pct_exp), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min_expr, 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))
max_expr <- ceiling(max(avg_exp$avg.exp)*10)/10
min_expr <- floor(min(avg_exp$avg.exp)*10)/10


# angpt dot plot
angpt_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(tmp_lr = factor(tmp_lr, levels = rev(c("Neutrophil", "Monocyte","Macrophage","Dendritic","Microglia","Div-Myeloid", "Fibroblast","Endothelial","Tip Cell","VSMC", "Pericyte","OPC", "Oligodendrocyte", "Astrocyte","Ependymal","Lymphocyte", "Neuron")))) %>%
  mutate(variable = plyr::mapvalues(variable, from = c('Tek'), to = c('Tie2'))) %>%
  ggplot(mapping = aes(x = variable, y = tmp_lr)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), pch = 21, color = 'black') +
  geom_vline(xintercept = 2.5, linetype = 'dashed', color = 'black', size = 1) +
  facet_grid(. ~ time) +
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
        strip.text.y = element_text(angle = 270, hjust = 0.5, face = 'bold', size = 14, color = 'black'),
        strip.text.x = element_text(size = 14, face = 'bold', color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 12, color = 'black', hjust = 0),
        legend.title = element_text(size = 12, color = 'black', hjust = 0.5),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.direction = 'horizontal',
        legend.box = 'horizontal',
        legend.position = 'bottom',
        legend.justification = 'left',
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
         size = guide_legend(title = '% detected', 
                             override.aes = list(fill = 'black'), 
                             title.position = 'top')); angpt_dotplot
ggsave(filename = paste0(results_out, 'angiopoeitin_dotplot.tiff'),
       plot = angpt_dotplot, device = 'tiff', height = 5.75, width = 5.5)
  



# vegf genes
vegf_ligands <- c('Flt1','Kdr','Flt4','Nrp1','Nrp2')
vegf_receptors <- c('Vegfa','Vegfb','Vegfc','Vegfd','Pgf')

DefaultAssay(select_sci) <- 'RNA'
avg_exp <- ScaleData(select_sci[['RNA']]@data, features = c(vegf_ligands, vegf_receptors))
avg_exp <- cbind(t(avg_exp), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- select_sci[['RNA']]@counts[c(vegf_ligands, vegf_receptors),]
pct_exp <- cbind(Matrix::t(pct_exp), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
# max_expr <- ceiling(max(avg_exp$avg.exp)*10)/10
max_expr <- 3
min_expr <- floor(min(avg_exp$avg.exp)*10)/10
myColors <- rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))(100))
myBreaks <- c(seq(min_expr, 0, length.out = 50),
              seq(max_expr/100, max_expr, length.out = 50))

# vegf dot plot
vegf_dotplot <- merge(avg_exp, pct_exp) %>%
  mutate(tmp_lr = factor(
    x = tmp_lr, levels = rev(c("Neutrophil", "Monocyte","Macrophage","Dendritic","Microglia","Div-Myeloid", "Fibroblast","Endothelial","Tip Cell","VSMC", "Pericyte","OPC", "Oligodendrocyte", "Astrocyte","Ependymal","Lymphocyte", "Neuron")))
  ) %>%
  mutate('variable' = plyr::mapvalues(
    x = variable,
    from =  c('Tek','Flt1','Kdr','Flt4','Pgf'),
    to = c('Tie2','Vegfr1','Vegfr2','Vegfr3','Plgf')
  )) %>%
  ggplot(mapping = aes(x = variable, y = tmp_lr)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), pch = 21, color = 'black') +
  geom_vline(xintercept = 5.5, linetype = 'dashed', color = 'black', size = 1) +
  facet_grid(. ~ time) +
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
        strip.text.y = element_text(angle = 270, hjust = 0.5, face = 'bold', size = 14, color = 'black'),
        strip.text.x = element_text(size = 14, face = 'bold', color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 12, color = 'black', hjust = 0),
        legend.title = element_text(size = 12, color = 'black', hjust = 0.5),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.direction = 'horizontal',
        legend.box = 'horizontal',
        legend.position = 'bottom',
        legend.justification = 'right',
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
         size = guide_legend(title = '% detected', 
                             override.aes = list(fill = 'black'), 
                             title.position = 'top'))
vegf_dotplot
ggsave(filename = paste0(results_out, 'vegf_dotplot.tiff'),
       plot = vegf_dotplot, device = 'tiff', height = 6, width = 11)



