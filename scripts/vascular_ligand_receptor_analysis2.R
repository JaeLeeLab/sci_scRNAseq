

########## Ligand-receptor interaction analysis of angiogenesis ###########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
devtools::install_github(
  repo = "https://github.com/JamesChoi94/SingleCellTools.git", 
  INSTALL_opts = '--no-multiarch'
)
require('Seurat')
require('dplyr')
require('ggplot2')
require('SingleCellTools')

results_in <- './results/ligand_receptor_analysis/'
results_out <- './results/vascular_ligand_receptor_analysis/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)


# data
sci <- readRDS(file = './data/sci.rds')
lr_results <- read.table(file = paste0(results_in, 'sci_cluster_LigandReceptorResults_manuscript_interactions.tsv'))



# Unbiased vascular LR analysis (Compute score)  -------------------------------

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
  from = c('C-Endothelial','A-Endothelial','V-Endothelial'),
  to = c('Endothelial', 'Endothelial','Endothelial')
)
sci$subcluster <- plyr::mapvalues(
  x = sci$subcluster, from = 'Tip Cell', to = 'Tip-Cell'
)
sci <- sci[, !sci$subcluster %in% c('U-Vascular','Unassigned')]
sci$subcluster <- gsub(pattern = ' ', replacement = '-', x = sci$subcluster)
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
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0, size = 10),
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
ggsave(filename = paste0(results_out, 'endothelial_1dpi_myeloid.tiff'),
       plot = endo_mye_1dpi, device = 'tiff', height = 9, width = 6)
# ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_myeloid.png',
#        plot = endo_mye_1dpi, device = 'png', height = 9, width = 6)
endo_vasc_1dpi <- filter_lr(endo_lr, comp = 'Vascular', time = '1dpi') %>% 
  plot_lr(split = TRUE, max_score = 2.7)
ggsave(filename = paste0(results_out, 'endothelial_1dpi_vascular.tiff'),
       plot = endo_vasc_1dpi, device = 'tiff', height = 16, width = 5)
# ggsave(filename = './results/revision_figures/LR_results/endothelial_1dpi_vascular.png',
#        plot = endo_vasc_1dpi, device = 'png', height = 16, width = 5)
endo_macro_1dpi <- filter_lr(endo_lr, comp = 'Macroglia', time = '1dpi') %>% 
  plot_lr(split = TRUE, max_score = 2.7)
ggsave(filename = paste0(results_out, 'endothelial_1dpi_macroglia.tiff'),
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




# Validate scores with expression dot plot (p-values) --------------------------


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
                            'C-Endothelial',
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
select_sci_subclusters <- plyr::mapvalues(
  x = select_sci_subclusters,
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
  levels = c("Neutrophil", "Monocyte","Macrophage","Dendritic","Microglia","Div-Myeloid", "Fibroblast","Endothelial","Tip Cell","VSMC", "Pericyte","OPC", "Oligodendrocyte", "Astrocyte","Ependymal")
)
Idents(select_sci) <- 'tmp_lr'



# angiopoeitin genes
angpt_ligands <- c('Angpt1','Angpt2')
angpt_receptors <- c('Tie1','Tek')

DefaultAssay(select_sci) <- 'RNA'
avg_exp_angpt <- ScaleData(select_sci[['RNA']]@data, features = c(angpt_ligands, angpt_receptors))
avg_exp_angpt <- cbind(t(avg_exp_angpt), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp_angpt <- select_sci[['RNA']]@counts[c(angpt_ligands, angpt_receptors),]
pct_exp_angpt <- cbind(Matrix::t(pct_exp_angpt), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)


# vegf genes
vegf_receptors <- c('Flt1','Kdr','Flt4','Nrp1','Nrp2')
vegf_ligands <- c('Vegfa','Vegfb','Vegfc','Vegfd','Pgf')

DefaultAssay(select_sci) <- 'RNA'
avg_exp_vegf <- ScaleData(select_sci[['RNA']]@data, features = c(vegf_ligands, vegf_receptors))
avg_exp_vegf <- cbind(t(avg_exp_vegf), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp_vegf <- select_sci[['RNA']]@counts[c(vegf_ligands, vegf_receptors),]
pct_exp_vegf <- cbind(Matrix::t(pct_exp_vegf), select_sci@meta.data[,c('tmp_lr','time')]) %>%
  reshape2::melt(id.vars = c('tmp_lr','time')) %>%
  group_by(tmp_lr, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)

angpt_dat <- merge(avg_exp_angpt, pct_exp_angpt) %>%
  mutate('class' = 'angpt')
vegf_dat <- merge(avg_exp_vegf, pct_exp_vegf) %>%
  mutate('class' = 'vegf')
vasc_dat <- rbind(angpt_dat, vegf_dat) %>%
  mutate(tmp_lr = factor(x = tmp_lr, 
                         levels = rev(c("Neutrophil",
                                        "Monocyte",
                                        "Macrophage",
                                        "Dendritic",
                                        "Microglia",
                                        "Div-Myeloid", 
                                        "Fibroblast",
                                        "Endothelial",
                                        "Tip Cell",
                                        "VSMC", 
                                        "Pericyte",
                                        "OPC",
                                        "Oligodendrocyte",
                                        "Astrocyte",
                                        "Ependymal",
                                        "Lymphocyte",
                                        "Neuron")))) %>%
  mutate(variable = plyr::mapvalues(x = variable,
                                    from =  c('Tek','Flt1','Kdr','Flt4','Pgf'),
                                    to = c('Tie2','Vegfr1','Vegfr2','Vegfr3','Plgf')))

max_z <- ifelse(
  test = max(vasc_dat$avg.exp) > min(vasc_dat$avg.exp),
  yes = max(vasc_dat$avg.exp),
  no = max(abs(vasc_dat$avg.exp))
)
max_z <- ceiling(max_z*10)/10
myColors <- rev(RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))


select_sci$tmp_de <- paste(select_sci$tmp_lr, select_sci$time, sep = '_')
Idents(select_sci) <- 'tmp_de'
tmp_de <- FindAllMarkers(
  object = select_sci,
  assay = 'RNA',
  slot = 'data',
  only.pos = TRUE,
  features = unique(c(angpt_ligands, angpt_receptors, vegf_ligands, vegf_receptors))
)
tmp_ident <- strsplit(x = as.character(tmp_de$cluster), split = '_')
tmp_celltype <- sapply(tmp_ident, `[`, 1)
tmp_time <- sapply(tmp_ident, `[`, 2)
tmp_de$celltype <- tmp_celltype
tmp_de$time <- tmp_time
tmp_de$gene <- plyr::mapvalues(
  x = tmp_de$gene,
  from =  c('Tek','Flt1','Kdr','Flt4','Pgf'),
  to = c('Tie2','Vegfr1','Vegfr2','Vegfr3','Plgf')
)

vasc_dat$significant <- 1
for (i in 1:nrow(vasc_dat)) {
  which_row <- which(as.character(tmp_de$celltype) == as.character(vasc_dat$tmp_lr[i]) & 
                       as.character(tmp_de$time) == as.character(vasc_dat$time[i]) &
                       as.character(tmp_de$gene) == as.character(vasc_dat$variable[i]))
  if (length(which_row) == 1) {
    vasc_dat$significant[i] <- tmp_de$p_val_adj[which_row]
  }
}
vasc_dat$significant <- ifelse(test = vasc_dat$significant < 1e-10,
                               yes = '*',
                               no = ' ')

make_plot <- function(x, vline = NULL) {
  tmp <- x %>%
    ggplot(mapping = aes(x = variable, y = tmp_lr)) +
    geom_point(mapping = aes(size = pct.exp, fill = avg.exp), pch = 21, color = 'black') +
    geom_vline(xintercept = vline, linetype = 'dashed', color = 'black', size = 1) +
    geom_text(mapping = aes(label = significant), size = 6, nudge_y = -0.15) +
    facet_grid(. ~ time) +
    scale_size(range = c(0,8), 
               limits = c(0,100),
               breaks = seq(25,100,25)) +
    scale_fill_gradientn(colors = myColors,
                         limits = c(-max_z, max_z),
                         breaks = c(-max_z, 0, max_z),
                         na.value = myColors[9]) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = 'black', face = 'italic'),
          axis.text.y = element_text(size = 12, color = 'black'),
          strip.background = element_rect(fill = NA, colour = 'black'),
          strip.text.y = element_text(angle = 270, hjust = 0.5, size = 14, color = 'black'),
          strip.text.x = element_text(size = 12, color = 'black'),
          strip.placement = 'outside',
          legend.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size = 12, color = 'black', hjust = 0.5),
          legend.title = element_text(size = 12, color = 'black', hjust = 0.5),
          legend.box.margin = margin(-10, 0, -10, -10),
          legend.direction = 'horizontal',
          legend.box = 'horizontal',
          legend.position = 'bottom',
          legend.justification = 'left',
          panel.border = element_rect(fill = NA, size = 1, color = 'grey60'),
          panel.background = element_rect(fill = NA)) +
    guides(fill = guide_colorbar(title = 'z-score', 
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
  return(tmp)
}

angpt_dotplot <- vasc_dat %>%
  filter(class == 'angpt') %>%
  make_plot(vline = 2.50)
vegf_dotplot <- vasc_dat %>%
  filter(class == 'vegf') %>%
  make_plot(vline = 5.5)
vasc_dotplot <- (angpt_dotplot + NoLegend()) + vegf_dotplot + patchwork::plot_layout(widths = c(1,2))
vasc_dotplot
ggsave(filename = paste0(results_out, 'vascular_angpt_vegf_dotplot_pval.tiff'),
       plot = vasc_dotplot, device = 'tiff', height = 5.5, width = 16.5)

