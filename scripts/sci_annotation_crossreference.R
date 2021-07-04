
######## Cell-type annotation using reference data ########


# Use SingleR package to annotate clusters using a reference dataset. Following
# tutorial here: https://bioconductor.org/packages/3.11/bioc/vignettes/SingleR/inst/doc/SingleR.html

# We implement the pseudo-bulk aggregation method for shorter compute times.


# Data import and libraries ------------------------------------------------

require('SingleR')
require('Matrix')
require('Seurat')
require('dplyr')
require('ggplot2')
# data_in <- './data/sci.rds'
results_out <- './results/sci_annotation_crossReference/'
dir.create(path = results_out)

# scseq dataset
sci <- readRDS(file = './data/sci.rds')



# Rosenberg, Roco, et al. comparison via SingleR --------------------------

# Rosenberg, Roco, et al. Single-cell profiling of the developing mouse brain 
# and spinal cord with split-pool barcoding. April 2018. (doi:10.1126/science.aam8999).
# 
# Authors isolated cells from P2 and P11 mice and used single-nucleus RNAseq via
# SPLiT-seq.
# 
# Abbreviations:
#   - OPC = oligodendrocyte precursor cell
#   - VLMC = vascular leptomeningeal cell


ref_in <- './ref/GSM3017261_150000_CNS_nuclei.mat'
splitseq <- R.matlab::readMat(con = ref_in)


# Note: .mat file imported as list. "DGE" slot has expression values. Genes run
# along columns. Cells run along rows. ~150000 total cells in dataset (brain + 
# spinal cord).

is_spine <- c(splitseq$sample.type == 'p11_spine')
rosenberg_expr <- splitseq$DGE[is_spine,]
colnames(rosenberg_expr) <- gsub(pattern = ' ', replacement = '', x = c(splitseq$genes))
rosenberg_anno <- c(splitseq$spinal.cluster.assignment[is_spine,])


# Tidy the labels and combine similar cell-types (specifically excitatory 
# neurons and inhibitory neurons). Remove unresolved/un-annotated clusters.
rosenberg_anno <- gsub(pattern = ' ', replacement = '', x = rosenberg_anno)
rosenberg_anno <- gsub(pattern = '^[0-9]+', replacement = '', x = rosenberg_anno)
rosenberg_anno[grepl(pattern = 'Inhibitory', x = rosenberg_anno)] <- 'Inhibitory'
rosenberg_anno[grepl(pattern = 'Excitatory', x = rosenberg_anno)] <- 'Excitatory'
# names(table(rosenberg_anno))
remove_unlabeled <- (rosenberg_anno == 'NA' | 
                       rosenberg_anno == 'Unresolved' | 
                       rosenberg_anno == 'Unassigned')
rosenberg_anno <- rosenberg_anno[!remove_unlabeled]
rosenberg_expr <- rosenberg_expr[!remove_unlabeled,]
rosenberg_expr <- t(rosenberg_expr)
colnames(rosenberg_expr) <- paste('rosenberg', 1:ncol(rosenberg_expr), sep = '_') # make barcode

# log-normalize the reference dataset
rosenberg_expr <- LogNormalize(data = rosenberg_expr)


# Take only genes that are shared between datasets (by name)
sci_genes <- rownames(sci@assays$RNA@data)
sci_only_genes <- sci_genes[!sci_genes %in% rownames(rosenberg_expr)]
shared_genes <- intersect(sci_genes, rownames(rosenberg_expr))
test_set <- sci@assays$RNA@data[shared_genes, ]
ref_set <- rosenberg_expr[shared_genes,]


# run singler
t1 <- Sys.time()
rosenberg_singler <- SingleR(test = test_set, 
                             ref = ref_set, 
                             labels = rosenberg_anno,
                             de.method = 'wilcox')
print(Sys.time() - t1)
sci[['Rosenberg_SingleR']] <- rosenberg_singler[['labels']]
write.table(x = data.frame(rosenberg_singler), 
            file = paste0(results_out, 'rosenberg_singleR_results.tsv'),
            sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)


# clear workspace
rm(test_set, ref_set, sci_genes, sci_only_genes, shared_genes, rosenberg_expr, rosenberg_anno, remove_unlabeled, is_spine, splitseq)




# Sathyamurthy et al. comparison via SingleR  -----------------------------


# Sathyamurthy et al. Massively Parallel Single Nucleus Transcriptional 
# Profiling Defines Spinal Cord Neurons and Their Activity during Behavior. Cell
# Reports Feb. 2018. (https://doi.org/10.1016/j.celrep.2018.02.003).
# 
# Authors isolated nuclei from adult mouse spinal (control, formalin, rotarod) 
# and used Drop-seq based method.
# 
# Control animal annotations in barcodes:
# GSM2785675	f1
# GSM2785676	m1
# GSM2785677	M4
# GSM2785678	M5
# GSM2785679	F3
# GSM2785680	F4
# GSM2943542	facx
# GSM2943543	fbcx
# GSM2943544	Macx


# Read in data and annotations. Note: Input matrix is not sparse. Gene names are
# stored in first column as "Gene". Genes along rows, cells along cols. 
ref_in <- './ref/GSE103892_Expression_Count_Matrix.txt'
sathy_expr <- read.table(file = ref_in, header = TRUE, sep = '\t')
sathy_expr <- sathy_expr %>%
  tibble::column_to_rownames(var = 'Gene')
sathy_expr <- Matrix::Matrix(data = data.matrix(sathy_expr), sparse = TRUE)
sathy_expr


# load cell-type annotations
ref_in <- './ref/GSE103892_Sample_Cell_Cluster_Information.txt'
sathy_anno <- read.table(file = ref_in, header = TRUE, sep = '\t')
sathy_anno <- sathy_anno[,c('sample_cellbarcode', 'cell.type')]

# Tidy the labels and combine similar cell-types (mostly neurons). Remove 
# "discarded".
names(table(sathy_anno$cell.type))
neuron <- c(paste('DE', 1:16, sep = '-'),
            paste('DI', 1:9, sep = '-'),
            paste('ME', 1:1, sep = '-'),
            paste('MI', 1:4, sep = '-'),
            paste('VE', 1:4, sep = '-'),
            paste('VI', 1:5, sep = '-'),
            paste('VC', 1:2, sep = '-'),
            paste('VM', 1:2, sep = '-'),
            'DE12') # stray label
sathy_anno$cell.type[sathy_anno$cell.type %in% neuron] <- 'neurons'
sathy_anno <- sathy_anno[sathy_anno$cell.type != '' & sathy_anno$cell.type != 'discarded',]


# Tidy expression matrix. Keep only cells with retained annotations. Keep only
# genes that are shared between datasets. Then log-transform data.
tmp <- match(x = sathy_anno$sample_cellbarcode, table = colnames(sathy_expr))
sathy_expr <- sathy_expr[,tmp]
sci_genes <- rownames(sci@assays$RNA@data)
sci_only_genes <- sci_genes[!sci_genes %in% rownames(sathy_expr)]
shared_genes <- intersect(sci_genes, rownames(sathy_expr))


# Log-transform reference data. Prep matrices with shared genes.
sathy_expr <- LogNormalize(sathy_expr)
test_set <- sci@assays$RNA@data[shared_genes,]
ref_set <- sathy_expr[shared_genes,]


# run singler
t1 <- Sys.time()
sathy_singler <- SingleR(test = test_set, 
                         ref = ref_set, 
                         labels = sathy_anno$cell.type,
                         de.method = 'wilcox')
print(Sys.time() - t1)
sci[['Sathyamurthy_SingleR']] <- sathy_singler[['labels']]
write.table(x = data.frame(sathy_singler), 
            file = paste0(results_out, 'Sathyamurthy_singleR_results.tsv'),
            sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)


# clear workspace
rm(test_set, ref_set, sci_genes, sci_only_genes, shared_genes, sathy_expr, sathy_anno)




# Zeisel et al. comparison via SingleR  --------------------------------------


# TO DO: Need to figure out best way to import data for SingleR. Will put on 
# backburner for now.

# This is the big one: Zeisel et al. Molecular Architecture of the Mouse Nervous
# System. Cell Aug. 2018. (https://doi.org/10.1016/j.cell.2018.06.021).

# Majority of cells collected via 10X. 




# Visualize predicted cell-type annotations ----------------------------------


# import and set seurat object metadata
rosenberg_singler <- read.table(file = paste0(results_out, 'Rosenberg_singleR_results.tsv'), sep = '\t', header = TRUE, row.names = 1)
sathy_singler <- read.table(file = paste0(results_out, 'Sathyamurthy_singleR_results.tsv'), sep = '\t', header = TRUE, row.names = 1)
sci@meta.data[['Rosenberg_SingleR']] <- rosenberg_singler[['labels']]
sci@meta.data[['Sathyamurthy_SingleR']] <- sathy_singler[['labels']]


# some text adjustments
sci@meta.data[['Rosenberg_SingleR']] <- plyr::mapvalues(
  x = sci@meta.data[['Rosenberg_SingleR']],
  from = c('CommittedOligodendroctyePrecursorCells', 
           'CerebrospinalFluid-ContactingNeurons(CSF-cNs)',
           'OligoMature',
           'OligodendroctyeMyelinating',
           'Gammamotorneurons',
           'Alphamotorneurons',
           'Excitatory',
           'Inhibitory',
           'VLMC'),
  to = c('Committed OPC',
         'CSF-contacting Neuron',
         'Mature Oligo',
         'Myelinating Oligo',
         'Gamma Motor Neuron',
         'Alpha Motor Neuron',
         'Excitatory Neuron',
         'Inhibitory Neuron',
         'Vascular Leptomeningeal')
)
sci@meta.data[['Sathyamurthy_SingleR']] <- plyr::mapvalues(
  x = sci@meta.data[['Sathyamurthy_SingleR']],
  from = c('astrocytes', 'endo/vascular', 'meninges/schwann','microglia','neurons','Oligos'),
  to = c('Astrocyte','Endo/Vascular','Meninges/Schwann','Microglia','Neuron','Oligo')
)


# visualize predictions with umaps
# theme settings
umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 16, color = 'black'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 12, color = 'black'))

# Rosenberg predicted labels
Idents(sci) <- 'Rosenberg_SingleR'
tmp <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','ident')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = ident), size = 0.1, alpha = 0.4) +
  xlab(label = 'UMAP 1') +
  ylab(label = 'UMAP 2') +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black'),
        legend.box.margin = margin(0,0,0,-1, unit = 'cm')) +
  guides(color = guide_legend(title = 'Cell-type predictions\n(Rosenberg, Roco, et al. 2018)', override.aes = list(size = 5, alpha = 1))) 
rosenberg_umap <- LabelClusters(plot = tmp, id = 'ident', repel = TRUE, size = 4)
ggsave(filename = paste0(results_out, 'Rosenberg_predicted_annotation.tiff'),
       plot = rosenberg_umap, device = 'tiff', height = 5, width = 7.5)


# Sathyamurthy predicted labels
Idents(sci) <- 'Sathyamurthy_SingleR'
tmp <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','ident')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = ident), size = 0.1, alpha = 0.4) +
  xlab(label = 'UMAP 1') +
  ylab(label = 'UMAP 2') +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black'),
        legend.box.margin = margin(0,0,0,-1, unit = 'cm')) +
  guides(color = guide_legend(title = 'Cell-type predictions\n(Sathyamurthy et al. 2018)', override.aes = list(size = 5, alpha = 1)))
sathyamurthy_umap <- LabelClusters(plot = tmp, id = 'ident', repel = TRUE, size = 5)
ggsave(filename = paste0(results_out, 'Sathyamurthy_predicted_annotation.tiff'),
       plot = sathyamurthy_umap, device = 'tiff', height = 5, width = 7.5)



# heatmap
rosenberg_heatmap <- table('default_cluster' = sci$default_cluster, 
                           'prediction' = sci$Rosenberg_SingleR) %>%
  proportions(margin = 1) %>%
  round(digits = 1) %>%
  reshape2::melt() %>%
  ggplot(mapping = aes(x = as.factor(default_cluster), y = prediction)) +
  geom_raster(mapping = aes(fill = value)) +
  scale_x_discrete(expand = c(0,0), breaks = as.character(levels(sci$default_cluster))) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = 'A') +
  xlab(label = 'SCI cluster') +
  ylab(label = 'Cell-type prediction\n(Rosenberg, Roco, et al. 2018)') +
  theme(axis.text.x = element_text(size = 11, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black')) +
  guides(fill = guide_colorbar(title = 'Proportion', frame.colour = 'black', ticks.colour = 'black'))
ggsave(filename = paste0(results_out, 'Rosenberg_predicted_annotation_heatmap.tiff'),
       plot = rosenberg_heatmap, device = 'tiff', height = 4.5, width = 12)



sathy_heatmap <- table('default_cluster' = sci$default_cluster, 
                           'prediction' = sci$Sathyamurthy_SingleR) %>%
  proportions(margin = 1) %>%
  round(digits = 1) %>%
  reshape2::melt() %>%
  ggplot(mapping = aes(x = as.factor(default_cluster), y = prediction)) +
  geom_raster(mapping = aes(fill = value)) +
  scale_x_discrete(expand = c(0,0), breaks = as.character(levels(sci$default_cluster))) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = 'A') +
  xlab(label = 'SCI cluster') +
  ylab(label = 'Cell-type prediction\n(Sathyamurthy et al. 2018)') +
  theme(axis.text.x = element_text(size = 11, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black')) +
  guides(fill = guide_colorbar(title = 'Proportion', frame.colour = 'black', ticks.colour = 'black'))
ggsave(filename = paste0(results_out, 'Sathyamurthy_predicted_annotation_heatmap.tiff'),
       plot = sathy_heatmap, device = 'tiff', height = 2.75, width = 11)



# rm(list = ls()); gc()




# Extended data figure 1 --------------------------------------------------


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
                    axis.title = element_text(size = 16, color = 'black'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 12, color = 'black'))


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
  theme(legend.text = element_text(size = 16, color = 'black'),
        legend.position = 'none') +
  guides(color = guide_legend(title = 'Cell-type (#)', override.aes = list(size = 8, alpha = 1)))

exfig1 <- cowplot::plot_grid(rosenberg_umap, sathyamurthy_umap, ncol = 1, align = 'hv', axis = 'ltbr')
ggsave(filename = paste0(results_out, 'extendedFigure1_predictions.tiff'),
       plot = exfig1, device = 'tiff', height = 10, width = 7.5)



# # Take UMAP coordinates and for each cell calculate median distance to 50 
# # closest cells. Take sample of cells from each time-point so that each group
# # has same number of cells.
# minCells_perTime <- min(table(sci@meta.data[['time']]))
# umap_sample <- split(x = sci@meta.data, f = sci@meta.data$time)
# umap_sample <- lapply(
#   X = umap_sample,
#   FUN = function(x) {
#     rownames(x)[sample(x = 1:nrow(x), size = minCells_perTime)]
#   }
# )
# umap_coord <- sci[['umap']]@cell.embeddings[unlist(umap_sample, use.names = FALSE),]
# sci_k_dist <- BiocNeighbors::findKNN(X = umap_coord, k = 50)[['distance']]
# median_k_dist <- median(sci_k_dist[, ncol(sci_k_dist)])
# 
# # For each time point, calculate the neighborhood density of each cell by
# # using a tri-cube kernel on UMAP coordinates. Median distance of all cells to 
# # 50 nearest neighbors used to define maximum distance for the kernel. Method
# # inspired by: https://www.nature.com/articles/s41586-019-0933-9#Sec7
# times <- levels(sci@meta.data[['time']])
# dens <- c()
# for (t in times) {
#   tmp_coord <- sci[['umap']]@cell.embeddings[umap_sample[[t]],]
#   dists <- BiocNeighbors::findNeighbors(X = tmp_coord,
#                                         threshold = median_k_dist)[['distance']]
#   tmp_dens <- sapply(
#     X = dists,
#     FUN = function(x) {
#       sum((1 - (x / median_k_dist)^3)^3)
#     }
#   )
#   names(tmp_dens) <- rownames(tmp_coord)
#   dens <- c(dens, tmp_dens)
# }
# umap_data <- expand.grid(unlist(umap_sample, use.names = FALSE),
#                          levels(sci@meta.data[['time']]),
#                          stringsAsFactors = FALSE)
# colnames(umap_data) <- c('barcode','time')
# umap_data <- cbind(umap_data,sci[['umap']]@cell.embeddings[umap_data$barcode,])
# umap_data[['dens']] <- NA
# umap_data[['dens']] <- ifelse(
#   test = umap_data$time == sci@meta.data[umap_data$barcode,'time'],
#   yes = dens[umap_data$barcode[umap_data$time == sci@meta.data[umap_data$barcode,'time']]],
#   no = NA
# )
# umap_data[['dens']] <- log2(umap_data[['dens']])
# umap_data[['time']] <- factor(
#   x = umap_data[['time']],
#   levels = levels(sci@meta.data[['time']])
# )
# 
# tmp_cols <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')[3:9])(50)
# umap_density <- umap_data[order(umap_data$dens, na.last = FALSE, decreasing = FALSE),] %>%
#   ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
#   geom_point(mapping = aes(color = dens), size = 0.5, alpha = 0.4) +
#   scale_color_gradientn(colors = tmp_cols,
#                         na.value = 'grey70',
#                         breaks = c(0, max(umap_data$dens, na.rm = TRUE)),
#                         labels = c('Low','High')) +
#   facet_wrap(. ~ time, ncol = 4) +
#   xlab(label = 'UMAP 1') + 
#   ylab(label = 'UMAP 2') +
#   theme(panel.background = element_rect(fill = NA, colour = NA),
#         panel.border = element_rect(fill = NA, color = NA),
#         axis.line = element_line(size = 1),
#         axis.title = element_text(size = 12, color = 'black'),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         strip.text = element_text(size = 18, color = 'black'),
#         legend.title = element_text(size = 12, angle = 90),
#         legend.text = element_text(size = 12, color = 'black')) +
#   guides(color = guide_colorbar(title = 'log2(density)',
#                                 title.position = 'left',
#                                 title.hjust = 0.5,
#                                 frame.colour = 'black',
#                                 ticks.colour = 'black'))
# 
# # joint heatmaps
# prediction_umaps <- (rosenberg_umap / sathyamurthy_umap)
# ggsave(filename = paste0(results_out, 'singleR_predictions_umaps.tiff'),
#        plot = prediction_umaps, device = 'tiff', height = 8, width = 9)
# fig1_umaps <- cowplot::plot_grid(celltype_umap, umap_density, ncol = 1, rel_heights = c(1,0.4))
# ggsave(filename = paste0(results_out, 'extendedDataFig1_celltypeUMAP.tiff'),
#        plot = fig1_umaps, device = 'tiff', height = 8, width = 9)
