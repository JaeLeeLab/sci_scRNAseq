

######### SCI cluster analysis #########

## !!!!!!!! Needs to be modified for batchelor results ####

# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
# data_in <- './data/'
# data_out <- './data/data_integration/'
results_out <- './results/sci_clustering/'
ref_in <- './ref/'
ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)

sci <- readRDS(file = './data/sci.rds')

# utils
firstup <- function(x) {x <- tolower(x); substr(x,1,1) <- toupper(substr(x,1,1)); return(x)}

# Set factors for metadata variables
sci$time <- factor(x = sci$time, levels = c('Uninjured','1dpi','3dpi','7dpi'))
sci$sample_id <- factor(sci$sample_id, levels = c('uninj_sample1', 'uninj_sample2', 'uninj_sample3', '1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2'))



# Dimensional reduction and initial clustering ------------------------------

# Score cell cycle as per Seurat vignette (https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html)
DefaultAssay(sci) <- 'RNAcorrected'
s_genes <- firstup(cc.genes$s.genes)
g2m_genes <- firstup(cc.genes$g2m.genes)
sci <- CellCycleScoring(
  object = sci,
  s.features = s_genes,
  g2m.features = g2m_genes,
  assay = 'RNAcorrected'
)
sci$CC.Difference <- sci$S.Score - sci$G2M.Score

# Dim Reduce and Cluster (default Seurat resolution)
DefaultAssay(sci) <- 'integrated'
sci <- ScaleData(sci, vars.to.regress = 'CC.Difference')
sci <- RunPCA(sci, npcs = 50)
npcs <- 1:20 # 20 PCs based on elbow plot
sci <- RunUMAP(sci, dims = npcs, min.dist = 0.3, n.neighbors = 30L, umap.method = 'uwot')
sci <- FindNeighbors(sci, dims = npcs)
sci <- FindClusters(sci, resolution = 0.8) 


# Identify marker genes per cluster to investigate potential low quality cells.
# We check the top 2 marker genes per cluster and generate a heatmap. We expect
# low quality clusters to have overlapping expression of markes and possibly
# mitochondrial enrichment.
DefaultAssay(sci) <- 'RNA'
markers <- FindAllMarkers(
  object = sci,
  only.pos = TRUE,
  logfc.threshold = 0.5,
  test.use = 'wilcox',
  assay = 'RNA',
  slot = 'data'
)
marker_genes <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC) %>%
  .[['gene']] %>%
  unique()
heatmap_data <- DietSeurat(
  object = sci,
  counts = FALSE,
  data = TRUE,
  assays = 'RNA'
)
heatmap_data <- ScaleData(object = heatmap_data, features = marker_genes)
marker_heatmap <- DoHeatmap(
  object = heatmap_data,
  features = marker_genes,
  angle = 0
) +
  labs(title = 'DE markers - all SCI cells (pre-cluster filter)') +
  theme(legend.position = 'none',
        plot.title = element_text(size = 12))

# Build summary result figure
umap <- DimPlot(sci, group.by = 'seurat_clusters', label = TRUE, label.size = 5)
initial_cluster_result <- cowplot::plot_grid(umap, marker_heatmap, ncol = 1, rel_heights = c(0.4, 1))

# save initial cluster result
ggsave(filename = paste0(results_out, 'initial_cluster_results.tiff'),
       plot = initial_cluster_result, device = 'tiff',
       height = 18, width = 18)


# Cluster 28 shares expression of Csf1r, Slc1a2, Hexb, mt-Co3, and Pou2f2. These
# are likely low quality cells. 
bad_cluster <- c(28)
bad_cells <- colnames(sci)[sci$seurat_clusters == 28] # 63 bad cells
write.table(x = bad_cells, file = paste0(results_out, 'low_quality_cluster_28.tsv'),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
keep_cells <- unique(sci$seurat_clusters)[!unique(sci$seurat_clusters) %in% bad_cluster]
sci <- subset(sci, idents = keep_cells)


# Repeat Dim reduce and clustering without low quality cells
DefaultAssay(sci) <- 'integrated'
sci <- RunPCA(sci, npcs = 50)
npcs <- 1:20
elbowplot <- ElbowPlot(sci, ndims = 50) +
  geom_vline(mapping = aes_string(xintercept = max(npcs)), 
             size = 1, 
             color = 'red',
             linetype = 'dashed') +
  labs(title = paste0('PCA of full SCI dataset (', ncol(sci), ' cells)'),
       subtitle = 'PCs up to dashed line used for downstream analysis')
ggsave(plot = elbowplot, filename = paste0(results_out, 'pca_elbowplot.tiff'),
       device = 'tiff', height = 3.5, width = 5)

# 20 PCs based on elbow plot
sci <- RunUMAP(sci, dims = npcs, min.dist = 0.3, n.neighbors = 30L, umap.method = 'uwot')
sci <- FindNeighbors(sci, dims = npcs)
sci <- FindClusters(sci, resolution = 0.8) # default resolution
sci$default_cluster <- sci$seurat_clusters




# Preliminary UMAPs -------------------------------------------------------

# Global options
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

sample_cols <- RColorBrewer::brewer.pal(n = 10, name = 'Spectral')
sample_cols <- c('firebrick','dodgerblue','goldenrod','firebrick','dodgerblue','goldenrod','firebrick','dodgerblue','firebrick','dodgerblue')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')

dissociationMethod_cols <- c('firebrick','dodgerblue')
names(dissociationMethod_cols) <- c('Standard','Enriched')

chemistry_cols <- c('firebrick','dodgerblue')
names(chemistry_cols) <- c('v2','v3')


# Visualize
umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 12, color = 'black'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 12, color = 'black'))

# Default clusters UMAP
tmp <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','default_cluster')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = default_cluster), size = 0.1, alpha = 0.25) +
  umap_theme +
  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
umap_defaultcluster <- LabelClusters(plot = tmp, id = 'default_cluster', size = 5, repel = FALSE)
ggsave(filename = paste0(results_out, 'umap_defaultcluster.tiff'),
       plot = umap_defaultcluster, device = 'tiff', height = 4, width = 4.25)

# Injury time-point UMAP
umap_time <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = time), size = 0.1, alpha = 0.25) +
  scale_color_manual(values = time_cols) +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Time after\ninjury', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_time.tiff'),
       plot = umap_time, device = 'tiff', height = 4, width = 5.5)

# Original sample UMAP
umap_sample <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 0.2, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 1) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_sample.tiff'),
       plot = umap_sample, device = 'tiff', height = 3.5, width = 15)

# dissociationMethod UMAP
umap_dissociationMethod <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','dissociationMethod')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = dissociationMethod), size = 0.1, alpha = 0.3) +
  scale_color_manual(values = dissociationMethod_cols) +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Dissociation Method', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_dissociationMethod.tiff'),
       plot = umap_dissociationMethod, device = 'tiff', height = 4, width = 6.25)

# Chemistry UMAP
umap_chemistry <- FetchData(object = sci, vars = c('UMAP_1','UMAP_2','chemistry')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = chemistry), size = 0.1, alpha = 0.15) +
  scale_color_manual(values = chemistry_cols) +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = '10X Chromium\nchemistry', override.aes = list(size = 7, alpha = 1)))
ggsave(filename = paste0(results_out, 'umap_chemistry.tiff'),
       plot = umap_chemistry, device = 'tiff', height = 4, width = 5.75)


