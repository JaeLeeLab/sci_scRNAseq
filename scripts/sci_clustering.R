

######### SCI cluster analysis #########



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

# Set factors for metadata variablesx
sci$time <- factor(x = sci$time, levels = c('Uninjured','1dpi','3dpi','7dpi'))
sci$sample_id <- factor(sci$sample_id, levels = c('uninj_sample1', 'uninj_sample2', 'uninj_sample3', '1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2'))


# Dimensional reduction and initial clustering ------------------------------


# dim reduce and cluster
DefaultAssay(sci) <- 'integrated'
sci <- RunPCA(sci, npcs = 50)

# 20 PCs based on elbow plot
npcs <- 1:20
sci <- RunUMAP(sci, dims = npcs, min.dist = 0.3, n.neighbors = 30L, umap.method = 'uwot')
sci <- FindNeighbors(sci, dims = npcs)
sci <- FindClusters(sci, resolution = 0.8) # default resolution
bad_cells <- colnames(sci)[sci$seurat_clusters == 28] # 252 bad cells
write.table(x = bad_cells, file = paste0(results_out, 'low_quality_cluster_28.tsv'),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
sci <- subset(sci, idents = c(0:27, 29)) # low quality cluster removal. high mito + overlapping markers.


# repeat dim reduce and cluster 
DefaultAssay(sci) <- 'integrated'
sci <- RunPCA(sci, npcs = 50)
elbowplot <- ElbowPlot(sci, ndims = 50) +
  geom_vline(mapping = aes(xintercept = 20), size = 1, color = 'red', linetype = 'dashed') +
  labs(title = paste0('PCA of full SCI dataset (', ncol(sci), ' cells)'),
       subtitle = 'PCs up to dashed line used for downstream analysis')
ggsave(plot = elbowplot, filename = paste0(results_out, 'pca_elbowplot.tiff'),
       device = 'tiff', height = 3.5, width = 5)

# 20 PCs based on elbow plot
npcs <- 1:20
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




# Density UMAP over time --------------------------------------------------


# Take UMAP coordinates and for each cell calculate median distance to 50 
# closest cells. Take sample of cells from each time-point so that each group
# has same number of cells.
minCells_perTime <- min(table(sci@meta.data[['time']]))
umap_sample <- split(x = sci@meta.data, f = sci@meta.data$time)
umap_sample <- lapply(
  X = umap_sample,
  FUN = function(x) {
    rownames(x)[sample(x = 1:nrow(x), size = minCells_perTime)]
  }
)
umap_coord <- sci[['umap']]@cell.embeddings[unlist(umap_sample, use.names = FALSE),]
sci_k_dist <- BiocNeighbors::findKNN(X = umap_coord, k = 50)[['distance']]
median_k_dist <- median(sci_k_dist[, ncol(sci_k_dist)])

# For each time point, calculate the neighborhood density of each cell by
# using a tri-cube kernel on UMAP coordinates. Median distance of all cells to 
# 50 nearest neighbors used to define maximum distance for the kernel. Method
# inspired by: https://www.nature.com/articles/s41586-019-0933-9#Sec7
times <- levels(sci@meta.data[['time']])
dens <- c()
for (t in times) {
  tmp_coord <- sci[['umap']]@cell.embeddings[umap_sample[[t]],]
  dists <- BiocNeighbors::findNeighbors(X = tmp_coord,
                                        threshold = median_k_dist)[['distance']]
  tmp_dens <- sapply(
    X = dists,
    FUN = function(x) {
      sum((1 - (x / median_k_dist)^3)^3)
    }
  )
  names(tmp_dens) <- rownames(tmp_coord)
  dens <- c(dens, tmp_dens)
}
umap_data <- expand.grid(unlist(umap_sample, use.names = FALSE),
                         levels(sci@meta.data[['time']]),
                         stringsAsFactors = FALSE)
colnames(umap_data) <- c('barcode','time')
umap_data <- cbind(umap_data,sci[['umap']]@cell.embeddings[umap_data$barcode,])
umap_data[['dens']] <- NA
umap_data[['dens']] <- ifelse(
  test = umap_data$time == sci@meta.data[umap_data$barcode,'time'],
  yes = dens[umap_data$barcode[umap_data$time == sci@meta.data[umap_data$barcode,'time']]],
  no = NA
)
umap_data[['dens']] <- log2(umap_data[['dens']])
umap_data[['time']] <- factor(
  x = umap_data[['time']],
  levels = levels(sci@meta.data[['time']])
)

tmp_cols <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')[3:9])(50)
umap_density <- umap_data[order(umap_data$dens, na.last = FALSE, decreasing = FALSE),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = dens), size = 1, alpha = 0.4) +
  scale_color_gradientn(colors = tmp_cols,
                        na.value = 'grey70',
                        breaks = c(0, max(umap_data$dens, na.rm = TRUE)),
                        labels = c('Low','High')) +
  facet_wrap(. ~ time) +
  xlab(label = 'UMAP 1') + 
  ylab(label = 'UMAP 2') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.line = element_line(size = 0),
        axis.title = element_text(size = 22, color = 'black', face = 'bold'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.text = element_text(size = 22, color = 'black', face = 'bold'),
        legend.title = element_text(size = 12, angle = 90, face = 'bold'),
        legend.text = element_text(size = 12, color = 'black', face = 'bold')) +
  guides(color = guide_colorbar(title = 'log2(density)',
                                title.position = 'left',
                                title.hjust = 0.5,
                                frame.colour = 'black',
                                ticks.colour = 'black'))
# umap_density
# saveRDS(umap_density, file = paste0(results_out, 'sci_density_umap.rds'))
ggsave(filename = paste0(results_out, 'sci_density_umap.tiff'),
       plot = umap_density, device = 'tiff', height = 8.25, width = 9.5)



# save
saveRDS(sci, file = './data/sci.rds')