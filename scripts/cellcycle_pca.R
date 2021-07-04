
require('Seurat')
require('dplyr')
require('ggplot2')


# Cell cycle PCA -------------------------------------------------------------

# Myeloid
myeloid <- readRDS(file = './data/myeloid.rds')
DefaultAssay(myeloid) <- 'integrated'
# myeloid <- RunPCA(myeloid, npcs = 20)
p2 <- PCAPlot(
  object = myeloid,
  group.by = 'Phase',
  dims = c(1,2),
  shuffle = TRUE,
  pt.size = 0.5
) + 
  theme_bw() +
  labs(title = 'Myeloids after cell-cycle regression') +
  scale_color_manual(values = c('G1' = '#e6194B',
                                'G2M' = '#000075',
                                'S' = 'goldenrod')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
myeloid <- ScaleData(object = myeloid)
myeloid <- RunPCA(myeloid, npcs = 20)
p1 <- PCAPlot(
  object = myeloid, 
  group.by = 'Phase', 
  dims = c(1,2),
  shuffle = TRUE, 
  pt.size = 0.5
) +
  theme_bw() +
  labs(title = 'Myeloids before cell-cycle regression') +
  scale_color_manual(values = c('G1' = '#e6194B',
                                'G2M' = '#000075',
                                'S' = 'goldenrod')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
p3 <- p1 | p2
# ggsave(filename = './results/revision_figures_round2/myeloid_pca_phase.tiff',
#        plot = p3, height = 3.75, width = 10, device = 'tiff')

# Vascular
vascular <- readRDS(file = './data/vascular.rds')
DefaultAssay(vascular) <- 'integrated'
p5 <- PCAPlot(
  object = vascular,
  group.by = 'Phase',
  dims = c(1,2),
  shuffle = TRUE,
  pt.size = 0.5
) + 
  theme_bw() +
  labs(title = 'Vascular after cell-cycle regression') +
  scale_color_manual(values = c('G1' = '#e6194B',
                                'G2M' = '#000075',
                                'S' = 'goldenrod')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
vascular <- ScaleData(object = vascular)
vascular <- RunPCA(object = vascular, npcs = 20)
p4 <- PCAPlot(
  object = vascular, 
  group.by = 'Phase', 
  dims = c(1,2),
  shuffle = TRUE, 
  pt.size = 0.5
) +
  theme_bw() +
  labs(title = 'Vascular before cell-cycle regression') +
  scale_color_manual(values = c('G1' = '#e6194B',
                                'G2M' = '#000075',
                                'S' = 'goldenrod')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
p6 <- p4 | p5
# ggsave(filename = './results/revision_figures_round2/vascular_pca_phase.tiff',
#        plot = p6, height = 3.75, width = 10, device = 'tiff')

# Macroglia
macroglia <- readRDS(file = './data/macroglia.rds')
DefaultAssay(macroglia) <- 'integrated'
p8 <- PCAPlot(
  object = macroglia,
  group.by = 'Phase',
  dims = c(1,2),
  shuffle = TRUE,
  pt.size = 0.5
) +
  theme_bw() +
  labs(title = 'Macroglia after cell-cycle regression')  +
  scale_color_manual(values = c('G1' = '#e6194B',
                                'G2M' = '#000075',
                                'S' = 'goldenrod')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
macroglia <- ScaleData(object = macroglia)
macroglia <- RunPCA(object = macroglia, npcs = 20)
p7 <- PCAPlot(
  object = macroglia,
  group.by = 'Phase',
  dims = c(1,2),
  shuffle = TRUE,
  pt.size = 0.5
) +
  theme_bw() + 
  labs(title = 'Macroglia before cell-cycle regression') +
  scale_color_manual(values = c('G1' = '#e6194B',
                                'G2M' = '#000075',
                                'S' = 'goldenrod')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
p9 <- p7 | p8
# ggsave(filename = './results/revision_figures_round2/macroglia_pca_phase.tiff',
#        plot = p9, height = 3.75, width = 10, device = 'tiff')

p10 <- p3 / p6 / p9
p10
ggsave(filename = './results/revision_figures_round2/pca_phase.tiff',
       plot = p10, device = 'tiff', height = 8, width = 8)

# SCI
sci <- readRDS(file = './data/sci.rds')
DefaultAssay(sci) <- 'integrated'
p10 <- DimPlot(
  object = sci, 
  group.by = 'celltype', 
  label = TRUE, 
  label.size = 4,
  repel = TRUE
) +
  theme_bw() +
  NoLegend()
p11 <- DimPlot(
  object = sci, 
  group.by = 'Phase',
  shuffle = TRUE
) +
  theme_bw()
p12 <- p10 + p11
ggsave(filename = './results/revision_figures_round2/sci_pca_phase.tiff',
       plot = p12, height = 4, width = 10, device = 'tiff')



# If need to repeat without regression ------------------------------------

myeloid <- FindNeighbors(myeloid, dims = 1:12)
myeloid <- RunUMAP(object = myeloid, dims = 1:12)
myeloid <- FindClusters(myeloid, resolution = 0.25)
# c12_markers <- FindMarkers(myeloid, ident.1 = 12, only.pos = TRUE, assay = 'RNA')

myeloid <- myeloid[,myeloid$integrated_snn_res.0.25 != 12]
p1 <- DimPlot(myeloid, label = TRUE, label.size = 5)
p2 <- DimPlot(myeloid, group.by = 'myeloid_subcluster', label = TRUE)
p1 + p2

myeloid$myeloid_subcluster <- plyr::mapvalues(
  x = myeloid$integrated_snn_res.0.25,
  from = 0:11,
  to = c('Homeostatic Microglia',
         'Inflammatory Microglia',
         'Chemotaxis-Inducing Mac',
         'Monocyte',
         'Dividing Microglia',
         'Inflammatory Mac',
         'Neutrophil',
         'Migrating Microglia',
         'Dendritic',
         'Interferon Myeloid',
         'Dividing Myeloid',
         'Border-Associated Mac')
)
p1 <- DimPlot(myeloid, label = TRUE, label.size = 5)
p2 <- DimPlot(myeloid, group.by = 'myeloid_subcluster', label = TRUE,
              repel = TRUE, label.size = 5)
p1 + p2




vascular <- readRDS(file = './data/vascular.rds')
DefaultAssay(vascular) <- 'integrated'
vascular <- ScaleData(object = vascular)
vascular <- RunPCA(vascular, npcs = 20)
vascular <- FindNeighbors(vascular, dims = 1:9)
vascular <- RunUMAP(object = vascular, dims = 1:9)
vascular <- FindClusters(vascular, resolution = 0.3)
p1 <- DimPlot(vascular, label = TRUE, label.size = 5) +
  theme_bw()
p2 <- DimPlot(vascular, group.by = 'vascular_subcluster', label = TRUE) +
  theme_bw()
p3 <- DotPlot(vascular, assay = 'RNA', features = c('Gkn3','Stmn2','Cldn5','Ly6a','Slc38a5','Icam1','Apln','Kcnj8','Tagln','Col1a1'), scale.by = 'size') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4 <- p3 + p1 + p2
ggsave(filename = './results/revision_figures_round2/vascular_noCCregression_reclustered.tiff', plot = p4, device = 'tiff', height = 4, width = 14)

macroglia <- readRDS(file = './data/macroglia.rds')
DefaultAssay(macroglia) <- 'integrated'
macroglia <- ScaleData(object = macroglia)
macroglia <- RunPCA(macroglia, npcs = 20)
macroglia <- FindNeighbors(macroglia, dims = 1:8)
macroglia <- RunUMAP(object = macroglia, dims = 1:8)
macroglia <- FindClusters(macroglia, resolution = 0.3)
p1 <- DimPlot(macroglia, label = TRUE, label.size = 5)
p2 <- DimPlot(macroglia, group.by = 'macroglia_subcluster', label = TRUE)
p1 + p2



# Setting C1-C2-Endothelial Resolution parameter --------------------------

vascular <- readRDS(file = './data/vascular.rds')
DefaultAssay(vascular) <- 'integrated'
res <- seq(0.1, 0.8, 0.05)
for (i in 1:length(res)) {
  vascular <- FindClusters(vascular, resolution = res[i])
}
res_umap <- vector(mode = 'list', length = length(res))
for (i in 1:length(res_umap)) {
  res_umap[[i]] <- DimPlot(
    object = vascular, 
    group.by = paste0('integrated_snn_res.', res[i]),
    label = TRUE,
    label.size = 5
  ) +
    theme_bw() +
    NoLegend() + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(title = paste('Resolution:', res[i]))
}
res_umap <- cowplot::plot_grid(plotlist = res_umap, ncol = 5)
res_umap

# Set resolution and regenerate figures:
vascular <- readRDS(file = './data/vascular.rds')
DefaultAssay(vascular) <- 'integrated'
vascular <- FindClusters(object = vascular, resolution = 0.2)
vascular$vascular_subcluster <- plyr::mapvalues(
  x = vascular$integrated_snn_res.0.2,
  from = 0:7,
  to = c('C-Endothelial',
         'Tip Cell',
         'U-Vascular',
         'A-Endothelial',
         'V-Endothelial',
         'Fibroblast',
         'Pericyte',
         'VSMC')
)
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
