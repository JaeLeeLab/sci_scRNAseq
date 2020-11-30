
################# U-Vascular doublet analysis ##################


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
results_out <- './results/vascular_doublet_analysis/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)

sci <- readRDS(file = './data/sci.rds')
vascular <- readRDS(file = './data/vascular.rds')

doublet_rates_10X <- read.table(
  file = paste0(ref_in, 'DoubletRates_10X.tsv'),
  header = TRUE
)

DefaultAssay(sci) <- 'RNA'
DefaultAssay(vascular) <- 'RNA'

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

vascular_cols <- c('A-Endothelial' = '#800000',
                   'C1-Endothelial' = '#e6194b',
                   'C2-Endothelial' = '#f58231',
                   'V-Endothelial' = '#808000',
                   'Tip Cell' = '#3cb44b',
                   'Pericyte' = '#008080',
                   'VSMC' = '#000075',
                   'Fibroblast' = '#4363d8',
                   'U-Vascular' = '#f032e6')



# Doublet score plots -----------------------------------------------------

doublet_vln_sci <- FetchData(
  object = sci,
  vars = c('doublet_scores','sample_id','celltype')
) %>%
  mutate(doublet_scores = log2(doublet_scores*1000)) %>%
  arrange(doublet_scores) %>%
  ggplot(mapping = aes(x = sample_id, y = doublet_scores)) +
  geom_violin(scale = 'width') +
  ggbeeswarm::geom_quasirandom(
    mapping = aes(color = celltype),
    size = 0.2,
    alpha = 0.3,
    width = 0.5
) +
  ylab(label = 'Doublet Score') +
  scale_color_manual(values = vascular_cols) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 14, color = 'black', angle = 45, hjust = 1),
        legend.key = element_rect(fill = NA),
        legend.title = element_rect(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(
    title = 'Cell-type',
    override.aes = list(
      size = 5,
      alpha = 1
    )
  ))

doublet_vln_vascular <- FetchData(
  object = vascular,
  vars = c('doublet_scores','sample_id','vascular_subcluster')
) %>%
  mutate(doublet_scores = log2(doublet_scores*1000)) %>%
  arrange(doublet_scores) %>%
  ggplot(mapping = aes(x = sample_id, y = doublet_scores)) +
  geom_violin(scale = 'width') +
  ggbeeswarm::geom_quasirandom(
    mapping = aes(color = vascular_subcluster),
    size = 1,
    alpha = 1,
    width = 0.5
  ) +
  ylab(label = 'Doublet Score') +
  scale_color_manual(values = vascular_cols) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 14, color = 'black', angle = 45, hjust = 1),
        legend.key = element_rect(fill = NA),
        # legend.title = element_rect(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(
    title = 'Vascular\nsubcluster',
    override.aes = list(
      size = 5,
      alpha = 1
    )
  ))
doublet_vln_vascular


doublet_vln_vascular_byType <- FetchData(
  object = vascular,
  vars = c('doublet_scores','sample_id','vascular_subcluster')
) %>%
  # mutate(doublet_scores = log2(doublet_scores*1000)) %>%
  arrange(doublet_scores) %>%
  ggplot(mapping = aes(x = vascular_subcluster, y = doublet_scores, fill = vascular_subcluster)) +
  geom_violin(scale = 'width') +
  # ggbeeswarm::geom_quasirandom(
  #   mapping = aes(color = vascular_subcluster),
  #   size = 1,
  #   alpha = 1,
  #   width = 0.5
  # ) +
  ylab(label = 'Doublet Score') +
  scale_fill_manual(values = vascular_cols) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 14, color = 'black', angle = 45, hjust = 1),
        legend.key = element_rect(fill = NA),
        # legend.title = element_rect(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black')) +
  guides(fill = guide_legend(title = 'Vascular\nsubcluster'))
doublet_vln_vascular_byType
ggsave(filename = paste0(results_out, 'vascular_doublet_score_violin.tiff'),
       plot = doublet_vln_vascular_byType, device = 'tiff', height = 3, width = 6)


# Proportion analysis -----------------------------------------------------


# Load counts for doublet calls from QC to estimate expected multiplet rates
doublet_results <- read.table(
  file = './results/quality_control/doublet_call.tsv',
  header = TRUE
)
# Set sample_id column
doublet_results[['sample_id']] <- gsub(
  pattern = '[^*]+_matrix_',
  replacement = '',
  x = rownames(doublet_results)
)

# Load counts for cell calling from preprocessing to estimate total loaded cells
# as an upper bound on multiplets.
cellcall_results <- read.csv(
  file = './results/preprocessing/cellcallresults_table.csv',
  header = TRUE,
  row.names = 1
)
# Set sample_id column
cellcall_results[['sample_id']] <- gsub(
  pattern = '[^*]+_matrix_',
  replacement = '',
  x = rownames(cellcall_results)
)
cellcall_results <- cellcall_results[match(cellcall_results$sample_id, doublet_results$sample_id),]


# Sum retained and discarded cells from cell calling to estimate total # loaded
doublet_results$nCells_Recovered <- apply(
  X = cellcall_results, 
  MARGIN = 1, 
  FUN = function(x) return(as.numeric(x[1]) + as.numeric(x[2]))
)
doublet_results$discarded <- cellcall_results$Discard

# Estimate multiplet counts by multiplying rate and # cells loaded
doublet_results$expected_multiplet <- 
  doublet_results$nCells_Recovered * doublet_results$Multiplet_rate_10X

# Estimate vascular-vascular multiplet rate by multiplying the estimated total
# multiplet count by the proportion of vascular/sci cells ^ 2. Dont bother
# subsetting U-vascular from other vascular cells because this will be an upper
# bound estimate of total multiplets.
doublet_results$expected_vascular2_multiplet <- 
  doublet_results$expected_multiplet * (ncol(vascular)/ncol(sci))^2*2
doublet_results$expected_vascular2_multiplet_rate <- 
  doublet_results$expected_vascular2_multiplet / doublet_results$nCells_Recovered * 100
U_vascular_count <- table(vascular$sample_id, vascular$vascular_subcluster)
doublet_results$U_vascular <- 
  U_vascular_count[match(doublet_results$sample_id, rownames(U_vascular_count)), 'U-Vascular']
doublet_results$U_vascular_rate <- 
  doublet_results$U_vascular / doublet_results$nCells_Recovered * 100

# Bar graph of summary
uvascular_doublet_bar <- doublet_results[c('sample_id','expected_vascular2_multiplet_rate','U_vascular_rate')] %>%
  reshape2::melt(id.vars = c('sample_id')) %>%
  mutate(sample_id = factor(
    x = sample_id,
    levels = levels(sci$sample_id)
  )) %>%
  mutate(variable = plyr::mapvalues(
    x = variable,
    from = c('expected_vascular2_multiplet_rate','U_vascular_rate'),
    to = c('Expected vascular\ndoublet rate', 'Observed\nU-Vascular rate'),
  )) %>%
  ggplot(mapping = aes(x = sample_id, y = value)) +
  geom_bar(mapping = aes(fill = variable),
           color = 'black',
           stat = 'identity',
           position = 'dodge') +
  ylab(label = '% of loaded cells') +
  scale_y_continuous(breaks = seq(0,10,0.5)) +
  theme(plot.title = element_text(face = 'bold'),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12, color = 'black'),
        legend.key.height = unit(10, units = 'mm'),
        legend.position = 'right',
        legend.direction = 'vertical')
uvascular_doublet_bar
ggsave(filename = paste0(results_out, 'UVascular_doublet_rates.tiff'),
       plot = uvascular_doublet_bar, device = 'tiff', height = 3, width = 6)





# # Previous analysis -------------------------------------------------------
# 
# total_cells <- table(sci$sample_id)
# celltype_counts <- table(sci$celltype, sci$sample_id)
# 
# 
# # 10x estimates
# rates <- c('1dpi_1' = 0.061, '1dpi_2' = 0.046, '3dpi_1' = 0.076, '3dpi_2' = 0.069, '7dpi_1' = 0.076, '7dpi_2' = 0.069, 'uninj_1' = 0.023, 'uninj_2' = 0.008)
# # Scrublet estimates
# rates <- c('1dpi_1' = 0.054, '1dpi_2' = 0.048, '3dpi_1' = 0.083, '3dpi_2' = 0.087, '7dpi_1' = 0.105, '7dpi_2' = 0.080, 'uninj_1' = 0.03, 'uninj_2' = 0.011)
# estimated_doublet_count <- total_cells * rates
# filtered_count <- c('1dpi_1' = 128, '1dpi_2' = 127, '3dpi_1' = 384, '3dpi_2' = 422, '7dpi_1' = 397, '7dpi_2' = 278, 'uninj_1' = 48, 'uninj_2' = 6)
# 
# umyeloid_count <- table(vascular$subcluster, vascular$orig.ident)['U-Myeloid',]
# uvascular_count <- table(vascular$subcluster, vascular$orig.ident)['U-Vascular',]
# myeloid_rate <- colSums(celltype_counts[c('Neutrophil','Monocyte','Macrophage','Dendritic','Microglia','Div-Myeloid'),])/colSums(celltype_counts)
# vascular_rate <- colSums(celltype_counts[c('Fibroblast','Endothelial','Pericyte'),])/colSums(celltype_counts)
# 
# myeloid_vascular_rate <- myeloid_rate * vascular_rate
# vascular_vascular_rate <- vascular_rate * vascular_rate
# estimated_umyeloid_count <- myeloid_vascular_rate * (estimated_doublet_count - filtered_count)
# estimated_uvascular_count <- vascular_vascular_rate * (estimated_doublet_count - filtered_count)
# 
# umyeloid_rate <- umyeloid_count/colSums(celltype_counts[c('Fibroblast','Endothelial','Pericyte'),])
# uvascular_rate <- uvascular_count/colSums(celltype_counts[c('Fibroblast','Endothelial','Pericyte'),])
# estimated_umyeloid_rate <- estimated_umyeloid_count/colSums(celltype_counts[c('Fibroblast','Endothelial','Pericyte'),])
# estimated_uvascular_rate <- estimated_uvascular_count/colSums(celltype_counts[c('Fibroblast','Endothelial','Pericyte'),])
# 
# umyeloid_doublet <- data.frame(cbind('Expected' = estimated_umyeloid_rate, 'Observed' = umyeloid_rate)) %>%
#   tibble::rownames_to_column(var = 'Sample') %>%
#   mutate('Sample' = factor(Sample, levels = c('uninj_1','uninj_2','1dpi_1','1dpi_2','3dpi_1','3dpi_2','7dpi_1','7dpi_2'))) %>%
#   reshape2::melt() %>%
#   ggplot(mapping = aes(x = Sample, y = value)) + 
#   geom_bar(mapping = aes(fill = variable), position = 'dodge', stat = 'identity') +
#   labs(title = 'U-Myeloid doublet likelihood') +
#   # xlab(label = 'Sample') +
#   ylab(label = 'Myeloid-vascular\ndoublet rates') +
#   scale_y_continuous(expand = c(0.025,0), breaks = seq(0,1,0.02), limits = c(0, 0.12)) +
#   theme(plot.title = element_text(face = 'bold'),
#         panel.background = element_rect(fill = NA, color = 'black'),
#         panel.border = element_rect(fill = NA, color = 'black'),
#         legend.title = element_blank(),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12),
#         axis.title.x = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.title.y = element_text(size = 14))
# 
# uvascular_doublet <- data.frame(cbind('Expected' = estimated_uvascular_rate, 'Observed' = uvascular_rate)) %>%
#   tibble::rownames_to_column(var = 'Sample') %>%
#   mutate('Sample' = factor(Sample, levels = c('uninj_1','uninj_2','1dpi_1','1dpi_2','3dpi_1','3dpi_2','7dpi_1','7dpi_2'))) %>%
#   reshape2::melt() %>%
#   ggplot(mapping = aes(x = Sample, y = value)) + 
#   geom_bar(mapping = aes(fill = variable), position = 'dodge', stat = 'identity') +
#   labs(title = 'U-Vascular doublet likelihood') +
#   # xlab(label = 'Sample') +
#   ylab(label = 'Vascular-vascular\ndoublet rates') +
#   scale_y_continuous(expand = c(0.025,0), breaks = seq(0,1,0.02), limits = c(0,0.12)) +
#   theme(plot.title = element_text(face = 'bold'),
#         panel.background = element_rect(fill = NA, color = 'black'),
#         panel.border = element_rect(fill = NA, color = 'black'),
#         legend.title = element_blank(),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 14))
# 
# doublet_likelihood <- cowplot::plot_grid(umyeloid_doublet + theme(legend.position = 'none'),
#                                          uvascular_doublet + theme(legend.position = 'none'),
#                                          ncol = 1)
# doublet_likelihood <- cowplot::plot_grid(doublet_likelihood, 
#                                          cowplot::get_legend(umyeloid_doublet),
#                                          ncol = 2, rel_widths = c(1,0.2))
# ggsave(filename = paste0(results_outpath, 'vascular_doublet_rates.tiff'), plot =  doublet_likelihood, device = 'tiff', height = 6, width = 7.5)