
####### Count-based quality control #######

# 1. library size
# 2. # unique genes
# 3. doublet detection
# 4. mitochondrial %
# 5. hemoglobin %
# 6. manual cluster + inspection


# Data Import -------------------------------------------------------------


# Load libraries and set directories
require('ggplot2')
require('Matrix')
require('dplyr')
require('reticulate')
data_in <- './data/filtered_feature_bc_matrix/'
data_out <- './data/QC_filtered_feature_bc_matrix/'
results_out <- './results/quality_control/'
ref_out <- './ref/'
dir.create(path = data_out)
dir.create(path = results_out)


# Import data
counts_in <- paste0(data_in, list.files(path = data_in))
counts <- vector(mode = 'list', length = length(counts_in))
for(ii in 1:length(counts_in)) {
  counts[[ii]] <- readRDS(file = counts_in[ii])
}
names(counts) <- gsub(pattern = paste(data_in, '.rds', sep = "|"), replacement = '', counts_in)


# For storing QC results for filtering downstream
metadata <- list()
metadata[['sample_id']] <- rep(x = gsub(pattern = 'filtered_feature_bc_matrix_',
                                        replacement = '', 
                                        x = names(counts)), 
                               times = sapply(X = counts, FUN = ncol))


# Experimental metadata ---------------------------------------------------


# time after injury
metadata[['time']] <- substr(x = metadata[['sample_id']],
                             start = 1,
                             stop = regexpr(pattern = '_', text = metadata[['sample_id']]) - 1)
metadata[['time']] <- plyr::mapvalues(x = metadata[['time']], from = 'uninj', to = 'Uninjured')
metadata[['time']] <- factor(x = metadata[['time']], levels = c('Uninjured','1dpi','3dpi','7dpi'))


# dissociation method
dissociationMethod <- data.frame(
  'sample_id' = c('uninj_sample1', 
                  'uninj_sample2', 
                  'uninj_sample3', 
                  '1dpi_sample1', 
                  '1dpi_sample2', 
                  '1dpi_sample3', 
                  '3dpi_sample1', 
                  '3dpi_sample2', 
                  '7dpi_sample1', 
                  '7dpi_sample2'),
  'dissociationMethod' = c('Standard',
                           'Standard',
                           'Enriched',
                           'Standard',
                           'Enriched',
                           'Enriched',
                           'Standard',
                           'Enriched',
                           'Standard',
                           'Enriched')
)
metadata[['dissociationMethod']] <- plyr::mapvalues(
  x = metadata[['sample_id']],
  from = dissociationMethod[['sample_id']],
  to = dissociationMethod[['dissociationMethod']]
)


# 10X chemistry
chemistry <- data.frame(
  'sample_id' = c('1dpi_sample1',
                  '1dpi_sample2',
                  '1dpi_sample3',
                  '3dpi_sample1',
                  '3dpi_sample2',
                  '7dpi_sample1',
                  '7dpi_sample2',
                  'uninj_sample1',
                  'uninj_sample2',
                  'uninj_sample3'),
  'chemistry' = c('v2',
                  'v2',
                  'v3',
                  'v2',
                  'v2',
                  'v2',
                  'v2',
                  'v2',
                  'v2',
                  'v3')
)
metadata[['chemistry']] <- plyr::mapvalues(
  x = metadata[['sample_id']],
  from = chemistry[['sample_id']],
  to = chemistry[['chemistry']]
)


# Library size QC ---------------------------------------------------------


# Assumption: Most cells of a sample are of acceptable quality.

# Compute size
library_size <- lapply(X = counts, FUN = Matrix::colSums)
metadata[['library_size']] <- do.call(c, c(library_size, use.names = FALSE)) # store results


# Compute MAD-based thresholds on log-library size
log_mad <- function(x, mad_dev = 3) {
  xmad <- mad(log10(x), constant = 1)
  lo <- 10^(median(log10(x)) - mad_dev * xmad)
  hi <- 10^(median(log10(x)) + mad_dev * xmad)
  return(c('low' = lo, 'high' = hi))
}
umi_cutoffs <- t(sapply(library_size, FUN = log_mad, mad_dev = 3))


# Filter by thresholds
pass_umi <- vector(mode = 'list', length = length(library_size))
names(pass_umi) <- names(library_size)
for(ii in 1:length(library_size)) {
  sample_id <- names(library_size)[ii]
  pass <- 
    library_size[[sample_id]] < umi_cutoffs[sample_id, 'high'] &
    library_size[[sample_id]] > umi_cutoffs[sample_id, 'low']
  pass_umi[[sample_id]] <- pass
}
metadata[['pass_umi']] <- do.call(c, c(pass_umi, use.names = FALSE)) # store results
pass_umi_table <- t(sapply(pass_umi, table))
colnames(pass_umi_table) <- c('fail','pass')
umi_cutoffs <- cbind(umi_cutoffs, pass_umi_table)
write.table(x = umi_cutoffs, file = paste0(results_out, 'UMI_madfilter.tsv'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE) # save results



# Unique genes QC ---------------------------------------------------------

# Compute size and thresholds
unique_genes <- function(x) return(Matrix::colSums(x > 0))
n_genes <- lapply(X = counts, FUN = unique_genes)
metadata[['n_genes']] <- do.call(c, c(n_genes, use.names = FALSE)) # store results
n_genes_cutoffs <- t(sapply(n_genes, FUN = log_mad, mad_dev = 3))


# Filter cells by threshold
pass_n_genes <- vector(mode = 'list', length = length(n_genes))
names(pass_n_genes) <- names(n_genes)
for(ii in 1:length(n_genes)) {
  sample_id <- names(n_genes)[ii]
  pass <- 
    n_genes[[sample_id]] < n_genes_cutoffs[sample_id, 'high'] &
    n_genes[[sample_id]] > n_genes_cutoffs[sample_id, 'low']
  pass_n_genes[[sample_id]] <- pass
}
metadata[['pass_n_genes']] <- do.call(c, c(pass_n_genes, use.names = FALSE)) # store results
pass_n_genes_table <- t(sapply(pass_n_genes, table))
colnames(pass_n_genes_table) <- c('fail','pass')
n_genes_cutoffs <- cbind(n_genes_cutoffs, pass_n_genes_table)
write.table(x = n_genes_cutoffs, file = paste0(results_out, 'n_genes_madfilter.tsv'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE) # save results



# Mitochondrial, ribosomal, and hemoglobin gene percentage -------------------

# Function to calculate % counts in feature set
percent_feature <- function(x, feature_pattern) {
  features <- grep(pattern = feature_pattern, x = rownames(x), value = TRUE)
  feature_subset <- Matrix::colSums(x[features,])
  library_size <- Matrix::colSums(x)
  feature_pct <- feature_subset / library_size
  return(feature_pct)
}


# Calculate mitochondrial gene %
percent_mt <- lapply(X = counts, FUN = percent_feature, feature_pattern = '^mt-')
metadata[['percent_mt']] <- do.call(c, c(percent_mt, use.names = FALSE)) # store results
percent_mt_cutoffs <- 0.2 # 20%
pass_percent_mt <- lapply(X = percent_mt, FUN = `<`, percent_mt_cutoffs)
metadata[['pass_percent_mt']] <- do.call(c, c(pass_percent_mt, use.names = FALSE)) # store results
mod_table <- function(x) table(factor(x, levels = c(FALSE, TRUE)))
pass_percent_mt_table <- t(sapply(X = pass_percent_mt, FUN = mod_table))
colnames(pass_percent_mt_table) <- c('fail','pass')
percent_mt_cutoffs <- cbind(percent_mt_cutoffs, pass_percent_mt_table)
write.table(x = percent_mt_cutoffs, file = paste0(results_out, 'percent_mitochondrial_flatfilter.tsv'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE) # save results


# Calculate ribosomal gene %
percent_rp <- lapply(X = counts, FUN = percent_feature, feature_pattern = '^Rp[ls]')
metadata[['percent_rp']] <- do.call(c, c(percent_rp, use.names = FALSE)) # store results
percent_rp_cutoffs <- t(sapply(X = percent_rp, FUN = log_mad, mad_dev = 4))
pass_percent_rp <- vector(mode = 'list', length = length(percent_rp))
names(pass_percent_rp) <- names(percent_rp)
for(ii in 1:length(percent_rp)) {
  sample_id <- names(percent_rp)[ii]
  pass <- 
    percent_rp[[sample_id]] < percent_rp_cutoffs[sample_id, 'high'] &
    percent_rp[[sample_id]] > percent_rp_cutoffs[sample_id, 'low']
  pass_percent_rp[[sample_id]] <- pass
}
metadata[['pass_percent_rp']] <- do.call(c, c(pass_percent_rp, use.names = FALSE)) # store results
pass_percent_rp_table <- t(sapply(pass_percent_rp, table))
colnames(pass_percent_rp_table) <- c('fail','pass')
percent_rp_cutoffs <- cbind(percent_rp_cutoffs, pass_percent_rp_table)
write.table(x = percent_rp_cutoffs, file = paste0(results_out, 'percent_riboprotein_madfilter.tsv'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE) # save results


# Calculate hemoglobin %
percent_hbb <- lapply(X = counts, FUN = percent_feature, feature_pattern = '^Hbb[-]')
metadata[['percent_hbb']] <- do.call(c, c(percent_hbb, use.names = FALSE))  # store results
percent_hbb_cutoffs <- 1/5000 # 1 in 5000 counts is Hbb
pass_percent_hbb <- lapply(X = percent_hbb, FUN = `<`, percent_hbb_cutoffs)
metadata[['pass_percent_hbb']] <- do.call(c, c(pass_percent_hbb, use.names = FALSE))  # store results
pass_percent_hbb_table <- t(sapply(X = pass_percent_hbb, FUN = table))
colnames(pass_percent_hbb_table) <- c('fail','pass')
percent_hbb_cutoffs <- cbind(percent_hbb_cutoffs, pass_percent_hbb_table)
write.table(x = percent_hbb_cutoffs, file = paste0(results_out, 'percent_hemoglobin_flatfilter.tsv'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE) # save results




# Doublet removal ---------------------------------------------------------

# Load multiplet rate data. Approximate doublet rates derived from 10X website:
# https://support.10xgenomics.com/permalink/3vzDu3zQjY0o2AqkkkI4CC

# DoubletRate_10x <- data.frame(
#   'Multiplet_rate' = c(0.004,
#                        0.008,
#                        0.016,
#                        0.023,
#                        0.031,
#                        0.039,
#                        0.046,
#                        0.054,
#                        0.061,
#                        0.069,
#                        0.076),
#   'nCells_Loaded' = c(870,
#                       1700,
#                       3500,
#                       5300,
#                       7000,
#                       8700,
#                       10500,
#                       12200,
#                       14000,
#                       15700,
#                       17400),
#   'nCells_Recovered' = c(500,
#                          1000,
#                          2000,
#                          3000,
#                          4000,
#                          5000,
#                          6000,
#                          7000,
#                          8000,
#                          9000,
#                          10000)
# )
# write.table(x = DoubletRate_10x, file = paste0(ref_out, 'DoubletRates_10X.tsv'),
#             sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

DoubletRate_10x <- read.table(file = paste0(ref_out, 'DoubletRates_10x.tsv'), header = TRUE)
# Build quick linear model to extrapolate to 30k cells
rate_prediction <- lm(formula = DoubletRate_10x$Multiplet_rate ~ DoubletRate_10x$nCells_Recovered)
x <- sort(c(seq(from = 0, by = 1000, length.out = 30),500))
y <- rate_prediction$coefficients[2]*x + rate_prediction$coefficients[1]
DoubletRate_10x_pred <- data.frame(
  'Multiplet_rate_pred' = y,
  'nCells_Recovered' = x
)


# Configs
scrub <- import(module = 'scrublet', convert = FALSE)
writeLines(text = str(py_config()), con = paste0(results_out, 'py_config.txt'))


# Run Scrublet (https://doi.org/10.1016/j.cels.2018.11.005)
doublet_results <- vector(mode = 'list', length = length(counts))
names(doublet_results) <- names(counts)
dir.create(path = paste0(results_out, 'scrublet_outs'))
for(ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  ncells_est <- round(ncol(counts[[sample_id]]), digits = -3)
  rate <- DoubletRate_10x_pred[['Multiplet_rate_pred']][which(DoubletRate_10x_pred[['nCells_Recovered']] == ncells_est)]
  scrublet_obj <- scrub$Scrublet(r_to_py(t(counts[[sample_id]]))$tocsc(),
                                 expected_doublet_rate = rate)
  scrublet_result <- py_capture_output(
    scrublet_obj$scrub_doublets(min_counts = 2,
                                min_cells = 3,
                                min_gene_variability_pctl = 85,
                                verbose = TRUE)
  )
  writeLines(text = c(sample_id, scrublet_result),
             con = paste0(results_out, 'scrublet_outs/', sample_id, '.txt'))
  message(paste('Done with', sample_id))
  doublet_results[[sample_id]][['Doublet_score']] <- py_to_r(scrublet_obj$doublet_scores_obs_)
  doublet_results[[sample_id]][['is_doublet']] <- 
    doublet_results[[sample_id]][['Doublet_score']] > py_to_r(scrublet_obj$threshold_)
  names(doublet_results[[sample_id]]) <- c('doublet_score', 'is_doublet')
  doublet_results[[sample_id]][['nCells_estimated']] <- ncells_est
  doublet_results[[sample_id]][['Multiplet_rate_10X']] <- rate
  doublet_results[[sample_id]][['Threshold_score']] <- py_to_r(scrublet_obj$threshold_)
  
  rm(scrublet_obj); gc()
}
doublet_scores <- lapply(doublet_results, FUN = `[[`, 'doublet_score')
metadata[['doublet_scores']] <- do.call(c, c(doublet_scores, use.names = FALSE)) # store results
is_doublet <- lapply(doublet_results, FUN = `[[`, 'is_doublet')
metadata[['is_doublet']]  <- do.call(c, c(is_doublet, use.names = FALSE)) # store results
doublet_table <- t(sapply(X = is_doublet, FUN = table))
colnames(doublet_table) <- c('singlet', 'multiplet')
doublet_table <- cbind('nCells_estimated' = sapply(doublet_results, FUN = `[[`, 'nCells_estimated'),
                       'Multiplet_rate_10X' = sapply(doublet_results, FUN = `[[`, 'Multiplet_rate_10X'),
                       'Threshold_score' = sapply(doublet_results, FUN = `[[`, 'Threshold_score'),
                       doublet_table)
write.table(x = doublet_table, file = paste0(results_out, 'doublet_call.tsv'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE) # save results




# Visualizations -----------------------------------------------------------


# ggplot params and setup
metadata_df <- data.frame(metadata)
sample_id <- gsub(pattern = 'filtered_feature_bc_matrix_',
                  replacement = '',
                  x = names(counts))
qc_theme <- theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = 'black'),
                  axis.text.y = element_text(size = 12, color = 'black'),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = 14, color = 'black'),
                  panel.background = element_rect(fill = NA, color = 'black'),
                  panel.border = element_rect(fill = NA, color = 'black'),
                  # panel.grid.major = element_line(color = 'grey', linetype = 'dashed'),
                  legend.position = 'none')


# Plot UMI quality control
umi_cutoffs_gg <- data.frame(umi_cutoffs[,c('low','high')])
umi_cutoffs_gg <- reshape2::melt(cbind(umi_cutoffs_gg, 'sample_id' = sample_id), id.vars = 'sample_id')
umi_plot <- metadata_df %>%
  ggplot(mapping = aes(x = sample_id, y = library_size)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  geom_boxplot(mapping = aes(fill = sample_id), outlier.size = 0, width = 0.2) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(10^seq(1,10,1), 3*10^seq(1,10,1)),
                     labels = scales::comma) +
  ylab(label = 'Library size (# UMI)') +
  labs(subtitle = 'Red lines label MAD thresholds.') +
  geom_point(data = umi_cutoffs_gg, mapping = aes(x = sample_id, y = value), 
             size = 10, pch = '-', color = 'indianred') +
  qc_theme
ggsave(filename = paste0(results_out, 'UMI_thresholds.tiff'), plot = umi_plot,
       device = 'tiff', height = 4, width = 1 + length(sample_id)/2.5)


# Plot unique gene quality control
n_genes_cutoffs_gg <- data.frame(n_genes_cutoffs[,c('low','high')])
n_genes_cutoffs_gg <- reshape2::melt(cbind(n_genes_cutoffs_gg, 'sample_id' = sample_id), id.vars = 'sample_id')
n_genes_plot <- metadata_df %>%
  ggplot(mapping = aes(x = sample_id, y = n_genes)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  geom_boxplot(mapping = aes(fill = sample_id), outlier.size = 0, width = 0.2) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(10^seq(1,10,1), 3*10^seq(1,10,1)),
                     labels = scales::comma) +
  ylab(label = '# Unique Genes') +
  labs(subtitle = 'Red lines label MAD thresholds.') +
  geom_point(data = n_genes_cutoffs_gg, mapping = aes(x = sample_id, y = value), 
             size = 10, pch = '-', color = 'indianred') +
  qc_theme
ggsave(filename = paste0(results_out, 'n_genes_thresholds.tiff'), plot = n_genes_plot,
       device = 'tiff', height = 4, width = 1 + length(sample_id)/2.5)


# Plot mitochondrial gene % quality control
percent_mt_plot <- metadata_df %>%
  ggplot(mapping = aes(x = sample_id, y = percent_mt)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  geom_boxplot(mapping = aes(fill = sample_id), outlier.size = 0, width = 0.2) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.001, 0.01, 0.03, 0.1, 0.2, 1, unique(percent_mt_cutoffs[,'percent_mt_cutoffs']))) +
  ylab(label = 'Mitochondrial [mt-] %') +
  geom_hline(mapping = aes(yintercept = unique(percent_mt_cutoffs[,1])),
             linetype = 'dashed',
             size = 1,
             color = 'indianred') +
  qc_theme
ggsave(filename = paste0(results_out, 'mitochondrial_percent_thresholds.tiff'),
       plot = percent_mt_plot, device = 'tiff', height = 4, width = 1 + length(sample_id)/2.5)


# Plot ribosomal gene % quality control
percent_rp_cutoffs_gg <- data.frame(percent_rp_cutoffs[,c('low','high')])
percent_rp_cutoffs_gg <- reshape2::melt(cbind(percent_rp_cutoffs_gg, 'sample_id' = sample_id), id.vars = 'sample_id')
percent_rp_plot <- metadata_df %>%
  ggplot(mapping = aes(x = sample_id, y = percent_rp)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  geom_boxplot(mapping = aes(fill = sample_id), outlier.size = 0, width = 0.2) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  ylab(label = 'Ribosomal [Rp] gene %') +
  labs(subtitle = 'Red lines label MAD thresholds.') +
  geom_point(data = percent_rp_cutoffs_gg, mapping = aes(x = sample_id, y = value), 
             size = 10, pch = '-', color = 'indianred') +
  qc_theme
ggsave(filename = paste0(results_out, 'ribosomal_percent_thresholds.tiff'),
       plot = percent_rp_plot, device = 'tiff', height = 4, width = 1 + length(sample_id)/2.5)


# Plot hemoglobin gene % quality control
percent_hbb_plot <- cbind(data.frame(percent_hbb_cutoffs), 'sample_id' = sample_id) %>%
  ggplot(mapping = aes(x = sample_id, y = fail)) +
  geom_bar(mapping = aes(fill = sample_id), stat = 'identity', color = 'black') +
  ylab(label = '# Hemoglobin(hi) cells') +
  labs(subtitle = 'Hemoglobin(hi) cell defined as cell\nwith greater than 1 in 5000 [Hbb-] counts ') +
  scale_y_continuous(breaks = seq(0, 10000, 50)) +
  qc_theme
ggsave(filename = paste0(results_out, 'hemoglobin_percent_thresholds.tiff'),
       plot = percent_hbb_plot, device = 'tiff', height = 4, width = 1 + length(sample_id)/2.5)


# Plot doublet score quality control
doublet_thresholds <- data.frame('Threshold_score' = sapply(X = doublet_results, FUN = `[[`, 'Threshold_score'),
                                 'sample_id' = sample_id)
doublet_plot <- metadata_df %>%
  ggplot(mapping = aes(x = sample_id, y = doublet_scores)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  scale_y_continuous(trans = 'log', breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1)) +
  ylab(label = 'Doublet score') + 
  geom_point(data = doublet_thresholds,
             mapping = aes(x = sample_id, y = Threshold_score),
             size = 10, pch = '-', color = 'indianred') +
  qc_theme
ggsave(filename = paste0(results_out, 'doublet_thresholds.tiff'), plot = doublet_plot,
       device = 'tiff', height = 4, width = 1 + length(sample_id)/2.5)


# Assemble figures
qc_plots <- cowplot::plot_grid(umi_plot, n_genes_plot, percent_mt_plot, percent_rp_plot, percent_hbb_plot, doublet_plot, ncol = 3)
ggsave(filename = paste0(results_out, 'QC_metrics_all.tiff'), plot = qc_plots,
       device = 'tiff', height = 6.5, width = 11)



# Filter out low quality cells --------------------------------------------

# Note: we do not filter on ribosomal gene % simply because it has not been
# performed regularly in other single-cell studies.
barcode_filter <- 
  metadata[['pass_umi']] &
  metadata[['pass_n_genes']] &
  metadata[['pass_percent_mt']] &
  # metadata[['pass_percent_rp']] &
  metadata[['pass_percent_hbb']] &
  !metadata[['is_doublet']]
table(barcode_filter)
good_barcodes <- names(barcode_filter[which(barcode_filter)])


filtered_counts <- vector(mode = 'list', length = length(counts))
names(filtered_counts) <- names(counts)
for(ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  good_sample_barcodes <- colnames(counts[[sample_id]]) %in% good_barcodes
  filtered_counts[[sample_id]] <- counts[[sample_id]][,good_sample_barcodes]
}
sum(sapply(filtered_counts, ncol))



# Save filtered count matrices --------------------------------------------

for(ii in 1:length(filtered_counts)) {
  sample_id <- names(filtered_counts)[ii]
  saveRDS(object = filtered_counts[[sample_id]],
          file = paste0(data_out, 'qc_', sample_id, '.rds'))
}
saveRDS(metadata, file = paste0(results_out, 'metadata.rds'))
rm(list = ls())




# For resubmission --------------------------------------------------------

qc_theme <- theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = 'black'),
                  axis.text.y = element_text(size = 12, color = 'black'),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = 14, color = 'black'),
                  panel.background = element_rect(fill = NA, color = 'black'),
                  panel.border = element_rect(fill = NA, color = 'black'),
                  # panel.grid.major = element_line(color = 'grey', linetype = 'dashed'),
                  legend.position = 'none')

umi_plot <- sci@meta.data %>%
  ggplot(mapping = aes(x = sample_id, y = library_size)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  geom_boxplot(mapping = aes(fill = sample_id), outlier.size = 0, width = 0.2) +
  scale_y_continuous(trans = 'log10',
                     limits = c(1000, 100000),
                     breaks = c(10^seq(1,10,1), 3*10^seq(1,10,1)),
                     labels = scales::comma) +
  ylab(label = 'Library size (# UMI)') +
  qc_theme

n_genes_plot <- sci@meta.data %>%
  ggplot(mapping = aes(x = sample_id, y = n_genes)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  geom_boxplot(mapping = aes(fill = sample_id), outlier.size = 0, width = 0.2) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(10^seq(1,10,1), 3*10^seq(1,10,1), 7*10^seq(1,10,1), 2*10^seq(1,10,1)),
                     labels = scales::comma) +
  ylab(label = '# Unique Genes') +
  qc_theme

percent_mt_plot <- sci@meta.data %>%
  ggplot(mapping = aes(x = sample_id, y = percent_mt)) +
  geom_violin(mapping = aes(fill = sample_id), scale = 'width') +
  geom_boxplot(mapping = aes(fill = sample_id), outlier.size = 0, width = 0.2) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.001, 0.01, 0.03, 0.1, 0.2, 1)) +
  ylab(label = 'Mitochondrial [mt-] %') +
  qc_theme

qc_plot <- umi_plot + n_genes_plot + percent_mt_plot
ggsave(filename = './results/revision_figures/qc_plot.tiff',
       plot = qc_plot, height = 4, width = 14)
