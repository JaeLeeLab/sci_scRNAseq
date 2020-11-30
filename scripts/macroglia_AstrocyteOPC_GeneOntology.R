
######## Pathway/GO comparison b/w Astrocytes and OPCs ########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('ggplot2')
require('scran')
require('dplyr')
require('org.Mm.eg.db')
require('topGO')
results_out <- './results/macroglia_AstrocyteOPC_GeneOntology/'
dir.create(path = results_out)

# Import seurat data
macroglia <- readRDS(file = './data/macroglia.rds')

# Import gene name conversion table
gene_convert <- read.table(file = './ref/gene_name_conversion.tsv',
                           header = TRUE, sep = '\t')

# Import wrapper functions for pathway analysis
source('./scripts/PathwayAnalysis_functions.R')


# Preset values for viz
macroglia_cols <- c('Ependymal-A' = '#800000',
                    'Ependymal-B' = '#e6194b',
                    'Astroependymal' = '#f58231',
                    'Astrocyte' = 'goldenrod',
                    'OPC-A' = '#3cb44b',
                    'OPC-B' = '#008080',
                    'Div-OPC' = '#4363d8',
                    'Pre-Oligo' = '#911eb4',
                    'Oligodendrocyte' = '#f032e6')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')









#  DE tests to identify changes in Astrocyte expression  ------------------------

DefaultAssay(macroglia) <- 'RNA'
Idents(macroglia) <- 'celltype'

astrocyte <- subset(macroglia, ident = 'Astrocyte')
Idents(astrocyte) <- 'time'
astrocyte[['RNAcorrected']] <- NULL
astrocyte[['integrated']] <- NULL

times <- levels(macroglia@meta.data[['time']])



# Use wilcox test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
astro_de_wilcox <- vector(mode = 'list', length = length(times) - 1)
names(astro_de_wilcox) <- paste(times[1:(length(times)-1)],
                                times[2:length(times)],
                                sep = '_')
for (t in 2:length(times)) {
  astro_de_wilcox[[t-1]] <- FindMarkers(
    object = astrocyte,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'wilcox'
  )
}
saveRDS(astro_de_wilcox,
        file = paste0(results_out, 'astrocyte_time_DE_wilcox.rds'))



# Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
astro_de_mast <- vector(mode = 'list', length = length(times) - 1)
names(astro_de_mast) <- paste(times[1:(length(times)-1)],
                              times[2:length(times)],
                              sep = '_')
for (t in 2:length(times)) {
  astro_de_mast[[t-1]] <- FindMarkers(
    object = astrocyte,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'MAST',
    latent.vars = 'nFeature_RNA'
  )
}
saveRDS(astro_de_mast,
        file = paste0(results_out, 'astrocyte_time_DE_mast.rds'))



# Use MAST with percent_rp regressed %%%%%%%%%%%%%%%%%%%%%%%%%%%%
astro_de_mast <- vector(mode = 'list', length = length(times) - 1)
names(astro_de_mast) <- paste(times[1:(length(times)-1)],
                              times[2:length(times)],
                              sep = '_')
for (t in 2:length(times)) {
  astro_de_mast[[t-1]] <- FindMarkers(
    object = astrocyte,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'MAST',
    latent.vars = c('nFeature_RNA', 'percent_rp')
  )
}
saveRDS(astro_de_mast,
        file = paste0(results_out, 'astrocyte_time_DE_mast_rp.rds'))





#  DE tests to identify changes in OPC expression  --------------------------

DefaultAssay(macroglia) <- 'RNA'
Idents(macroglia) <- 'celltype'

opc <- subset(macroglia, ident = 'OPC')
Idents(opc) <- 'time'
opc[['RNAcorrected']] <- NULL
opc[['integrated']] <- NULL

times <- levels(macroglia@meta.data[['time']])
# table(opc$sample_id, opc$time)
# Table of OPC counts by sample shows that data from uninj_sample1 and 
# uninj_sample2 are low. All other samples have sufficient number of cells (at 
# least 100). Because of unbalanced design, we do not use blocking factors.


# Use wilcox test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_wilcox <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_wilcox) <- paste(times[1:(length(times)-1)],
                                times[2:length(times)],
                                sep = '_')
for (t in 2:length(times)) {
  opc_de_wilcox[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'wilcox'
  )
}
saveRDS(opc_de_wilcox, file = paste0(results_out, 'opc_time_DE_wilcox.rds'))




# Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_mast <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_mast) <- paste(times[1:(length(times)-1)],
                              times[2:length(times)],
                              sep = '_')
for (t in 2:length(times)) {
  opc_de_mast[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'MAST',
    latent.vars = 'nFeature_RNA'
  )
}
saveRDS(opc_de_mast, file = paste0(results_out, 'opc_time_DE_mast.rds'))



# Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_mast <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_mast) <- paste(times[1:(length(times)-1)],
                            times[2:length(times)],
                            sep = '_')
for (t in 2:length(times)) {
  opc_de_mast[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'MAST',
    latent.vars = c('nFeature_RNA', 'percent_rp')
  )
}
saveRDS(opc_de_mast, file = paste0(results_out, 'opc_time_DE_mast_rp.rds'))





#  DE tests to identify changes in OPC-A expression  --------------------------

DefaultAssay(macroglia) <- 'RNA'
Idents(macroglia) <- 'macroglia_subcluster'

opc <- subset(macroglia, ident = 'OPC-A')
Idents(opc) <- 'time'
opc[['RNAcorrected']] <- NULL
opc[['integrated']] <- NULL

times <- levels(macroglia@meta.data[['time']])
# table(opc$sample_id, opc$time)
# Table of OPC counts by sample shows that data from uninj_sample1 and 
# uninj_sample2 are low. All other samples have sufficient number of cells (at 
# least 100). Because of unbalanced design, we do not use blocking factors.



# Use wilcox test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_wilcox <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_wilcox) <- paste(times[1:(length(times)-1)],
                              times[2:length(times)],
                              sep = '_')
for (t in 2:length(times)) {
  opc_de_wilcox[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'wilcox'
  )
}
saveRDS(opc_de_wilcox, file = paste0(results_out, 'opcA_time_DE_wilcox.rds'))



# Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_mast <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_mast) <- paste(times[1:(length(times)-1)],
                            times[2:length(times)],
                            sep = '_')
for (t in 2:length(times)) {
  opc_de_mast[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'MAST',
    latent.vars = 'nFeature_RNA'
  )
}
saveRDS(opc_de_mast, file = paste0(results_out, 'opcA_time_DE_mast.rds'))








#  DE tests to identify changes in non-OPC-B expression  --------------------------

DefaultAssay(macroglia) <- 'RNA'
Idents(macroglia) <- 'macroglia_subcluster'

opc <- subset(macroglia, ident = c('OPC-A','Div-OPC','Pre-Oligo'))
Idents(opc) <- 'time'
opc[['RNAcorrected']] <- NULL
opc[['integrated']] <- NULL

times <- levels(macroglia@meta.data[['time']])
# table(opc$sample_id, opc$time)
# Table of OPC counts by sample shows that data from uninj_sample1 and 
# uninj_sample2 are low. All other samples have sufficient number of cells (at 
# least 100). Because of unbalanced design, we do not use blocking factors.



# Use wilcox test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_wilcox <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_wilcox) <- paste(times[1:(length(times)-1)],
                              times[2:length(times)],
                              sep = '_')
for (t in 2:length(times)) {
  opc_de_wilcox[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'wilcox'
  )
}
saveRDS(opc_de_wilcox, file = paste0(results_out, 'opcNotB_time_DE_wilcox.rds'))



# Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_mast <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_mast) <- paste(times[1:(length(times)-1)],
                            times[2:length(times)],
                            sep = '_')
for (t in 2:length(times)) {
  opc_de_mast[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'MAST',
    latent.vars = 'nFeature_RNA'
  )
}
saveRDS(opc_de_mast, file = paste0(results_out, 'opcNotB_time_DE_mast.rds'))







#  DE tests to identify changes in OPC-B expression  --------------------------

DefaultAssay(macroglia) <- 'RNA'
Idents(macroglia) <- 'macroglia_subcluster'

opc <- subset(macroglia, ident = 'OPC-B')
Idents(opc) <- 'time'
opc[['RNAcorrected']] <- NULL
opc[['integrated']] <- NULL

times <- levels(macroglia@meta.data[['time']])[2:4]
# table(opc$sample_id, opc$time)
# Table of OPC counts by sample shows that data from uninj_sample1 and 
# uninj_sample2 are low. All other samples have sufficient number of cells (at 
# least 100). Because of unbalanced design, we do not use blocking factors.



# Use wilcox test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_wilcox <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_wilcox) <- paste(times[1:(length(times)-1)],
                              times[2:length(times)],
                              sep = '_')
for (t in 2:length(times)) {
  opc_de_wilcox[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'wilcox'
  )
}
saveRDS(opc_de_wilcox, file = paste0(results_out, 'opcB_time_DE_wilcox.rds'))



# Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_de_mast <- vector(mode = 'list', length = length(times) - 1)
names(opc_de_mast) <- paste(times[1:(length(times)-1)],
                            times[2:length(times)],
                            sep = '_')
for (t in 2:length(times)) {
  opc_de_mast[[t-1]] <- FindMarkers(
    object = opc,
    ident.1 = times[t-1],
    ident.2 = times[t],
    assay = 'RNA',
    slot = 'data',
    test.use = 'MAST',
    latent.vars = 'nFeature_RNA'
  )
}
saveRDS(opc_de_mast, file = paste0(results_out, 'opcB_time_DE_mast.rds'))







#  DE tests to identify differnces between OPC-A/B  --------------------------

DefaultAssay(macroglia) <- 'RNA'
Idents(macroglia) <- 'time'


# Use wilcox test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_opcB_de_srat <- FindMarkers(
  object = macroglia,
  subset.ident = '1dpi',
  group.by = 'macroglia_subcluster',
  ident.1 = 'OPC-A',
  ident.2 = 'OPC-B',
  assay = 'RNA',
  slot = 'data',
  logfc.threshold = 0
)
saveRDS(opcA_opcB_de_srat, 
        file = paste0(results_out, 'opcA_opcB_DE_wilcox_seurat_1dpi.rds'))


# Using scran blocking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_sce <- as.SingleCellExperiment(x = macroglia[,macroglia$macroglia_subcluster %in% c('OPC-A','OPC-B') & macroglia$time %in% c('1dpi')])
opcA_opcB_de_scran <- findMarkers(
  x = opc_sce,
  groups = colData(opc_sce)[['macroglia_subcluster']],
  test.type = 'wilcox',
  pval.type = 'all',
  block = colData(opc_sce)[['time']]
)
opcA_opcB_de_scran <- opcA_opcB_de_scran[c('OPC-A','OPC-B')]
saveRDS(opcA_opcB_de_scran,
        file = paste0(results_out, 'opcA_opcB_DE_wilcox_scranTimeBlock_1dpi.rds'))



# Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Idents(macroglia) <- 'time'
opcA_opcB_de_mast <- FindMarkers(
  object = macroglia,
  subset.ident = c('1dpi'),
  group.by = 'macroglia_subcluster',
  ident.1 = 'OPC-A',
  ident.2 = 'OPC-B',
  assay = 'RNA',
  slot = 'data',
  logfc.threshold = 0,
  test.use = 'MAST',
  latent.vars = c('nFeature_RNA')
)
saveRDS(opcA_opcB_de_mast, 
        file = paste0(results_out, 'opcA_opcB_DE_MAST_1dpi.rds'))



# Use MAST with percent_rp regressed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Idents(macroglia) <- 'time'
opcA_opcB_de_mast <- FindMarkers(
  object = macroglia,
  subset.ident = c('1dpi'),
  group.by = 'macroglia_subcluster',
  ident.1 = 'OPC-A',
  ident.2 = 'OPC-B',
  assay = 'RNA',
  slot = 'data',
  logfc.threshold = 0,
  test.use = 'MAST',
  latent.vars = c('nFeature_RNA','percent_rp')
)
saveRDS(opcA_opcB_de_mast, 
        file = paste0(results_out, 'opcA_opcB_DE_MAST_1dpi_rp.rds'))



# # Use MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Idents(macroglia) <- 'time'
# opcA_opcB_de_mast <- FindMarkers(
#   object = macroglia,
#   subset.ident = c('3dpi'),
#   group.by = 'macroglia_subcluster',
#   ident.1 = 'OPC-A',
#   ident.2 = 'OPC-B',
#   assay = 'RNA',
#   slot = 'data',
#   logfc.threshold = 0.25,
#   test.use = 'MAST',
#   latent.vars = c('nFeature_RNA')
# )
# saveRDS(opcA_opcB_de_mast, 
#         file = paste0(results_out, 'opcA_opcB_DE_MAST_3dpi.rds'))







# Astrocyte Gene Ontology analysis w/ "count" based enrichment tests -----------


# Set criteria for "interesting" genes vs not. Uses tests like the %%%%%%%%%%%%%%%%%%%%%%%
# hypergeometric test, Fisher's exact test, and binomial test.
astro_de_wilcox <- readRDS(file = paste0(results_out, 'astrocyte_time_DE_wilcox.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster 
sig_genes <- lapply(
  X = astro_de_wilcox,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
astro_go_wilcox <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = astro_go_wilcox, 
        file = paste0(results_out, 'Astrocyte_time_GOresults_wilcox.rds'))





# Alternatively, take significant results with from MAST DE %%%%%%%%%%%%%%%%%%%%%%%
# Take significant results with positive FC per macroglia subcluster
astro_de_mast <- readRDS(file = paste0(results_out, 'astrocyte_time_DE_mast.rds'))
# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)
sig_genes <- lapply(
  X = astro_de_mast,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < -0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
astro_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = astro_go_mast, 
        file = paste0(results_out, 'Astrocyte_time_GOresults_MAST.rds'))








# OPC Gene Ontology analysis w/ "count" based enrichment tests ----------------

# Set criteria for "interesting" genes vs not. Uses tests like the %%%%%%%%%%%%%%%%%%%%%%%
# hypergeometric test, Fisher's exact test, and binomial test.
opc_de_wilcox <- readRDS(file = paste0(results_out, 'opc_time_DE_wilcox.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- lapply(
  X = opc_de_wilcox,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opc_go_wilcox <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opc_go_wilcox, 
        file = paste0(results_out, 'OPC_time_GOresults_wilcox.rds'))





# Alternatively, take significant results with from MAST DE %%%%%%%%%%%%%%%%%%%%%%%
# Take significant results with positive FC per macroglia subcluster
opc_de_mast <- readRDS(file = paste0(results_out, 'opc_time_DE_mast.rds'))

sig_genes <- lapply(
  X = opc_de_mast,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opc_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opc_go_mast, 
        file = paste0(results_out, 'OPC_time_GOresults_MAST.rds'))








# OPC-A Gene Ontology analysis w/ "count" based enrichment tests ----------------

# Set criteria for "interesting" genes vs not. Uses tests like the %%%%%%%%%%%%%%%%%%%%%%%
# hypergeometric test, Fisher's exact test, and binomial test.
opcA_de_wilcox <- readRDS(file = paste0(results_out, 'opcA_time_DE_wilcox.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- lapply(
  X = opcA_de_wilcox,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcA_go_wilcox <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcA_go_wilcox, 
        file = paste0(results_out, 'opcA_time_GOresults_wilcox.rds'))





# Alternatively, take significant results with from MAST DE %%%%%%%%%%%%%%%%%%%%%%%
# Take significant results with positive FC per macroglia subcluster
opcA_de_mast <- readRDS(file = paste0(results_out, 'opcA_time_DE_mast.rds'))

sig_genes <- lapply(
  X = opcA_de_mast,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcA_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcA_go_mast, 
        file = paste0(results_out, 'opcA_time_GOresults_MAST.rds'))







# Not OPC-B Gene Ontology analysis w/ "count" based enrichment tests ----------------

# Set criteria for "interesting" genes vs not. Uses tests like the %%%%%%%%%%%%%%%%%%%%%%%
# hypergeometric test, Fisher's exact test, and binomial test.
opcNotB_de_wilcox <- readRDS(file = paste0(results_out, 'opcNotB_time_DE_wilcox.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- lapply(
  X = opcNotB_de_wilcox,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < -0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcNotB_go_wilcox <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcNotB_go_wilcox, 
        file = paste0(results_out, 'opcNotB_time_GOresults_wilcox.rds'))




# Alternatively, take significant results with from MAST DE %%%%%%%%%%%%%%%%%%%%%%%
# Take significant results with positive FC per macroglia subcluster
opcNotB_de_mast <- readRDS(file = paste0(results_out, 'opcNotB_time_DE_mast.rds'))

sig_genes <- lapply(
  X = opcNotB_de_mast,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < -0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcNotB_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcNotB_go_mast, 
        file = paste0(results_out, 'opcNotB_time_GOresults_MAST.rds'))






# OPC-B Gene Ontology analysis w/ "count" based enrichment tests ----------------

# Set criteria for "interesting" genes vs not. Uses tests like the
# hypergeometric test, Fisher's exact test, and binomial test.
opcB_de_wilcox <- readRDS(file = paste0(results_out, 'opcB_time_DE_wilcox.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- lapply(
  X = opcB_de_wilcox,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcB_go_wilcox <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcB_go_wilcox, 
        file = paste0(results_out, 'opcB_time_GOresults_wilcox.rds'))



# Alternatively, take significant results with from MAST DE
# Take significant results with positive FC per macroglia subcluster
opcB_de_mast <- readRDS(file = paste0(results_out, 'opcB_time_DE_mast.rds'))

sig_genes <- lapply(
  X = opcB_de_mast,
  FUN = function(x) {
    up <- x %>%
      dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
      dplyr::filter(p_val_adj < 1e-03) %>%
      dplyr::select(avg_logFC, p_val_adj)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcB_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcB_go_mast, 
        file = paste0(results_out, 'opcB_time_GOresults_MAST.rds'))





# OPC-A vs OPC-B Gene Ontology analysis w/ "count" based enrichment tests ---------------

# Set criteria for "interesting" genes vs not. Uses tests like the
# hypergeometric test, Fisher's exact test, and binomial test.
opcA_opcB_de_wilcox <- readRDS(file = paste0(results_out, 'opcA_opcB_DE_wilcox_seurat_1dpi.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- {
  x <- opcA_opcB_de_wilcox
  up <- x %>%
    dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  up <- factor(as.integer(gene_superset %in% rownames(up)))
  names(up) <- gene_superset
  names(up) <- plyr::mapvalues(
    x = names(up),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  down <- x %>%
    dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  down <- factor(as.integer(gene_superset %in% rownames(down)))
  names(down) <- gene_superset
  names(down) <- plyr::mapvalues(
    x = names(down),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
    all(grepl(pattern = 'ENSMUS', x = names(down)))
  print(paste('All genes mapped to Ensembl ID?:', sanity_check))
  up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
  down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
  list('up' = up, 'down' = down)
}

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcA_opcB_go_wilcox <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcA_opcB_go_wilcox,
        file = paste0(results_out, 'opcA_opcB_GOresults_wilcox_seurat_1dpi.rds'))



# Alternatively, take significant results with from scran blocked DE
# Take significant results with positive FC per macroglia subcluster
opcA_opcB_de_scran <- readRDS(file = paste0(results_out, 'opcA_opcB_DE_wilcox_scranTimeBlock.rds'))

sig_genes <- lapply(
  X = opcA_opcB_de_scran,
  FUN = function(x) {
    up <- x %>%
      as.data.frame() %>%
      dplyr::filter(summary.AUC > 0.5) %>% # log fold-change is positive
      dplyr::filter(FDR < 1e-03) %>%
      dplyr::select(summary.AUC, FDR)
    up <- factor(as.integer(gene_superset %in% rownames(up)))
    names(up) <- gene_superset
    names(up) <- plyr::mapvalues(
      x = names(up),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    down <- x %>%
      as.data.frame() %>%
      dplyr::filter(summary.AUC < 0.5) %>% # log fold-change is positive
      dplyr::filter(FDR < 1e-03) %>%
      dplyr::select(summary.AUC, FDR)
    down <- factor(as.integer(gene_superset %in% rownames(down)))
    names(down) <- gene_superset
    names(down) <- plyr::mapvalues(
      x = names(down),
      from = gene_convert[['mgi_symbol']],
      to = gene_convert[['ensembl_gene_id']],
      warn_missing = FALSE
    )
    sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
      all(grepl(pattern = 'ENSMUS', x = names(down)))
    print(paste('All genes mapped to Ensembl ID?:', sanity_check))
    up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
    down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
    return(list('up' = up, 'down' = down))
  }
)
sig_genes <- unlist(x = sig_genes, recursive = FALSE, use.names = TRUE)

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcA_opcB_go_scran <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcA_opcB_go_scran,
        file = paste0(results_out, 'opcA_opcB_GOresults_wilcox_scranTimeBlock.rds'))



# Set criteria for "interesting" genes vs not. Uses tests like the
# hypergeometric test, Fisher's exact test, and binomial test.
opcA_opcB_de_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_DE_mast_1dpi.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- {
  x <- opcA_opcB_de_mast
  up <- x %>%
    dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  up <- factor(as.integer(gene_superset %in% rownames(up)))
  names(up) <- gene_superset
  names(up) <- plyr::mapvalues(
    x = names(up),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  down <- x %>%
    dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  down <- factor(as.integer(gene_superset %in% rownames(down)))
  names(down) <- gene_superset
  names(down) <- plyr::mapvalues(
    x = names(down),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
    all(grepl(pattern = 'ENSMUS', x = names(down)))
  print(paste('All genes mapped to Ensembl ID?:', sanity_check))
  up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
  down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
  list('up' = up, 'down' = down)
}

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcA_opcB_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcA_opcB_go_mast,
        file = paste0(results_out, 'opcA_opcB_GOresults_mast_1dpi.rds'))


# Set criteria for "interesting" genes vs not. Uses tests like the
# hypergeometric test, Fisher's exact test, and binomial test.
opcA_opcB_de_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_DE_MAST_1dpi_rp.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- {
  x <- opcA_opcB_de_mast
  up <- x %>%
    dplyr::filter(avg_logFC > 0.5) %>% # log fold-change is positive
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  up <- factor(as.integer(gene_superset %in% rownames(up)))
  names(up) <- gene_superset
  names(up) <- plyr::mapvalues(
    x = names(up),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  down <- x %>%
    dplyr::filter(avg_logFC < -0.5) %>% # log fold-change is negative
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  down <- factor(as.integer(gene_superset %in% rownames(down)))
  names(down) <- gene_superset
  names(down) <- plyr::mapvalues(
    x = names(down),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
    all(grepl(pattern = 'ENSMUS', x = names(down)))
  print(paste('All genes mapped to Ensembl ID?:', sanity_check))
  up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
  down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
  list('up' = up, 'down' = down)
}

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcA_opcB_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcA_opcB_go_mast,
        file = paste0(results_out, 'opcA_opcB_GOresults_mast_1dpi_rp.rds'))


# Set criteria for "interesting" genes vs not. Uses tests like the
# hypergeometric test, Fisher's exact test, and binomial test.
opcA_opcB_de_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_DE_wilcox_mast_3dpi.rds'))

# Set of all tested genes used as universe for Gene Ontology
gene_superset <- rownames(macroglia)

# Take significant results with positive FC per macroglia subcluster
sig_genes <- {
  x <- opcA_opcB_de_mast
  up <- x %>%
    dplyr::filter(avg_logFC > 0) %>% # log fold-change is positive
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  up <- factor(as.integer(gene_superset %in% rownames(up)))
  names(up) <- gene_superset
  names(up) <- plyr::mapvalues(
    x = names(up),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  down <- x %>%
    dplyr::filter(avg_logFC < 0) %>% # log fold-change is negative
    dplyr::filter(p_val_adj < 1e-03) %>%
    dplyr::select(avg_logFC, p_val_adj)
  down <- factor(as.integer(gene_superset %in% rownames(down)))
  names(down) <- gene_superset
  names(down) <- plyr::mapvalues(
    x = names(down),
    from = gene_convert[['mgi_symbol']],
    to = gene_convert[['ensembl_gene_id']],
    warn_missing = FALSE
  )
  sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(up))) &&
    all(grepl(pattern = 'ENSMUS', x = names(down)))
  print(paste('All genes mapped to Ensembl ID?:', sanity_check))
  up <- up[grepl(pattern = 'ENSMUS', x = names(up))]
  down <- down[grepl(pattern = 'ENSMUS', x = names(down))]
  list('up' = up, 'down' = down)
}

# Run GO analysis
go_bp <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'BP'
)
go_mf <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'MF'
)
go_cc <- lapply(
  X = sig_genes,
  FUN = runGO,
  ontology = 'CC'
)
opcA_opcB_go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = opcA_opcB_go_mast,
        file = paste0(results_out, 'opcA_opcB_GOresults_mast_seurat_3dpi.rds'))




# Visualize Astrocyte GO results ------------------------------------------------


# Plot results by wilcox DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
astro_go_wilcox <- readRDS(file = paste0(results_out, 'Astrocyte_time_GOresults_wilcox.rds')) 

# Gather astrocyte GO results
astro_top_odds <- lapply(
  X = astro_go_wilcox,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 3) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        return(tmp[1:15,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(astro_top_odds)) {
  comp <- astro_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(astro_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  astro_top_odds[[p]] <- comp
}
astro_top_odds <- unlist(astro_top_odds, recursive = FALSE, use.names = FALSE)
astro_top_odds <- Reduce(f = rbind, x = astro_top_odds)
astro_top_odds[['comparison']] <- factor(
  x = astro_top_odds[['comparison']],
  levels = unique(astro_top_odds[['comparison']])
)
astro_top_odds[['tmp_id']] <- paste(
  astro_top_odds[['comparison']],
  astro_top_odds[['Term']],
  sep = '_'
)

astro_go_wilcox_df <- astro_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

astro_GO_wilcox_plot <- astro_go_wilcox_df  %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
astro_GO_wilcox_plot
ggsave(filename = paste0(results_out, 'astrocyte_time_GO_wilcox_plot.tiff'),
       plot = astro_GO_wilcox_plot, device = 'tiff', height = 9, width = 7.5)



# Plot results by mast DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
astro_go_mast <- readRDS(file = paste0(results_out, 'Astrocyte_time_GOresults_mast.rds')) 

# Gather astrocyte GO results
astro_top_odds <- lapply(
  X = astro_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          # dplyr::filter(Significant > 3) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        return(tmp[1:15,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(astro_top_odds)) {
  comp <- astro_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(astro_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  astro_top_odds[[p]] <- comp
}
astro_top_odds <- unlist(astro_top_odds, recursive = FALSE, use.names = FALSE)
astro_top_odds <- Reduce(f = rbind, x = astro_top_odds)
astro_top_odds[['comparison']] <- factor(
  x = astro_top_odds[['comparison']],
  levels = unique(astro_top_odds[['comparison']])
)
astro_top_odds[['tmp_id']] <- paste(
  astro_top_odds[['comparison']],
  astro_top_odds[['Term']],
  sep = '_'
)

astro_go_mast_df <- astro_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

astro_GO_mast_plot <- astro_go_mast_df %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
astro_GO_mast_plot
ggsave(filename = paste0(results_out, 'astrocyte_time_GO_mast_plot.tiff'),
       plot = astro_GO_mast_plot, device = 'tiff', height = 9, width = 7.5)







# Visualize OPC GO results ------------------------------------------------

# Plot results by wilcox DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_go_wilcox <- readRDS(file = paste0(results_out, 'opc_time_GOresults_wilcox.rds')) 

# Gather opc GO results
opc_top_odds <- lapply(
  X = opc_go_wilcox,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 3) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        return(tmp[1:15,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opc_top_odds)) {
  comp <- opc_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opc_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opc_top_odds[[p]] <- comp
}
opc_top_odds <- unlist(opc_top_odds, recursive = FALSE, use.names = FALSE)
opc_top_odds <- Reduce(f = rbind, x = opc_top_odds)
opc_top_odds[['comparison']] <- factor(
  x = opc_top_odds[['comparison']],
  levels = unique(opc_top_odds[['comparison']])
)
opc_top_odds[['tmp_id']] <- paste(
  opc_top_odds[['comparison']],
  opc_top_odds[['Term']],
  sep = '_'
)

opc_go_wilcox_df <- opc_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

opc_GO_wilcox_plot <- opc_go_wilcox_df  %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opc_GO_wilcox_plot
ggsave(filename = paste0(results_out, 'opc_time_GO_wilcox_plot.tiff'),
       plot = opc_GO_wilcox_plot, device = 'tiff', height = 9, width = 7.5)



# Plot results by mast DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opc_go_mast <- readRDS(file = paste0(results_out, 'opc_time_GOresults_mast.rds')) 

# Gather opc GO results
opc_top_odds <- lapply(
  X = opc_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 3) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        return(tmp[1:15,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opc_top_odds)) {
  comp <- opc_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opc_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opc_top_odds[[p]] <- comp
}
opc_top_odds <- unlist(opc_top_odds, recursive = FALSE, use.names = FALSE)
opc_top_odds <- Reduce(f = rbind, x = opc_top_odds)
opc_top_odds[['comparison']] <- factor(
  x = opc_top_odds[['comparison']],
  levels = unique(opc_top_odds[['comparison']])
)
opc_top_odds[['tmp_id']] <- paste(
  opc_top_odds[['comparison']],
  opc_top_odds[['Term']],
  sep = '_'
)

opc_go_mast_df <- opc_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(comparison, desc(log_pval))

opc_GO_mast_plot <- opc_go_mast_df  %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opc_GO_mast_plot
ggsave(filename = paste0(results_out, 'opc_time_GO_mast_plot.tiff'),
       plot = opc_GO_mast_plot, device = 'tiff', height = 9, width = 7.5)





# Visualize OPC-A GO results ------------------------------------------------

# Plot results by wilcox DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_go_wilcox <- readRDS(file = paste0(results_out, 'opcA_time_GOresults_wilcox.rds')) 

# Gather opcA GO results
opcA_top_odds <- lapply(
  X = opcA_go_wilcox,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          # dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:15,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcA_top_odds)) {
  comp <- opcA_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_top_odds[[p]] <- comp
}
opcA_top_odds <- unlist(opcA_top_odds, recursive = FALSE, use.names = FALSE)
opcA_top_odds <- Reduce(f = rbind, x = opcA_top_odds)
opcA_top_odds[['comparison']] <- factor(
  x = opcA_top_odds[['comparison']],
  levels = unique(opcA_top_odds[['comparison']])
)
opcA_top_odds[['tmp_id']] <- paste(
  opcA_top_odds[['comparison']],
  opcA_top_odds[['Term']],
  sep = '_'
)

opcA_go_wilcox_df <- opcA_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))
# opcA_go_wilcox_df$direction <- c('up','down')[as.numeric(grepl(pattern = 'up', x = as.character(opcA_go_wilcox_df$comparison)))]

opcA_GO_wilcox_plot <- opcA_go_wilcox_df %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_GO_wilcox_plot
ggsave(filename = paste0(results_out, 'opcA_time_GO_wilcox_plot.tiff'),
       plot = opcA_GO_wilcox_plot, device = 'tiff', height = 9, width = 7)



# Plot results by mast DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_go_mast <- readRDS(file = paste0(results_out, 'opcA_time_GOresults_mast.rds'))

# Gather opcA GO results
opcA_top_odds <- lapply(
  X = opcA_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          # dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:15,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcA_top_odds)) {
  comp <- opcA_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_top_odds[[p]] <- comp
}
opcA_top_odds <- unlist(opcA_top_odds, recursive = FALSE, use.names = FALSE)
opcA_top_odds <- Reduce(f = rbind, x = opcA_top_odds)
opcA_top_odds[['comparison']] <- factor(
  x = opcA_top_odds[['comparison']],
  levels = unique(opcA_top_odds[['comparison']])
)
opcA_top_odds[['tmp_id']] <- paste(
  opcA_top_odds[['comparison']],
  opcA_top_odds[['Term']],
  sep = '_'
)

opcA_go_mast_df <- opcA_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))
# opcA_go_mast_df$direction <- c('up','down')[as.numeric(grepl(pattern = 'up', x = as.character(opcA_go_mast_df$comparison)))]

opcA_GO_mast_plot <- opcA_go_mast_df %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_GO_mast_plot
ggsave(filename = paste0(results_out, 'opcA_time_GO_mast_plot.tiff'),
       plot = opcA_GO_mast_plot, device = 'tiff', height = 9, width = 7)





# Visualize Not OPC-B GO results ------------------------------------------------

# Plot results by wilcox DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcNotB_go_wilcox <- readRDS(file = paste0(results_out, 'opcNotB_time_GOresults_wilcox.rds')) 

# Gather opcNotB GO results
opcNotB_top_odds <- lapply(
  X = opcNotB_go_wilcox,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          # dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:15,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcNotB_top_odds)) {
  comp <- opcNotB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcNotB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcNotB_top_odds[[p]] <- comp
}
opcNotB_top_odds <- unlist(opcNotB_top_odds, recursive = FALSE, use.names = FALSE)
opcNotB_top_odds <- Reduce(f = rbind, x = opcNotB_top_odds)
opcNotB_top_odds[['comparison']] <- factor(
  x = opcNotB_top_odds[['comparison']],
  levels = unique(opcNotB_top_odds[['comparison']])
)
opcNotB_top_odds[['tmp_id']] <- paste(
  opcNotB_top_odds[['comparison']],
  opcNotB_top_odds[['Term']],
  sep = '_'
)

opcNotB_go_wilcox_df <- opcNotB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))
# opcNotB_go_wilcox_df$direction <- c('up','down')[as.numeric(grepl(pattern = 'up', x = as.character(opcNotB_go_wilcox_df$comparison)))]

opcNotB_GO_wilcox_plot <- opcNotB_go_wilcox_df %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcNotB_GO_wilcox_plot
ggsave(filename = paste0(results_out, 'opcNotB_time_GO_wilcox_plot.tiff'),
       plot = opcNotB_GO_wilcox_plot, device = 'tiff', height = 9, width = 7)



# Plot results by mast DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcNotB_go_mast <- readRDS(file = paste0(results_out, 'opcNotB_time_GOresults_mast.rds'))

# Gather opcNotB GO results
opcNotB_top_odds <- lapply(
  X = opcNotB_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          # dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:15,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcNotB_top_odds)) {
  comp <- opcNotB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcNotB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcNotB_top_odds[[p]] <- comp
}
opcNotB_top_odds <- unlist(opcNotB_top_odds, recursive = FALSE, use.names = FALSE)
opcNotB_top_odds <- Reduce(f = rbind, x = opcNotB_top_odds)
opcNotB_top_odds[['comparison']] <- factor(
  x = opcNotB_top_odds[['comparison']],
  levels = unique(opcNotB_top_odds[['comparison']])
)
opcNotB_top_odds[['tmp_id']] <- paste(
  opcNotB_top_odds[['comparison']],
  opcNotB_top_odds[['Term']],
  sep = '_'
)

opcNotB_go_mast_df <- opcNotB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))
# opcNotB_go_mast_df$direction <- c('up','down')[as.numeric(grepl(pattern = 'up', x = as.character(opcNotB_go_mast_df$comparison)))]

opcNotB_GO_mast_plot <- opcNotB_go_mast_df %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcNotB_GO_mast_plot
ggsave(filename = paste0(results_out, 'opcNotB_time_GO_mast_plot.tiff'),
       plot = opcNotB_GO_mast_plot, device = 'tiff', height = 9, width = 7)





# # Visualize OPC-B GO results ------------------------------------------------
# 
# # Plot results by wilcox DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# opcB_go_wilcox <- readRDS(file = paste0(results_out, 'opcB_time_GOresults_wilcox.rds')) 
# 
# # Gather opcB GO results
# opcB_top_odds <- lapply(
#   X = opcB_go_wilcox,
#   FUN = function(xx) {
#     tmp <- lapply(
#       X = xx,
#       FUN = function(yy) {
#         too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
#         yy[['GO_table']][['pvalue']][too_low] <- 1e-30
#         tmp <- yy[['GO_table']] %>%
#           dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
#           # dplyr::filter(Significant > 5) %>%
#           dplyr::filter(pvalue < 1e-03) %>%
#           dplyr::arrange(pvalue) %>%
#           dplyr::select(Term, pvalue, odds_ratio)
#         tmp <- tmp[1:15,]
#         tmp <- tmp[!is.na(tmp$Term),]
#         return(tmp)
#       }
#     )
#     names(tmp) <- names(xx)
#     return(tmp)
#   }
# )
# for (p in 1:length(opcB_top_odds)) {
#   comp <- opcB_top_odds[[p]]
#   for (c in 1:length(comp)) {
#     go_table <- comp[[c]]
#     go_table[['ontology']] <- names(opcB_top_odds)[p]
#     go_table[['comparison']] <- names(comp)[c]
#     comp[[c]] <- go_table
#   }
#   opcB_top_odds[[p]] <- comp
# }
# opcB_top_odds <- unlist(opcB_top_odds, recursive = FALSE, use.names = FALSE)
# opcB_top_odds <- Reduce(f = rbind, x = opcB_top_odds)
# opcB_top_odds[['comparison']] <- factor(
#   x = opcB_top_odds[['comparison']],
#   levels = unique(opcB_top_odds[['comparison']])
# )
# opcB_top_odds[['tmp_id']] <- paste(
#   opcB_top_odds[['comparison']],
#   opcB_top_odds[['Term']],
#   sep = '_'
# )
# 
# opcB_go_wilcox_df <- opcB_top_odds %>%
#   mutate('log_odds' = log2(odds_ratio)) %>%
#   mutate('log_pval' = -log10(pvalue)) %>%
#   arrange(desc(log_pval))
# # opcB_go_wilcox_df$direction <- c('up','down')[as.numeric(grepl(pattern = 'up', x = as.character(opcB_go_wilcox_df$comparison)))]
# 
# opcB_GO_wilcox_plot <- opcB_go_wilcox_df %>%
#   mutate('comparison_tmp' = plyr::mapvalues(
#     x = comparison,
#     from = c('1dpi_3dpi.down',
#              '3dpi_7dpi.down'),
#     to = c('Upreg. from 1dpi to 3dpi',
#            'Upreg. from 3dpi to 7dpi')
#   )) %>%
#   filter(grepl(pattern = 'down', x = comparison)) %>%
#   filter(ontology == 'BP') %>%
#   ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
#   geom_bar(fill = 'grey80', 
#            color = 'black', 
#            stat = 'identity') +
#   facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
#   scale_x_continuous(breaks = seq(0, 100, 5)) +
#   scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
#   ylab(label = 'GO Term') +
#   xlab(label = 'log10(p-value)') +
#   theme(panel.background = element_rect(fill = NA, color = 'black'),
#         panel.border = element_rect(fill = NA, color = 'black'),
#         plot.title = element_text(size = 16, color = 'black'),
#         strip.text = element_text(size = 12, color = 'black'),
#         strip.background = element_rect(color = 'black'),
#         axis.title = element_text(size = 14, color = 'black'),
#         axis.title.y = element_blank(),
#         axis.text = element_text(size = 12, color = 'black'),
#         legend.title = element_text(size = 14, color = 'black'),
#         legend.text = element_text(size = 12, color = 'black'))
# opcB_GO_wilcox_plot
# ggsave(filename = paste0(results_out, 'opcB_time_GO_wilcox_plot.tiff'),
#        plot = opcB_GO_wilcox_plot, device = 'tiff', height = 7, width = 7)
# 
# 
# 
# 
# 
# # Plot results by mast DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# opcB_go_mast <- readRDS(file = paste0(results_out, 'opcB_time_GOresults_mast.rds'))
# 
# # Gather opcB GO results
# opcB_top_odds <- lapply(
#   X = opcB_go_mast,
#   FUN = function(xx) {
#     tmp <- lapply(
#       X = xx,
#       FUN = function(yy) {
#         too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
#         yy[['GO_table']][['pvalue']][too_low] <- 1e-30
#         tmp <- yy[['GO_table']] %>%
#           dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
#           # dplyr::filter(Significant > 5) %>%
#           dplyr::filter(pvalue < 1e-03) %>%
#           dplyr::arrange(pvalue) %>%
#           dplyr::select(Term, pvalue, odds_ratio)
#         tmp <- tmp[1:15,]
#         tmp <- tmp[!is.na(tmp$Term),]
#         return(tmp)
#       }
#     )
#     names(tmp) <- names(xx)
#     return(tmp)
#   }
# )
# for (p in 1:length(opcB_top_odds)) {
#   comp <- opcB_top_odds[[p]]
#   for (c in 1:length(comp)) {
#     go_table <- comp[[c]]
#     go_table[['ontology']] <- names(opcB_top_odds)[p]
#     go_table[['comparison']] <- names(comp)[c]
#     comp[[c]] <- go_table
#   }
#   opcB_top_odds[[p]] <- comp
# }
# opcB_top_odds <- unlist(opcB_top_odds, recursive = FALSE, use.names = FALSE)
# opcB_top_odds <- Reduce(f = rbind, x = opcB_top_odds)
# opcB_top_odds[['comparison']] <- factor(
#   x = opcB_top_odds[['comparison']],
#   levels = unique(opcB_top_odds[['comparison']])
# )
# opcB_top_odds[['tmp_id']] <- paste(
#   opcB_top_odds[['comparison']],
#   opcB_top_odds[['Term']],
#   sep = '_'
# )
# 
# opcB_go_mast_df <- opcB_top_odds %>%
#   mutate('log_odds' = log2(odds_ratio)) %>%
#   mutate('log_pval' = -log10(pvalue)) %>%
#   arrange(desc(log_pval))
# # opcB_go_mast_df$direction <- c('up','down')[as.numeric(grepl(pattern = 'up', x = as.character(opcB_go_mast_df$comparison)))]
# 
# opcB_GO_mast_plot <- opcB_go_mast_df %>%
#   mutate('comparison_tmp' = plyr::mapvalues(
#     x = comparison,
#     from = c('1dpi_3dpi.down',
#              '3dpi_7dpi.down'),
#     to = c('Upreg. from 1dpi to 3dpi',
#            'Upreg. from 3dpi to 7dpi')
#   )) %>%
#   filter(grepl(pattern = 'down', x = comparison)) %>%
#   filter(ontology == 'BP') %>%
#   ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
#   geom_bar(fill = 'grey80', 
#            color = 'black', 
#            stat = 'identity') +
#   facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
#   scale_x_continuous(breaks = seq(0, 100, 5)) +
#   scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
#   ylab(label = 'GO Term') +
#   xlab(label = 'log10(p-value)') +
#   theme(panel.background = element_rect(fill = NA, color = 'black'),
#         panel.border = element_rect(fill = NA, color = 'black'),
#         plot.title = element_text(size = 16, color = 'black'),
#         strip.text = element_text(size = 12, color = 'black'),
#         strip.background = element_rect(color = 'black'),
#         axis.title = element_text(size = 14, color = 'black'),
#         axis.title.y = element_blank(),
#         axis.text = element_text(size = 12, color = 'black'),
#         legend.title = element_text(size = 14, color = 'black'),
#         legend.text = element_text(size = 12, color = 'black'))
# opcB_GO_mast_plot
# ggsave(filename = paste0(results_out, 'opcB_time_GO_mast_plot.tiff'),
#        plot = opcB_GO_mast_plot, device = 'tiff', height = 9, width = 7)





# Visualize OPC-A vs OPC-B GO results -----------------------------------------

# Plot results by wilcox DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_opcB_go_wilcox <- readRDS(file = paste0(results_out, 'opcA_opcB_GOresults_wilcox_seurat_1dpi.rds')) 

# Gather opcA_opcB GO results
opcA_opcB_top_odds <- lapply(
  X = opcA_opcB_go_wilcox,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          # dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        return(tmp[1:30,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcA_opcB_top_odds)) {
  comp <- opcA_opcB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_opcB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_opcB_top_odds[[p]] <- comp
}
opcA_opcB_top_odds <- unlist(opcA_opcB_top_odds, recursive = FALSE, use.names = FALSE)
opcA_opcB_top_odds <- Reduce(f = rbind, x = opcA_opcB_top_odds)
opcA_opcB_top_odds[['comparison']] <- factor(
  x = opcA_opcB_top_odds[['comparison']],
  levels = unique(opcA_opcB_top_odds[['comparison']])
)
opcA_opcB_top_odds[['tmp_id']] <- paste(
  opcA_opcB_top_odds[['comparison']],
  opcA_opcB_top_odds[['Term']],
  sep = '_'
)

opcA_opcB_go_wilcox_df <- opcA_opcB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

opcA_opcB_GO_wilcox_plot <- opcA_opcB_go_wilcox_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('OPC-A','OPC-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison + ontology, scales = 'free', drop = TRUE, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_opcB_GO_wilcox_plot
ggsave(filename = paste0(results_out, 'opcA_opcB_GO_wilcox_seurat_plot.tiff'),
       plot = opcA_opcB_GO_wilcox_plot, device = 'tiff', height = 10, width = 22)



# Plot results by scran DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_opcB_go_scran <- readRDS(file = paste0(results_out, 'opcA_opcB_GOresults_wilcox_scranTimeBlock.rds')) 

# Gather opcA_opcB GO results
opcA_opcB_top_odds <- lapply(
  X = opcA_opcB_go_scran,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          # dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        return(tmp[1:30,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
opcA_opcB_top_odds <- lapply(
  X = opcA_opcB_top_odds,
  FUN = function(xx) {
    return(xx[1:2])
  }
)
for (p in 1:length(opcA_opcB_top_odds)) {
  comp <- opcA_opcB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_opcB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_opcB_top_odds[[p]] <- comp
}
opcA_opcB_top_odds <- unlist(opcA_opcB_top_odds, recursive = FALSE, use.names = FALSE)
opcA_opcB_top_odds <- Reduce(f = rbind, x = opcA_opcB_top_odds)
opcA_opcB_top_odds[['comparison']] <- factor(
  x = opcA_opcB_top_odds[['comparison']],
  levels = unique(opcA_opcB_top_odds[['comparison']])
)
opcA_opcB_top_odds[['tmp_id']] <- paste(
  opcA_opcB_top_odds[['comparison']],
  opcA_opcB_top_odds[['Term']],
  sep = '_'
)

opcA_opcB_go_scran_df <- opcA_opcB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(comparison, desc(log_pval))

opcA_opcB_go_scran_plot <- opcA_opcB_go_scran_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('OPC-A.up','OPC-A.down'),
    to = c('OPC-A','OPC-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison + ontology, scales = 'free', drop = TRUE, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_opcB_go_scran_plot
ggsave(filename = paste0(results_out, 'opcA_opcB_GO_wilcox_scranTimeBlock_plot.tiff'),
       plot = opcA_opcB_go_scran_plot, device = 'tiff', height = 9, width = 7.5)



# Plot results by MAST DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_opcB_go_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_GOresults_mast_1dpi.rds')) 

# Gather opcA_opcB GO results
opcA_opcB_top_odds <- lapply(
  X = opcA_opcB_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          # dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:20,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcA_opcB_top_odds)) {
  comp <- opcA_opcB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_opcB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_opcB_top_odds[[p]] <- comp
}
opcA_opcB_top_odds <- unlist(opcA_opcB_top_odds, recursive = FALSE, use.names = FALSE)
opcA_opcB_top_odds <- Reduce(f = rbind, x = opcA_opcB_top_odds)
opcA_opcB_top_odds[['comparison']] <- factor(
  x = opcA_opcB_top_odds[['comparison']],
  levels = unique(opcA_opcB_top_odds[['comparison']])
)
opcA_opcB_top_odds[['tmp_id']] <- paste(
  opcA_opcB_top_odds[['comparison']],
  opcA_opcB_top_odds[['Term']],
  sep = '_'
)

opcA_opcB_go_mast_df <- opcA_opcB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

opcA_opcB_GO_mast_plot <- opcA_opcB_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('OPC-A','OPC-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison + ontology, scales = 'free_y', drop = TRUE, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_opcB_GO_mast_plot
ggsave(filename = paste0(results_out, 'opcA_opcB_GO_mast_seurat_plot_1dpi.tiff'),
       plot = opcA_opcB_GO_mast_plot, device = 'tiff', height = 10, width = 22)



# Plot results by MAST DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_opcB_go_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_GOresults_mast_1dpi_rp.rds'))

# Gather opcA_opcB GO results
opcA_opcB_top_odds <- lapply(
  X = opcA_opcB_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          # dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:20,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcA_opcB_top_odds)) {
  comp <- opcA_opcB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_opcB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_opcB_top_odds[[p]] <- comp
}
opcA_opcB_top_odds <- unlist(opcA_opcB_top_odds, recursive = FALSE, use.names = FALSE)
opcA_opcB_top_odds <- Reduce(f = rbind, x = opcA_opcB_top_odds)
opcA_opcB_top_odds[['comparison']] <- factor(
  x = opcA_opcB_top_odds[['comparison']],
  levels = unique(opcA_opcB_top_odds[['comparison']])
)
opcA_opcB_top_odds[['tmp_id']] <- paste(
  opcA_opcB_top_odds[['comparison']],
  opcA_opcB_top_odds[['Term']],
  sep = '_'
)

opcA_opcB_go_mast_df <- opcA_opcB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

opcA_opcB_GO_mast_plot <- opcA_opcB_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('OPC-A','OPC-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison + ontology, scales = 'free_y', drop = TRUE, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_opcB_GO_mast_plot
ggsave(filename = paste0(results_out, 'opcA_opcB_GO_mast_plot_1dpi_rp.tiff'),
       plot = opcA_opcB_GO_mast_plot, device = 'tiff', height = 10, width = 22)


# Plot results by MAST DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opcA_opcB_go_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_GOresults_mast_seurat_3dpi.rds')) 

# Gather opcA_opcB GO results
opcA_opcB_top_odds <- lapply(
  X = opcA_opcB_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          # dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:20,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcA_opcB_top_odds)) {
  comp <- opcA_opcB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_opcB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_opcB_top_odds[[p]] <- comp
}
opcA_opcB_top_odds <- unlist(opcA_opcB_top_odds, recursive = FALSE, use.names = FALSE)
opcA_opcB_top_odds <- Reduce(f = rbind, x = opcA_opcB_top_odds)
opcA_opcB_top_odds[['comparison']] <- factor(
  x = opcA_opcB_top_odds[['comparison']],
  levels = unique(opcA_opcB_top_odds[['comparison']])
)
opcA_opcB_top_odds[['tmp_id']] <- paste(
  opcA_opcB_top_odds[['comparison']],
  opcA_opcB_top_odds[['Term']],
  sep = '_'
)

opcA_opcB_go_mast_df <- opcA_opcB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

opcA_opcB_GO_mast_plot <- opcA_opcB_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('OPC-A','OPC-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison + ontology, scales = 'free_y', drop = TRUE, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_opcB_GO_mast_plot
ggsave(filename = paste0(results_out, 'opcA_opcB_GO_mast_seurat_plot_3dpi.tiff'),
       plot = opcA_opcB_GO_mast_plot, device = 'tiff', height = 10, width = 22)





# Extended Figure 13 formatting -------------------------------------------

astro_go_mast <- readRDS(file = paste0(results_out, 'Astrocyte_time_GOresults_mast.rds')) 

# Gather astrocyte GO results
astro_top_odds <- lapply(
  X = astro_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        return(tmp[1:15,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(astro_top_odds)) {
  comp <- astro_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(astro_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  astro_top_odds[[p]] <- comp
}
astro_top_odds <- unlist(astro_top_odds, recursive = FALSE, use.names = FALSE)
astro_top_odds <- Reduce(f = rbind, x = astro_top_odds)
astro_top_odds[['comparison']] <- factor(
  x = astro_top_odds[['comparison']],
  levels = unique(astro_top_odds[['comparison']])
)
astro_top_odds[['tmp_id']] <- paste(
  astro_top_odds[['comparison']],
  astro_top_odds[['Term']],
  sep = '_'
)

astro_go_mast_df <- astro_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

astro_GO_mast_plot <- astro_go_mast_df  %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 30)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))


opc_go_mast <- readRDS(file = paste0(results_out, 'OPC_time_GOresults_mast.rds')) 

# Gather opc GO results
opc_top_odds <- lapply(
  X = opc_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:15,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opc_top_odds)) {
  comp <- opc_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opc_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opc_top_odds[[p]] <- comp
}
opc_top_odds <- unlist(opc_top_odds, recursive = FALSE, use.names = FALSE)
opc_top_odds <- Reduce(f = rbind, x = opc_top_odds)
opc_top_odds[['comparison']] <- factor(
  x = opc_top_odds[['comparison']],
  levels = unique(opc_top_odds[['comparison']])
)
opc_top_odds[['tmp_id']] <- paste(
  opc_top_odds[['comparison']],
  opc_top_odds[['Term']],
  sep = '_'
)

opc_go_mast_df <- opc_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))
# opc_go_mast_df$direction <- c('up','down')[as.numeric(grepl(pattern = 'up', x = as.character(opc_go_mast_df$comparison)))]

opc_GO_mast_plot <- opc_go_mast_df %>%
  mutate('comparison_tmp' = plyr::mapvalues(
    x = comparison,
    from = c('Uninjured_1dpi.down',
             '1dpi_3dpi.down',
             '3dpi_7dpi.down'),
    to = c('Upreg. from Uninj to 1dpi',
           'Upreg. from 1dpi to 3dpi',
           'Upreg. from 3dpi to 7dpi')
  )) %>%
  filter(grepl(pattern = 'down', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_grid(comparison_tmp ~ ., scales = 'free_y', drop = TRUE, switch = 'y') +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 30)) +
  scale_y_discrete(labels = function(x) sub("[^_]*_[^_]*_", "", x),
                   position = 'right') +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black')) +
  scale_x_reverse(breaks = seq(0,30,5))

# Compile
tmp1 <- astro_GO_mast_plot + 
  labs(title = 'Astrocyte Gene Ontology') + 
  theme(plot.title = element_text(size = 18),
        plot.margin = margin(0, 5, 0, 0))
tmp2 <- opc_GO_mast_plot + 
  labs(title = 'OPC Gene Ontology') +
  theme(plot.title = element_text(size = 18, hjust = 1),
        plot.margin = margin(0, 0, 0, 5))
spacer <- ggplot() + theme_void()
exfig12 <- tmp1 + spacer + tmp2 + 
  patchwork::plot_layout(widths = c(1,0.15,1))
exfig12
ggsave(filename = paste0(results_out, 'Astrocyte_vs_OPC_GO_overTime.tiff'),
       plot = exfig12, device = 'tiff', height = 10, width = 15)





# Extended data figure 12 formatting -------------------------------------------

opcA_opcB_de_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_DE_mast_1dpi_rp.rds'))
opcA_opcB_go_mast <- readRDS(file = paste0(results_out, 'opcA_opcB_GOresults_mast_1dpi_rp.rds'))

# Gather opcA_opcB GO results
opcA_opcB_top_odds <- lapply(
  X = opcA_opcB_go_mast,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          dplyr::filter(Significant > 5) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          dplyr::arrange(pvalue) %>%
          # dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:20,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(opcA_opcB_top_odds)) {
  comp <- opcA_opcB_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(opcA_opcB_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  opcA_opcB_top_odds[[p]] <- comp
}
opcA_opcB_top_odds <- unlist(opcA_opcB_top_odds, recursive = FALSE, use.names = FALSE)
opcA_opcB_top_odds <- Reduce(f = rbind, x = opcA_opcB_top_odds)
opcA_opcB_top_odds[['comparison']] <- factor(
  x = opcA_opcB_top_odds[['comparison']],
  levels = unique(opcA_opcB_top_odds[['comparison']])
)
opcA_opcB_top_odds[['tmp_id']] <- paste(
  opcA_opcB_top_odds[['comparison']],
  opcA_opcB_top_odds[['Term']],
  sep = '_'
)

opcA_opcB_go_mast_df <- opcA_opcB_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

opcA_opcB_GO_mast_plot1 <- opcA_opcB_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('OPC-A (1dpi)','OPC-B (1dpi)')
  )) %>%
  filter(grepl(pattern = 'OPC-A \\(1dpi\\)', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison, scales = 'free_y', drop = TRUE, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 15)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 16, color = 'black', face = 'bold'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_opcB_GO_mast_plot2 <- opcA_opcB_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('OPC-A (1dpi)','OPC-B (1dpi)')
  )) %>%
  filter(grepl(pattern = 'OPC-B \\(1dpi\\)', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison, scales = 'free_y', drop = TRUE, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 15)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 16, color = 'black', face = 'bold'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
opcA_opcB_GO_mast_plot <- opcA_opcB_GO_mast_plot1 + opcA_opcB_GO_mast_plot2


opcA_opcB_volcano <- opcA_opcB_de_mast %>%
  mutate('gene' = rownames(.)) %>%
  mutate('avg_log2FC' = -log2(exp(x = avg_logFC))) %>%
  mutate('log_pval' = -log10(p_val_adj)) %>%
  mutate('label' = gene) %>%
  mutate('label' = ifelse(abs(avg_log2FC) <= log2(2), yes = NA, no = label)) %>%
  mutate('label' = ifelse(p_val_adj >= 1e-03, yes = NA, no = label)) %>%
  mutate('pt_size' = ifelse(is.na(label), yes = 1, no = 2)) %>%
  mutate('pt_color' = ifelse(pt_size == 1, yes = 'grey80', no = 'red')) %>%
  ggplot(mapping = aes(x = avg_log2FC, y = log_pval)) +
  geom_point(mapping = aes(color = pt_color, size = pt_size), alpha = 0.5) + 
  geom_vline(xintercept = log2(2), linetype = 'dashed', size = 1, color = 'indianred') +
  geom_vline(xintercept = -log2(2), linetype = 'dashed', size = 1, color = 'indianred') +
  geom_hline(yintercept = 3, linetype = 'dashed', size = 1, color = 'dodgerblue') +
  ggrepel::geom_text_repel(mapping = aes(label = label), 
                           size = 4, 
                           segment.alpha = 0.6,
                           force = 2) +
  scale_x_continuous(limits = c(-3,3),
                     labels = seq(-5, 5, 1),
                     breaks = seq(-5, 5, 1)) +
  xlab(label = 'Log2(fold-change)') +
  ylab(label = '-Log10(adjusted p-value)') +
  scale_size(range = c(1,3)) +
  scale_color_manual(values = c('grey80' = 'grey50', 'red' = 'red')) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black', size = 1),
        axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 16, color = 'black'),
        panel.grid = element_line(color = 'grey60', linetype = 'dotted'),
        legend.position = 'none')
opcA_opcB_volcano <- cowplot::plot_grid(NULL, opcA_opcB_volcano, NULL, ncol = 3, rel_widths = c(0.1, 1, 0.1))
opcAB_figure <- cowplot::plot_grid(opcA_opcB_GO_mast_plot, NULL, opcA_opcB_volcano, nrow = 3, rel_heights = c(1, 0.1, 1))
opcAB_figure
ggsave(filename = paste0(results_out, 'opcA_opcB_GO_mast_plot_1dpi_rp_summary.tiff'),
       plot = opcAB_figure, device = 'tiff', height = 10, width = 12)
