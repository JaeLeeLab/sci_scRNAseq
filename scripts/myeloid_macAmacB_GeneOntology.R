
######## Pathway/GO comparison b/w Macrophage-A and Macrophage-B ########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('ggplot2')
require('scran')
require('dplyr')
# require('biomaRt')
require('org.Mm.eg.db')
require('topGO')
# require('grid')
# require('ComplexHeatmap')
# data_in <- './data/'
# data_out <- './data/data_integration/'
results_out <- './results/myeloid_macAmacB_GeneOntology/'
# ref_in <- './ref/'
# ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)

# Import seurat data
myeloid <- readRDS(file = './data/myeloid.rds')

# Import gene name conversion table
gene_convert <- read.table(file = './ref/gene_name_conversion.tsv',
                           header = TRUE, sep = '\t')

# Import wrapper functions for pathway analysis
source('./scripts/PathwayAnalysis_functions.R')


# Preset values for viz
myeloid_cols <- c('Neutrophil' = '#800000',    
                  'Monocyte' = '#9a6324',
                  'Macrophage-A' = '#e6194b',
                  'Macrophage-B' = '#f58231',
                  'BA-Macrophage' = '#CCCC00',
                  'Dendritic' = '#808000',
                  'Div-Myeloid' = '#3cb44b',
                  'Microglia-A' = '#008080',
                  'Microglia-B' = 'cyan3',
                  'Microglia-C' = '#000075',
                  'Div-Microglia' = '#4363d8',
                  'IFN-Myeloid' = '#911eb4')
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')






#  DE tests to determine enriched genes  ------------------------------------

DefaultAssay(myeloid) <- 'RNA'
Idents(myeloid) <- 'myeloid_subcluster'

mac_de_srat <- Seurat::FindMarkers(
  object = myeloid,
  ident.1 = 'Macrophage-A',
  ident.2 = 'Macrophage-B',
  assay = 'RNA',
  slot = 'data',
  test.use = 'wilcox',
  logfc.threshold = 0.25,
  min.pct = 0.1
)
saveRDS(mac_de_srat, file = paste0(results_out, 'macA_macB_DE_wilcox_seurat.rds'))


# # Subset data to test DE using scran::findMarkers. We use this function because
# # it allows for blocking factors. We use the wilcox test while blocking for 
# # dissociation protocol as a source of variation.
# mac <- subset(myeloid, idents = c('Macrophage-A','Macrophage-B'))
# mac_sce <- as.SingleCellExperiment(x = mac, assay = 'RNA')
# mac_sce$myeloid_subcluster <- factor(
#   x = mac_sce$myeloid_subcluster,
#   levels = c('Macrophage-A','Macrophage-B')
# )
# mac_de_wilcox <- findMarkers(
#   x = mac_sce,
#   groups = mac_sce$myeloid_subcluster,
#   # block = mac_sce$dissociationMethod,
#   pval.type = 'all',
#   test.type = 'wilcox'
# )
# saveRDS(
#   object = mac_de_wilcox,
#   file = paste0(results_out, 'macAmacB_wilcox_scranResults.rds')
# )



# Alternatively, use MAST as implemented in Seurat. Include the gene detection
# rate in MAST's hurdle model (equivalnet to cellular detection rate in MAST
# vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html)
mac_de_mast <- Seurat::FindMarkers(
  object = myeloid,
  ident.1 = 'Macrophage-A',
  ident.2 = 'Macrophage-B',
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = 'nFeature_RNA',
  logfc.threshold = 0.25,
  min.pct = 0.1
)
saveRDS(object = mac_de_mast, file = paste0(results_out, 'macA_macB_DE_mast.rds'))


# Given the enrichment of ribosomal protein gene expression, which dominates
# downstream GO analysis, we also try regressing out percent_rp.
mac_de_mast_rp <- Seurat::FindMarkers(
  object = myeloid,
  ident.1 = 'Macrophage-A',
  ident.2 = 'Macrophage-B',
  assay = 'RNA',
  slot = 'data',
  test.use = 'MAST',
  latent.vars = c('nFeature_RNA', 'percent_rp'),
  logfc.threshold = 0.25,
  min.pct = 0.1
)
saveRDS(object = mac_de_mast_rp, file = paste0(results_out, 'macA_macB_DE_mast_rp.rds'))


# Subset data to test DE using scran::findMarkers. We use this function because
# it allows for blocking factors. We use the wilcox test while blocking for 
# dissociation protocol as a source of variation.
mac <- subset(myeloid, idents = c('Macrophage-A','Macrophage-B'))
mac_sce <- as.SingleCellExperiment(x = mac, assay = 'RNA')
mac_sce$myeloid_subcluster <- factor(
  x = mac_sce$myeloid_subcluster,
  levels = c('Macrophage-A','Macrophage-B')
)
mac_de_scran <- findMarkers(
  x = mac_sce,
  groups = mac_sce$myeloid_subcluster,
  block = mac_sce$dissociationMethod,
  pval.type = 'all',
  test.type = 'wilcox'
)
saveRDS(
  object = mac_de_scran,
  file = paste0(results_out, 'macA_macB_DE_wilcox_scran_methodBlock.rds')
)


# Gene Ontology analysis w/ "count" based enrichment tests -------------------

# Set criteria for "interesting" genes vs not. Uses tests like the
# hypergeometric test, Fisher's exact test, and binomial test.
mac_de_wilcox <- readRDS(file = paste0(results_out, 'macA_macB_DE_wilcox_seurat.rds'))

# Set of all tested genes used as universe for Gene Ontology
# table(rownames(mac_de_wilcox[['Macrophage-A']]) %in% rownames(myeloid[['RNA']]@data))
gene_superset <- rownames(myeloid)

# Take significant results with positive FC per macrophage subcluster
sig_genes <- {
  x <- mac_de_wilcox
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
go_wilcox <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = go_wilcox, 
        file = paste0(results_out, 'Macrophage_GOresults_wilcox.rds'))



# Alternatively, take significant results with from MAST DE
mac_de_mast <- readRDS(file = paste0(results_out, 'macA_macB_DE_mast.rds'))
# Take significant results with positive FC per macrophage subcluster
sig_genes <- {
  x <- mac_de_mast
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
go_mast <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = go_mast, 
        file = paste0(results_out, 'Macrophage_GOresults_MAST.rds'))



# Alternatively, take significant results with from MAST DE with percent_rp 
# regression.
mac_de_mast_rp <- readRDS(file = paste0(results_out, 'macA_macB_DE_mast_rp.rds'))
gene_superset <- rownames(myeloid)

# Take significant results with positive FC per macrophage subcluster
sig_genes <- {
  x <- mac_de_mast_rp
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
go_mast_rp <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = go_mast_rp, 
        file = paste0(results_out, 'Macrophage_GOresults_MAST_rp.rds'))


# Alternatively, take significant results with from scran::findMarkers() with
# dissociation method blocked
mac_de_scran <- readRDS(file = paste0(results_out, 'macA_macB_DE_wilcox_scran_methodBlock.rds'))

sig_genes <- lapply(
  X = mac_de_scran,
  FUN = function(x) {
    up <- x %>%
      as.data.frame() %>%
      dplyr::filter(summary.AUC > 0.7) %>% # log fold-change is positive
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
      dplyr::filter(summary.AUC < 0.3) %>% # log fold-change is positive
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
mac_go_scran <- list(
  'BP' = go_bp,
  'MF' = go_mf,
  'CC' = go_cc
)
saveRDS(object = mac_go_scran,
        file = paste0(results_out, 'Macrophage_GOresults_scran_methodBlock.rds'))






# # Gene Ontology analysis w/ "rank" based enrichment tests -------------------
# 
# !!! NOT WORKING PROPERLY - LOW P-VALUES DESPITE NO SIGNIFICANT GENES COMPARED
# TO EXPECTED FREQUENCY !!!
# 
# # Method 2: Use output values of DE tests to weight genes. Uses enrichment tests
# # like Kolmogorov-Smirnov like tests (also known as GSEA).
# 
# # Set of all tested genes used as universe for Gene Ontology
# # table(rownames(mac_de_wilcox[['Macrophage-A']]) %in% rownames(myeloid[['RNA']]@data))
# gene_superset <- rownames(myeloid)
# 
# # Take significant results with good AUC per myeloid subcluster
# sig_genes <- lapply(
#   X = mac_de_wilcox,
#   FUN = function(x) {
#     tmp_df <- as.data.frame(x) %>%
#       dplyr::arrange(dplyr::desc(summary.AUC)) %>%
#       # dplyr::filter(FDR < 1e-03) %>%
#       # dplyr::filter(summary.AUC > 0.7) %>% # genes w/ greater concordance prob in macA
#       dplyr::select(FDR, summary.AUC)
#     all_sig <- all(tmp_df[['FDR']][tmp_df[['summary.AUC']] >= 0.7] < 1e-03)
#     print(paste('All genes with AUC > 0.7 have FDR < 1e-03:', all_sig))
#     tmp <- tmp_df[['summary.AUC']]
#     names(tmp) <- rownames(tmp_df)
#     names(tmp) <- plyr::mapvalues(
#       x = names(tmp),
#       from = gene_convert[['mgi_symbol']],
#       to = gene_convert[['ensembl_gene_id']],
#       warn_missing = FALSE
#     )
#     sanity_check <- all(grepl(pattern = 'ENSMUS', x = names(tmp)))
#     print(paste('All genes mapped to Ensembl ID?:', sanity_check))
#     tmp <- tmp[grepl(pattern = 'ENSMUS', x = names(tmp))]
#     return(tmp)
#   }
# )
# 
# sig_fxn <- function(x) return(x >= 0.7)
# go_bp <- lapply(
#   X = sig_genes,
#   FUN = runGO,
#   ontology = 'BP',
#   use_scores = TRUE,
#   score_fxn = sig_fxn
# )
# go_mf <- lapply(
#   X = sig_genes,
#   FUN = runGO,
#   ontology = 'MF',
#   use_scores = TRUE,
#   score_fxn = function(x) return(x >= 0.7)
# )
# go_cc <- lapply(
#   X = sig_genes,
#   FUN = runGO,
#   ontology = 'CC',
#   use_scores = TRUE,
#   score_fxn = function(x) return(x >= 0.7)
# )





# Visualize GO results ----------------------------------------------------

# Plot results by wilcox DE
mac_go_wilcox <- readRDS(file = paste0(results_out, 'Macrophage_GOresults_wilcox.rds'))

# Gather GO results
mac_top_odds <- lapply(
  X = mac_go_wilcox,
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
          # dplyr::arrange(pvalue) %>%
          dplyr::arrange(desc(odds_ratio)) %>%
          # dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[1:30,]
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp)
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
for (p in 1:length(mac_top_odds)) {
  comp <- mac_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(mac_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  mac_top_odds[[p]] <- comp
}
mac_top_odds <- unlist(mac_top_odds, recursive = FALSE, use.names = FALSE)
mac_top_odds <- Reduce(f = rbind, x = mac_top_odds)
mac_top_odds[['comparison']] <- factor(
  x = mac_top_odds[['comparison']],
  levels = unique(mac_top_odds[['comparison']])
)
mac_top_odds[['tmp_id']] <- paste(
  mac_top_odds[['comparison']],
  mac_top_odds[['Term']],
  sep = '_'
)

mac_go_wilcox_df <- mac_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

mac_GO_wilcox_plot <- mac_go_wilcox_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('Macrophage-A','Macrophage-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ ontology + comparison, scales = 'free', drop = TRUE, ncol = 2) +
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
mac_GO_wilcox_plot
ggsave(filename = paste0(results_out, 'mac_GO_wilcox_seurat_plot.tiff'),
       plot = mac_GO_wilcox_plot, device = 'tiff', height = 16, width = 15)



# Plot results by MAST DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mac_go_mast <- readRDS(file = paste0(results_out, 'Macrophage_GOresults_MAST.rds')) 

# Gather mac GO results
mac_top_odds <- lapply(
  X = mac_go_mast,
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
          # dplyr::arrange(pvalue) %>%
          dplyr::arrange(desc(odds_ratio)) %>%
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
for (p in 1:length(mac_top_odds)) {
  comp <- mac_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(mac_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  mac_top_odds[[p]] <- comp
}
mac_top_odds <- unlist(mac_top_odds, recursive = FALSE, use.names = FALSE)
mac_top_odds <- Reduce(f = rbind, x = mac_top_odds)
mac_top_odds[['comparison']] <- factor(
  x = mac_top_odds[['comparison']],
  levels = unique(mac_top_odds[['comparison']])
)
mac_top_odds[['tmp_id']] <- paste(
  mac_top_odds[['comparison']],
  mac_top_odds[['Term']],
  sep = '_'
)

mac_go_mast_df <- mac_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

mac_GO_mast_plot <- mac_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('Macrophage-A','Macrophage-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ ontology + comparison, scales = 'free_y', drop = TRUE, ncol = 2) +
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
mac_GO_mast_plot
ggsave(filename = paste0(results_out, 'mac_GO_mast_plot.tiff'),
       plot = mac_GO_mast_plot, device = 'tiff', height = 12, width = 15)




# Plot results by MAST and precent_rp regression DE %%%%%%%%%%%%%%%%%%%%%%%%%
mac_go_mast_rp <- readRDS(file = paste0(results_out, 'Macrophage_GOresults_MAST_rp.rds')) 

# Gather mac GO results
mac_top_odds <- lapply(
  X = mac_go_mast_rp,
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
for (p in 1:length(mac_top_odds)) {
  comp <- mac_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(mac_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  mac_top_odds[[p]] <- comp
}
mac_top_odds <- unlist(mac_top_odds, recursive = FALSE, use.names = FALSE)
mac_top_odds <- Reduce(f = rbind, x = mac_top_odds)
mac_top_odds[['comparison']] <- factor(
  x = mac_top_odds[['comparison']],
  levels = unique(mac_top_odds[['comparison']])
)
mac_top_odds[['tmp_id']] <- paste(
  mac_top_odds[['comparison']],
  mac_top_odds[['Term']],
  sep = '_'
)

mac_go_mast_rp_df <- mac_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

mac_go_mast_rp_plot <- mac_go_mast_rp_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('Macrophage-A','Macrophage-B')
  )) %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ ontology + comparison, scales = 'free', drop = TRUE, ncol = 2) +
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
mac_go_mast_rp_plot
ggsave(filename = paste0(results_out, 'mac_go_mast_rp_plot.tiff'),
       plot = mac_go_mast_rp_plot, device = 'tiff', height = 16, width = 18)




# Plot results by scran DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mac_go_scran <- readRDS(file = paste0(results_out, 'Macrophage_GOresults_scran_methodBlock.rds')) 

# Gather mac GO results
mac_top_odds <- lapply(
  X = mac_go_scran,
  FUN = function(xx) {
    tmp <- lapply(
      X = xx,
      FUN = function(yy) {
        too_low <- yy[['GO_table']][['pvalue']] == '< 1e-30'
        yy[['GO_table']][['pvalue']][too_low] <- 1e-30
        tmp <- yy[['GO_table']] %>%
          dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
          # dplyr::filter(Significant > 7) %>%
          dplyr::filter(pvalue < 1e-03) %>%
          # dplyr::arrange(pvalue) %>%
          dplyr::arrange(desc(odds_ratio)) %>%
          dplyr::select(Term, pvalue, odds_ratio)
        tmp <- tmp[!is.na(tmp$Term),]
        return(tmp[1:30,])
      }
    )
    names(tmp) <- names(xx)
    return(tmp)
  }
)
mac_top_odds <- lapply(
  X = mac_top_odds,
  FUN = function(xx) {
    return(xx[1:2])
  }
)
for (p in 1:length(mac_top_odds)) {
  comp <- mac_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(mac_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  mac_top_odds[[p]] <- comp
}
mac_top_odds <- unlist(mac_top_odds, recursive = FALSE, use.names = FALSE)
mac_top_odds <- Reduce(f = rbind, x = mac_top_odds)
mac_top_odds[['comparison']] <- factor(
  x = mac_top_odds[['comparison']],
  levels = unique(mac_top_odds[['comparison']])
)
mac_top_odds[['tmp_id']] <- paste(
  mac_top_odds[['comparison']],
  mac_top_odds[['Term']],
  sep = '_'
)

mac_go_scran_df <- mac_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(comparison, desc(log_pval))

mac_go_scran_plot <- mac_go_scran_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('Macrophage-A.up','Macrophage-A.down'),
    to = c('Macrophage-A','Macrophage-B')
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
mac_go_scran_plot


mac_go_scran_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('Macrophage-A.up','Macrophage-A.down'),
    to = c('Macrophage-A','Macrophage-B')
  )) %>%
  group_by(comparison, ontology) %>%
  top_n(n = 10, wt = log_pval) %>%
  ungroup() %>%
  # filter(grepl(pattern = 'down', x = comparison)) %>%
  # filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(mapping = aes(fill = ontology), 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison, scales = 'free', drop = TRUE, ncol = 3) +
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




# Figure for paper --------------------------------------------------------


# Plot results by MAST DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mac_go_mast <- readRDS(file = paste0(results_out, 'Macrophage_GOresults_MAST.rds')) 

# Gather mac GO results
mac_top_odds <- lapply(
  X = mac_go_mast,
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
          # dplyr::arrange(pvalue) %>%
          dplyr::arrange(desc(odds_ratio)) %>%
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
for (p in 1:length(mac_top_odds)) {
  comp <- mac_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(mac_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  mac_top_odds[[p]] <- comp
}
mac_top_odds <- unlist(mac_top_odds, recursive = FALSE, use.names = FALSE)
mac_top_odds <- Reduce(f = rbind, x = mac_top_odds)
mac_top_odds[['comparison']] <- factor(
  x = mac_top_odds[['comparison']],
  levels = unique(mac_top_odds[['comparison']])
)
mac_top_odds[['tmp_id']] <- paste(
  mac_top_odds[['comparison']],
  mac_top_odds[['Term']],
  sep = '_'
)

mac_go_mast_df <- mac_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

mac_GO_mast_plot1 <- mac_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('Macrophage-A','Macrophage-B')
  )) %>%
  filter(grepl(pattern = 'Macrophage-A', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison, scales = 'free_y', drop = TRUE, ncol = 2) +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 15)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 16, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))

mac_GO_mast_plot2 <- mac_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('Macrophage-A','Macrophage-B')
  )) %>%
  filter(grepl(pattern = 'Macrophage-B', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison, scales = 'free_y', drop = TRUE, ncol = 2) +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 15)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = 'log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        strip.text = element_text(size = 16, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
mac_GO_mast_plot <- mac_GO_mast_plot1 + mac_GO_mast_plot2
mac_GO_mast_plot
ggsave(filename = paste0(results_out, 'macrophage_GO_result.tiff'),
       plot = mac_GO_mast_plot, device = 'tiff', height = 4, width = 13)



# Figure for paper (functional annotation) ------------------------------------


# Plot results by MAST DE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mac_go_mast <- readRDS(file = paste0(results_out, 'Macrophage_GOresults_MAST.rds')) 

# Gather mac GO results
mac_top_odds <- lapply(
  X = mac_go_mast,
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
          # dplyr::arrange(pvalue) %>%
          dplyr::arrange(desc(odds_ratio)) %>%
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
for (p in 1:length(mac_top_odds)) {
  comp <- mac_top_odds[[p]]
  for (c in 1:length(comp)) {
    go_table <- comp[[c]]
    go_table[['ontology']] <- names(mac_top_odds)[p]
    go_table[['comparison']] <- names(comp)[c]
    comp[[c]] <- go_table
  }
  mac_top_odds[[p]] <- comp
}
mac_top_odds <- unlist(mac_top_odds, recursive = FALSE, use.names = FALSE)
mac_top_odds <- Reduce(f = rbind, x = mac_top_odds)
mac_top_odds[['comparison']] <- factor(
  x = mac_top_odds[['comparison']],
  levels = unique(mac_top_odds[['comparison']])
)
mac_top_odds[['tmp_id']] <- paste(
  mac_top_odds[['comparison']],
  mac_top_odds[['Term']],
  sep = '_'
)

mac_go_mast_df <- mac_top_odds %>%
  mutate('log_odds' = log2(odds_ratio)) %>%
  mutate('log_pval' = -log10(pvalue)) %>%
  arrange(desc(log_pval))

mac_GO_mast_plot1 <- mac_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('Chemotaxis-Inducing Mac',
           'Inflammatory Mac')
  )) %>%
  filter(grepl(pattern = 'Chemotaxis-Inducing Mac', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison, scales = 'free_y', drop = TRUE, ncol = 2,
             labeller = label_wrap_gen(width = 15)) +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 15)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = '-log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 14, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 14, color = 'black'))

mac_GO_mast_plot2 <- mac_go_mast_df %>%
  mutate('comparison' = plyr::mapvalues(
    x = comparison,
    from = c('up','down'),
    to = c('Chemotaxis-Inducing Mac',
           'Inflammatory Mac')
  )) %>%
  filter(grepl(pattern = 'Inflammatory Mac', x = comparison)) %>%
  filter(ontology == 'BP') %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ comparison, scales = 'free_y', drop = TRUE, ncol = 2,
             labeller = label_wrap_gen(width = 15)) +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(0, 15)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = '-log10(p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        plot.title = element_text(size = 14, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 14, color = 'black'))
mac_GO_mast_plot <- mac_GO_mast_plot1 + mac_GO_mast_plot2
mac_GO_mast_plot
ggsave(filename = paste0(results_out, 'macrophage_GO_result.tiff'),
       plot = mac_GO_mast_plot, device = 'tiff', height = 4, width = 14)



# Microglia GO ------------------------------------------------------------


mg_markers <- read.table(
  file = 'results/myeloid_annotation_markers/microglia_DE_wilcox.tsv', 
  sep = '\t',
  header = TRUE,
  row.names = 1
)
go_genes <- mg_markers %>%
  mutate(pct.diff = pct.1 - pct.2) %>%
  filter(pct.diff > 0) %>%
  group_by(cluster)
ensembl_convert <- read.table('ref/gene_name_conversion.tsv', sep = '\t', header = TRUE)
go_genes$ensembl <- plyr::mapvalues(
  x = go_genes$gene,
  from = ensembl_convert$mgi_symbol,
  to = ensembl_convert$ensembl_gene_id,
  warn_missing = FALSE
)
write.csv(
  x = go_genes,
  file = paste0(results_out, 'microglia_go_genes.csv')
)

# Use gProfiler, use output
mg_go_results <- readxl::read_xlsx(path = paste0(results_out, 'microglia_GO_terms.xlsx'))
h_mg <- mg_go_results[c(1,2)]
a_mg <- mg_go_results[c(4,5)]
b_mg <- mg_go_results[c(7,8)]
c_mg <- mg_go_results[c(10,11)]
mg_go_results <- list('Homeostatic Microglia' = h_mg, 
                      'Inflammatory Microglia' = a_mg, 
                      'Dividing Microglia' = b_mg,
                      'Migrating Microglia' = c_mg)
for (i in 1:length(mg_go_results)) {
  colnames(mg_go_results[[i]]) <- c('GO term', 'adj. p-value')
  mg_go_results[[i]]$Cluster = names(mg_go_results)[i]
  mg_go_results[[i]] <- mg_go_results[[i]][!apply(is.na(mg_go_results[[i]]), any, MARGIN = 1),]
}
mg_go_results <- Reduce(rbind, mg_go_results)
mg_go_results$log_pval <- -log10(mg_go_results$`adj. p-value`)
mg_go_results$tmp_id <- paste(
  mg_go_results$Cluster,
  mg_go_results$`GO term`,
  sep = '_'
)
mg_go_results$Cluster <- factor(
  x = mg_go_results$Cluster,
  levels = c('Homeostatic Microglia','Inflammatory Microglia','Dividing Microglia','Migrating Microglia')
)

mg_go_plot <- mg_go_results %>%
  group_by(Cluster) %>%
  top_n(n = 15, wt = log_pval) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = log_pval, y = reorder(tmp_id, log_pval))) +
  geom_bar(fill = 'grey80', 
           color = 'black', 
           stat = 'identity') +
  facet_wrap(. ~ Cluster, 
             scales = 'free_y', 
             drop = TRUE, 
             ncol = 2, 
             labeller = label_wrap_gen()) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_discrete(labels = function(x) sub("[^*_]+_", "", x)) +
  ylab(label = 'GO Term') +
  xlab(label = '-log10(adj. p-value)') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        strip.background = element_rect(color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'))
mg_go_plot
ggsave(filename = paste0(results_out, 'microglia_GO_plot.tiff'),
       plot = mg_go_plot, device = 'tiff', height = 8, width = 14)

