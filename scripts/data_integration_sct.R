

####### Data Integration and Batch Correction #######


# Setup -------------------------------------------------------------

# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
data_in <- './data/QC_filtered_feature_bc_matrix/'
data_out <- './data/data_integration/'
results_in <- './results/quality_control/'
results_out <- './results/data_integration/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = data_out)
dir.create(path = results_out)


# Import data ------------------------------------------------------------

# counts
counts_in <- paste0(data_in, list.files(path = data_in))
counts_in <- counts_in[grepl(pattern = 'qc_filtered', x = counts_in)]
counts <- vector(mode = 'list', length = length(counts_in))
for(ii in 1:length(counts_in)) {
  counts[[ii]] <- readRDS(file = counts_in[ii])
}
names(counts) <- gsub(pattern = paste(data_in, '.rds', sep = "|"), replacement = '', counts_in)
metadata <- readRDS(file = paste0(results_in, 'metadata.rds'))
metadata <- data.frame(metadata)

# Metadata prep
sample_metadata <- vector(mode = 'list', length = length(counts))
names(sample_metadata) <- names(counts)
for(ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  filtered_barcodes <- colnames(counts[[sample_id]])
  sample_metadata[[sample_id]] <- metadata[rownames(metadata) %in% filtered_barcodes,]
}

# Cell-cycle gene setup
firstup <- function(x) {x <- tolower(x); substr(x,1,1) <- toupper(substr(x,1,1)); return(x)}
s_genes <- firstup(cc.genes.updated.2019$s.genes)
g2m_genes <- firstup(cc.genes.updated.2019$g2m.genes)



# Data normalization ------------------------------------------------------

# Set unique gene names (i.e. remove duplicates)
remove_duplicated_genes <- function(x) {
  bad_gene <- duplicated(rownames(x)) | (nchar(rownames(x)) == 0)
  counts <- x[!bad_gene,]
  return(counts)
}
counts <- lapply(X = counts, FUN = remove_duplicated_genes)


# Initialize Seurat objects
for(ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  sample_name <- gsub(pattern = 'qc_filtered_feature_bc_matrix_',
                      replacement = '',
                      x = sample_id)
  counts[[sample_id]] <- counts[[sample_id]] %>%
    CreateSeuratObject(project = sample_name) %>%
    AddMetaData(metadata = sample_metadata[[sample_id]]) %>%
    NormalizeData() %>%
    CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
}

# Cell-cycle module scoring
counts <- lapply(X = counts, 
                 FUN = function(x) {
                   x$CC.Difference <- x$S.Score - x$G2M.Score
                   return(x)
                  })

# Variance-stabilizing normalization
counts <- lapply(X = counts, FUN = SCTransform, vars.to.regress = 'CC.Difference')



# Integration / Batch correction ------------------------------------------


# Note: Integration is memory-intensive. I ran these on Miami Project's Imaris
# workstation. Approx compute time ~1hr with 128Gb RAM and ~3.5GHz clock speed.

# Integration
sci_features <- SelectIntegrationFeatures(object.list = counts, nfeatures = 3000)
mito_hemo <- grepl(pattern = '^mt-|^Hbb-', x = sci_features)
sci_features <- sci_features[!mito_hemo]
counts <- PrepSCTIntegration(object.list = counts, anchor.features = sci_features, assay = 'SCT')
sci_anchors <- FindIntegrationAnchors(object.list = counts, 
                                      anchor.features = sci_features,
                                      normalization.method = 'SCT')
sci <- IntegrateData(anchorset = sci_anchors, normalization.method = 'SCT')


# # Less memory-intensive alternative
# memory.limit(size = 72000)
# sci_features <- SelectIntegrationFeatures(object.list = counts, nfeatures = 3000)
# mito_hemo <- grepl(pattern = '^mt-|^Hbb-', x = sci_features)
# sci_features <- sci_features[!mito_hemo]
# counts <- PrepSCTIntegration(object.list = counts, anchor.features = sci_features, assay = 'SCT')
# rna <- lapply(counts, FUN = GetAssay, assay = 'RNA') # extract/store RNA assay
# saveRDS(rna, paste0(data_out, 'rna_data.rds'))
# counts <- lapply(counts, FUN = function(x) {x[['RNA']] <- NULL; return(x)})
# counts <- lapply(counts, FUN = function(x) {
#   DietSeurat(x, counts = FALSE, data = TRUE, scale.data = TRUE, assays = 'SCT')
#   }
# )
# rm(rna);gc()
# library('future')
# plan("multiprocess", workers = 4)
# options(future.globals.maxSize = 4000 * 1024^2)
# sci_anchors <- FindIntegrationAnchors(object.list = counts, 
#                                       anchor.features = sci_features,
#                                       normalization.method = 'SCT')
# saveRDS(sci_anchors, file = paste0(data_out, 'sci_anchors.rds'))
# rm(counts); gc()
# sci <- IntegrateData(anchorset = sci_anchors, normalization.method = 'SCT')
# rm(sci_anchors); gc()
# rna <- readRDS(file = paste0(data_out, 'rna_data.rds'))
# rna <- merge(x = rna[[1]], rna[2:length(rna)])
# rna <- CreateAssayObject(counts = rna)
# sci[['RNA']] <- NULL
# sci[['RNA']] <- rna


# save integration features for reference
write.table(x = sci@assays$integrated@var.features,
            file = paste0(results_out, 'integration_variable_features.tsv'),
            sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)


# the final product
saveRDS(object = sci, file = './data/sci.rds')


