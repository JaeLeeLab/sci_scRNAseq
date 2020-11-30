


####### Data Integration and Batch Correction #######


# Setup -------------------------------------------------------------

# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
require('SingleCellExperiment')
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
names(counts) <- gsub(pattern = 'uninj', replacement = 'Uninjured', x = names(counts)) # Important because metadata contains "Uninjured"
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





# Batch Correction between replicates -------------------------------------

# We can correct for batch effects (primarily due to sequencing depth) using
# the "batchelor" package. The function below is a wrapper for 
# batchelor::rescaleBatches(), which performs a linear regression-based 
# correction (ie downsampling counts). STRONG assumption: composition of cells
# in datasets are the same. Correction can fail if composition varies.
# Inputs:
#   obj_list: list of dgCMatrix objects. These should be technical replicates of a 
#     given condition/group.
#   metadata_ls: list of data.frames containing metadata corresponding to each 
#     object in obj_list.
correct_replicates <- function(obj_list, metadata_ls) {
  
  # Identify shared genes across datasets and take common subset
  gene_universe <- Reduce(intersect, lapply(X = obj_list, FUN = rownames))
  gene_universe <- gene_universe[nchar(gene_universe) > 0] # to fix weird '' rowname error
  obj_list <- lapply(
    X = obj_list,
    FUN = function(x) {
      x <- x[gene_universe,]
      return(x)
    }
  )
  
  # Init SingleCellExperiment object
  obj_list_sce <- lapply(
    X = obj_list,
    FUN = function(x) {
      x <- SingleCellExperiment(assays = list('counts' = x))
      logcounts(x) <- Seurat::NormalizeData(object = counts(x), verbose = FALSE)
      return(x)
    }
  )
  
  # Gather metadata
  group_barcodes <- unlist(sapply(X = obj_list_sce, FUN = colnames), use.names = FALSE)
  group_metadata <- do.call(rbind, metadata_ls)
  
  # Extract raw count and cell-level meta data. Merge all into
  # single matrix. These values derived directly from Seurat objects.
  raw_counts <- do.call(cbind, lapply(X = obj_list_sce, FUN = counts))
  log_raw_counts <- do.call(cbind, lapply(X = obj_list_sce, FUN = logcounts))
  # cell_metadata <- do.call(rbind, lapply(X = obj_list_sce, FUN = colData))
  # cell_metadata <- data.frame(cell_metadata)
  
  # Perform linear regression-based batch correction (ie scale down counts).
  # Automatically generates single, merged matrix. Recalculate the corrected 
  # "raw counts".
  corrected_log <- batchelor::rescaleBatches(obj_list_sce, 
                                             log.base = exp(1),
                                             pseudo.count = 1,
                                             assay.type = 'logcounts')
  corrected_log <- Matrix::Matrix(data = assays(corrected_log)[['corrected']],
                                  sparse = TRUE)

  # Set corrected values as new Seurat Assay slot
  corrected_log <- CreateAssayObject(counts = corrected_log)
  out_seurat <- CreateSeuratObject(counts = raw_counts,
                                   assay = 'RNA')
  out_seurat[['RNAcorrected']] <- corrected_log
  
  # Import log counts and metadata. Return assembled Seurat object.
  slot(out_seurat[['RNA']], 'data') <- log_raw_counts
  out_seurat@meta.data <- cbind(out_seurat@meta.data, group_metadata)
  return(out_seurat)
}


# Setup structures by condition
inj_groups <- levels(metadata[['time']])
sci_corrected <- vector(mode = 'list', length = length(inj_groups))
names(sci_corrected) <- inj_groups

# Perform the correction
for(ii in 1:length(inj_groups)) {
  inj_time <- inj_groups[ii]
  inj_samples <- counts[grep(pattern = inj_time, x = names(counts))]
  inj_samples <- correct_replicates(
    obj_list = inj_samples,
    metadata_ls = sample_metadata[grep(pattern = inj_time, x = names(sample_metadata))]
  )
  sci_corrected[[inj_time]] <- inj_samples
  message(paste('Done with:', inj_time))
}

# clean up
rm(inj_time, inj_samples, counts); gc()



# Integration across conditions (ie time) -----------------------------------


memory.limit(64000)

# To identify shared cell-types/states across conditions, we perform data 
# integration according to Seurat's pipeline. We use RNA assay to identify 
# variable genes, but identify anchors with the RNAcorrected assay. 
sci_corrected <- lapply(
  X = sci_corrected,
  FUN = FindVariableFeatures,
  assay = 'RNA',
  nfeatures = 2000
)
sci_corrected <- lapply(
  X = sci_corrected,
  FUN = function(x) {
    slot(x[['RNAcorrected']], 'var.features') <- slot(x[['RNA']], 'var.features')
    return(x)
  }
)
sci_anchors <- FindIntegrationAnchors(
  object.list = sci_corrected, 
  assay = rep('RNAcorrected', length(sci_corrected)),
  normalization.method = 'LogNormalize'
)
rm(sci_corrected); gc()

# Integrate into single dataset
sci <- IntegrateData(sci_anchors)
DefaultAssay(sci) <- 'integrated'


# save integration features for reference
write.table(x = sci@assays$integrated@var.features,
            file = paste0(data_out, 'integration_variable_features.tsv'),
            sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)


# the final product
saveRDS(object = sci, file = paste0(data_out, 'sci_batchelor.rds'))


rm(list = ls()); gc()