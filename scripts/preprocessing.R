
###### Preprocessing #####


# Data Import -------------------------------------------------------------

# Load libraries and set directories
data_in <- './data/raw_feature_bc_matrix/'
data_out <- './data/filtered_feature_bc_matrix/'
results_out <- './results/preprocessing/'
ref_out <- './ref/'
dir.create(path = data_out)
dir.create(path = results_out)


# Import data
counts_in <- paste0(data_in, list.files(path = data_in))
# counts_in <- counts_in[grep(pattern = 'raw.+h5$', x = counts_in)]
counts <- vector(mode = 'list', length = length(counts_in))
for(ii in 1:length(counts_in)) {
  counts[[ii]] <- Seurat::Read10X_h5(file = counts_in[ii],
                                     use.names = FALSE,
                                     unique.features = TRUE)
}
names(counts) <- substr(x = counts_in, 
                        start = nchar(data_in)+1, 
                        stop = regexpr(pattern = '.h5$', text = counts_in)-1)


# Give cell barcodes unique names across samples (remove "-1" (flow cell #)).
for(ii in 1:length(counts)) {
  barcode <- substr(x = colnames(counts[[ii]]),
                    start = 1,
                    stop = 16) # Barcodes are 16nt long
  sample_id <- gsub(pattern = 'raw_feature_bc_matrix_', 
                    x = names(counts)[ii],
                    replacement = '')
  barcode <- paste(barcode, sample_id, sep = '_')
  colnames(counts[[ii]]) <- barcode
  rm(barcode, sample_id)
}



# Cell Calling ------------------------------------------------------------

# Given raw gene-count matrix, determine which barcodes (droplets) contain cells
# and which do not using total UMI rank.


# Determine barcode ranks by total UMI counts using DropletUtils package.
# Params:
#   min_UMI: minimum UMI count below which presumed to be empty droplet.
#   fit_bounds: boundary UMI values to fit curve for estimating inflection/knee points
min_UMI <- 1000
fit_bounds <- c(min_UMI, 1e6)
knees <- vector(mode = 'numeric', length = length(counts))
names(knees) <- names(counts)
inflections <- vector(mode = 'numeric', length = length(counts))
names(inflections) <- names(counts)
for(ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  ranks <- DropletUtils::barcodeRanks(m = counts[[sample_id]],
                                      lower = min_UMI,
                                      fit.bounds = fit_bounds)
  knees[sample_id] <- ranks@metadata$knee # UMI value for knee 
  inflections[sample_id] <- ranks@metadata$inflection # UMI value for inflection
  rm(ranks)
}


# Determine which barcodes (droplets) contain cells and which are empty using 
# emptyDrops() function from DropletUtils package.
# Params:
#   drops_fdr: false-discovery rate for statistical test of significant deviation against ambient mRNA profile
drops_fdr <- 0.001
drops <- vector(mode = 'list', length = length(counts))
names(drops) <- names(counts)
for(ii in 1:length(drops)) {
  sample_id <- names(drops)[ii]
  drops[[ii]] <- DropletUtils::emptyDrops(m = counts[[sample_id]],
                                          retain = knees[sample_id],
                                          lower = inflections[sample_id])
  message(paste('Finished with', sample_id))
}
for(ii in 1:length(drops)) {
  sample_id <- names(drops)[ii]
  write.table(x = drops[[sample_id]],
              file = paste0(results_out, 'EmptyDropsOutput_', sample_id, '.tsv'),
              sep = '\t', quote = FALSE, col.names = NA, row.names = TRUE)
}
# note: v3 chemistry attempts ~6*10e6 drops, while v2 attempts ~7*10e5 drops


# Plot barcode-ranks and label which droplets pass cell call
nc <- ceiling(sqrt(length(drops)))
nr <- ceiling(length(drops)/nc)
tiff(filename = paste0(results_out, 'barcoderank.tiff'), height = 4*nr, width = 4*nc, units = 'in', res = 720)
par(mfrow = c(nr, nc))
for(ii in 1:length(drops)) {
  sample_id <- names(drops)[ii]
  drops[[sample_id]]$Rank <- rank(-drops[[sample_id]]$Total)
  drops[[sample_id]]$Retain <- (drops[[sample_id]]$FDR < drops_fdr)
  drops[[sample_id]]$Col <- 'black'
  drops[[sample_id]]$Col[drops[[sample_id]]$Retain] <- 'red'
  tmp_uniq <- !duplicated(drops[[sample_id]]$Rank)
  
  # plot points
  plot(x = drops[[sample_id]][['Rank']][tmp_uniq],
       y = drops[[sample_id]][['Total']][tmp_uniq],
       col = drops[[ii]]$Col[tmp_uniq],
       log = 'xy',
       main = sample_id, xlab = 'Rank', ylab = 'Total UMI',
       cex.lab = 1.4, cex.main = 1.4, cex.axis = 1.2)
  # horizontal lines for cutoffs
  abline(h = knees[sample_id], col = 'indianred', lwd = 1.5, lty = 'dashed')
  abline(h = inflections[sample_id], col = 'dodgerblue', lwd = 1.5, lty = 'dashed')
  # line and point color legends
  legend(x = 1, y = 10, cex = 1.2, 
         legend = c('Knee','Inflection'), col = c('indianred','dodgerblue'),
         lty = 'dashed', lwd = 1.5, text.font = 20)
  legend(legend = c('Retain', 'Discard'), col = c('red','black'),
         x = 1, y = 100, cex = 1.2, pch = c(1,1))
  # Label exact values
  knee_label <- c('x' = 10^(0.8*max(log10(nrow(drops[[sample_id]])))),
                  'y' = unname(knees[sample_id]))
  inflection_label <- c('x' = 10^(0.8*max(log10(nrow(drops[[sample_id]])))),
                        'y' = unname(inflections[sample_id]))
  text(x = knee_label['x'], y = knee_label['y'], labels = knees[sample_id],
       col = 'indianred', cex = 1.4, adj = c(0, -1.25))
  text(x = inflection_label['x'], y = inflection_label['y'], labels = inflections[sample_id],
       col = 'dodgerblue', cex = 1.4, adj = c(0, 1.25))
}
dev.off()


# Summary table of cell calling results
retain <- do.call(rbind, lapply(X = drops, FUN = function(xx) table(xx$Retain)))
values <- matrix(data = c(unlist(x = knees), unlist(x = inflections)), ncol = 2)
results <- cbind(retain, values)
colnames(x = results) <- c('Discard', 'Retain', 'Knee', 'Inflection')
write.table(x = results, file = paste0(results_out, 'cellcallresults_table.tsv'),
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)




# Gene nomenclature conversion --------------------------------------------


# Use latest Ensembl version (101 as of 2020/10/08) for all genes: 
# In CellRanger v4, uses Ensembl 98 and not 101. But there is an error with 
# most querying archived versions using biomaRt. 
# https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_2020A
ensembl_ids <- unique(unlist(sapply(counts, rownames)))
ensembl <- biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
                               # version = 98,
                               dataset = 'mmusculus_gene_ensembl')
fields <- c('mgi_symbol','ensembl_gene_id', 'external_gene_name')
convert_table <- biomaRt::getBM(attributes = fields,
                                filters = 'ensembl_gene_id',
                                values = ensembl_ids,
                                mart = ensembl)
# Not all 10X Ensembl IDs found in biomaRt query.
length(ensembl_ids[!ensembl_ids %in% convert_table$ensembl_gene_id]) # 579

# Save
write.table(file = paste0(ref_out, 'gene_name_conversion.tsv'), x = convert_table,
            quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# Convert names (will add convert_table data into Seurat metadata slots)
convert_names <- function(x) {
  rownames(x) <- plyr::mapvalues(x = rownames(x),
                                 from = convert_table$ensembl_gene_id,
                                 to = convert_table$mgi_symbol,
                                 warn_missing = FALSE)
  return(x)
}
counts <- lapply(X = counts, FUN = convert_names)




# Filter called cells -----------------------------------------------------

filtered_counts <- vector(mode = 'list', length = length(counts))
names(filtered_counts) <- names(counts)
for(ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  retain_col <- which(drops[[sample_id]][['Retain']])
  filtered_counts[[sample_id]] <- counts[[sample_id]][,retain_col]
}
names(filtered_counts) <- gsub(pattern = 'raw', 
                               replacement = 'filtered', 
                               x = names(counts))

# save filtered count matrices
for(ii in 1:length(filtered_counts)) {
  sample_id <- names(filtered_counts)[ii]
  saveRDS(object = filtered_counts[[sample_id]],
          file = paste0(data_out, sample_id, '.rds'))
}
sapply(filtered_counts, ncol)

# cleanup
rm(list = ls())
gc()
