

######### Data preparation for GEO submission ########

library('Matrix')
library('Seurat')



files_dir <- './data/QC_filtered_feature_bc_matrix/'
outs_dir <- './data/GEO_submission/'
dir.create(path = outs_dir)
gene_convert <- read.table(file = './ref/gene_name_conversion.tsv',
                           header = TRUE, sep = '\t')


my_files <- list.files(path = files_dir)
my_files <- my_files[grepl(pattern = '.rds', x = my_files)]
for (f in my_files) {
  dat <- readRDS(file = paste0(files_dir, f))
  dat <- Matrix::Matrix(data = dat, sparse = FALSE)
  tmp_filename <- gsub(
    pattern = '.rds',
    replacement = '.txt',
    x = f
  )
  write.table(
    x = dat,
    file = paste0(outs_dir, tmp_filename),
    sep = '\t',
    row.names = TRUE,
    col.names = TRUE
  )
}



sci <- readRDS(file = './data/sci.rds')
sci[['RNA']]@meta.features$mgi_symbol <- rownames(sci[['RNA']]@meta.features)
sci[['RNA']]@meta.features$ensembl_gene_id <- plyr::mapvalues(
  x = sci[['RNA']]@meta.features$mgi_symbol,
  from = gene_convert$mgi_symbol,
  to = gene_convert$ensembl_gene_id
)
sci_mat <- sci[['RNA']]@counts
genes_tsv <- rownames(sci_mat)
barcodes_tsv <- colnames(sci_mat)

writeMM(sci_mat, file = paste0(outs_dir, 'sci_mat.mtx'))
write(x = genes_tsv, file = paste0(outs_dir, 'genes.tsv'))
write(x = barcodes_tsv, file = paste0(outs_dir, 'barcodes.tsv'))

gene_metadata <- sci[['RNA']]@meta.features
barcode_metadata <- sci@meta.data

write.table(x = gene_metadata, file = paste0(outs_dir, 'gene_metadata.tsv'),
            sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(x = barcode_metadata, file = paste0(outs_dir, 'barcode_metadata.tsv'),
            sep = '\t', row.names = TRUE, col.names = TRUE)
