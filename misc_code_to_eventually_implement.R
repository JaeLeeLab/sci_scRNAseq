

# testing blocking factors

library('scran')
rna_sce <- as.SingleCellExperiment(sci, assay = 'RNA')
dim(logcounts(rna_sce))

# no blocking
t <- Sys.time()
markers_ctrl <- findMarkers(x = rna_sce,
                            groups = sci$default_cluster,
                            test.type = 'wilcox',
                            direction = 'up',
                            pval.type = 'some')
print(Sys.time() - t) # ~ 4 mins

# blocking time
t <- Sys.time()
markers_time <- findMarkers(x = rna_sce,
                            groups = sci$default_cluster,
                            test.type = 'wilcox',
                            direction = 'up',
                            block = sci$time,
                            pval.type = 'some')
print(Sys.time() - t) # ~6.5 mins



min_pval <- 1e-03
pvals_ctrl <- sapply(markers_ctrl, FUN = function(x) {sum(x[['FDR']] < min_pval)})
pvals_time <- sapply(markers_time, FUN = function(x) {sum(x[['FDR']] < min_pval)})

data.frame('ctrl' = pvals_ctrl, 
           'time' = pvals_time) %>%
  tibble::rownames_to_column(var = 'default_cluster') %>%
  reshape2::melt(id.vars = 'default_cluster') %>%
  ggplot(mapping = aes(x = default_cluster, y = value, fill = variable)) +
  geom_bar(position = 'dodge', stat = 'identity')


# pval.type = 'any' yields thousands of genes per cluster
# pval.type = 'all' yields less than 200 genes for some clusters (probably most abundant cell-types)
# pval.type = 'some' ranges from between several hundred to 3000 genes per cluster (more for unique cells like OPCs)




# pseudobulk methods for comparing groups (time) --------------------------

summed <- aggregateAcrossCells(merged, 
                               id=colData(merged)[,c("celltype.mapped", "sample")])
summed


Idents(sci) <- 'tmp_cell'
den <- AverageExpression(object = sci, 
                         assay = 'RNA',
                         slot = 'data', 
                         add.ident = 'sample_id',
                         features = sci@assays$integrated@var.features)
de

dat.dist <- dist(x = t(avg))
dat.tree <- hclust(d = dat.dist, method = 'ward.D2')
dat.dend <- dendsort::dendsort(as.dendrogram(dat.tree, hang = 0.1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
dat.cols <- gg_color_hue(n = length(levels(dat$seurat_clusters)))
dendextend::labels_colors(dat.dend) <- dat.cols[dat.tree$labels][order.dendrogram(dat.dend)]
dendextend::labels_colors(dat.dend) <- dat.cols
dendextend::labels_cex(dat.dend) <- 0.75
tiff(filename = paste0(results_outpath, 'Cluster_dendrogram.tiff'), height = 630, width = 630)
{par(mar = c(2,2,2,2), cex = 3, lwd = 2); dendextend::plot_horiz.dendrogram(dat.dend, side = FALSE, main = 'Cluster dendrogram')}
dev.off()
dendextend::plot_horiz.dendrogram(dat.dend, side = FALSE, main = 'Cluster dendrogram')






# -------------------------------------------------------------------------


ReadsPerCell <- function(
  mol_info, 
  cell_barcodes, 
  barcode_length = 16
) {
  set.seed(1)
  # Read in molecule info (.h5 file, CellRanger output)
  if (class(mol_info) == 'character'){
    mol_info <- DropletUtils::read10xMolInfo(sample = mol_info, get.gem = FALSE)
    
  } 
  # Cut cell barcode names to leave only barcode, depending on its length
  search <- substr(cell_barcodes, 1, barcode_length)
  total_reads <- sum(mol_info$data$reads[mol_info$data$cell %in% search])
  out <- c(total_reads, total_reads/length(search))
  names(out) <- c('total_mapped_reads','reads_per_cell')
  gc()
  h5closeAll()
  return(out)
}


scran::bootstrapCluster()


# Pairwise comparison for identifyign marker genes: https://osca.bioconductor.org/marker-detection.html

# The t-test also allows us to specify a non-zero log-fold change as the null hypothesis. This allows us to consider the magnitude of the log-fold change in our p-value calculations, in a manner that is more rigorous than simply filtering directly on the log-fold changes (McCarthy and Smyth 2009). (Specifically, a simple threshold does not consider the variance and can enrich for genes that have both large log-fold changes and large variances.)
markers.pbmc.up2 <- scran::findMarkers(sce.pbmc, direction="up", lfc=1)




# GUNZIP reference data ----------------------------------------------------------------

ref_path <- './ref/GSE121654_RAW_Hammond2019_Immunity/'
to_unzip <- list.files(path = ref_path)
for(ii in 1:length(to_unzip)) {
  out_name <- gsub(pattern = '\\.gz',
                   replacement = '',
                   x = to_unzip[ii])
  out_name <- gsub(pattern = '\\.dge',
                   replacement = '',
                   x = out_name)
  GEOquery::gunzip(filename = paste0(ref_path, to_unzip[ii]),
                   destname = paste0(ref_path, out_name))
}

tmp <- read.table(file = './ref/GSE121654_RAW_Hammond_Immunity2019/GSM3442006_E14_F_B10.txt',
                  sep = '\t', row.names = 1, header = TRUE)
str(tmp)
