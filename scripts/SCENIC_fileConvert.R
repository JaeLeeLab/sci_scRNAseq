
# Directory management. Write all files needed for SCENIC into results/SCENIC
setwd('D:/MiamiProject/sci_scRNAseq/')
data_out <- './results/SCENIC/'

# SeuratDisk package (currently in development 2021-03-04) used for saving 
# Seurat objects as .h5 files, which can be read by scanpy for SCENIC.
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
require('SeuratDisk')
require('Seurat')

# Load myeloid data, subset peripheral myeloid cells
myeloid <- readRDS(file = 'data/myeloid.rds')
DefaultAssay(myeloid) <- 'RNA'
macrophage <- c('Monocyte','Macrophage-A','Macrophage-B')
macrophage <- myeloid[,myeloid$myeloid_subcluster %in% macrophage]
macrophage$myeloid_subcluster <- as.character(macrophage$myeloid_subcluster)

cells_expr <- Matrix::rowSums(macrophage[['RNA']]@counts)
cells_expr <- cells_expr > 1/2*sqrt(ncol(macrophage))
macrophage <- macrophage[cells_expr,]


# Two-step process: save Seurat object as h5Seurat, which is an intermediate 
# on-disk file-type. Then convert h5Seurat to h5ad on-disk.
SaveH5Seurat(object = macrophage, 
             filename = paste0(data_out, 'macrophage.h5Seurat'),
             overwrite = TRUE)
Convert(source = paste0(data_out, 'macrophage.h5Seurat'),
        dest = 'h5ad',
        overwrite = TRUE)


test <- macrophage[,1:500]
SaveH5Seurat(object = test, 
             filename = paste0(data_out, 'test.h5Seurat'),
             overwrite = TRUE)
Convert(source = paste0(data_out, 'test.h5Seurat'),
        dest = 'h5ad',
        overwrite = TRUE)
