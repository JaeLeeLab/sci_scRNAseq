
tmp <- readRDS(file = './data/filtered_feature_bc_matrix/filtered_feature_bc_matrix_uninj_sample3.rds')
x <- Matrix::colSums(tmp)
tmp <- tmp[,x > 1935 & x < 26504]
tmp <- tmp[!duplicated(rownames(tmp)) & nchar(rownames(tmp)) != 0,]
tmp <- tmp %>%
  CreateSeuratObject() %>%
  NormalizeData() %>%
  PercentageFeatureSet(pattern = '^mt-', col.name = 'percent.mt') %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 20) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:15)

p1 <- VlnPlot(tmp, 'percent.mt', pt.size = 0) + NoLegend()
p2 <- VlnPlot(tmp, 'Aqp4', pt.size = 0) + NoLegend()
p3 <- VlnPlot(tmp, 'Mbp', pt.size = 0) + NoLegend()
p4 <- VlnPlot(tmp, 'Olig1', pt.size = 0) + NoLegend()
p5 <- VlnPlot(tmp, 'Bmp4', pt.size = 0) + NoLegend()
p6 <- p1/p2/p3/p4/p5
ggsave(filename = './results/quality_control/__uninj_sample3_astrocyte_percentmt.tiff', plot = p6, device = 'tiff', height = 10, width = 8)
