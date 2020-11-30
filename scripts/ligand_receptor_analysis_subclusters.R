
########## Ligand receptor co-expression analysis ###########


# Data import -------------------------------------------------------------


# For stochastic methods
set.seed(123)

# libraries and directories
require('Seurat')
require('dplyr')
require('ggplot2')
results_out <- './results/ligand_receptor_analysis/'
data_out <- './data/ligand_receptor_analysis/'
ref_in <- './ref/'
ref_out <- './ref/'
dir.create(path = results_out)
dir.create(path = data_out)

source('./scripts/ligand_receptor_analysis_functions.R')

sci <- readRDS(file = './data/sci.rds')
# sci_rna_1 <- readRDS(file = paste0(data_out, 'sci_rna_1.rds'))
# sci_rna_2 <- readRDS(file = paste0(data_out, 'sci_rna_2.rds'))



# Data setup --------------------------------------------------------------


# We first have to map subcluster identities back to full SCI dataset.
subcluster_paths <- list(
  'myeloid_subclusters' = './results/myeloid_annotation_markers/myeloid_subcluster.tsv',
  'vascular_subclusters' = './results/vascular_annotation_markers/vascular_subcluster.tsv',
  'macroglia_subclusters' = './results/macroglia_annotation_markers/macroglia_subcluster.tsv'
)
# Load
subclusters_bySet <- lapply(
  X = subcluster_paths,
  FUN = read.table,
  sep = '\t',
  header = TRUE,
  row.names = 1
)
# Extract labels and barcodes
subclusters <- Reduce(
  f = c, 
  x = sapply(
    X = subclusters_bySet,
    FUN = `[[`,
    1
  )
)
barcodes <- Reduce(
  f = c,
  x = sapply(
    X = subclusters_bySet,
    FUN = rownames
  )
)
names(subclusters) <- barcodes
rm(barcodes, subcluster_paths)

# Add to sci metadata
sci@meta.data[['subcluster']] <- subclusters[match(x = rownames(sci@meta.data), 
                                                   table = names(subclusters))]
# label unassigned cells
sci@meta.data[['subcluster']][is.na(sci@meta.data[['subcluster']])] <- 'Unassigned'
Idents(sci) <- 'subcluster'


# Clusters to merge
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('OPC-A','OPC-B'),
  to = c('OPC','OPC')
)
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('Ependymal-A','Ependymal-B'),
  to = c('Ependymal','Ependymal')
)
sci@meta.data[['subcluster']] <- plyr::mapvalues(
  x = sci@meta.data[['subcluster']],
  from = c('C1-Endothelial','C2-Endothelial'),
  to = c('C-Endothelial', 'C-Endothelial')
)


# Interaction analysis at subcluster level --------------------------------


# Genes for myeloid, vascular, and macroglia clusters identified from previous.
previous_genes <- c('Tie1','Tek','Angpt1','Angpt2','Flt1','Kdr','Flt4','Nrp1','Nrp2','Vegfa', 'Vegfb', 'Vegfc', 'Vegfd', 'Pgf', 'Il6st', 'Il6ra','Ilfra','Cntfr', 'Osmr', 'Apoe', 'Spp1', 'Il6', 'Cntf', 'Lif', 'Ctf1', 'Osm', 'Agt', 'Angpt4', 'Angptl1', 'Apln', 'Csf1', 'Csf2', 'Csf3', 'Cxcl10', 'Cxcl12', 'Dlk1', 'Dll1', 'Dll4', 'Edn1', 'Edn2', 'Edn3', 'Efna1', 'Efna2', 'Efna3', 'Efna4', 'Efna5', 'Efnb1', 'Efnb2' , 'Efnb3', 'Egf', 'Fn1', 'Hbegf', 'Il10', 'Il11', 'Il31', 'Pdgfa', 'Pdgfb', 'Pdgfc', 'Pdgfd', 'Sema3a', 'Sema3b', 'Sema3c', 'Sema3e', 'Sema3f', 'Sema4d', 'Sema5a', 'Sema6d', 'Sema7a', 'Tgfb1', 'Tgfb2', 'Tgfb3')

lr_ref <- read.csv(file = './ref/fantom_PairsLigRec_mouse.csv')
keep_rows <- lr_ref$Ligand.ApprovedSymbol %in% previous_genes | 
  lr_ref$Receptor.ApprovedSymbol %in% previous_genes
lr_ref <- lr_ref[keep_rows,]


# Subset subcluster identities of interest
select_sci_subclusters <- c('Neutrophil',
                            'Monocyte',
                            'Macrophage-A',
                            'Macrophage-B',
                            'Div-Myeloid',
                            'Microglia-A',
                            'Microglia-B',
                            'Microglia-C',
                            'Div-Microglia',
                            'BA-Macrophage',
                            'A-Endothelial',
                            'C-Endothelial',
                            'V-Endothelial',
                            'Tip Cell',
                            'Fibroblast',
                            'Pericyte',
                            'VSMC',
                            'Ependymal',
                            'Astroependymal',
                            'Astrocyte',
                            'OPC',
                            'Div-OPC',
                            'Oligodendrocyte')
Idents(sci) <- 'subcluster'
select_sci <- subset(sci, idents = select_sci_subclusters)
rm(sci); gc()

# Make lightweight data structures
DefaultAssay(select_sci) <- 'RNA'
select_sci[['SCT']] <- NULL
select_sci[['integrated']] <- NULL

# Set identities and calculate setup values (avg_exp, pct_exp, cell counts)
Idents(select_sci) <- 'subcluster'
sci_lr_setup <- setupLR(seurat_object = select_sci,
                        lr_ref = lr_ref_1,
                        split_by = 'time',
                        assay = 'RNA',
                        slot = 'data')
gc()

# Permutation test
sci_result <- calculateLR(setup = sci_lr_setup,
                          resample = 1000,
                          adjust_pval = FALSE)

# Save results
write.table(x = sci_result,
            file = paste0(results_out, 'sci_cluster_LigandReceptorResults_manuscript_interactions.tsv'),
            sep = '\t',
            col.names = TRUE)

# Sample plot
plotLR(results = lr_results, ligands = 'Apoe', l_cells = c('Macrophage-A','Macrophage-B','Microglia-A'))
