
runGO <- function(
  genes,
  use_scores = FALSE,
  score_fxn = NULL,
  ontology = 'BP'
) {
  require('org.Mm.eg.db')
  require('topGO')
  require('dplyr')
  
  if (use_scores) {
    if (is.null(score_fxn)) {
      stop('If use_scores, provide filtering function.')
    }
    go_input <- new(
      Class = 'topGOdata',
      ontology = ontology,
      allGenes = genes,
      geneSel = score_fxn,
      annot = annFUN.org,
      mapping = 'org.Mm.eg.db',
      ID = 'ensembl',
      nodeSize = 10
    )
    go_result <- runTest(
      object = go_input,
      statistic = 'ks',
      algorithm = 'classic'
    )
  } else {
    go_input <- new(
      Class = 'topGOdata',
      ontology = ontology,
      allGenes = genes,
      annot = annFUN.org,
      mapping = 'org.Mm.eg.db',
      ID = 'ensembl',
      nodeSize = 10
    )
    go_result <- runTest(
      object = go_input,
      statistic = 'fisher',
      algorithm = 'classic'
    )
  }
  
  go_table <- GenTable(
    object = go_input,
    pvalue = go_result,
    topNodes = 250,
    numChar = 75
  )
  go_table[['odds_ratio']] <- go_table[['Significant']] / go_table[['Expected']]
  go_terms <- genesInTerm(
    object = go_input,
    whichGO = go_table[['GO.ID']]
  )
  genes_by_GO <- lapply(
    X = go_terms,
    genes = genes,
    FUN = function(x, genes) {
      sig_genes <- names(genes)[which(genes == 1)]
      return(intersect(sig_genes, x))
    }
  )
  outs <- list(
    'GO_table' = go_table,
    'gene_by_GO' = genes_by_GO
  )
}
