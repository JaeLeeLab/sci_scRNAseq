
runGO <- function(
  genes,
  use_scores = FALSE,
  score_fxn = NULL,
  statistic = 'fisher',
  ontology = 'BP'
) {
  require('org.Mm.eg.db')
  require('topGO')
  require('dplyr')
  
  if (use_scores) {
    if (is.null(score_fnx)) {
      stop('If use_scores, provide filtering function.')
    }
    go_input <- new(
      Class = 'topGOdata',
      # description = 'Upregulated in Mac-A over Mac-B (wilcox DE)',
      ontology = ontology,
      allGenes = genes,
      geneSel = score_fxn,
      annot = annFUN.org,
      mapping = 'org.Mm.eg.db',
      ID = 'ensembl'
    )
  } else {
    go_input <- new(
      Class = 'topGOdata',
      # description = 'Upregulated in Mac-A over Mac-B (wilcox DE)',
      ontology = ontology,
      allGenes = genes,
      annot = annFUN.org,
      mapping = 'org.Mm.eg.db',
      ID = 'ensembl'
    )
  }
  go_result <- runTest(
    object = go_input,
    statistic = statistic,
    algorithm = 'classic'
  )
  go_table <- GenTable(
    object = go_input,
    pvalue = go_result,
    topNodes = 250,
    numChar = 75
  )
  go_table[['odds_ratio']] <- go_table[['Significant']] / go_table[['Expected']]
  genes_byGO <- genesInTerm(
    object = go_input,
    whichGO = go_table[['GO.ID']]
  )
  genes_byGO <- lapply(
    X = genes_byGO,
    genes = genes,
    FUN = function(x, genes) {
      return(intersect(names(genes), x))
    }
  )
  outs <- list(
    'GO_table' = go_table,
    'gene_by_GO' = genes_byGO
  )
}
