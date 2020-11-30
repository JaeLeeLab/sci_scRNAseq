
######### Ligand-receptor analysis functions ########

# Written by: JSC
# Date: 2020/07/27
# Method inspired by CellPhoneDB

# Packages
require('Seurat')
require('ggplot2')
require('dplyr')
require('tibble')





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setupLR <- function(
  seurat_object,
  ref_path = NULL,
  lr_ref = NULL,
  split_by = NULL,
  min_pct = 0.1,
  assay = "RNA",
  slot = "data",
  source_filter = NULL
) {
  
  # Make sure seurat data is present.
  if (class(seurat_object) != 'Seurat') {
    stop('\"seurat_object\" must be of class Seurat')
  }
  
  # Make sure ligand-receptor reference list is provided.
  if (is.null(ref_path) && is.null(lr_ref)) {
    stop('Must provide either \"ref_path\" (path of LR reference) or \"lr_ref\" (data.frame of imported LR reference).')
  }
  
  # Load reference csv.
  if (!is.null(ref_path)) {
    lr_ref <- read.csv(file = ref_path, stringsAsFactors = FALSE)
  }
  
  if (!is.null(source_filter)) {
    lr_ref <- lr_ref %>%
      filter(Pair.Source == source_filter)
  }
  
  # Check for Pair.Name column. Error if not present.
  if (!any(colnames(lr_ref) == 'Pair.Name')) {
    stop('Ligand-receptor reference list requires a column titled \"Pair.Name\".
         Entries in \"Pair.Name\" should be formatted as: [Ligand gene]_[Receptor gene]. 
         E.g. Apoe_Lrp1"')
  } else {
    pair_column <- which(colnames(lr_ref) == 'Pair.Name')
  }
  
  # Split Pair.Name into ligand and receptor names
  tmp <- strsplit(x = lr_ref[['Pair.Name']], split = '_')
  ligand_names <- sapply(X = tmp, FUN = `[`, 1)
  receptor_names <- sapply(X = tmp, FUN = `[`, 2)
  
  # Retain all ligand-receptor pairs where both gene names are present in the gene expression dataset
  all_genes <- rownames(slot(object = seurat_object[[assay]], 'counts'))
  lr_ref <- lr_ref[ligand_names %in% all_genes & receptor_names %in% all_genes,]
  if (nrow(lr_ref) == 0) {
    stop('No LR pairs were detected in provided Seurat data.')
  }

  # Get ligands and receptors that were detected in data
  tmp <- strsplit(x = lr_ref[['Pair.Name']], split = '_')
  ligand_names <- sapply(X = tmp, FUN = `[`, 1)
  receptor_names <- sapply(X = tmp, FUN = `[`, 2)
  
  # Message regarding use of seurat identities
  if (!length(unique(seurat_object@active.ident)) > 0) {
    stop('Only one identity set. Cannot calculate scores with one identity')
  }
  tmp <- paste0('Using active identities: ', paste(unique(seurat_object@active.ident), collapse = ', '))
  message(tmp)
  
  # New vector of all genes to retrieve data
  retrieve_genes <- union(ligand_names, receptor_names)
  active_idents <- slot(object = seurat_object, name = 'active.ident')
  
  # Extract data
  DefaultAssay(seurat_object) <- assay
  lr_data <- FetchData(object = seurat_object, vars = c(split_by, retrieve_genes), slot = slot)
  active_idents <- slot(object = seurat_object, name = 'active.ident')
  lr_data <- cbind(active_idents, lr_data)
  
  # Calculate average/percent expression for each gene, by "split_by" if provided.
  exp_avg <- lr_data %>% group_by(active_idents, .add = TRUE)
  exp_pct <- lr_data %>% group_by(active_idents, .add = TRUE)
  if (!is.null(split_by)) {
    split_by_name <- as.name(split_by)
    exp_avg <- exp_avg %>% group_by(!!split_by_name, .add = TRUE)
    exp_pct <- exp_pct %>% group_by(!!split_by_name, .add = TRUE)
  }
  exp_avg <- exp_avg %>% summarise(across(where(is.numeric), .fns = mean))
  exp_pct <- exp_pct %>% summarise(across(where(is.numeric), .fns = function(x) round(mean(x > 0), 3)))
  
  # Determine which genes meet minimum percent detection threshold
  minpct_genes <- sapply(X = exp_pct[sapply(exp_pct, is.numeric)], FUN = function(x) any(x > 0))
  minpct_genes <- names(minpct_genes)[minpct_genes]
  lr_ref_out <- lr_ref[ligand_names %in% minpct_genes & receptor_names %in% minpct_genes,]
  lr_genes_out <- sort(unique(minpct_genes))
  
  # Table of cell counts, further split if "split_by" provided.
  if (!is.null(split_by)) {
    cell_counts <- table(slot(object = seurat_object, name = 'active.ident'),
                         slot(object = seurat_object, name = 'meta.data')[[split_by]])
  } else {
    cell_counts <- table(slot(object = seurat_object, name = 'active.ident'))
  }
  
  outs <- list(
    lr_data = lr_data,
    lr_ref = lr_ref_out,
    lr_genes = lr_genes_out,
    exp_avg = exp_avg,
    exp_pct = exp_pct,
    cell_counts = cell_counts,
    split_by = split_by,
    assay = assay,
    slot = slot
  )
  
  return(outs)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calculateLR <- function(
  setup,
  resample = 1000,
  adjust_pval = FALSE,
  tmp_out_dir = NULL,
  ligand_subsets = NULL,
  receptor_subsets = NULL
) {
  # Import data
  lr_data <- setup[['lr_data']]
  lr_ref <- setup[['lr_ref']]
  lr_genes <- setup[['lr_genes']]
  exp_avg <- setup[['exp_avg']]
  exp_pct <- setup[['exp_pct']]
  split_by <- setup[['split_by']]
  cell_counts <- setup[['cell_counts']]
  assay <- setup[['assay']]
  slot <- setup[['slot']]
  
  # Cell-level expression matrix for genes in LR reference
  exp_mat <- as.matrix(lr_data[sapply(lr_data, is.numeric)])
  
  # Extract ligand/receptor gene names + expression matrices
  tmp <- strsplit(x = lr_ref[['Pair.Name']], split = '_')
  ligand_names <- sapply(X = tmp, FUN = `[`, 1)
  receptor_names <- sapply(X = tmp, FUN = `[`, 2)
  exp_mat_ligands <- exp_mat[,match(ligand_names, colnames(exp_mat))]
  exp_mat_receptors <- exp_mat[,match(receptor_names, colnames(exp_mat))]
  
  
  # Data.frame with results and null distribution score values
  if (class(lr_data[['active_idents']]) != 'factor') {
    stop('Seurat identities must be of class factor.')
  }
  if (!is.null(split_by)) {
    var_set <- expand.grid(levels(lr_data[['active_idents']]), 
                           levels(lr_data[['active_idents']]), 
                           levels(lr_data[[split_by]]))
    colnames(var_set) <- c('Ligand_cell', 'Receptor_cell', 'split_by')
  } else {
    var_set <- expand.grid(levels(lr_data[['active_idents']]), 
                           levels(lr_data[['active_idents']]))
    colnames(var_set) <- c('Ligand_cell', 'Receptor_cell')
  }
  
  # If subsets of cells are specified, calculate only scores for those.
  if (!is.null(ligand_subsets)) {
    if (!all(ligand_subsets %in% levels(lr_data[['active_idents']]))) {
      stop('Not all ligand_subsets are in active identities.')
    }
    var_set <- var_set[var_set[['Ligand_cell']] %in% ligand_subsets,]
  }
  if (!is.null(receptor_subsets)) {
    if (!all(receptor_subsets %in% levels(lr_data[['active_idents']]))) {
      stop('Not all receptor_subsets are in active identities.')
    }
    var_set <- var_set[var_set[['Receptor_cell']] %in% receptor_subsets,]
  }
  
  null_scores <- vector(mode = 'list', length = nrow(var_set))
  names(null_scores) <- apply(X = var_set, MARGIN = 1, FUN = paste, collapse = '_')
  
  # Use maxT method for multiple hypotheses p-value adjustment (Benjamini-Hochberg appears too harsh for 1000x permutations)
  if (adjust_pval) {
    maxT_scores <- matrix(0, nrow = resample, ncol = nrow(lr_ref))
    colnames(maxT_scores) <- paste(colnames(exp_mat_ligands), colnames(exp_mat_receptors), sep = '_')
  }

  # progress bar
  pb = txtProgressBar(min = 0, max = nrow(var_set), initial = 0, style = 3)
  
  # calculate null distribution of randomly permuted ligand-receptor score values
  for (i in 1:nrow(var_set)) {
    cellx <- as.character(var_set[['Ligand_cell']][i]) # these are factors, but cell_counts is sorted by level
    celly <- as.character(var_set[['Receptor_cell']][i])
    if (!is.null(split_by)) {
      split_byz <- as.character(var_set[['split_by']][i])
      countx <- cell_counts[cellx, split_byz]
      county <- cell_counts[celly, split_byz]
    } else {
      countx <- cell_counts[cellx]
      county <- cell_counts[celly]
    }
    if (countx == 0 | county == 0) {
      next
    }
    
    # Util function to randomly select n_c cells (n = #, c = cell-type) and extract element positions.
    cell_sample <- function(x) {
      tmp <- rep(FALSE, nrow(exp_mat))
      tmp[sample(nrow(exp_mat), size = x, replace = FALSE)] <- TRUE
      return(tmp)
    }
    null_index_x <- t(replicate(n = resample, expr = cell_sample(countx)))
    null_index_y <- t(replicate(n = resample, expr = cell_sample(county)))
    
    # Calculate average expression of all ligands and receptors for the n_c cells
    null_avg_lig <- (null_index_x %*% exp_mat_ligands) / countx
    null_avg_rec <- (null_index_y %*% exp_mat_receptors) / county
    
    # Calculate null LR scores
    tmp_scores <- 1/2 * (null_avg_lig + null_avg_rec)
    colnames(tmp_scores) <- paste(colnames(null_avg_lig), colnames(null_avg_rec), sep = '_')
    null_scores[[i]] <- tmp_scores
    
    # For max-T method, take element-wise maximum values (this iterates nrow(var_set) times)
    if (adjust_pval) {
      maxT_scores[tmp_scores > maxT_scores] <- tmp_scores[tmp_scores > maxT_scores]
    }
    
    # progress bar
    setTxtProgressBar(pb, i); gc()
    
  }
  
  # Calculate average expression matrix
  exp_names <- exp_avg[['active_idents']]
  if (!is.null(split_by)) {
    exp_names <- paste(exp_names, exp_avg[[split_by]], sep = '_')
  }
  
  # Get indices from expression matrix for ligands/receptors
  index_x <- var_set[['Ligand_cell']]
  index_y <- var_set[['Receptor_cell']]
  if (!is.null(split_by)) {
    index_x <- paste(index_x, var_set[['split_by']], sep = '_')
    index_y <- paste(index_y, var_set[['split_by']], sep = '_')
  }
  index_x <- match(x = index_x, table = exp_names)
  index_y <- match(x = index_y, table = exp_names)
  index_l <- match(ligand_names, colnames(exp_avg))
  index_r <- match(receptor_names, colnames(exp_avg))
  exp_avg_lig <- as.matrix(exp_avg[index_x, index_l])
  exp_avg_rec <- as.matrix(exp_avg[index_y, index_r])
  
  # Sanity check - do all column LR pairs match for null maxT matrix and actual expression matrix?
  sanity_check <- all(paste(colnames(exp_avg_lig), colnames(exp_avg_rec), sep = '_') == colnames(tmp_scores))
  if (!sanity_check) {
    stop('LR column indexing error: null score LR column names do not match sample score LR column names. Refer to source code. Sorry!')
  }
  
  # Calculate LR scores (for data, not nulls)
  exp_scores <- 1/2 * (exp_avg_lig + exp_avg_rec)
  
  # NOTE: Percent threshold step is removed here and left to the visualization step for removal of data points.
  # # Determine which cells do not have at least 10% expression of their gene
  # exp_pct_l <- exp_pct[index_x, index_l] < 0.1
  # exp_pct_r <- exp_pct[index_y, index_r] < 0.1
  # exp_pct_lr <- exp_pct_l | exp_pct_r
  
  # Replace values that don't meet 10% threshold
  # exp_scores[exp_pct_lr] <- NA
  rownames(exp_scores) <- apply(X = var_set, MARGIN = 1, FUN = paste, collapse = '_')
  colnames(exp_scores) <- paste(colnames(exp_avg_lig), colnames(exp_avg_rec), sep = '_')
  
  # Calculate p-values from ecdf using null scores for each LR-pair. 
  # NOTE: A true permutation test will test all possible permutations of cell-
  # sampling. Since this is computationally prohibitive, we estimate with 1000 
  # permutations. Thus if LR-scores ("effects") are large, it is possible that 
  # p-values can be zero.  
  message('\nCalculating p-values...')
  pvals <- matrix(NA, nrow = nrow(var_set), ncol = ncol(exp_scores))
  
  # For raw p-values:
  for (i in 1:nrow(pvals)) {
    tmp_nulls <- null_scores[[i]]
    tmp_scores <- exp_scores[i,]
    for (j in 1:length(tmp_scores)) {
      # Skip when cell counts were zero
      if (is.null(tmp_nulls)) {
        next
      } else {
        pvals[i,j] <- 1-ecdf(tmp_nulls[,j])(tmp_scores[j])
      }
    }
  }
  colnames(pvals) <- paste(colnames(exp_avg_lig), colnames(exp_avg_rec), sep = '_')
  rownames(pvals) <- apply(X = var_set, MARGIN = 1, FUN = paste, collapse = '_')
  
  # For adjusted p-values:
  if (adjust_pval) {
    adj_pvals <- matrix(1, nrow = nrow(var_set), ncol = ncol(exp_scores))
    for (i in 1:ncol(adj_pvals)) {
      adj_pvals[,i] <- sapply(X = exp_scores[,i],
                              nulls = maxT_scores[,i],
                              FUN = function(x, nulls) {1-ecdf(nulls)(x)})
    }
    colnames(adj_pvals) <- paste(colnames(exp_avg_lig), colnames(exp_avg_rec), sep = '_')
    rownames(adj_pvals) <- apply(X = var_set, MARGIN = 1, FUN = paste, collapse = '_')
  }
  
  # Long-form data.table of results
  tmp <- expand.grid(rownames(exp_scores), colnames(exp_scores), stringsAsFactors = FALSE)
  tmp_idents <- strsplit(x = tmp[[1]], split = '_')
  tmp_pairs <- strsplit(x = tmp[[2]], split = '_')
  ligand_cell <- sapply(X = tmp_idents, FUN = `[[`, 1)
  receptor_cell <- sapply(X = tmp_idents, FUN = `[[`, 2)
  ligand <- sapply(X = tmp_pairs, FUN = `[[`, 1)
  receptor <- sapply(X = tmp_pairs, FUN = `[[`, 2)
  cell_pair <- paste(ligand_cell, receptor_cell, sep = '_')
  lr_pair <- paste(ligand, receptor, sep = '_')
  tmp_scores <- exp_scores %>% reshape2::melt() %>% .[['value']]
  tmp_avg_l <- exp_avg_lig %>% reshape2::melt() %>% .[['value']]
  tmp_avg_r <- exp_avg_rec %>% reshape2::melt() %>% .[['value']]
  tmp_pvals <- pvals %>% reshape2::melt() %>% .[['value']]
  if (adjust_pval) {
    tmp_adj_pvals <- adj_pvals %>% reshape2::melt() %>% .[['value']]
  } else {
    tmp_adj_pvals <- tmp_pvals
  }
  tmp_pct <- exp_pct %>% as.matrix()
  tmp_pct_l <- c(tmp_pct[index_x, index_l])
  tmp_pct_r <- c(tmp_pct[index_y, index_r])
  split_var <- NA
  if (!is.null(split_by)) {
    split_var <- sapply(X = tmp_idents, FUN = `[[`, 3)
    split_var <- factor(x = split_var, levels = levels(lr_data[[split_by]]))
    cell_prop <- round(prop.table(cell_counts, margin = 1), 3)*100
    tmp_count_l <- mapply(FUN = function(x,y) cell_prop[x,y], ligand_cell, split_var)
    tmp_count_r <- mapply(FUN = function(x,y) cell_prop[x,y], receptor_cell, split_var)
  }
  
  # Final result compilation
  results <- data.frame(
    'Pair_name' = lr_pair,
    'Score' = tmp_scores,
    'pval' = tmp_pvals,
    'adj_pval' = tmp_adj_pvals,
    'Ligand_cell' = ligand_cell,
    'Receptor_cell' = receptor_cell,
    'split_by' = split_var,
    'Ligand' = ligand,
    'Receptor' = receptor,
    'Ligand_avgExp' = tmp_avg_l,
    'Receptor_avgExp' = tmp_avg_r,
    'Ligand_pct' = tmp_pct_l,
    'Receptor_pct' = tmp_pct_r,
    'Cell_pair' = cell_pair,
    'LR_pair' = lr_pair,
    stringsAsFactors = FALSE
  )
  if (!is.null(split_by)) {
    results[['Ligand_cell_pct']] <- tmp_count_l
    results[['Receptor_cell_pct']] <- tmp_count_r
  }
  message('Done!')
  return(results)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotLR <- function(
  results,
  ligands = NULL,
  receptors = NULL,
  l_cells = NULL,
  r_cells = NULL,
  subset_split = NULL,
  split_along_y = FALSE,
  use_adj_pval = FALSE,
  min_pval = NULL,
  min_exp_percent = 0.1,
  min_cell_percent = 0,
  resample = 1000
) {
  # Import locally
  tmp_results <- results
  
  # Remove p-values for which expression percents do not meet provided threshold
  tmp_results[['pval']] <- ifelse(test = results[['Ligand_pct']] < min_exp_percent | 
                                    results[['Receptor_pct']] < min_exp_percent,
                                  yes = NA,
                                  no = results[['pval']])
  if (use_adj_pval) {
    if (all(is.na(tmp_results[['adj_pval']]))) {
      stop('Adjusted p-values not calculated')
    }
    tmp_results[['adj_pval']] <- ifelse(test = results[['Ligand_pct']] < min_exp_percent | 
                                          results[['Receptor_pct']] < min_exp_percent,
                                        yes = NA,
                                        no = results[['adj_pval']])
  }
  
  
  # Series of filters based on those provided by user.
  if (!is.null(ligands)) {
    tmp_results <- tmp_results %>% filter(Ligand %in% ligands)
    if (nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  if (!is.null(receptors)) {
    tmp_results <- tmp_results %>% filter(Receptor %in% receptors)
    if (nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  if (!is.null(l_cells)) {
    tmp_results <- tmp_results %>% filter(Ligand_cell %in% l_cells)
    if (nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  if (!is.null(r_cells)) {
    tmp_results <- tmp_results %>% filter(Receptor_cell %in% r_cells)
    if (nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  # Remove rows where at least "min_cell_percent" of a cell-type is not present in a given subset of "split_by".
  if (!is.null(subset_split)) {
    if (!all(is.na(tmp_results[['split_by']]))) {
      stop('Provided \"subset_split\" but no values present in \"split_by\" column of \"results\".')
    }
    tmp_results <- tmp_results %>% filter(split_by %in% subset_split)
    if (nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  # Remove high p-values, if min_pval provided.
  if (!is.null(min_pval)) {
    tmp_results <- tmp_results %>% filter(pval < min_pval)
    if (nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  tmp_results <- tmp_results %>% ungroup()
  tmp_pairs <- unique(tmp_results[['Pair_name']][!is.na(tmp_results[['pval']])])
  tmp_results <- tmp_results[tmp_results[['Pair_name']] %in% tmp_pairs,]
  
  # Calculate -log10 of p-values for visualization. 
  min_col <- floor(min(tmp_results[['Score']])/0.5) * 0.5
  tmp_results[['log_pval']] <- -log10(tmp_results[['pval']])
  tmp_results[['log_pval']][is.infinite(tmp_results[['log_pval']])] <- log10(resample)
  if (use_adj_pval) {
    tmp_results[['log_adj_pval']] <- -log10(tmp_results[['adj_pval']])
    tmp_results[['log_adj_pval']][is.infinite(tmp_results[['log_adj_pval']])] <- log10(resample)
  }
  if (use_adj_pval) {
    tmp_pval <- 'log_adj_pval'
  } else {
    tmp_pval <- 'log_pval'
  }
  tmp_plot <- tmp_results %>%
    ggplot() + 
    geom_point(mapping = aes_string(x = 'Pair_name', y = 'Receptor_cell', size = tmp_pval, fill = 'Score'), 
               color = 'black', pch = 21)
  if (!all(is.na(tmp_results[['split_by']]))) {
    if (split_along_y) {
      tmp_plot <- tmp_plot + facet_grid(Ligand_cell + split_by ~ ., switch = 'y', drop = TRUE)
    } else {
      tmp_plot <- tmp_plot + facet_grid(Ligand_cell ~ split_by, switch = 'y', drop = TRUE)
    }
  } else {
    tmp_plot <- tmp_plot + facet_grid(Ligand_cell ~ ., switch = 'y')
  }
  tmp_plot <- tmp_plot + 
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = 'Spectral')),
                         limits = c(min_col, NA)) +
    scale_radius(limits = c(0,NA), range = c(1,6)) +
    scale_y_discrete(position = 'right') +
    theme(strip.text = element_text(size = 12, color = 'black', face = 'bold'),
          strip.background = element_rect(fill = NA, color = 'black'),
          axis.title = element_blank(),
          axis.text = element_text(size = 12, color = 'black'),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.text = element_text(size = 12, color = 'black'),
          legend.title = element_text(size = 12, color = 'black', face = 'bold'),
          legend.key = element_rect(fill = NA),
          panel.background = element_rect(fill = NA, color = 'black'),
          panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70')) +
    guides(fill = guide_colorbar(title = 'Score',
                                 frame.linewidth = 1,
                                 ticks.linewidth = 1,
                                 frame.colour = 'black',
                                 ticks.colour = 'black'),
           size = guide_legend(title = '-log10(p-value)', 
                               override.aes = list(fill = 'black')))
  return(tmp_plot)
}
