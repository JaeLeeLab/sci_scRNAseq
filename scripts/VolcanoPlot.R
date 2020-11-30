
# Classic DE test visualization is a volcano plot.

VolcanoPlot <- function(
  de_results,
  label_by_fc = NULL,
  label_by_pval = NULL,
  label_size = 4,
  use_adj_pval = TRUE,
  return_result = FALSE
) {
  require('ggrepel')
  de_results <- de_results %>%
    mutate('gene' = rownames(.)) %>%
    mutate('avg_log2FC' = log2(exp(x = avg_logFC)))
  if(use_adj_pval) {
    de_results <- de_results %>% 
      mutate('log_pval' = -log10(x = p_val_adj))
  } else {
    de_results <- de_results %>%
      mutate('log_pval' = -log10(x = p_val))
  }
  if(!is.null(label_by_fc) | !is.null(label_by_pval)) {
    de_results[['label']] <- de_results[['gene']]
    if(!is.null(label_by_fc)) {
      if(!any(abs(de_results[['avg_logFC']]) > label_by_fc)) {
        stop('No genes above provided avg_logFC threshold.')
      }
      de_results <- de_results %>%
        mutate('label' = ifelse(abs(avg_logFC) <= label_by_fc, NA, label))
    }
    if(!is.null(label_by_pval)) {
      if(use_adj_pval) {
        if(!any(de_results[['p_val_adj']] < label_by_pval)) {
          stop('No genes below provided adjusted p-value threshold.')
        }
        de_results <- de_results %>%
          mutate('label' = ifelse(p_val_adj >= label_by_pval, NA, label))
        } else {
          if(!any(de_results[['p_val']] < label_by_pval)) {
            stop('No genes below provided adjusted p-value threshold.')
          }
          de_results <- de_results %>%
            mutate('label' = ifelse(p_val >= label_by_pval, NA, label))
        }
    }
  }
  # return(de_results)
  de_volcano <- de_results %>%
    ggplot(mapping = aes(x = avg_logFC, y = log_pval)) +
    geom_point(color = 'grey40') + 
    geom_vline(xintercept = 0, linetype = 'dashed', size = 1, color = 'black') +
    scale_x_continuous() +
    xlab(label = 'Log2(fold-change)') +
    ylab(label = '-Log10(adjusted p-value)') +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black', size = 1),
          axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 16, color = 'black'),
          panel.grid = element_line(color = 'grey60', linetype = 'dotted'))
  if(!is.null(label_by_fc) | !is.null(label_by_pval)) {
    de_volcano <- de_volcano +
      geom_text_repel(mapping = aes(label = label), size = label_size)
    if(!is.null(label_by_fc)) {
      de_volcano <- de_volcano +
        geom_vline(xintercept = label_by_fc, color = 'red', linetype = 'dashed', size = 1) +
        geom_vline(xintercept = -label_by_fc, color = 'red', linetype = 'dashed', size = 1)
    }
    if(!is.null(label_by_pval)) {
      de_volcano <- de_volcano +
        geom_hline(yintercept = -log10(label_by_pval), color = 'blue', linetype = 'dashed', size = 1)
    }
  }
  if(return_result) {
    de_results <- de_results[!is.na(de_results[['label']]),]
    de_results[c('log_pval', 'label')] <- NULL
    return(de_results)
  }
  return(de_volcano)
}
