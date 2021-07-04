
# HELPER FUNCTIONS

load(file = 'SCI_portal_data.RData')

is.gene <- function(x) {x %in% gene}
is.categorical <- function(x) {x %in% categorical_data}
is.numerical <- function(x) {x %in% numerical_data}
cols_fxn <- colorRamp(c("grey", "lightblue", "dodgerblue4", "royalblue4"))
exp_scale <- function(x) {
  if ((max(x) - min(x)) == 0) {
    return('grey')
  }
  x <- (x - min(x)) / (max(x) - min(x))
  tmp <- colorRamp(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(x)
  mapply(
    FUN = rgb,
    red = tmp[,1],
    green = tmp[,2],
    blue = tmp[,3],
    maxColorValue = 255
  )
}


# Gene expression tabPanel ----

## expression_gene1_UMAPplot ----
expression_gene1_UMAPplot <- function(
  gene1 = "Cx3cr1",
  dataset = 'macroglia'
) {
  df <- data.frame(
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_1', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_2', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    as.numeric(log_x_sci[gene1,])
  )
  colnames(df) <- c('UMAP_1','UMAP_2','gene1')
  df <- df[order(df$gene1, decreasing = FALSE),]
  # p1 <- ggplot(data = df,
  #        mapping = aes(x = UMAP_1, y = UMAP_2)) +
  #   geom_point(mapping = aes(color = gene1)) +
  #   labs(title = paste('Gene1:', gene1), subtitle = paste('Dataset:', dataset)) +
  #   theme(plot.title = element_text(size = 14),
  #         plot.subtitle = element_text(size = 12),
  #         panel.background = element_rect(colour = NA, fill = NA),
  #         panel.border = element_rect(colour = NA, fill = NA),
  #         axis.line = element_line(color = 'black', size = 0.5),
  #         axis.title = element_text(hjust = 0, size = 12),
  #         axis.ticks = element_blank(),
  #         axis.text = element_blank(),
  #         legend.key = element_rect(fill = NA),
  #         legend.text = element_text(size = 12),
  #         legend.title = element_text(size = 12)) +
  #   guides(color = guide_legend(title = 'log-norm\nexpression', 
  #                               override.aes = list(size = 4)))
  # p1 <- ggplotly(p = p1, width = 425, height = 325) %>%
  #   layout(autosize = FALSE)
  # return(p1)
  plot(x = df$UMAP_1, y = df$UMAP_2, col = exp_scale(df$gene1),
       pch = 21, cex = 0.2, xlab = 'UMAP_1', ylab = 'UMAP_2',
       main = paste('Gene 1:', gene1))
}

## expression_gene2_UMAPplot ----
expression_gene2_UMAPplot <- function(
  gene2 = "Cx3cr1",
  dataset = 'macroglia'
) {
  df <- data.frame(
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_1', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_2', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    as.numeric(log_x_sci[gene2,])
  )
  colnames(df) <- c('UMAP_1','UMAP_2','gene2')
  df <- df[order(df$gene2, decreasing = FALSE),]
  # p2 <- ggplot(data = df,
  #              mapping = aes(x = UMAP_1, y = UMAP_2)) +
  #   geom_point(mapping = aes(color = gene2)) +
  #   labs(title = paste('Gene2:', gene2), subtitle = paste('Dataset:', dataset)) +
  #   theme(plot.title = element_text(size = 14),
  #         plot.subtitle = element_text(size = 12),
  #         panel.background = element_rect(colour = NA, fill = NA),
  #         panel.border = element_rect(colour = NA, fill = NA),
  #         axis.line = element_line(color = 'black', size = 0.5),
  #         axis.title = element_text(hjust = 0, size = 12),
  #         axis.ticks = element_blank(),
  #         axis.text = element_blank(),
  #         legend.key = element_rect(fill = NA),
  #         legend.text = element_text(size = 12),
  #         legend.title = element_text(size = 12)) +
  #   guides(color = guide_legend(title = 'log-norm\nexpression', 
  #                               override.aes = list(size = 4)))
  # p2 <- ggplotly(p = p2, width = 425, height = 325) %>% 
  #   layout(autosize = FALSE)
  # return(p2)
  plot(x = df$UMAP_1, y = df$UMAP_2, col = exp_scale(df$gene2),
       pch = 21, cex = 0.2, xlab = 'UMAP_1', ylab = 'UMAP_2',
       main = paste('Gene 2:', gene2))
}

## expression_gene1_VlnPlot ----
expression_gene1_VlnPlot <- function(
  gene1 = "Cx3cr1", 
  dataset = 'macroglia', 
  group.by = 'L1_taxon'
) {
  df <- data.frame(
    obs_sci[grepl(pattern = group.by, x = colnames(obs_sci))],
    as.numeric(log_x_sci[gene1,])
  )
  colnames(df) <- c('Group','gene1')
  exp_table <- df %>%
    group_by(Group) %>%
    summarise(ncells = n(),
              pct_pos = sum(gene1 > 0)/n() * 100)
  exp_table$height <- max(df$gene1)
  exp_table$pct_pos <- round(exp_table$pct_pos, digits = 2)
  p1 <- df[!is.na(df$Group),] %>%
    ggplot(mapping = aes(x = Group, y = gene1)) +
    geom_violin(scale = 'width') +
    labs(title = paste('Gene 1:', gene1)) +
    xlab(label = group.by) +
    ylab(label = 'normalized exp.') +
    geom_text(data = exp_table, 
              mapping = aes(x = Group,
                            y = height - 0.1*height, 
                            label = paste('%+:', pct_pos))) +
    geom_text(data = exp_table,
              mapping = aes(x = Group,
                            y = height + 0.1*height,
                            label = paste0('# cells:', ncells))) +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12))
  # p1 <- ggplotly(p = p1, width = 425, height = 50) %>%
  #   layout(autosize = FALSE)
  return(p1)
}


## expression_gene2_VlnPlot ----
expression_gene2_VlnPlot <- function(
  gene2 = "Cx3cr1", 
  dataset = 'macroglia', 
  group.by = 'L1_taxon'
) {
  df <- data.frame(
    obs_sci[grepl(pattern = group.by, x = colnames(obs_sci))],
    as.numeric(log_x_sci[gene2,])
  )
  colnames(df) <- c('Group','gene2')
  exp_table <- df %>%
    group_by(Group) %>%
    summarise(ncells = n(),
              pct_pos = sum(gene2 > 0)/n() * 100)
  exp_table$height <- max(df$gene2)
  exp_table$pct_pos <- round(exp_table$pct_pos, digits = 2)
  p2 <- df[!is.na(df$Group),] %>%
    ggplot(mapping = aes(x = Group, y = gene2)) +
    geom_violin(scale = 'width') +
    labs(title = paste('Gene 2:', gene2)) +
    xlab(label = group.by) +
    ylab(label = 'normalized exp.') +
    geom_text(data = exp_table, 
              mapping = aes(x = Group,
                            y = height - 0.1*height, 
                            label = paste('%+:', pct_pos))) +
    geom_text(data = exp_table,
              mapping = aes(x = Group,
                            y = height + 0.1*height,
                            label = paste0('# cells:', ncells))) +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12))
  # p2 <- ggplotly(p = p2, width = 425, height = 50) %>%
  #   layout(autosize = FALSE)
  return(p2)
}
