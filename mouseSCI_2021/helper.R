
# HELPER FUNCTIONS

load(file = 'SCI_portal_data.RData')

is.gene <- function(x) {x %in% gene}
is.categorical <- function(x) {x %in% categorical_data}
is.numerical <- function(x) {x %in% numerical_data}
exp_scale <- function(x, cols) {
  if ((max(x) - min(x)) == 0) {
    return('grey')
  }
  x <- (x - min(x)) / (max(x) - min(x))
  tmp <- colorRamp(cols)(x)
  mapply(
    FUN = rgb,
    red = tmp[,1],
    green = tmp[,2],
    blue = tmp[,3],
    maxColorValue = 255
  )
}

obs_sci$celltype <- factor(
  x = obs_sci$celltype,
  levels = c('Neutrophil',
             'Monocyte',
             'Macrophage',
             'Dendritic',
             'Microglia',
             'Div-Myeloid',
             'Fibroblast',
             'Endothelial',
             'Pericyte',
             'OPC',
             'Oligodendrocyte',
             'Astrocyte',
             'Ependymal',
             'Lymphocyte',
             'Neuron')
)
obs_sci$myeloid_subcluster <- factor(
  x = obs_sci$myeloid_subcluster,
  levels = c('Neutrophil',
             'Monocyte',
             'Chemotaxis-Inducing Mac',
             'Inflammatory Mac',
             'Border-Associated Mac',
             'Dendritic',
             'Dividing Myeloid',
             'Homeostatic Microglia',
             'Inflammatory Microglia',
             'Dividing Microglia',
             'Migrating Microglia',
             'Interferon Myeloid')
)
obs_sci$vascular_subcluster <- factor(
  x = obs_sci$vascular_subcluster,
  levels = c('A-Endothelial',
             'C-Endothelial',
             'V-Endothelial',
             'Tip Cell',
             'Pericyte',
             'VSMC',
             'Fibroblast',
             'U-Vascular')
)
obs_sci$macroglia_subcluster <- factor(
  x = obs_sci$macroglia_subcluster,
  levels = c('Ependymal-A',
             'Ependymal-B',
             'Astroependymal',
             'Astrocyte',
             'OPC-A',
             'OPC-B',
             'Div-OPC',
             'Pre-Oligo',
             'Oligodendrocyte')
)
obs_sci$time <- factor(
  x = obs_sci$time,
  levels = c('Uninjured','1dpi','3dpi','7dpi')
)
obs_sci$L1_taxon <- factor(
  x = obs_sci$L1_taxon,
  levels = c('Myeloid','Vascular','Macroglia','Neural')
)
obs_sci$L2_taxon <- factor(
  x = obs_sci$L2_taxon,
  levels = levels(obs_sci$celltype)
)
obs_sci$L3_taxon <- factor(
  x = obs_sci$L3_taxon,
  levels = c(levels(obs_sci$myeloid_subcluster),
             levels(obs_sci$vascular_subcluster),
             levels(obs_sci$macroglia_subcluster))
)


# Gene expression tabPanel ----

## expression_gene_UMAPplot ----
expression_gene_UMAPplot <- function(
  gene = "Cx3cr1",
  dataset = 'sci'
) {
  df <- data.frame(
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_1', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_2', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    as.numeric(log_x_sci[gene,])
  )
  colnames(df) <- c('UMAP_1','UMAP_2','gene')
  df <- df[order(df$gene, decreasing = FALSE),]
  par(mai = c(0.1, 0.1, 1, 0.1), mar = c(0.2, 0.2, 2, 0.2))
  plot(x = df$UMAP_1, y = df$UMAP_2, 
       col = exp_scale(x = df$gene, cols = c("grey", "lightblue", "dodgerblue4", "royalblue4")),
       pch = 16, cex = 0.25, xlab = 'UMAP_1', ylab = 'UMAP_2', 
       xaxt = 'n', yaxt = 'n',
       main = paste('Gene:', gene))
}

## expression_gene_splitUMAPplot ----
expression_gene_splitUMAPplot <- function(
  gene = "Cx3cr1",
  dataset = 'sci'
) {
  df <- data.frame(
    obs_sci[c(paste(dataset, 'UMAP_1', sep = '_'),
              paste(dataset, 'UMAP_2', sep = '_'),
              'time')],
    as.numeric(log_x_sci[gene,])
  )
  colnames(df) <- c('UMAP_1','UMAP_2','time','gene')
  df <- df[order(df$gene, decreasing = FALSE),]
  subplot <- function(df_in, time_in) {
    tmp <- df_in[df_in$time == time_in,]
    tmp_col <- exp_scale(x = df_in$gene, 
                         cols=c("grey","lightblue","dodgerblue4","royalblue4"))
    tmp_col <- tmp_col[df_in$time == time_in]
    plot(x = tmp$UMAP_1, y = tmp$UMAP_2,
         xlab = '', ylab = '',
         col = tmp_col, pch = 21, cex = 0.2575, main = time_in,
         cex.main = 1.5, xaxt='n', yaxt = 'n')
  }
  par(mfrow = c(2,2), mai = c(0.1, 0.1, 1, 0.1), mar = c(0.2,0.2,2,0.2))
  {
    subplot(df = df, time = levels(obs_sci$time)[1])
    subplot(df = df, time = levels(obs_sci$time)[2])
    subplot(df = df, time = levels(obs_sci$time)[3])
    subplot(df = df, time = levels(obs_sci$time)[4])
  }
  par(mfrow = c(1,1))
}

## cluster_UMAPplot ----
cluster_UMAPplot <- function(
  dataset = 'sci',
  groupby = 'celltype'
) {
  df <- data.frame(
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_1', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    obs_sci[grepl(pattern = paste(dataset, 'UMAP_2', sep = '_'), ignore.case = TRUE, x = colnames(obs_sci))],
    obs_sci[grepl(pattern = groupby, x = colnames(obs_sci))]
  )
  df <- df[sample(x = 1:nrow(df), size = nrow(df)),]
  colnames(df) <- c('UMAP_1','UMAP_2','group')
  if (!is.factor(df$group)) {
    df$group <- factor(df$group, levels = sort(unique(df$group)))
  }
  labels_df <- data.frame(
    'UMAP_1' = sapply(split(x = df$UMAP_1, f = df$group), median),
    'UMAP_2' = sapply(split(x = df$UMAP_2, f = df$group), median),
    'label' = levels(df$group)
  )
  y_range <- range(df$UMAP_2, na.rm = TRUE)
  umap_px_ratio <- signif((y_range[2] - y_range[1]), 1)
  y_range[1] <- y_range[1] - 0.9*(y_range[2] - y_range[1])
  legend_nrow <- ceiling(length(unique(df$group))/3)
  new_y_intersp <- function(n) {
    1.15 * signif(7/n, digits = 2)
  }
  {
    par(mai = c(0.1, 0.1, 1, 0.1), mar = c(0.2,0.2,2,0.2), xpd = NA)
    plot(x = df$UMAP_1, y = df$UMAP_2, 
       col = rainbow(n = length(unique(df$group)), v = 0.9)[df$group],
       pch = 16, cex = 0.25, xlab = 'UMAP_1', ylab = 'UMAP_2', # checkCEX
       main = paste('Cells grouped by:', groupby),
       xaxt = 'n', yaxt = 'n', cex.main = 1.1,
       ylim = y_range)
    text(x = labels_df$UMAP_1,
         y = labels_df$UMAP_2,
         labels = labels_df$label,
         cex = 1)
    tmp <- (max(df$UMAP_1, na.rm = TRUE) - min(df$UMAP_1, na.rm = TRUE))/3.5
    legend(x = 'bottom',
           x.intersp = 0.6,
           y.intersp = new_y_intersp(legend_nrow),
           text.width = tmp,
           legend = labels_df$label,
           col = rainbow(n = length(unique(df$group)), v = 0.9),
           pch = 16, xpd = TRUE, ncol = 3, plot = TRUE,
           pt.cex = 1.5, cex = 0.95)
  }
}

## expression_gene_DotPlot ----
expression_gene_DotPlot <- function(
  gene = "Cx3cr1", 
  dataset = 'sci', 
  groupby = 'celltype'
) {
  df <- data.frame(
    obs_sci[c(groupby, 'time')],
    as.numeric(log_x_sci[gene,])
  )
  colnames(df) <- c('Group','time','gene')
  pct_exp <- sapply(
    X = split(x = df$gene, f = paste(df$Group, df$time, sep ='__')), 
    FUN = function(x) sum(x > 0)/length(x) * 100
  )
  avg_exp <- sapply(
    X = split(x = df$gene, f = paste(df$Group, df$time, sep ='__')), 
    FUN = mean
  )
  group <- sapply(strsplit(names(pct_exp), '__'), `[`, 1)
  time <- sapply(strsplit(names(pct_exp), '__'), `[`, 2)
  ncells <- table(paste(df$Group, df$time, sep = '__'))
  class(ncells) <- 'numeric'
  exp_table <- data.frame(
    'ncells' = ncells,
    'pct_exp' = pct_exp,
    'avg_exp' = avg_exp,
    'group' = group,
    'time' = time
  )
  if (!is.factor(exp_table$group)) {
    exp_table$group <- factor(x = exp_table$group, 
                              levels = sort(unique(df$Group)))
  }
  # exp_table <- exp_table[!exp_table$group %in% "NA",]
  exp_table$time <- factor(exp_table$time,
                           levels = c('Uninjured','1dpi','3dpi','7dpi'))
  max_exp <- ceiling(max(avg_exp) * 10) / 10
  p1 <- ggplot(data = exp_table, 
         mapping = aes(x = time, y = group)) +
    geom_point(mapping = aes(size = pct_exp, fill = avg_exp), pch = 21)+
    scale_fill_viridis_c(option = 'A',
                         breaks = c(0, max_exp),
                         limits = c(0, max_exp),
                         labels = c(0, max_exp)) +
    scale_size(range = c(0, 8), limits = c(0, 100)) +
    theme_bw() +
    xlab(label = 'Time after SCI') +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                      size = 12),
          legend.title.align = 1,
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 14),) +
    guides(fill = guide_colorbar(title = 'Log(normalized counts)',
                                 title.position = 'left',
                                 frame.colour = 'black',
                                 ticks = FALSE),
           size = guide_legend(title = 'Percent detected',
                               title.position = 'left',
                               override.aes = list(color = 'black',
                                                   fill = 'black')))
  return(p1)
}

expression_summary_table <- function(
  gene = 'Cx3cr1',
  dataset = 'sci',
  groupby = 'celltype'
) {
  df <- data.frame(
    obs_sci[c(groupby, 'time', 'sample_id')],
    as.numeric(log_x_sci[gene,])
  )
  colnames(df) <- c('Group','time','sample_id','gene')
  pct_exp <- sapply(
    X = split(x = df$gene, f = paste(df$Group, df$time, sep = '__')), 
    FUN = function(x) sum(x > 0)/length(x) * 100
  )
  avg_exp <- sapply(
    X = split(x = df$gene, f = paste(df$Group, df$time, sep = '__')), 
    FUN = mean
  )
  group <- sapply(strsplit(names(pct_exp), '__'), `[`, 1)
  time <- sapply(strsplit(names(pct_exp), '__'), `[`, 2)
  ncells <- table(paste(df$Group, df$time, sep = '__'))
  class(ncells) <- 'numeric'
  exp_table <- data.frame(
    'ncells' = ncells,
    'pct_exp' = pct_exp,
    'avg_exp' = avg_exp,
    'group' = group,
    'time' = time
  )
}