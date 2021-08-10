
library('shiny')
library('ggplot2')

source('helper.R')

# UI: for app to query SCI scRNAseq -------------------------------------

ui <- fluidPage(
  
  ## App title ----
  titlePanel(h3('Single-cell atlas of mouse spinal cord injury')),
  
  ## Main tabset panel with query modes ----
  tabsetPanel(
    
    ### About tabPanel ----
    tabPanel(
      title = 'About',
      
      #### Publication reference ----
      fluidRow(
        column(
          width = 12,
          br(),
          HTML(
            text = "<p style='font-size:15px'>This website accompanies the 
            paper: Lindsay M. Milich, James S. Choi, Christine Ryan, Susana R. 
            Cerqueira, Sofia Benavides, Stephanie L. Yahn, Pantelis Tsoulfas, 
            Jae K. Lee; Single-cell analysis of the cellular heterogeneity and 
            interactions in the injured mouse spinal cord. <i>J Exp Med</i> 2 
            August 2021; 218 (8): e20210040. 
            <a href='https://doi.org/10.1084/jem.20210040'>DOI: https://doi.org/10.1084/jem.20210040</a></p><br>" 
          ),
          hr()
        )
      ),
      
      #### Description of tabPanels ----
      fluidRow(
        column(
          width = 8,
          HTML(
            text = '<p>Each of the tabs allows exploring the data presented in the paper.</p>'
          ),
          HTML(
            text = '<p><li><b>Gene expression:</b> type in a gene of interest to plot its expression in the UMAP low dimensional space.</li><p>'
          ),
          HTML(
            text = '<p><li><b>Gene expression time course:</b> type in a gene of interest to plot its expression in the UMAP low dimensional space, with cells from different injury time-points plotted separately.</li></p>'
          ),
          HTML(
          text = "<p><li><b>Cluster marker genes:</b> select a cluster of interest to display the genes that show preferential expression in that cluster over others. Clicking on a row in the table will show the corresponding expression plot. On the left, you can select the type of test, which controls how many clusters (all or 75%) have lower expression of the gene, compared to the selected cluster.</li></p>" 
          )
        ),
        column(
          width = 4,
          HTML(
            text = '<b>Insert UMAP summaries here.</b>'
          )
        ),
        hr()
      ),
      
      #### Data availability and code ----
      fluidRow(
        column(
          width = 12,
          HTML(
            text = "<p>Raw data is available from the Sequence Read Archive database under study <a href='https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP295673'>SRP295673</a>. Gene-count matrices are available from the Gene Expression Omnibus under accession <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162610'>GSE162610</a>. Relevant sample-level and cell-level metadata are available in under the GEO accession.</p><br>
            <p>Code used to analyze the single-cell RNA-seq data from <i>Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord</i> are available on <a href='https://github.com/JamesChoi94/sci_scRNAseq'>Github</a></p>"
          )
        ),
        hr()
      )
    ),
    
    ### Gene expresion tabPanel ----
    tabPanel(
      title = 'Gene expression',
      
      #### Gene and dataset selection ----
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          width = 2,
          div(style = 'height:10px'),
          selectInput(
            inputId = 'dataset_in',
            label = 'Select dataset',
            choices = list('sci','myeloid','vascular','macroglia'),
            selected = 'sci',
            multiple = FALSE
          ),
          selectInput(
            inputId = 'groupby_in',
            label = 'Group cells by:',
            choices = categorical_data,
            selected = 'celltype',
            multiple = FALSE
          ),
          selectInput(
            inputId = 'gene_in',
            label = 'Select gene of interest',
            choices = genes,
            selected = 'Cx3cr1',
            multiple = FALSE
          )
        ),
        
        #### Visualize gene expression ----
        mainPanel = mainPanel(
          width = 10,
          
          ##### Query gene and cluster ----
          fluidRow(
            column(
              width = 4,
              plotOutput(
                outputId = 'cluster_UMAPplot',
                height = '500px',
              )
            ),
            column(
              width = 8,
              plotOutput(
                outputId = 'expression_gene_splitUMAPplot',
                height = '500px'
              )
            )
          ),
          br(),
          fluidRow(
            column(
              width = 4,
              plotOutput(
                outputId = 'expression_gene_DotPlot',
                height = '550px',
              )
            ),
            column(
              width = 1
              # insert data table
            )
          )
        )
      )
    )
  )
)



# Server: logic required to generate figures ---------------------------

server <- function(input, output, session) {
  
  ## populate inputs ----
#   ## electInput(
#   inputId = 'dataset_in',
#   label = 'Select dataset_in',
#   choices = list('sci','myeloid','vascular','macroglia'),
#   selected = 'sci',
#   multiple = FALSE
#   ),
# selectInput(
#   inputId = 'groupby_in',
#   label = 'Group cells by:',
#   choices = categorical_data,
#   selected = 'celltype',
#   multiple = FALSE
# ),
# selectInput(
#   inputId = 'gene',
#   label = 'Select gene of interest',
#   choices = gene,
#   selected = 'Cx3cr1',
#   multiple = FALSE
# ),
# actionButton(
#   inputId = 'geneplot',
#   label = 'Plot expression'
# )
# ),
  updateSelectInput(
    session = session,
    label = 'Select dataset',
    inputId = 'dataset_in',
    choices = list('sci','myeloid','vascular','macroglia'),
    selected = 'sci'
  )
  updateSelectInput(
    session = session,
    label = 'Group cells by:',
    inputId = 'groupby_in',
    selected = 'celltype',
    choices = categorical_data
  )
  updateSelectizeInput(
    session = session,
    label = 'Select gene of interest',
    inputId = 'gene_in',
    choices = genes,
    selected = 'Cx3cr1',
    server = TRUE
  )
  
  
  ## cluster_UMAPplot ----
  output$cluster_UMAPplot <- renderPlot(expr = {
    cluster_UMAPplot(dataset = input$dataset_in, groupby = input$groupby_in)
  })
  
  ## expression_gene_splitUMAPplot ----
  output$expression_gene_splitUMAPplot <- renderPlot(expr = {
    expression_gene_splitUMAPplot(gene = input$gene_in, 
                                  dataset = input$dataset_in)
  })
  
  # ## expression_gene_UMAPplot ----
  # output$expression_gene_UMAPplot <- renderPlot(expr = {
  #   expression_gene_UMAPplot(gene = input$gene_in, dataset = input$dataset_in)
  # })
  
  ## expresssion_gene_DotPlot ----
  output$expression_gene_DotPlot <- renderPlot(expr = {
    expression_gene_DotPlot(gene = input$gene_in, 
                            dataset = input$dataset_in,
                            groupby = input$groupby_in)
  })
  
  # observeEvent(input$geneplot, {
  #   renderPlot(expr = {
  #     cluster_UMAPplot(dataset = input$dataset_in, groupby = input$groupby_in)
  #   })
  # })
  
  ## expression_gene2_UMAPplot ----
  # output$expression_gene2_UMAPplot <- renderPlot(expr = {
  #   expression_gene2_UMAPplot(gene2 = input$gene2, dataset = input$dataset_in)
  # })
  
  ## expression_gene2_VlnPlot ----
  # output$expression_gene2_VlnPlot <- renderPlot(expr = {
  #   expression_gene2_VlnPlot(
  #     gene2 = input$gene2,
  #     dataset = input$dataset_in,
  #     groupby = input$vln.groupby_in
  #   )
  # })
}



# RunApp ------------------------------------------------------------------

shiny::shinyApp(ui = ui, server = server)

# setwd('D:/MiamiProject/sci_scRNAseq/')
# rsconnect::deployApp('mouseSCI_2021')
# rsconnect::accounts()
# rsconnect::accountInfo()
# rsconnect::