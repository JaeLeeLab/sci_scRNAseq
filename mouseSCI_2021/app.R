
library('shiny')
library('ggplot2')
library('dplyr')

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
            inputId = 'gene1',
            label = 'Type 1st gene of interest',
            choices = gene,
            selected = 'Cx3cr1',
            multiple = FALSE
          ),
          selectInput(
            inputId = 'gene2',
            label = 'Type 2nd gene of interest',
            choices = gene,
            selected = 'Cx3cr1',
            multiple = FALSE
          ),
          selectInput(
            inputId = 'dataset',
            label = 'Select dataset',
            choices = list('sci','myeloid','vascular','macroglia'),
            selected = 'macroglia',
            multiple = FALSE
          ),
          selectInput(
            inputId = 'vln.group.by',
            label = 'Group violin by:',
            choices = categorical_data,
            selected = 'L2_taxon',
            multiple = FALSE
          )
        ),
        
        #### Visualize gene expression ----
        mainPanel = mainPanel(
          width = 10,
          
          ##### Query two genes at once ----
          fluidRow(
            column(
              width = 5,
              # div(style = 'height:70px'),
              plotOutput(
                outputId = 'expression_gene1_UMAPplot',
                height = '550px'
              )
            ),
            column(
              width = 5,
              # div(style = 'height:70px'),
              plotOutput(
                outputId = 'expression_gene2_UMAPplot',
                height = '550px'
              )
            )
          ),
          
          #### Expression violin plot ----
          fluidRow(
            column(
              width = 10,
              # div(style = 'margin-top:50px'),
              plotOutput(
                outputId = 'expression_gene1_VlnPlot',
                height = '250px'
              ),
              plotOutput(
                outputId = 'expression_gene2_VlnPlot',
                height = '250px'
              ),
              br()
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
  updateSelectizeInput(
    session = session,
    label = 'Type 1st gene of interest',
    inputId = 'gene1',
    choices = gene,
    selected = 'Cx3cr1',
    server = TRUE
  )
  updateSelectizeInput(
    session = session,
    label = 'Type 2nd gene of interest',
    inputId = 'gene2',
    choices = gene,
    selected = 'Cx3cr1',
    server = TRUE
  )
  updateSelectInput(
    session = session,
    label = 'Select dataset',
    inputId = 'dataset',
    choices = list('sci','myeloid','vascular','macroglia'),
    selected = 'macroglia'
  )
  updateSelectInput(
    session = session,
    label = 'Group violin by:',
    inputId = 'vln.group.by',
    selected = 'L2_taxon',
    choices = categorical_data
  )
  
  ## expression_gene1_UMAPplot ----
  output$expression_gene1_UMAPplot <- renderPlot(expr = {
    expression_gene1_UMAPplot(gene1 = input$gene1, dataset = input$dataset)
  })
  
  # ## expression_gene2_UMAPplot ----
  output$expression_gene2_UMAPplot <- renderPlot(expr = {
    expression_gene2_UMAPplot(gene2 = input$gene2, dataset = input$dataset)
  })
  
  # ## expression_gene1_VlnPlot ----
  output$expression_gene1_VlnPlot <- renderPlot(expr = {
    expression_gene1_VlnPlot(
      gene1 = input$gene1,
      dataset = input$dataset,
      group.by = input$vln.group.by
    )
  })

  ## expression_gene2_VlnPlot ----
  output$expression_gene2_VlnPlot <- renderPlot(expr = {
    expression_gene2_VlnPlot(
      gene2 = input$gene2,
      dataset = input$dataset,
      group.by = input$vln.group.by
    )
  })
}



# RunApp ------------------------------------------------------------------

shiny::shinyApp(ui = ui, server = server)

# rsconnect::deployApp('shinyApp')
# rsconnect::accounts()
# rsconnect::accountInfo()
# rsconnect::