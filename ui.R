library(shinythemes)

# Define UI for application that displays GEO results based on user choices
ui <- shinyUI(fluidPage(theme = shinytheme("cosmo"), 
                        h1(strong("REALGAR"), align="center", style = "color: #9E443A;"),
                        h3(strong("Reducing Associations by Linking Genes And transcriptomic Results"), align="center", style = "color: #9E443A;"), hr(),
                        p("REALGAR integrates publically available microarray gene expression data from the", 
                          a("Gene Expression Omnibus (GEO),", href="http://www.ncbi.nlm.nih.gov/geo/"), 
                          " allowing researchers to easily access a breadth of information about tissue-specific gene expression", 
                          "related to asthma and glucocorticoid responsiveness.",
                          " The table below summarizes the datasets relevant to the conditions and tissue types selected,",
                          "with links to the original study and the full set of microarray expression data.",
                          " The plots aggregate data from these studies to display fold changes and p-values across all",
                          "relevant studies.  For the gene selected, the tracks show exons, obtained from the ",
                          a("UCSC Genome Browser,", href="https://genome.ucsc.edu/"),
                          " as well as genome-wide association study (GWAS) data from ",
                          a("GRASP,", href="https://grasp.nhlbi.nih.gov/Search.aspx"), 
                          "and glucocorticoid receptor NR3C1 binding sites from ",
                          a("UCSC's ENCODE", href="https://genome.ucsc.edu/ENCODE/"),
                          " where applicable. "), br(),
                        
                        wellPanel(fluidRow(align = "left",
                                           column(5,  
                                                  
                                                  fluidRow(column(5, h4("Select options:"),align = "left")), 
                                                  fluidRow(column(5,checkboxGroupInput(inputId="Tissue", label="Tissue", choices=c("Bronchial epithelium"="BE","Lens epithelial cell" = "LEC",
                                                                                                                              "Nasal epithelium"="NE","CD4"="CD4","CD8"="CD8","PBMC"="PBMC","White blood cell"="WBC", "Airway smooth muscle"="ASM",
                                                                                                                              "BAL"="BAL", "Whole lung"="Lung","Lymphoblastic leukemia cell" = "chALL","MCF10A-Myc" = "MCF10A-Myc",
                                                                                                                              "Macrophage" = "MDM","Osteosarcoma U2OS cell" = "U2O", "Lymphoblastoid cell" = "LCL"),selected="BE"),
                                                                  actionButton("selectall","Select/unselect all tissues", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")), 
                                                      
                                                      column(5,
                                                             fluidRow(checkboxGroupInput(inputId="Asthma", label="Asthma Type", choices=c("Allergic asthma"="allergic_asthma",
                                                                                                                                          "Severe asthma"="severe_asthma","Asthma"="asthma","Mild asthma"="mild_asthma","Non-allergic asthma"=
                                                                                                                                              "non_allergic_asthma", "Asthma and rhinitis"="asthma_and_rhinitis"), selected="asthma")),
                                                             fluidRow(radioButtons(inputId="GC_included", label="Treatment", 
                                                                                   choice = c("Glucocorticoid treatment" = "GC", "No treatment" = ""), select = "")),
                                                             fluidRow(textInput(inputId="curr_gene",label="Type the official gene symbol:", value= "GAPDH"),
                                                                      tags$head(tags$style(type="text/css", "#curr_gene {width: 190px}"))),
                                                             fluidRow(img(src="http://shiny.rstudio.com/tutorial/lesson2/www/bigorb.png", height=35, width=38),
                                                                      "Created with RStudio's ",
                                                                      a("Shiny", href="http://www.rstudio.com/shiny"))))),
                                           column(7,align="left",
                                                  fluidRow(h4("Datasets loaded:")),
                                                  fluidRow(DT::dataTableOutput(outputId="GEO_table"))))),
                        
                        
                        mainPanel(hr(),
                                  fluidRow(br(),
                                           column(4,plotOutput(outputId="forestplot_asthma",width="auto", height="650px"),align="left"),
                                           column(3,plotOutput(outputId="forestplot_GC",width="auto", height="650px"),align="left"),
                                           column(5,plotOutput(outputId="pval_plot_outp",width="auto", height="650px"),align="left"),
                                           column(4, downloadButton(outputId="asthma_fc_download",label="Download asthma forest plot"), align="center"),
                                           column(4, downloadButton(outputId="GC_fc_download",label="Download GC forest plot"), align="center"),
                                           column(4, downloadButton(outputId="pval_download", label="Download p-value plot"), align="center")),
                                  
                                  fluidRow(br(),
                                      column(12, h5(strong("Data shown in plots above")), align = "center")),
                                  
                                  fluidRow(
                                           column(10, offset= 1, DT::dataTableOutput(outputId="tableforgraph"), align="center"),
                                           column(12, downloadButton(outputId="table_download", label="Download data"), align="center")),br(), hr(), width = 12,
                                  
                                  fluidRow(column(10, align="center", plotOutput(outputId="gene_tracks_outp"), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br()),
                                           column(12, downloadButton(outputId="gene_tracks_download", label="Download gene tracks"), align="center"),br()))))
