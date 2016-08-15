library(shiny)

###abc 


# Define UI for application that displays GEO results based on user choices
shinyUI(fluidPage(theme = shinytheme("cosmo"),
                  h1(strong("Asthma Gene Explorer"), align="center", style = "color: #000080;"), hr(),
                  sidebarLayout(
                    sidebarPanel(width = 3, h4("Options:",align = "center"), hr(),
                                 checkboxGroupInput(inputId="Tissue", label="Tissue", choices=c("Bronchial Epithelium"="BE",
                                                                                                "Nasal Epithelium"="NE","CD4"="CD4","CD8"="CD8","PBMC"="PBMC","White Blood Cell"="WBC", "Airway Smooth Muscle"="ASM",
                                                                                                "BAL"="BAL", "Lung"="Lung"),selected="BE" ),
                                 checkboxGroupInput(inputId="Asthma", label="Condition", choices=c("Allergic Asthma"="allergic_asthma",
                                                                                                                                "Severe Asthma"="severe_asthma","Asthma"="asthma","Mild Asthma"="mild_asthma","Non Allergic Asthma"=
                                                                                                                                  "non_allergic_asthma", "Asthma and Rhinitis"="asthma_and_rhinitis", "Other"="NA"), selected="asthma"),
                                 textInput(inputId="curr_gene",label="Type the Official Gene Symbol", value= "GAPDH"),
                                 actionButton(inputId="updategraph",label="Update"),
                                 hr(),
                                 fixedRow(img(src="http://shiny.rstudio.com/tutorial/lesson2/www/bigorb.png", height=42, width=48),
                                          "Created with ",
                                          a("RStudio's shiny", href="http://www.rstudio.com/shiny"))),
                    mainPanel(
                      p("Asthma Gene Explorer integrates publically available gene expression data from the", 
                        a("Gene Expression Omnibus (GEO),", href="http://www.ncbi.nlm.nih.gov/geo/"), 
                        " allowing researchers to easily access a breadth of information about tissue-specific gene expression", 
                        "related to asthma and glucocorticoid responsiveness.",
                        " The table below summarizes the datasets relevant to the conditions and tissue types selected,",
                        "with links to the original study and the full set of microarray expression data.",
                        " The plots aggregate data from these studies to display fold changes and p-values across all",
                        "relevant studies, as well as genome-wide association study (GWAS) data from ",
                        a("GRASP,", href="https://grasp.nhlbi.nih.gov/Search.aspx"), " for the gene of interest. ",
                        "Note that fold changes were computed with a reference patten of asthma vs. [other type of asthma]."),
                      br(),
                      fluidRow(align="left",tableOutput(outputId="GEO_table")),
                      fluidRow(align="left",
                               column(6,plotOutput(outputId="fc_plot",width="95%", height="650px")),
                               column(6,plotOutput(outputId="pval_plot",width="95%", height="650px"))),
                      fluidRow(align="center", br(), column(12,tableOutput(outputId="tableforgraph"))),
                      fixedRow(align="center",br(),
                               downloadButton(outputId="geo_download",label="Download List of GEO Datasets"),
                               downloadButton(outputId="fc_download",label="Download Fold Change Plot"),
                               downloadButton(outputId="pval_download", label="Download P-Value Plot"),
                               downloadButton(outputId="table_download", label="Download heatmap data")),br()))))
