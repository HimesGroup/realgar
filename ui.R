library(shinythemes)

# Define UI for application that displays GEO results based on user choices
ui <- shinyUI(fluidPage(theme = shinytheme("cosmo"), 
                        #h1(strong("REALGAR"), align="center", style = "color: #9E443A;"),
                        h3(strong("Reducing Associations by Linking Genes And transcriptomic Results"), align="center", style = "color: #9E443A;"), hr(),
                        p("REALGAR is an integrated resource of tissue-specific results from expression studies. ",
                          "This app brings together microarray expression data from the ", 
                          a("Gene Expression Omnibus (GEO),", href="http://www.ncbi.nlm.nih.gov/geo/"),
                          " RNA-Seq results from ",
                          a("Himes Lab,", href="http://himeslab.org/asthmagenes/"),
                          " transcript data from ",
                          a("GENCODE,", href="http://www.gencodegenes.org/"), 
                          " genome-wide association study (GWAS) results from ",
                          a("EVE", href="http://eve.uchicago.edu/"), ", ",
                          a("GABRIEL", href="https://www.cng.fr/gabriel/index.html"), "and ",
                          a("GRASP", href="https://grasp.nhlbi.nih.gov/Search.aspx"), 
                          " and glucocorticoid receptor binding sites from ",
                          a("UCSC's ENCODE,", href="https://genome.ucsc.edu/ENCODE/"),
                          " allowing researchers to access a breadth of information with a click. ",  
                          "We hope REALGAR's disease-specific and tissue specific information ",
                          "will facilitate prioritization and experiment design, leading to clinically actionable insights."), br(),
                        
                        wellPanel(fluidRow(align = "left",
                                           column(12,  
                                                  
                                                  fluidRow(h4(strong("Select options:")),align = "left"), 
                                                  fluidRow(div(style="margin-right: 25px", column(5, fluidRow(checkboxGroupInput(inputId="Tissue", label="Tissue", choices=c("Airway smooth muscle"="ASM", "Bronchial epithelium"="BE", 
                                                                                                                                             "Bronchoalveolar lavage"="BAL", "CD4"="CD4", "CD8"="CD8",
                                                                                                                                             "Lens epithelium" = "LEC","Lymphoblastic leukemia cell" = "chALL", 
                                                                                                                                             "Lymphoblastoid cell" = "LCL","Macrophage" = "MACRO", "MCF10A-Myc" = "MCF10A-Myc",
                                                                                                                                             "Nasal epithelium"="NE","Osteosarcoma U2OS cell" = "U2O", 
                                                                                                                                             "Peripheral blood mononuclear cell"="PBMC","White blood cell"="WBC","Whole lung"="Lung"),selected="BE", inline = TRUE), align="left"),
                                                           fluidRow(actionButton("selectall","Select/unselect all tissues", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), align="left"))),
                                                          column(2, fluidRow(checkboxGroupInput(inputId="Asthma", label="Asthma Type", choices=c("Allergic asthma"="allergic_asthma",
                                                                                                                                          "Severe asthma"="severe_asthma","Asthma"="asthma","Mild asthma"="mild_asthma","Non-allergic asthma"=
                                                                                                                                              "non_allergic_asthma", "Asthma and rhinitis"="asthma_and_rhinitis"), selected="asthma"), inline = TRUE), align="left"),
                                                  column(2, fluidRow(radioButtons(inputId="GC_included", label="Treatment", 
                                                                                   choice = c("Glucocorticoid treatment" = "GC", "No treatment" = ""), select = "")),
                                                         fluidRow(checkboxGroupInput(inputId="which_SNPs", label="GWAS results to include", choices=c("EVE"="snp_eve_subs", "GABRIEL"="snp_gabriel_subs", "GRASP"="snp_subs"), selected=c("snp_eve_subs", "snp_gabriel_subs", "snp_subs"))), align="left"),
                                                  column(3, fluidRow(textInput(inputId="curr_gene",label="Type official gene symbol or SNP ID:", value= "GAPDH"),
                                                                  tags$head(tags$style(type="text/css", "#curr_gene {width: 190px}"))),
                                                         fluidRow(img(src="http://shiny.rstudio.com/tutorial/lesson2/www/bigorb.png", height=35, width=38), "Created with RStudio's ", a("Shiny", href="http://www.rstudio.com/shiny")), align="center")),
                                                  hr()), 
                                           fluidRow(column(12, h4(strong("Data download:")),align = "left")),
                                           fluidRow(column(12, h5("The results displayed in the forest plots and gene tracks below may also be downloaded directly here:"))), br(),
                                           fluidRow(column(6, downloadButton(outputId="table_download", label="Download fold change and p-value results displayed in forest plots below"), align="left"),
                                                    column(6, downloadButton(outputId="SNP_data_download", label="Download GWAS results displayed in gene tracks below"), align="left")))),
                                           #,
                                           # column(7,align="left",
                                           #        fluidRow(h4("Datasets loaded:")),
                                           #        fluidRow(p("Click on links to access the datasets and studies referenced in the table.")),
                                           #        fluidRow(DT::dataTableOutput(outputId="GEO_table")))

                        
                        
                        mainPanel(hr(),
                                  fluidRow(br(),
                                           column(10, plotOutput(outputId="forestplot_asthma",width="900px", height="650px"), align="left"), 
                                           div(style="margin-top: 45px", column(2, imageOutput("color_scale1"), align="right")), # margin-top needed to align color scale w/ forest plot
                                           column(6, downloadButton(outputId="asthma_fc_download",label="Download asthma forest plot"), align="center"),
                                           column(10, conditionalPanel(condition = "input.GC_included == 'GC'", plotOutput(outputId="forestplot_GC",width="900px", height="650px")),align="left"),
                                           column(2, div(style="margin-top: 45px", conditionalPanel(condition = "input.GC_included == 'GC'", imageOutput("color_scale2")), align="right")), # margin-top needed to align color scale w/ forest plot
                                           column(6, conditionalPanel(condition = "input.GC_included == 'GC'", downloadButton(outputId="GC_fc_download",label="Download GC forest plot")), align="center")),
                                           
                                  fluidRow(br()
                                      # column(12, h5(strong("Results Used to Create Plots Above")), align = "center")
                                      ),
                                  
                                  fluidRow(
                                           # column(10, offset= 1, DT::dataTableOutput(outputId="tableforgraph"), align="center")
                                           # column(12, downloadButton(outputId="table_download", label="Download table"), align="center")
                                           ),br(), hr(), width = 12,
                                  
                                  fluidRow(column(12, p("Transcripts for the selected gene are displayed here. ",
                                                        "Any SNPs and/or GR binding sites that fall within the gene ",
                                                        "or within +/- 10kb of the gene are also displayed, ",
                                                        "each in a separate track. GR binding sites are colored by the ",
                                                        "ENCODE binding score, ",
                                                        "with the highest binding scores corresponding to the darkest color. ",
                                                        "SNPs are colored by p-value, with the lowest p-values corresponding to the darkest color.",
                                                        "All SNP p-values are <=0.05 and are obtained directly",
                                                        "from the study in which the association was published.")),
                                           column(12, downloadButton(outputId="gene_tracks_download", label="Download gene tracks"), align="center"), br(),
                                           column(12, align="center", plotOutput(outputId="gene_tracks_outp2"), br(), br())))))
