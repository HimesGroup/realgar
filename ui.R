#Install ui libraries in path /opt/R/4.0.2/lib/R/library

#Footer
footstyle = "position:fixed;
       bottom:0;
right:0;
left:0;
background:#EBF4FA;
padding:5px;
box-sizing:border-box;"

ui <- shinyUI(fluidPage(theme = shinytheme("lumen"),
                        h1(strong("Reducing Associations by Linking Genes And omics Results (REALGAR)"), align="center", style = "color: #9E443A;"), hr(),
                        p("REALGAR is a tissue-specific, disease-focused resource that integrates results from omics studies", 
                          " to facilitate the design of experiments to validate genome-wide association study (GWAS) findings and generate novel hypotheses.",
                          " Omics results are currently available for gene expression microarrays,",
                          " RNA-Seq, and ChIP-Seq datasets.",
                          "To use REALGAR, input an official gene symbol or SNP ID, and select tissues, conditions, treatments,",
                          " and GWAS of interest. The 'Omics Results’ tab allows you to visualize and ", 
                          "download results for the studies matching your selection criteria. ",
                          "The ‘Datasets loaded’ tab provides more information about the datasets selected. ",
                          "The ‘About’ tab provides additional information about the app."),br(),                        
                        navbarPage( "", 
                                    tabPanel(h4(strong("Omics Results"), style = "color: #9E443A;"),
                                             wellPanel(fluidRow(align = "left",
                                                                column(12, # change mk
                                                                       fluidRow(
                                                                         tags$head(
                                                                           tags$style(type="text/css", "#inline label{ display: table-cell; text-align: center; vertical-align: middle; }
                                                                                      #inline .form-group { display: table-row;}")
                                                                           #tags$style(HTML(.shiny-input-container {margin-right: 0;}))
                                                                         ),
                                                                         tags$div(id = "inline", column(12,
                                                                                                         h4(textInput(inputId = "current", label = strong(paste0("Official Gene Symbol or SNP ID:", stri_dup(intToUtf8(160), 6))), value="FKBP5"))))
                                                                )), br(), br(), br(), br(),
                                                                column(12,
                                                                       fluidRow(div(style="margin-right: 25px", 
                                                                                    column(2,fluidRow(div(style="margin-left: 25px", h4(pickerInput("Tissue", label = strong("Tissue"),
                                                                                           choices = tissue_choices, selected=tissue_choices, multiple = TRUE,
                                                                                           options = list(`actions-Box` = TRUE, `none-Selected-Text` = "None",
                                                                                                          `selected-Text-Format`= "count")))))),  # change mk
                                                                                    #column(1, withSpinner(uiOutput("loadProxy"),color= "#9E443A",type=8)),
                                                                                    column(1, HTML("<br>")),
                                                                                    column(2, fluidRow(h4(pickerInput(inputId="Treatment", label=strong("Treatment"), choices = treatment_choices, selected=treatment_selected, multiple = TRUE,
                                                                                                                   options = list(`actions-Box` = TRUE, `none-Selected-Text` = "None",`selected-Text-Format`= "count")))), # change mk
                                                                                           br()),
                                                                                    column(1, HTML("<br>")),
                                                                                    column(2,fluidRow(h4(pickerInput(inputId="Asthma", label=strong("Condition"), choices=asthma_choices, selected=asthma_selected, multiple = TRUE, 
                                                                                                                     options = list(`actions-Box` = TRUE, `none-Selected-Text` = "None",`selected-Text-Format`= "count")))),
                                                                                           br()),
                                                                                    column(1, HTML("<br>")),
                                                                                    column(2,
                                                                                          fluidRow(h4(pickerInput(inputId="which_SNPs", label=strong("GWAS Results"), choices = gwas_choices, selected=gwas_selected, multiple = TRUE,
                                                                                                                  options = list(`actions-Box` = TRUE, `none-Selected-Text` = "None",`selected-Text-Format`= "count")))), br() # change mk
                                                                                          )
                                                                      )),
                                                                      ),)),
                                             fluidRow(
                                               column(12,
                                                      tags$p(HTML(paste0('Once you have made above selection, click Go. Data downloads and figure displays rely on the parameters specified ', em('before'), ' the most recent time Go was clicked. The intial figure display and data downloads correspond to the initial parameters shown.'))),
                                                      align = 'center')
                                             ),
                                             
                                             fluidRow(
                                               column(12, actionButton("go", strong("Go")), align = 'center'),
                                             ),
                                             
                                                                                                                                                          
                                                                column(width=12,
                                                                       
                                                                                  tabsetPanel(tabPanel(strong("Differential Gene Expression"),br(),
                                                                                            p("Differential expression results for individual microarray and RNA-Seq studies were obtained using ",
                                                                                              a("RAVED", href="https://github.com/HimesGroup/raved", target="_blank"),
                                                                                              ". Meta-analysis results displayed were obtained using three summary statistics-based approaches: ",
                                                                                              "(1) an effect size-based method that applies a random-effects model, ",
                                                                                              "(2) a p-value-based method that applies Fisher's sum-of-logs method, ",
                                                                                              "and  (3) a rank-based method that adopts the Rank Product. ",
                                                                                              "Note that p-values from the p-value-based and effect size-based methods are not adjusted for multiple comparisons corrections. ",
                                                                                              "We suggest that users apply a stringent multiple comparisons corrections threshold to assign significance when using those methods ",
                                                                                              HTML(paste0("(e.g., correct for tests in 25,000 genes, or a threshold of 2x10",tags$sup(-6),". ")),
                                                                                              "For the rank-based method, an analytic rank product is provided instead of a permutated empirical p-value, ",
                                                                                              "so we suggest that users refer to the rank score when prioritizing the genes for functional validation. ",
                                                                                              "More information can be found in the 'Datasets loaded' and 'About' tabs."),
                                                                                            fluidRow(column(12, downloadButton(outputId="table_download", label="Download gene expression results"), align="center")),
                                                                                            br(),
                                                                                            tabsetPanel(
                                                                                                    tabPanel(strong("Treatments"),br(),
                                                                                                      column(12, withSpinner(plotOutput(outputId="forestplot_GC",width="1200px", height="auto"),color= "#9E443A"),align="center"),
                                                                                                      column(12, textOutput("GC_pcomb_text"), align="center"), # output combined p-values
                                                                                                      column(12, downloadButton(outputId="GC_fc_download",label="Download GC forest plot"), align="center"),
                                                                                                      column(12, fluidRow(br(), br(), br()))),

                                                                                                    tabPanel(strong("Pollutants"),br(),
                                                                                                      column(12, withSpinner(plotOutput(outputId="forestplot_cig",width="1200px", height="auto"),color= "#9E443A"),align="center"),
                                                                                                      column(12, textOutput("cig_pcomb_text"), align="center"), # output combined p-values
                                                                                                      column(12, downloadButton(outputId="cig_fc_download",label="Download pollutant forest plot"), align="center"),
                                                                                                      column(12, br(), br(), br())),
                                                                                                      
                                                                                                      tabPanel(strong("Conditions"),br(),                    
                                                                                                        column(12, withSpinner(plotOutput(outputId="forestplot_asthma",width="1200px", height="auto"),color= "#9E443A"), align="center"), #1355px #1250px #1000px
                                                                                                        column(12, textOutput("asthma_pcomb_text"), align="center"), # output combined p-values
                                                                                                        column(12, downloadButton(outputId="asthma_fc_download",label="Download asthma forest plot"), align="center"),
                                                                                                        column(12, fluidRow(br(), br(), br()))))
                                                                                            ),

                                                                                                                               
                                                                                   tabPanel(strong("Gene Expression Levels"),
                                                                                            sidebarLayout(
                                                                                              sidebarPanel(
                                                                                                column(12, selectInput("tissue_te", label = "RNA-Seq Dataset (tissue)", choices = rnaseq_choices, selected="SRP033351")),
                                                                                                br(),
                                                                                                p("Normalized read counts of each gene were obtained with DESeq2 based on the HTSeq matrix for the selected dataset."),
                                                                                                br()),
                                                                                            mainPanel(  
                                                                                              plotOutput("GeneBoxPlot"),
                                                                                              br(),br(),br(),br(), br(), br(),
                                                                                              uiOutput("studyText"),
                                                                                              downloadButton('downloadPic', 'Download Figure'),br(),br(),
                                                                                              htmlOutput("table_title", align="center"),
                                                                                              dataTableOutput("diffResults"),
                                                                                            br()))),
                                                                                   
                                                                                   tabPanel(strong("Gene Tracks"), br(),
                                                                                            p("Reference transcript for the selected gene is displayed along with any SNPs and/or ",
                                                                                              "transcription factor binding sites that fall within 20kb of the gene boundaries. ",
                                                                                              "Transcription factor binding site tracks contain ChIP-Seq results obtained with ",
                                                                                              a("brocade", href="https://github.com/HimesGroup/brocade", target="_blank"),  "for the ",
                                                                                              "Glucocorticoid Receptor (GR) and RNA Polymerase 2 (RNAP2) under dexamethasone (dex) or vehicle control (EtOH) exposure conditions. ",
                                                                                              "SNPs from GWAS studies meeting the selected p-value threshold are colored by p-value reported by the original study, ",
                                                                                              "with the lowest p-values corresponding to the brightest colors. ",
                                                                                              "cis-eQTL results for whole blood, lung, and skeletal muscle tissues were obtained from the ",
                                                                                              a("GTEx V8 release", href="https://gtexportal.org/home/datasets", target="_blank"),
                                                                                              ". Glucocorticoid  response element (GRE) motifs within Glucocorticoid Receptor (GR) binding sites were identified with",
                                                                                              a("the FIMO tool", href="https://meme-suite.org/meme/doc/fimo.html", target="_blank"),
                                                                                              ". More information can be found in the 'Datasets loaded' and 'About' tabs."),br(),
                                                                                            column(12,
                                                                                                   fluidRow(
                                                                                                     tags$head(
                                                                                                       tags$style(type="text/css", "#inline label{ display: table-cell; text-align: center; vertical-align: middle; }
                                                                                                               #inline .form-group { display: table-row;}")
                                                                                                     ),
                                                                                                     tags$div(id = "inline", selectInput("pval_thr", label = paste0("GWAS p-value threshold:", stri_dup(intToUtf8(160), 6)), choices = pval_for_select, selected=pval_select, width = '400px')),
                                                                                                     tags$head(
                                                                                                       tags$style(HTML('
                                                                                                        .selectize-input {
                                                                                                        white-space: nowrap;
                                                                                                      }
                                                                                                      #pval_thr + div>.selectize-dropdown{
                                                                                                      width: 660px !important;
                                                                                                      }
                                                                                                      '
                                                                                                       )
                                                                                                       )
                                                                                                     )),
                                                                                                   align="left"),
                                                                                            column(12, withSpinner(plotOutput("karyoPlot",height = "1300px",width="1200px")),color= "#9E443A",align="center"), br(),
                                                                                            column(12, 
                                                                                                   fluidRow(
                                                                                                     #column(4, downloadButton(outputId="gene_tracks_download", label="Download gene track figure"), align="center"),
                                                                                                     column(4, downloadButton(outputId="SNP_data_download", label="Download GWAS results"), align="center"),
                                                                                                     column(4, downloadButton(outputId="GRbinding_data_download", label="Download GR-binding results"), align="center"),
                                                                                                     column(4, downloadButton(outputId="GRE_data_download", label="Download putative GRE motifs"), align="center"))),
                                                                                            column(12, fluidRow(br(), br(), br())))  #tabPanel #1200px
                                                                                   
                                                                                   
                                                                                   
                                                                                   )#tabsetPanel

                                                                                   )
                                             ),
                                                                              
                                                     tabPanel(h4(strong("Datasets loaded"), style = "color: #9E443A;"),
                                                                 column(12,align="left",
                                                                        fluidRow(h4("Gene expression")),
                                                                        fluidRow(p("The following information is provided in each column: "),
                                                                                 p("Dataset: GEO accession numbers that link to corresponding GEO entries"),
                                                                                 p("PMID: PubMed IDs for papers, when available, that link to corresponding PubMed entries"),
                                                                                 p("Report: Quality Control report generated by RAVED"),
                                                                                 p("Description: brief overview of study design"), br()),
                                                                        fluidRow(DT::dataTableOutput(outputId="GEO_table"), br()),
                                                                        
                                                                        fluidRow(h4("Transcription factor binding sites")),
                                                                        fluidRow(p("The following information is provided in each column: "),
                                                                                 p("Dataset: GEO accession numbers that link to corresponding GEO entries"),
                                                                                 p("PMID: PubMed IDs for papers, when available, that link to corresponding PubMed entries"),
                                                                                 p("Report: Quality Control report generated by brocade"),
                                                                                 p("Description: brief overview of study design"), br()),
                                                                        fluidRow(DT::dataTableOutput(outputId="chipseq_table"), br()),
                                                                        
                                                                        fluidRow(h4("GWAS")),
                                                                        fluidRow(p("The following information is provided in each column: "),
                                                                                 p("Dataset: name of the study, which links to the study website or publication"),
                                                                                 p("Description: brief overview of GWAS study design")),
                                                                        fluidRow(DT::dataTableOutput(outputId="GWAS_table"), br())),
                                                              tags$style(HTML("
                                                                p
                                                                  {
                                                                    margin: 0px;
                                                                    padding: 1px;
                                                                  }
                                                              "))),
                                                 
                                    tabPanel(h4(strong("About"), style = "color: #9E443A;"),
                                             p("This Rstudio ", a("Shiny", href="https://www.rstudio.com/shiny", target="_blank"), " app was created by members of the ",
                                               a("Himes Lab", href="https://himeslab.org/", target="_blank"), ".",
                                               "Full code is available in ", a("GitHub", href="https://github.com/HimesGroup/realgar", target="_blank"), "."), br(),
                                             h4("Contributors"),
                                             p("Coding and design by Maya Shumyatcher, Mengyuan Kan, Avantika Diwadkar and Blanca Himes"), 
                                             p("Analysis of individual datasets by Avantika Diwadkar, Jaehyun Joo, Mengyuan Kan, Nisha Narayanan, Supriya Saxena, Haoyue Shuai, Maya Shumyatcher, Gabriel Soliman"), br(),
                                             h4("References"),
                                             p("If you use REALGAR in your research, please cite the following papers:"),
                                             p("Kan M, Diwadkar AR, Saxena S, Shuai H, Joo J, Himes BE.", a("REALGAR: a web app of integrated respiratory omics data.", href="REALGAR_bioinformatics.pdf", target="blank"), "Bioinformatics. 2022 Sep 15;38(18):4442-4445.(PMID:",a("35863045",href="https://www.ncbi.nlm.nih.gov/pubmed/35863045", target="_blank"), ")."),
                                             p("Shumyatcher M, Hong R, Levin J, Himes BE.", a("Disease-Specific Integration of Omics Data to Guide Functional Validation of Genetic Associations.",href="REALGAR_resubmission_07.06.2017_FINAL_Corrected.pdf", target="_blank"),"AMIA Annu Symp Proc. 2018;2017:1589–1596. 
                        Published 2018 Apr 16. (PMID:",a("29854229).",href="https://www.ncbi.nlm.nih.gov/pubmed/29854229", target="_blank")), br(),
                                             p("Analysis of gene expression microarray and RNA-Seq data was performed as described in this paper:"),
                                             p("Kan M, Shumyatcher M, Diwadkar A, Soliman G, Himes BE.", a("Integration of Transcriptomic Data Identifies Global and 
                         Cell-Specific Asthma-Related Gene Expression Signatures.",href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6371257/pdf/2973744.pdf",target="_blank"),"AMIA Annu Symp Proc. 2018;2018:1338–1347. 
                         Published 2018 Dec 5. (PMID:", a("30815178)",href="https://www.ncbi.nlm.nih.gov/pubmed/30815178", target="_blank"),
                                               ". You can refer to the code",a("here.",href="https://github.com/HimesGroup/raved", target="_blank")), br(),
                                             p("Analysis of ChIP-Seq data was performed as described in this paper:"),
                                             p("Diwadkar AR, Kan M, Himes BE.", a("Facilitating Analysis of Publicly Available ChIP-Seq Data for Integrative Studies.",href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7153109/pdf/3200301.pdf",target="_blank"),"AMIA Annu Symp Proc. 2019;2019:371-379. 
                         Published 2020 Mar 4. (PMID:", a("32308830)",href="https://www.ncbi.nlm.nih.gov/pubmed/32308830", target="_blank"),
                                               ". You can refer to the code",a("here.",href="https://github.com/HimesGroup/brocade", target="_blank")), br(),
                                             h4("Updates"),
                                             p("2022-01 Added Gene-expression levels tab"),
                                             p("2021-12 Redesigned integrative gene tracks"),
                                             p("2021-08 Added initial PAH and PM transcriptomic datasets"),
                                             p("2020-09 Incorporated UKBB asthma and COPD GWAS results"),
                                             p("2020-08 Added initial cigarette and e-cig transcriptomic datasets"),
                                             p("2019-08 Added initial PDE inhibitor transcriptomic datasets"),
                                             p("2019-01 Incorporated glucocorticoid receptor ChIP-Seq datasets"),
                                             p("2018-08 Added initial Beta2-agonist transcriptomic datasets"),
                                             p("2018-02 Incorporated TAGC and Ferreira asthma GWAS results"),
                                             p("2018-01 Added reactive meta-analysis scores for transcriptomic datasets"),
                                             p("2017-08 Published initial version of REALGAR"),
                                             br(), br(), br()
                                             )
                                    ),
                                      tags$style(HTML("
                                        #.navbar-default .navbar-brand {color:white;}
                                        #.navbar-default .navbar-brand:hover {color:white;}
                                        #.navbar { background-color:red;} # background color for the whole panel
                                        #.navbar-default .navbar-nav > li > a {color:white;}
                                        .navbar-default .navbar-nav > .active > a,
                                        .navbar-default .navbar-nav > .active > a:focus,
                                        .navbar-default .navbar-nav > .active > a:hover {background-color:	#ebebeb;} # grey 92
                                        #.navbar-default .navbar-nav > li > a:hover {background-color:#e0e0e0;text-decoration}

                                    "))                        
                        
                        )) # change mk
                                                   #Footer
                                                    # tags$footer(img(src="RStudio-Logo-final.png", height=124*.65, width=110.4*.65), 
                                                    #             "Created with RStudio's ", a("Shiny", href="https://www.rstudio.com/shiny", target="_blank"),
                                                    #             p("This app was created by Maya Shumyatcher, Mengyuan Kan, Avantika Diwadkar and Blanca Himes at the ", a("Himes Lab.", href="https://himeslab.org/", target="_blank")),style = footstyle)))))