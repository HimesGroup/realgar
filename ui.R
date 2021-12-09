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
                        p("REALGAR is a tissue-specific, disease-focused resource for integrating omics results. ",
                          "This app brings together genome-wide association study (GWAS) results, "," transcript data from ",
                          a("GENCODE,", href="https://www.gencodegenes.org/", target="_blank"),
                          "ChIP-Seq data, microarray gene expression and gene-level RNA-Seq results from the ", 
                          a("Gene Expression Omnibus (GEO).", href="https://www.ncbi.nlm.nih.gov/geo/", target="_blank"),
                          " REALGAR facilitates prioritization of genes and experiment design of functional validation studies."),
                        p("To use REALGAR, input an official gene symbol or SNP ID, and select tissues, phenotypes, treatments/exposures",
                          " and GWAS of interest. The 'Results' tab allows you to visualize and download ", 
                          "results for the studies matching your selection criteria. ",
                          "The 'Datasets loaded' tab provides more information about the datasets selected. "),br(),                        
                        navbarPage( "", 
                                    tabPanel(h4(strong("Results"), style = "color: #9E443A;"),
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
                                                                                    column(2,fluidRow(h4(pickerInput(inputId="Asthma", label=strong("Condition"), choices=asthma_choices, selected=asthma_selected, multiple = TRUE, 
                                                                                                                  options = list(`actions-Box` = TRUE, `none-Selected-Text` = "None",`selected-Text-Format`= "count")))),
                                                                                           br()),
                                                                                    column(1, HTML("<br>")),
                                                                                    column(2, fluidRow(h4(pickerInput(inputId="Treatment", label=strong("Treatment"), choices = treatment_choices, selected=treatment_selected, multiple = TRUE,
                                                                                                                   options = list(`actions-Box` = TRUE, `none-Selected-Text` = "None",`selected-Text-Format`= "count")))), # change mk
                                                                                           br()),
                                                                                    column(1, HTML("<br>")),
                                                                                    column(2,
                                                                                          fluidRow(h4(pickerInput(inputId="which_SNPs", label=strong("GWAS Results"), choices = gwas_choices, selected=c("GRASP"="snp_subs", "UKBiobank Asthma"="snp_UKBB_asthma_subs"), multiple = TRUE,
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
                                             
                                                                                                                                                          
                                                                column(width=12,tabsetPanel(tabPanel(strong("Forestplots"),br(),
                                                                                            p("Differential expression results for individual microarray and RNA-Sequencing study were obtained using ",
                                                                                              a("RAVED", href="https://github.com/HimesGroup/raved", target="_blank"),
                                                                                              ". Integrated results were obtained using three summary statistics-based approches: ",
                                                                                              "(1) an effect size-based method that applied a random-effects model, ",
                                                                                              "(2) p-value-based method that applied Fisher's sum-of-logs method, ",
                                                                                              "and  (3) rank-based method that adpoted the Rank Product. ",
                                                                                              "Note that p-values from the p-value-based and effect size-based methods are not adjusted for multiple corrections in this app, ",
                                                                                              "so we suggest that users apply a stringent threshold of multiple corrections corresponding to 25,000 genes (i.e. 2x10-6). ",
                                                                                              "For the rank-based method, an analytic rank product is provided instead of the permutated empirical p-value, ",
                                                                                              "so we suggest users refer to the rank score when prioritizing the genes for functional validation. For more information, you can check the References tab."),
                                                                                            fluidRow(column(12, downloadButton(outputId="table_download", label="Download gene expression results"), align="center")),
                                                                                            br(),
                                                                                        tabsetPanel(tabPanel(strong("Asthma"),br(),                    
                                                                                                             # column(12, withSpinner(plotOutput(outputId="forestplot_asthma",width="1200px", height="auto"),color= "#9E443A"), align="center"), #1355px #1250px #1000px
                                                                                                             # column(12, textOutput("asthma_pcomb_text"), align="center"), # output combined p-values
                                                                                                             # hr(), hr(),
                                                                                                             # column(12, downloadButton(outputId="asthma_fc_download",label="Download asthma forest plot"), align="center"),
                                                                                                             # column(12, fluidRow(br(), br(), br()))),
                                                                                                    fluidRow(column(12, withSpinner(plotOutput(outputId="forestplot_asthma",width="1200px", height="auto"),color= "#9E443A"), align="center")), #1355px #1250px #1000px
                                                                                                    fluidRow(column(12, textOutput("asthma_pcomb_text"), align="center")), # output combined p-values
                                                                                                    fluidRow(
                                                                                                      #column(6, downloadButton(outputId="asthma_fc_download",label="Download asthma forest plot"), align="center"),
                                                                                                      column(12, downloadButton(outputId="asthma_fc_download",label="Download asthma forest plot"), align="center")),
                                                                                                      #column(12, downloadButton(outputId="table_download",label="Download gene expression results"), align="center")),
                                                                                                    fluidRow(column(12, br(), br(), br()))),                                                                                                    
                                                                                                    tabPanel(strong("Treatment"),br(),
                                                                                                      column(12, withSpinner(plotOutput(outputId="forestplot_GC",width="1200px", height="auto"),color= "#9E443A"),align="center"),
                                                                                                      column(12, textOutput("GC_pcomb_text"), align="center"), # output combined p-values
                                                                                                      column(12, downloadButton(outputId="GC_fc_download",label="Download GC forest plot"), align="center"),
                                                                                                      column(12, fluidRow(br(), br(), br()))),
                                                                                                    #tabPanel(strong("Pollutant"),br(),
                                                                                                    #         column(12, withSpinner(plotOutput(outputId="forestplot_cig",width="1200px", height="auto"),color= "#9E443A"),align="center"),
                                                                                                    #         column(12, textOutput("cig_pcomb_text"), align="center"), # output combined p-values
                                                                                                    #         column(12, downloadButton(outputId="cig_fc_download",label="Download smoking forest plot"), align="center"),
                                                                                                    #         column(12, fluidRow(br(), br(), br()))),
                                                                                                    tabPanel(strong("Pollutants"),br(),
                                                                                                      column(12, withSpinner(plotOutput(outputId="forestplot_cig",width="1200px", height="auto"),color= "#9E443A"),align="center"),
                                                                                                      column(12, textOutput("cig_pcomb_text"), align="center"), # output combined p-values
                                                                                                      column(12, downloadButton(outputId="cig_fc_download",label="Download smoking forest plot"), align="center"),
                                                                                                      column(12, fluidRow(br(), br(), br()))))),
                                                                                                                       
                                                                                   # tabPanel("Gene Tracks", br(),
                                                                                   #          p("Transcripts for the selected gene are displayed here. ",
                                                                                   #            "Any SNPs and/or transcription factor binding sites that fall within the gene ",
                                                                                   #            "or within +/- 10kb of the gene are also displayed, ",
                                                                                   #            "each in a separate track. Transcription factor binding sites are colored by the ",
                                                                                   #            "ENCODE binding score (shown below each binding site), ",
                                                                                   #            "with the highest binding scores corresponding to the brightest colors. ",
                                                                                   #            "Only those SNPs with p-value <= 0.05 are included. ",
                                                                                   #            "SNPs are colored by p-value, with the lowest p-values corresponding to the brightest colors. ",
                                                                                   #            "All SNP p-values are obtained directly from the study in which the association was published."),br(),
                                                                                   #               column(12, downloadButton(outputId="gene_tracks_download", label="Download gene tracks"), align="center"), br(),
                                                                                   #               column(12, HTML("<div style='height: 90px;'>"), imageOutput("color_scale3"), align="center", HTML("</div>")), 
                                                                                   #               column(12, align="center", plotOutput(outputId="gene_tracks_outp2")),
                                                                                   #               column(12, fluidRow(br(), br(), br(),br()))),
                                                                                                 # )#tabsetPanel
                                                                                                 # )),
                                                                                  tabPanel(strong("Gene Tracks"), br(),
                                                                                            p("Transcripts for the selected gene are displayed here. ",
                                                                                              "Any SNPs and/or transcription factor binding sites that fall within the gene ",
                                                                                              "or within +/- 20kb of the gene are also displayed, ",
                                                                                              "each in a separate track. Transcription factor binding sites are colored by the ",
                                                                                              "adjusted p-value from the analysis, ",
                                                                                              "with the lowest p-values corresponding to the brightest colors. ",
                                                                                              "Only those SNPs with p-value <= 0.05 are included. ",
                                                                                              "SNPs are colored by p-value, with the lowest p-values corresponding to the brightest colors. ",
                                                                                              "All SNP p-values are obtained directly from the study in which the association was published."),br(),
                                                                                            #column(12, downloadButton(outputId="gene_tracks_download", label="Download gene track figure"), align="center"), br(),
                                                                                            #column(12, HTML("<div style='height: 90px;'>"), imageOutput("color_scale3"), align="center", HTML("</div>")), 
                                                                                            #column(12, withSpinner(plotOutput("karyoPlot",height = "1300px",width="1200px")),color= "#9E443A"),align="center"), #tabPanel #1200px
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
                                                                                                    )
                                                                                                    ),
                                                                                                    align="left"),
                                                                                           column(12, withSpinner(plotOutput("karyoPlot",height = "1300px",width="1200px")),color= "#9E443A",align="center"), br(),
                                                                                           column(12, 
                                                                                                  fluidRow(
                                                                                                    column(4, downloadButton(outputId="gene_tracks_download", label="Download gene track figure"), align="center"),
                                                                                                    column(4, downloadButton(outputId="SNP_data_download", label="Download GWAS results"), align="center"),
                                                                                                    column(4, downloadButton(outputId="GRbinding_data_download", label="Download GR-binding results"), align="center"))),
                                                                                           column(12, fluidRow(br(), br(), br()))),  #tabPanel #1200px
                                                                                   tabPanel(strong("Gene Expression Levels"),
                                                                                            sidebarLayout(
                                                                                              sidebarPanel(
                                                                                                column(12, selectInput("tissue_te", label = "Dataset (tissue)", choices = rnaseq_choices, selected="SRP033351")),
                                                                                                br(),
                                                                                                p("Transcripts for the selected gene in all available RNA-Seq studies are displayed here."),
                                                                                                p("Read counts of each gene is obtained by HTSeq and normalized by",
                                                                                                  a("DESeq2.", href="https://www.ncbi.nlm.nih.gov/pubmed/25516281", target="_blank")),
                                                                                                #verbatimTextOutput('gene_name_out'), # change mk
                                                                                              br()),
                                                                                            mainPanel(  
                                                                                              plotOutput("GeneBoxPlot"),
                                                                                              br(),br(),br(),br(), br(), br(),
                                                                                              uiOutput("studyText"),
                                                                                              downloadButton('downloadPic', 'Download Figure'),br(),br(),
                                                                                              htmlOutput("table_title", align="center"),
                                                                                              dataTableOutput("diffResults"),
                                                                                            br())
                                                                                   ))
                                                                                   )#tabsetPanel
                                                                                   )),
                                                                              
                                                     tabPanel(h4(strong("Datasets loaded"),style = "color: #9E443A;"),
                                                                 column(12,align="left",
                                                                        fluidRow(h4("Gene expression datasets:")),
                                                                        fluidRow(p("For gene expression datasets, the following information is provided:",
                                                                                   "(1) GEO accession numbers that link directly to GEO entries, ", 
                                                                                   "(2) Quality control report generated by RAVED",
                                                                                   "(3) PMIDs for papers, when available, that link directly to PubMed entries, and ", 
                                                                                   "(4) brief descriptions for all gene expression datasets that match selected criteria. "), br()),
                                                                        fluidRow(DT::dataTableOutput(outputId="GEO_table"), br()),
                                                                        fluidRow(h4("GWAS datasets:")),
                                                                        fluidRow(p("For GWAS datasets, the following information is provided:",
                                                                                   "(1) names of the studies selected, which link directly to the study website or publication and",
                                                                                   "(2) brief descriptions for all GWAS studies selected.")),
                                                                        fluidRow(DT::dataTableOutput(outputId="GWAS_table"), br()))),
                                    
                                                # tabPanel("Transcriptomic Explorer",
                                                #     sidebarLayout(
                                                #     sidebarPanel(
                                                #     p("Expression results for individual genes using publicly available RNA-Seq datasets of interest."),
                                                #     br(),
                                                #     selectInput("tissue_te", label = "Tissue:",
                                                #                 choices = rnaseq_choices, 
                                                #                 selected="SRP033351"),
                                                #     #uiOutput("genesAvail_te"), tags$head(tags$style(type="text/css", "#curr_gene_te {width: 190px}")),
                                                #     br(),
                                                #     p("Read counts of each gene is obtained by HTSeq and normalized by",
                                                #       a("DESeq2.", href="https://www.ncbi.nlm.nih.gov/pubmed/25516281", target="_blank"),
                                                #       "Differential expression results were obtained by DESeq2. Genes with a total read count <10 were filtered out before differential analysis."),
                                                #     br(),br(),br(),br()),
                                                # 
                                                #     mainPanel(p(""),
                                                #       plotOutput("GeneBoxPlot"),
                                                #       br(),br(),br(),br(),br(),
                                                #       uiOutput("studyText"),
                                                #       downloadButton('downloadPic', 'Download Figure'),br(),br(),
                                                #       htmlOutput("table_title", align="center"),
                                                #       dataTableOutput("diffResults"),br()))),
                                                 
                                                 tabPanel(h4(strong("References"),style = "color: #9E443A;"),
                                                   p("If you use REALGAR in your research, please cite the following papers:"),
                                                   p("Shumyatcher M, Hong R, Levin J, Himes BE.", a("Disease-Specific Integration of Omics Data to Guide Functional Validation of Genetic Associations.",href="REALGAR_resubmission_07.06.2017_FINAL_Corrected.pdf", target="_blank"),"AMIA Annu Symp Proc. 2018;2017:1589–1596. 
                        Published 2018 Apr 16. (PMID:",a("29854229).",href="https://www.ncbi.nlm.nih.gov/pubmed/29854229", target="_blank"),"You can refer to the code",a("here.",href="https://github.com/HimesGroup/taffeta",target="_blank")),
                                                   p("Kan M, Shumyatcher M, Diwadkar A, Soliman G, Himes BE.", a("Integration of Transcriptomic Data Identifies Global and 
                         Cell-Specific Asthma-Related Gene Expression Signatures.",href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6371257/pdf/2973744.pdf",target="_blank"),"AMIA Annu Symp Proc. 2018;2018:1338–1347. 
                         Published 2018 Dec 5. (PMID:", a("30815178).",href="https://www.ncbi.nlm.nih.gov/pubmed/30815178", target="_blank"),
                                                     "You can refer to the code",a("here.",href="https://github.com/HimesGroup/raved", target="_blank")),
                                                   p("Diwadkar AR, Kan M, Himes BE.", a("Facilitating Analysis of Publicly Available ChIP-Seq Data for Integrative Studies.",href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7153109/pdf/3200301.pdf",target="_blank"),"AMIA Annu Symp Proc. 2019;2019:371-379. 
                         Published 2020 Mar 4. (PMID:", a("32308830).",href="https://www.ncbi.nlm.nih.gov/pubmed/32308830", target="_blank"),
                                                     "You can refer to the code",a("here.",href="https://github.com/HimesGroup/brocade", target="_blank")),
                                                   p("The app is created with RStudio's ", a("Shiny", href="https://www.rstudio.com/shiny", target="_blank"),
                                                      "by Maya Shumyatcher, Mengyuan Kan, Avantika Diwadkar and Blanca Himes at the ", a("Himes Lab.", href="https://himeslab.org/", target="_blank"),
                                                     "Full code for REALGAR is available ", a("here.", href="https://github.com/HimesGroup/realgar", target="_blank")))))) # change mk
                                                   #Footer
                                                    # tags$footer(img(src="RStudio-Logo-final.png", height=124*.65, width=110.4*.65), 
                                                    #             "Created with RStudio's ", a("Shiny", href="https://www.rstudio.com/shiny", target="_blank"),
                                                    #             p("This app was created by Maya Shumyatcher, Mengyuan Kan, Avantika Diwadkar and Blanca Himes at the ", a("Himes Lab.", href="https://himeslab.org/", target="_blank")),style = footstyle)))))