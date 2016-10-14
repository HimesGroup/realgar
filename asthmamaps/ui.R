library(shiny)
library(leaflet)

county_names <- read.csv("../databases/map_data/county_names.csv")
counties <- as.character(county_names$x)

shinyUI(navbarPage(id="nav",
                   
                   tabPanel("BRFSS Asthma Data 2007-2012","BRFSS Asthma Data 2007-2012",
                            div(class="outer",
                                tags$head(
                                  includeCSS("styles.css"),
                                  includeScript("gomap.js")),
                                
                                leafletOutput("map", width = "100%", height = "100%"),
                                
                                
                                absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                                              draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
                                              width = 500, height = "auto",
                                              
                                              h2("County data explorer"),
                                              
                                              #selectInput("type", "Asthma Defintion:", choices = c("Current", "Lifetime"), width ="100%"),
                                              selectizeInput("county_input", "County:", choices = counties, selected = "Philadelphia, Pennsylvania", 
                                                             multiple=FALSE, width="100%"),
                                              selectInput("var", label = "Year:", choices = c(as.character(c("2007":"2012")),"2007-2012 (all years)"),
                                                          selected = "2007-2012 (all years)", width="100%"),
                                              br(),
                                              h4(textOutput("info"), align="center"),
                                              plotOutput("income_plot", height = 140),
                                              plotOutput("race_plot", height = 140),
                                              plotOutput("bmi_plot", height = 140),
                                              plotOutput("gender_plot", height = 140),
                                              plotOutput("age_plot", height = 140),
                                              plotOutput("smoking_plot", height = 140)
                                )
                            ))
                   
))

