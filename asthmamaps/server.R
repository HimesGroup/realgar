library(shiny)

library(ggplot2)

library(maps)
library(mapproj)
library(dplyr)
data("county.fips")
source("helpers.R")
library(reshape2)

# simpleCap2 <- function(x) {
#   if(length(x) > 0 & !is.na(x)){
#     segments <- strsplit(x, ",")[[1]]
#     a <- segments[1]
#     b <- segments[2]
#     a <- strsplit(a, " ")[[1]]
#     b <- strsplit(b, " ")[[1]]
#     a <- paste(toupper(substring(a, 1, 1)), substring(a, 2), sep = "", collapse = " ")
#     b <- paste(toupper(substring(b, 1, 1)), substring(b, 2), sep = "", collapse = " ")
#     c <- paste(b, a, sep = ", ")
#     return(c)
#   } else return(NA)
# }

##getting names of counties as options
# weighted_asthma_prev <- read.csv("data/weighted_current_asthma_w_counts_leaflet.csv")
# counties <- filter(weighted_asthma_prev, YEAR == "all")
# matches <- match(counties$CntyFIPS, county.fips$fips)
# polynames <- as.character(county.fips[matches,]$polyname)
# polynames_ed <- sapply(polynames, function(i) simpleCap2(i))
# polynames_final <- polynames_ed[!is.na(polynames_ed)]
# write.csv(polynames_final, "data/county_names.csv")

county_names <- read.csv("../databases/map_data/county_names.csv")

weighted_current_asthma_prev <- read.csv("../databases/map_data/weighted_current_asthma_w_counts_leaflet.csv")
weighted_current_asthma_prev$asthma_percent <- weighted_current_asthma_prev$asthnow*100

weighted_current_vars <- read.csv("../databases/map_data/weighted_current_variables_leaflet.csv")
colnames(weighted_current_vars) <- c("id", "CntyFIPS", "YEAR", "ASTHMA", "<$25,000", "$25,000-$75,000", ">$75,000", 
                                     "Male", "Female", "White", "Asian/Pacific Islander",
                                     "Black", "Hispanic", "American Indian/Alaskan Native",
                                     "Not overweight or obese", "Overweight", "Grade 1 Obese", "Grades 2 & 3 Obese",
                                     "Never smoked", "Former smoker", "Current smoker",
                                     "Less than high school", "High school", "Some college or more",
                                     "22-34", "35-44", "45-54", "55-64", "65+")

current.asthma.2007 = filter(weighted_current_asthma_prev, YEAR == 2007)
current.asthma.2008 = filter(weighted_current_asthma_prev, YEAR == 2008)
current.asthma.2009 = filter(weighted_current_asthma_prev, YEAR == 2009)
current.asthma.2010 = filter(weighted_current_asthma_prev, YEAR == 2010)
current.asthma.2011 = filter(weighted_current_asthma_prev, YEAR == 2011)
current.asthma.2012 = filter(weighted_current_asthma_prev, YEAR == 2012)
current.asthma.all = filter(weighted_current_asthma_prev, YEAR == "all")
asthma_prev_data <- list(current.asthma.2007, current.asthma.2008, current.asthma.2009,
                         current.asthma.2010, current.asthma.2011, current.asthma.2012, current.asthma.all)

current.vars.2007 = filter(weighted_current_vars, YEAR == 2007)
current.vars.2008 = filter(weighted_current_vars, YEAR == 2008)
current.vars.2009 = filter(weighted_current_vars, YEAR == 2009)
current.vars.2010 = filter(weighted_current_vars, YEAR == 2010)
current.vars.2011 = filter(weighted_current_vars, YEAR == 2011)
current.vars.2012 = filter(weighted_current_vars, YEAR == 2012)
current.vars.all = filter(weighted_current_vars, YEAR == "all")
vars_data <- list(current.vars.2007, current.vars.2008, current.vars.2009,
                  current.vars.2010, current.vars.2011, current.vars.2012, current.vars.all)

simpleCap <- function(x) {
  if(length(x) > 0 & !is.na(x)){
    segments <- strsplit(x, ",")[[1]]
    a <- segments[1]
    b <- segments[2]
    a <- strsplit(a, " ")[[1]]
    b <- strsplit(b, " ")[[1]]
    a <- paste(toupper(substring(a, 1, 1)), substring(a, 2), sep = "", collapse = " ")
    b <- paste(toupper(substring(b, 1, 1)), substring(b, 2), sep = "", collapse = " ")
    c <- paste0(paste(b, a, sep = ", "), ": ")
    return(c)
  } else return(paste("or click on a county to view its data. "))
}

shinyServer(
  function(input, output, session) {
    
    brfss_year <- reactive({
      switch(input$var,
             "2007" = asthma_prev_data[[1]],
             "2008" = asthma_prev_data[[2]],
             "2009" = asthma_prev_data[[3]],
             "2010" = asthma_prev_data[[4]],
             "2011" = asthma_prev_data[[5]],
             "2012" = asthma_prev_data[[6]],
             "2007-2012 (all years)" = asthma_prev_data[[7]]
      )
    })
    
    full_dat <- reactive({
      switch(input$var,
             "2007" = vars_data[[1]],
             "2008" = vars_data[[2]],
             "2009" = vars_data[[3]],
             "2010" = vars_data[[4]],
             "2011" = vars_data[[5]],
             "2012" = vars_data[[6]],
             "2007-2012 (all years)" = vars_data[[7]]
      )
    })
    
    output$map <- renderLeaflet({
      percent_map(brfss_year())
    })
    
    county.click <- reactive ({
      as.character(county_names[match(input$county_input, 
                                      county_names$x), "X"]) })
    
    fips.click <- reactive ({ county.fips[match(county.click(), county.fips$polyname),]$fips })
    
    asthma_percent_raw <- reactive ({ if( !is.null(fips.click()) & !is.na(fips.click()) ) {
      return(brfss_year()[which(brfss_year()$CntyFIPS == fips.click()), "asthma_percent"]) } })
    
    asthma_percent <- reactive ({ if(length(asthma_percent_raw() == 1)) return(paste0(as.character(round(asthma_percent_raw(), 2)), "%"))
                                  else return ("") })
    
    sample_size <- reactive({ ifelse(!is.null(fips.click()) & !is.na(fips.click()), paste0(" (N = ", formatC(brfss_year()[match(fips.click(), brfss_year()$CntyFIPS), "count"], format = "d", big.mark=","), ")"), "") })
    
    output$info <- renderText({ if(length(county.click()) > 0 & nchar(asthma_percent()) > 0) { 
      paste0(
        "Weighted Asthma Prevalence: ", asthma_percent(), sample_size()) }
      else if (!is.na(county.click())) { paste0(simpleCap( county.click() ), " No data for this county/year") }
      else paste0("")
    }) 
    
    county_yearly_dat <- reactive ({ if(!is.null(fips.click())) return(full_dat()[which(full_dat()$CntyFIPS == fips.click()),]) else return(data.frame(0)) })
    
    bmi <- reactive({ melt(select(county_yearly_dat(), 4, 15, 16, 17, 18), id.vars=c("ASTHMA")) })
    race <- reactive({ melt(select(county_yearly_dat(), 4, 10, 11, 12, 13, 14), id.vars=c("ASTHMA")) })
    income <- reactive({ melt(select(county_yearly_dat(), 4, 5, 6, 7), id.vars=c("ASTHMA")) })
    smoking <- reactive({ melt(select(county_yearly_dat(), 4, 19, 20, 21), id.vars=c("ASTHMA")) })
    age_cat <- reactive({ melt(select(county_yearly_dat(), 4, 25:29), id.vars=c("ASTHMA")) })
    gender <- reactive({ melt(select(county_yearly_dat(), 4, 8, 9), id.vars=c("ASTHMA")) })
    
    output$bmi_plot <- renderPlot({ ggplot(bmi(), aes(x = factor(ASTHMA), y = value, fill = factor(variable))) + 
                                      geom_bar(stat="identity", position="fill") +
                                      scale_x_discrete(labels=c("No Asthma", "Asthma")) +
                                      theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) +
                                      guides(fill = guide_legend(title = "Body Mass Index (BMI)", title.position = "top", title.hjust = 0.5, nrow=2)) +
                                      theme(axis.line = element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank()) +
                                      coord_flip() })
    
    output$race_plot <- renderPlot({ ggplot(race(), aes(x = factor(ASTHMA), y = value, fill = factor(variable))) + 
                                       geom_bar(stat="identity", position="fill") +
                                       scale_x_discrete(labels=c("No Asthma", "Asthma")) +
                                       theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) +
                                       guides(fill = guide_legend(title = "Race/Ethnicity", title.position = "top", title.hjust = 0.5, nrow=3)) +
                                       theme(axis.line = element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank()) +
                                       coord_flip() })
    
    output$income_plot <- renderPlot({ ggplot(income(), aes(x = factor(ASTHMA), y = value, fill = factor(variable))) + 
                                         geom_bar(stat="identity", position="fill") +
                                         scale_x_discrete(labels=c("No Asthma", "Asthma")) +
                                         theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) +
                                         guides(fill = guide_legend(title = "Income", title.position = "top", title.hjust = 0.5)) +
                                         theme(axis.line = element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank()) +
                                         coord_flip() })
    
    output$smoking_plot <- renderPlot({ ggplot(smoking(), aes(x = factor(ASTHMA), y = value, fill = factor(variable))) + 
                                          geom_bar(stat="identity", position="fill") +
                                          scale_x_discrete(labels=c("No Asthma", "Asthma")) +
                                          theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) +
                                          guides(fill = guide_legend(title = "Smoking Status", title.position = "top", title.hjust = 0.5)) +
                                          theme(axis.line = element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank()) +
                                          coord_flip() })
    
    output$age_plot <- renderPlot({ ggplot(age_cat(), aes(x = factor(ASTHMA), y = value, fill = factor(variable))) + 
                                      geom_bar(stat="identity", position="fill") +
                                      scale_x_discrete(labels=c("No Asthma", "Asthma")) +
                                      theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) +
                                      guides(fill = guide_legend(title = "Age", title.position = "top", title.hjust = 0.5)) +
                                      theme(axis.line = element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank()) +
                                      coord_flip() })
    
    output$gender_plot <- renderPlot({ ggplot(gender(), aes(x = factor(ASTHMA), y = value, fill = factor(variable))) + 
                                         geom_bar(stat="identity", position="fill") +
                                         scale_x_discrete(labels=c("No Asthma", "Asthma")) +
                                         theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) +
                                         guides(fill = guide_legend(title = "Gender", title.position = "top", title.hjust = 0.5)) +
                                         theme(axis.line = element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank()) +
                                         coord_flip() })
    
    
  })










