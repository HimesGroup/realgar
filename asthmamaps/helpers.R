percent_map <- function(by_county) {
  
  simpleCap2 <- function(x) {
    if(length(x) > 0 & !is.na(x)){
      segments <- strsplit(x, ",")[[1]]
      a <- segments[1]
      b <- segments[2]
      a <- strsplit(a, " ")[[1]]
      b <- strsplit(b, " ")[[1]]
      a <- paste(toupper(substring(a, 1, 1)), substring(a, 2), sep = "", collapse = " ")
      b <- paste(toupper(substring(b, 1, 1)), substring(b, 2), sep = "", collapse = " ")
      c <- paste(b, a, sep = ", ")
      return(c)
    } else return(NA)
  }
  
  mapCounties <- map("county", fill = TRUE,
                     plot = FALSE)
  mapNames <- map("county", plot=FALSE, namesonly =TRUE)
  
  names_matches <- match(mapNames, county.fips$polyname)
  fips_to_match <- county.fips$fips[names_matches]
  matches <- match(fips_to_match, by_county$CntyFIPS)
  vals <- by_county$asthma_percent[matches]
  ns <- format(by_county$count,big.mark=",")[matches]
  
  binpal <- colorBin("Blues", vals, 6, pretty = FALSE)
  
  names_reformatted <- rep(0, 3085)
  for(i in c(1:3085)){
    names_reformatted[i] <- simpleCap2(mapNames[i])
  }
  percents_reformatted <- rep(0, 3085)
  for(i in c(1:3085)){
    percents_reformatted[i] <- ifelse(!is.na(vals[i]), paste0(round(vals[i], 2), "%"), "No data")
  }
  
  county_popup <- paste0("<strong>County: </strong>", 
                         names_reformatted, 
                         "<br><strong>Weighted Asthma Prevalence: </strong>", 
                         percents_reformatted, "<br><strong>N: </strong>", ns)
  
  leaflet() %>% 
    addTiles(options = tileOptions(opacity = 1)) %>%
    addPolygons(data=mapCounties,
                color=~binpal(vals),
                weight = 1,
                fillOpacity = 0.6,
                smoothFactor = 0.5,
                stroke = FALSE,
                popup=county_popup) %>%
    addPolylines(data=mapCounties,
                 color = "black",
                 weight = 0.2) %>%
    addLegend("bottomleft", pal = binpal, 
              title = "Asthma Prevalence",
              values = vals, 
              opacity = 1,
              na.label = "No data",
              labFormat = labelFormat(digits = 1, suffix = "%"))
  
}





