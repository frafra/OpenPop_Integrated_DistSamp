
setupMap_NorwayMunic <- function(shp.path, d_trans,
                                 areas, areaAggregation = TRUE){
  
  ## Check whether analysis uses areas (function does not support locality-level plotting)
  if(!areaAggregation){
    stop("Plotting on the map (including this function) does not presently support visualization of locality-level results. Please set areaAggregation = TRUE for the entire workflow.")
  }
  
  ## Set the ecoding (there are Norwegian letters in the data tables)
  Sys.setlocale(locale = 'no_NB.utf8')
  
  ## Reading in map of Norwegian municipalities
  mapNM <- maptools::readShapeSpatial(shp.path, 
                                      proj4string=CRS("+proj=longlat +ellps=WGS84"))
  
  ## Extract list of localities and areas with corresponding municipality names
  areaMunic <- d_trans %>%
    dplyr::filter(verbatimLocality %in% areas) %>%
    dplyr::select(locality, verbatimLocality, municipality) %>%
    dplyr::distinct()
  
  ## Resolve conflicts for plotting (i.e. cases where different areas include the same municipality)
  areaMunic <- areaMunic %>%
    dplyr::mutate(municipality_adj = dplyr::case_when(verbatimLocality == "Soknedal Fjellstyre" ~ NA, 
                                                      verbatimLocality == "Kongsvoll" ~ NA,
                                                      TRUE ~ municipality)) %>%
    dplyr::filter(!(locality %in% c("Harodalen", "Savngovann"))) %>%
    dplyr::select(verbatimLocality, municipality_adj) %>%
    dplyr::distinct() %>%
    dplyr::rename(navn = municipality_adj,
                  Area = verbatimLocality)
  
  ## Add area names to map data
  matchData <- as_tibble(mapNM@data) %>%
    dplyr::left_join(., areaMunic, by = "navn") 

  mapNM@data <- as.data.frame(matchData)
  
  ## Return set up map
  return(mapNM)
}
