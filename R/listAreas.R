#' List the areas that should be used in the analysis
#'
#' @return a character vector with UTF-8 encoding. 
#' @export
#'
#' @examples

listAreas <- function(){
  
  areas <- c(
     "Åfjord",
     "Ålen og Haltdalen Fjellstyre",
     "Budal Fjellstyre",
     "Dovre Fjellstyre",
     "Eidfjord Fjellstyre",
     "Engerdal Fjellstyre",
     "Fron Fjellstyre",
     "Gausdal Fjellstyre",
     "Helgeland 1",
     "Helgeland 2",
     "Helgeland 3",
     "Indre Finnmark",
     "Indre Salten",
     "Kongsvoll",
     "Kvikne Fjellstyre",
     "Lesja Fjellstyre",
     "Lierne Fjellstyre",
     "Midtre Salten",
     "Namdalseid Fjellstyre",
     "Namskogan Fjellstyre",
     "Njardarheim",
     "Nordre Salten",
     "Os Fjellstyre",
     "Osen Fjellstyre",
     "Øst Finnmark",
     "Øvre Numedal Fjellstyre",
     "Øystre Slidre Fjellstyre",
     "Regionfelt Østlandet",
     "Ringebu Fjellstyre",
     "Røyrvik Fjellstyre",
     "Selbu Fjellstyre",
     "Snåsa Fjellstyre",
     "Soknedal Fjellstyre",
     "Sollia Fjellstyre",
     "Statskog og Klinga utm.",
     "Statskog Røros",
     "Steinkjer Lirype",
     "Stjørdal Fjellstyre",
     "Vest Finnmark kyst",
     "Vestre Slidre Fjellstyre",
     "Ytre Salten"
  )
  
  Encoding(areas) <- "UTF-8"
  return(areas)
}

