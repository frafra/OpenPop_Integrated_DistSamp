

################################################################################
################################################################################

## Downloading rodent data - from all areas; encoding = UTF-8"

Sys.setlocale(locale = "no_NB.utf8")

library(tidyverse)
library(lubridate)
library(RODBC)

# myconn <- "THIS IS SET IN FILE 0_Rodent_credentials.R"

Kommune <- as_tibble(sqlFetch(myconn, "FYLKEKOMMUNE"))

Taks <- as_tibble(sqlQuery(myconn, paste("SELECT TakseringUUID,
                               TakseringID,
                               LengdeTaksert,
                               FK_LinjeID, 
                               Aar,
                               Dato,
                               SettSmagnager, 
                               SettMus,
                               SettLemen,
                               SettUkjentSmagnager
                               FROM Taksering")))

TaksLin <- as_tibble(sqlQuery(myconn, paste("SELECT LInjeUUID, LinjeID, FK_OmradeID, STAsText FROM Takseringslinje")))
TaksOmr <- as_tibble(sqlQuery(myconn, paste("SELECT OmradeID, FK_Fylkekomnr, OmradeNavn, FK_RapporteringsnivaaID FROM Takseringsomrade")))

RappOmr <- as_tibble(sqlQuery(myconn, paste("SELECT ID, Navn, DataPublishGBIF FROM Rapporteringsnivaa")))


close(myconn)

Kommune <- Kommune %>% mutate(Kommunenavn = str_trim(Kommunenavn), Fylkesnavn = str_trim(Fylkesnavn))


############################################################################################

Rodent <- Taks %>%
  left_join(., TaksLin, by = c("FK_LinjeID" = "LinjeID")) %>%
  left_join(., TaksOmr, by = c("FK_OmradeID" = "OmradeID")) %>%
  left_join(., Kommune, by = c("FK_Fylkekomnr" = "Fylkekomnr")) %>%
  left_join(., RappOmr, by = c("FK_RapporteringsnivaaID" = "ID")) %>%
  filter(DataPublishGBIF == 1) %>%
  filter(!is.na(SettSmagnager)) %>%
  rename(
    EventID = TakseringUUID,
    date = Dato,
    year = Aar,
    LineID = FK_LinjeID,
    locationID = LInjeUUID,
    municipality = Kommunenavn,
    locality = OmradeNavn,
    stateProvince = Fylkesnavn,
    verbatimLocality = Navn,
    sampleSizeValue = LengdeTaksert,
    rodentOcc = SettSmagnager
  ) %>%
  mutate(
    eventRemarks = paste0("LineTransect"),
    sampleSizeUnit = paste0("metre")
  ) %>%
  select(
    EventID,
    eventRemarks,
    date,
    year,
    locationID,
    stateProvince,
    municipality,
    locality,
    verbatimLocality,
    sampleSizeValue,
    sampleSizeUnit,
    rodentOcc
  )

#### Writing to file(s)

write.table(Rodent, "Rodent_Data/Rodent_data.csv", dec = ",", sep = ";", row.names = FALSE, fileEncoding = "windows-1252")
saveRDS(Rodent, "Rodent_Data/Rodent_data.rds")

# test <- read.csv("Rodent_Data/Rodent_data.csv",  dec=",", sep=";", fileEncoding = "windows-1252")
