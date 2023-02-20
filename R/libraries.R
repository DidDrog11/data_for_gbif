if (!require("pacman")) install.packages("pacman")

pkgs =
  c("eml",
    "here",
    "tidyverse",
    "bib2df",
    "googledrive",
    "geodata",
    "readxl",
    "countrycode",
    "sf",
    "rosm",
    "lubridate",
    "taxize"
    )
pacman::p_load(pkgs, character.only = T)

# if(!exists("google_api")) {
# google_api <- rstudioapi::askForSecret("Google API Key")
# }

if(!exists("ENTREZ_KEY")) {
  ENTREZ_KEY <- rstudioapi::askForSecret("Entrez API")
}

#source(here::here("R", "functions.R"))