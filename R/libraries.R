if (!require("pacman")) install.packages("pacman")

pkgs =
  c("EML",
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
    "taxize",
    "xtable"
    )
pacman::p_load(pkgs, character.only = T)

if(!exists("google_api")) {
google_api <- rstudioapi::askForSecret("Google API Key")
}
