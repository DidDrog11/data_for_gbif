# Load in required packages
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
    "parzer",
    "taxize",
    "xtable"
    )
pacman::p_load(pkgs, character.only = T)
