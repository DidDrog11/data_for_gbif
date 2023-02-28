# This is the main script to replicate the analysis in the manuscript

# Load libraries and read in functions

source(here::here("R", "libraries.R"))

# The following scripts ensure local data is up to date

## create download folders
dir.create("data_download")
dir.create(here("data_download", "admin_spatial"))

## spatial data is obtained from GADM and is saved in the data_download/admin_spatial folder

source(here("R", "1_download_spatial.R"))

dir.create("data_raw")
source(here("R", "1_update_data.R"))

# The following scripts clean and reshape the data

dir.create("data_clean")
source(here("R", "2_1_cleaning.R"))

source(here("R", "2_2_clean_species.R"))

# Produce the metadata file
source(here("R", "3_create_metadata.R"))

source(here("R", "2_4_clean_pathogen.R"))
