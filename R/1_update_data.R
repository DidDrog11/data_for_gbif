source(here::here("R", "libraries.R"))

# Download the data from the Google Drive
drive_download("https://docs.google.com/spreadsheets/d/1rQYjHhk6uk1PoKZZVsgFlmuqWGicU2tTisk9ddfAwTM/edit#gid=0", path = here("data_download", "included_studies.xlsx"), overwrite = T)

# Save the studies sheet as an RDS
studies <- read_xlsx(here("data_download", "included_studies.xlsx"), sheet = "study", col_types = "text")

write_rds(studies, here("data_raw", "studies.rds"))

# Save the rodent data as an RDS
rodent_data <- read_xlsx(here("data_download", "included_studies.xlsx"), sheet = "trapping", col_types = c("guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess",
                                                                                                           "guess", "guess", "text", "text", "guess", "guess"))
write_rds(rodent_data, here("data_raw", "rodent_data.rds"))

# Save the pathogen data as an RDS
pathogen_data <- read_xlsx(here("data_download", "included_studies.xlsx"), sheet = "pathogen", col_types = c("guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess",
                                                                                                             "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess",
                                                                                                             "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess",
                                                                                                             "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                                                                             "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
write_rds(pathogen_data, here("data_raw", "pathogen.rds"))
