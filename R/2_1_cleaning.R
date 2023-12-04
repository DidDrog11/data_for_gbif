source(here::here("R", "libraries.R"))

# Studies -----------------------------------------------------------------
# Convert the unique ID to a factor
studies <- read_rds(here("data_raw", "studies.rds")) %>%
  mutate(unique_id = as_factor(unique_id))

write_rds(studies, here("data_clean", "studies.rds"))