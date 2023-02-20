source(here::here("R", "libraries.R"))

if(!file.exists(here("data_download", "admin_spatial", "level_2_admin.rds"))) {
  
  # Level 2 -------------------------------------------------
  
  BEN_2 <- gadm(country ="BEN", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  BFA_2 <- gadm(country ="BFA", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  CIV_2 <- gadm(country ="CIV", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  GHA_2 <- gadm(country ="GHA", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  GIN_2 <- gadm(country ="GIN", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  GNB_2 <- gadm(country ="GNB", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  LBR_2 <- gadm(country ="LBR", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  MLI_2 <- gadm(country ="MLI", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  MRT_2 <- gadm(country ="MRT", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  NER_2 <- gadm(country ="NER", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  NGA_2 <- gadm(country ="NGA", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  SEN_2 <- gadm(country ="SEN", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  SLE_2 <- gadm(country ="SLE", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  CPV_2 <- gadm(country ="CPV", level = 1, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  GO_2 <- gadm(country ="TGO", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  GMB_2 <- gadm(country ="GMB", level = 2, path = here("data_download", "admin_spatial")) %>%
    st_as_sf()
  
  
  level_two <- ls(pattern = "_2")
  countries_level_two <- do.call("list",mget(level_two))
  write_rds(bind_rows(countries_level_two), here("data_download", "admin_spatial", "level_2_admin.rds"))
  
}

level_2_admin <- read_rds(here("data_download", "admin_spatial", "level_2_admin.rds"))
