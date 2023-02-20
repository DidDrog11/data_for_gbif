source(here::here("R", "libraries.R"))

rodent_data <- read_rds(here("data_raw", "rodent_data.rds")) %>%
  mutate(country = as_factor(country),
         record_id = 1:nrow(.))

study_titles = read_rds(here("data_clean", "studies.rds")) %>%
  select(unique_id, first_author, title, reference_uid) %>%
  mutate(title = paste0("Data from: ", title, " DOI/ISSN: ", reference_uid))

# Cleaning species names --------------------------------------------------

genus_synonym <-  as.list(c("praomys", "mus"))
names(genus_synonym) <- c("myomys", "nannomys")

# The accepted term to look up on GBIF
accepted_names <- as.list(c("mus musculus", "oenomys hypoxanthus", "suncus megalura", "dasymys incomtus", "atelerix albiventris", "gerbilliscus gambiana",
                            "gerbilliscus guinea", "gerbilliscus guinea", "gerbilliscus kempii", "dipodillus campestris", "gerbilliscus gambiana",
                            "grammomys poensis", "graphiurus nagtglasii", "graphiurus kelleni", "paraechinus aethiopicus", "lophuromys flavopunctatus",
                            "massoutiera mzabi", "mastomys natalensis", "steatomys caurinus", "gerbilliscus kempii", "gerbilliscus guineae",
                            "gerbilliscus guineae", "gerbilliscus kempii", "gerbilliscus robusta", "taterillus gracilis", "grammomys poensis",
                            "crocidura olivieri", "paraxerus poensis", "protoxerus aubinnii", "epixerus ebii", "funisciurus pyrropus", "gerbilliscus guineae"))
# The term in the source
names(accepted_names) <- c("mus domesticus", "aethomys hypoxanthus", "crocidura megalura", "dasymys bentleyae", "erinaceus albiventris", "gerbilliscus gambianus",
                           "tatera guineae", "gerbilliscus guineae", "gerbilliscus kempi", "gerbillus campestris", "gerbillus gambianus", "grammomys rutilans",
                           "graphiurus hueti", "graphiurus parvus", "hemiechinus aethiopicus", "lophuromys flavipunctatus", "massouteria mzabi",
                           "mastomys hildebrandtii", "steatomys caurianus", "tatera gambiana", "tatera guinea", "tatera guineae", "tatera kempi",
                           "tatera robusta", "taterillus gracillis", "thamnomys rutilans", "crocidura occidentalis", "aethosciurus poensis", "allosciurus aubinii",
                           "epixerus jonesii", "funisciurus leonis", "gerbilliscus guinea")

write_rds(genus_synonym, here("data_clean", "genus_dictionary.rds"))
write_rds(accepted_names, here("data_clean", "species_dictionary.rds"))

as.data.frame.gbifid <- function(x, ...){
  data.frame(ids = as.character(unclass(x)),
             class = "gbifid",
             match = attr(x, "match"),
             multiple_matches = attr(x, "multiple_matches"),
             pattern_match = attr(x, "pattern_match"),
             uri = attr(x, "uri"),
             stringsAsFactors = FALSE)
}

# Cleaning species names and matching to GBIF
rodent_data %<>%
  mutate(genus = recode(genus, !!!genus_synonym),
         classification = paste(genus, ifelse(species == "-", "sp.", species), sep = " "),
         classification = recode(classification, !!!accepted_names)) %>%
  separate(col = classification, into = c("genus", "species"), sep = " ", remove = F) %>%
  mutate(species = ifelse(species %in% c("sp.", "sp.1", "sp.2"), "", species))

# Link the genus of trapped rodents to a gbif id
if(!file.exists(here("data_clean", "genus_gbif.rds"))){
  genus <- tibble(rodent_data %>%
                    distinct(genus)) %>%
    arrange(genus)
  genus$gbif_id <- get_gbifid(snakecase::to_sentence_case(genus$genus), ask = T)
  write_rds(genus, here("data_clean", "genus_gbif.rds"))
}

# Pull the hierarchy for each genus from gbif
if(!file.exists(here("data_clean", "genus_hierarchy.rds"))){
  genera <- classification(snakecase::to_sentence_case(genus$genus), db = "gbif")
  genera <- as_tibble(do.call(rbind,c(genera, make.row.names = T)), rownames = "genera") %>%
    mutate(genera = str_sub(genera, 1, -3)) %>%
    pivot_wider(id_cols = genera, names_from = rank, values_from = name) %>%
    dplyr::select(-genera) %>%
    write_rds(here("data_clean", "genus_hierarchy.rds"))
}

# Link the species of trapped rodents to a gbif id
if(!file.exists(here("data_clean", "species_gbif.rds"))){
  species <- tibble(rodent_data %>%
                      filter(species != "-") %>%
                      distinct(classification)) %>%
    arrange(classification)
  species$gbif_id <- get_gbifid(species$classification, ask = T)
  write_rds(species, here("data_clean", "species_gbif.rds"))
}

genus <- read_rds(here("data_clean", "genus_gbif.rds"))
genus_hierarchy <- read_rds(here("data_clean", "genus_hierarchy.rds")) %>%
  mutate(across(.cols = everything(), .fns = str_to_lower))
species <- read_rds(here("data_clean", "species_gbif.rds")) %>%
  drop_na(gbif_id)

# Add GBIF ID for species and genus identification to the raw data
rodent_classifications <- rodent_data %>%
  full_join(., genus %>%
              rename("genus_gbif" = gbif_id), by = "genus") %>%
  full_join(., genus_hierarchy, by = "genus") %>%
  full_join(., species %>%
              rename("species_gbif" = gbif_id), by = "classification") %>%
  mutate(gbif_id = ifelse(is.na(species_gbif), genus_gbif, species_gbif),
         iso3c = countrycode(as.character(country), "country.name", "iso3c")) %>%
  distinct() %>%
  drop_na(number)


# Import GBIF template ----------------------------------------------------
gbif_template <- read_xlsx(path = here("gbif_template", "occurrence_ipt_template_v2_example_data.xlsx"))
gbif_col_names <- colnames(gbif_template)[!colnames(gbif_template) %in% colnames(rodent_classifications)]

# Handling coordinate formats ---------------------------------------------
# Cleaning coordinates
rodent_coords <- rodent_data %>%
  separate(col = longitude_DMS_W, into = c("long_degrees", "long_minutes", "long_seconds"), "_", remove = F) %>%
  separate(col = latitude_DMS_N, into = c("lat_degrees", "lat_minutes", "lat_seconds"), "_") %>%
  mutate(across(all_of(c("long_degrees", "long_minutes", "long_seconds","lat_degrees", "lat_minutes", "lat_seconds")), as.double),
         long_hemi = ifelse(long_degrees<0, "E", #assign -ve numbers to E
                            ifelse(substring(longitude_DMS_W, 1, 1) == "-", "E", "W")), #as 0 cannot be -ve we can check the sign on the text entry
         lat_hemi = ifelse(lat_degrees<0, "S", "N"),
         gps_dms = ifelse(is.na(lat_degrees|long_degrees), F, T),
         long_dms = ifelse(gps_dms == T,
                           paste(long_hemi, long_degrees, " ",
                                 ifelse(is.na(long_minutes), 0, long_minutes), '.',
                                 ifelse(is.na(long_seconds), 0, long_seconds),
                                 sep = ""), NA),
         lat_dms = ifelse(gps_dms == T,
                          paste(lat_hemi, lat_degrees, " ",
                                ifelse(is.na(lat_minutes), 0, lat_minutes), '.',
                                ifelse(is.na(lat_seconds), 0, lat_seconds),
                                sep = ""), NA),
         long_dms = gsub("-", "", long_dms),
         lat_dms = gsub("-", "", lat_dms),
         iso3c = countrycode(as.character(country), "country.name", "iso3c")) %>%
  dplyr::select(-longitude_DMS_W)

# Converting coordinate types into consistent decimal degrees
dms <- rodent_coords %>%
  drop_na(long_dms, lat_dms)

dms <- dms %>%
  mutate(lon_dd = parzer::parse_lon(long_dms),
         lat_dd = parzer::parse_lat(lat_dms)) %>%
  st_as_sf(coords = c("lon_dd", "lat_dd"), crs = "+proj=longlat +datum=WGS84")

dd <- rodent_coords %>%
  drop_na(longitude_D_E,latitude_D_N) %>%
  mutate(lon_dd = longitude_D_E,
         lat_dd = latitude_D_N) %>%
  st_as_sf(coords = c("lon_dd", "lat_dd"), crs = "+proj=longlat +datum=WGS84")

utm <- rodent_coords %>%
  drop_na(UTM_coordinates) %>%
  separate(col = UTM_coordinates, into = c("zone", "easting", "northing"), "_")

utm_q <- utm %>%
  filter(zone == "28Q") %>%
  st_as_sf(coords = c("easting", "northing")) %>%
  st_set_crs(value = "EPSG:2161") %>%
  st_transform(crs = "+proj=longlat +datum=WGS84")

utm_p <- utm %>%
  filter(zone == "28P") %>%
  st_as_sf(coords = c("easting", "northing"), crs = "EPSG:2161") %>%
  st_transform(crs = "+proj=longlat +datum=WGS84")

rodent_gps <- bind_rows(dms, dd) %>%
  bind_rows(utm_q, utm_p)

st_crs(rodent_gps) = 4326

no_gps <- rodent_coords %>%
  filter(!record_id %in% rodent_gps$record_id) %>%
  mutate(long_dms = as.double(long_dms),
         lat_dms = as.double(lat_dms))


# Create coords for lowest level of resolution ----------------------------
admin <- read_rds(here("data_download", "admin_spatial", "level_2_admin.rds")) %>%
  filter(COUNTRY %in% unique(no_gps$country))

region_conversion <- c("Mount Nimba" = "Lola", "Coastal Region" = "Boké|Kindia|Conakry", "Savannah Region" = "Labé|Mamou",
                       "Savannah/Forest Transition" = "Faranah|Kankan", "Forest Region" = "Nzérékoré")

no_gps_list <- no_gps %>%
  mutate(location = recode(region, !!!region_conversion),
         location = coalesce(location, country),
         accuracy = case_when(is.na(region) ~ "national",
                              TRUE ~ "regional")) %>% 
  group_by(country) %>%
  group_split()

no_gps_df <- lapply(no_gps_list, function(x) {
  
  country = unique(x$country)
  locations = unique(x$location)
  coords <- list()
  
  for(i in 1:length(locations)) {
    
    if(unique(x$accuracy) == "regional") {
      
      coords[[i]] <- admin %>%
        filter(COUNTRY == country) %>%
        filter(str_detect(NAME_1, paste(str_c(locations[i], "|"), collapse = ""))) %>%
        summarise() %>%
        st_centroid() %>%
        mutate(location = locations[i])
    
    } else {
    
      coords[[i]] <- admin %>%
        filter(COUNTRY == country) %>%
        st_make_valid() %>%
        summarise() %>%
        st_centroid() %>%
        mutate(location = locations[i])
      
  }}
  
  coords <- bind_rows(coords)
  
  imputed_coords <- left_join(x %>%
                                select(unique_id, record_id, location, accuracy), coords, by = "location")
  }) %>%
  bind_rows()  %>%
  mutate(occurrenceID = paste0(unique_id, "_", record_id)) %>%
  select(-unique_id, -record_id)
  

all_rodent_coords <- rodent_gps %>%
  mutate(occurrenceID = paste0(unique_id, "_", record_id),
         accuracy = case_when(!is.na(town_village) ~ "site",
                              TRUE ~ region)) %>%
  select(-unique_id, -record_id) %>%
  bind_rows(no_gps_df) %>%
  select(occurrenceID, accuracy, geometry)


# Build the GBIF dataframe ------------------------------------------------


final_columns <- c("occurrenceID", "basisOfRecord", "scientificName", "eventDate",
                   "countryCode", "taxonRank", "kingdom", "phylum", "order", "family", "genus",
                   "taxonRank", "decimalLatitude", "decimalLongitude",
                   "geodeticDatum", "coordinateUncertaintyInMeters", "individualCount", "dataGeneralizations",
                   "country", "locality", "verbatimLocality", "identifiedBy", "datasetName", "eventRemarks")

dataGeneralisations <- c("region" = "Centre of smallest associated administrative district used",
                         "regional" = "Centre of smallest associated administrative district used",
                         "Niakhar" = "Centre of smallest associated administrative district used",
                         "country" = "Centroid of the country used as no finer scale geographic information available",
                         "site" = "Coordinates reflect the location of the study site where the sample was collected")

prep_df <- rodent_classifications %>%
  bind_cols(as_tibble(matrix(nrow = nrow(rodent_classifications), ncol = length(gbif_col_names)), .name_repair = ~ gbif_col_names)) %>%
  relocate(colnames(gbif_template)) %>%
  mutate(occurrenceID = paste0(unique_id, "_", record_id),  # The occurrenceID submitted to GBIF will be a combination of the study identifier and record identifier
         basisOfRecord = "LivingSpecimen", # All rodents were trapped in either rodent traps or brought to the researchers by local hunters
         month_start = str_split(month_trapping, "-", simplify = TRUE)[ , 1],
         month_start = case_when(nchar(month_start) == 3 ~ match(month_start, month.abb),
                                 nchar(month_start) > 3 ~ match(month_start, month.name),
                                 month_start == "Nove" ~ match(str_sub(month_start, 1, 3), month.abb)), # Convert to numeric months for the starting date of collections
         month_end = str_split(month_trapping, "-", simplify = TRUE)[ , 2],
         month_end = case_when(nchar(month_end) == 3 ~ match(month_end, month.abb),
                                 nchar(month_end) > 3 ~ match(month_end, month.name),
                                 month_end == "Nove" ~ match(str_sub(month_end, 1, 3), month.abb)), # Convert to numeric months for the ending date of collections
         year_start = str_split(year_trapping, "-", simplify = TRUE)[, 1],
         year_start = case_when(is.na(year_trapping) ~ str_split(unique_id, "_", simplify = TRUE)[, 2],
                                TRUE ~ year_start), # For studies not specifying the year of trapping the year of publication is used
         year_end = str_split(year_trapping, "-", simplify = TRUE)[ , 2],
         year_end = case_when(year_end == "" ~ as.character(NA),
                              TRUE ~ year_end),
         year_end = coalesce(year_end, year_start), # Set the year end to equal the start if only one year is noted
         eventDate = case_when(is.na(year_trapping) ~ year_start,
                               year_start == year_end & is.na(month_start) ~ year_start,
                               year_start != year_end & is.na(month_start) ~ paste0(year_start, "/", year_end),
                               year_start == year_end & is.na(month_end) ~ paste0(year_start, "-", month_start),
                               year_start == year_end & month_start == month_end ~ paste0(year_start, "-", month_start),
                               year_start == year_end & month_start != month_end ~ paste0(year_start, "-", month_start, "/", year_end, "-", month_end),
                               year_start != year_end ~ paste0(year_start, "-", month_start, "/", year_end, "-", month_end))) %>%
  select(-endDayOfYear, -year, -month, -day, -verbatimEventDate) %>%
  mutate(eventRemarks = case_when(is.na(year_trapping) ~ "No dates supplied in publication, year of publication used",
                                  TRUE ~ as.character(NA))) %>%
  left_join(., species %>%
              rename(gbif_species = classification) %>%
              mutate(gbif_id = as.character(gbif_id)), by = "gbif_id") %>% # Merge GBIF IDs with accepted scientific name
  mutate(scientificName = case_when(!is.na(gbif_species) ~ gbif_species,
                                           !is.na(genus) ~ genus,
                                           TRUE ~ as.character(NA)),  # Use scientific name if available, otherwise genus association
         countryCode = countrycode(country, "country.name", "iso2c")) %>%
  select(-higherClassification, -specificEpithet, -infraspecificEpithet) %>%
  mutate(taxonRank = case_when(!is.na(gbif_species) ~ "species",
                               gbif_id == "1459" ~ "order",
                               TRUE ~ "genus")) %>%
  select(-dateIdentified,	-nomenclaturalCode, -year_trapping, -month_trapping,
         -region, -latitude_DMS_N, -longitude_DMS_W, -latitude_D_N, -longitude_D_E, -UTM_coordinates) %>%
  left_join(., all_rodent_coords %>%
              mutate(lat = st_coordinates(geometry)[ ,2],
                     lon = st_coordinates(geometry)[ ,1]), by = "occurrenceID") %>%
  mutate(decimalLatitude = lat,
         decimalLongitude = lon,
         geodeticDatum = "EPSG:4326",
         dataGeneralizations = recode(accuracy, !!!dataGeneralisations),
         coordinateUncertaintyInMeters = case_when(accuracy == "site" ~ 1000,
                                                   str_detect(accuracy, "regional|Niakhar") ~ 100000,
                                                   accuracy == "national" ~ 500000),
         genus = case_when(taxonRank == "order" ~ as.character(NA),
                           TRUE ~ genus)) %>%
  select(-any_of(matches("verbatimCoo|georef|higher|continent|island|state|county"))) %>%
  mutate(locality = town_village,
         verbatimLocality = habitat,
         individualCount = number) %>%
  left_join(., study_titles, by = "unique_id") %>%
  mutate(identifiedBy = first_author,
         datasetName = title,
         order = case_when(scientificName == "rodentia" ~ "rodentia",
                           TRUE ~ order),
         phylum = "chordata",
         kingdom = "animalia")

final_df <- prep_df %>%
  select(any_of(final_columns))

dir.create(here("final_data"))
write_rds(final_df, here("data_clean", "occurrence_df.rds"))
write_tsv(final_df, here("final_data", "taxa.txt"))
