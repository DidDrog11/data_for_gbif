source(here::here("R", "libraries.R"))

rodent_data <- read_rds(here("data_raw", "rodent_data.rds")) %>%
  mutate(country = as_factor(country),
         record_id = 1:nrow(.))

study_titles = read_rds(here("data_clean", "studies.rds")) %>%
  select(unique_id, first_author, title, reference_uid, rightsHolder, license)


# Clean citations ---------------------------------------------------------

citations <- bib2df(here("citations", "include_final.bib"))

clean_citations <- citations %>%
  distinct(TITLE, .keep_all = TRUE) %>%
  select(BIBTEXKEY, AUTHOR, TITLE, JOURNAL, YEAR, VOLUME, NUMBER, PAGES, DOI) %>%
  group_by(BIBTEXKEY) %>%
  unnest(AUTHOR) %>%
  summarise(authors = paste(AUTHOR, collapse = ", "),
            title = unique(TITLE),
            journal = unique(JOURNAL),
            year = unique(YEAR),
            volume = unique(VOLUME),
            number = unique(NUMBER),
            pages = unique(PAGES),
            doi = unique(DOI)
            ) %>%
  mutate(authors = str_remove_all(authors, "\\{|\\}")) %>%
  left_join(study_titles %>%
              select(unique_id, title),
            by = "title") %>%
  left_join(study_titles %>%
              select(unique_id, reference_uid),
            by = c("doi" = "reference_uid")) %>%
  mutate(unique_id = coalesce(unique_id.x, unique_id.y)) %>%
  select(-unique_id.x, -unique_id.y) %>%
  mutate(full_citation = paste0(authors, ". ", year, ", ", title, ". ", 
                                if_else(is.na(journal), "", journal),
                                ". ", 
                                if_else(is.na(volume), "", volume),
                                if_else(is.na(number), "", paste0("(", number, ")")),
                                if_else(is.na(pages), "", paste0(" p.", pages, ". ")),
                                if_else(is.na(doi), "", paste0("DOI/ISSN/ISBN: ", doi))))

# Manually link the articles that cannot be done based on titles or doi
clean_citations$unique_id[clean_citations$BIBTEXKEY == "barnett_ecology_2000"] <- "ab_2000_sierraleone"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "decher_rapid_2004"] <- "jd_2004_guinea"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "ebenezer_effects_2012"] <- "ae_2012_nigeria"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "granjon_small_2002"] <- "lg_2002_mauritania"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "fichet-calvet_diversity_2009"] <- "fc_2009_guinea"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "diallo_bacteriological_1982"] <- "ad_1982_nigeria"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "monadjem_rapid_2005"] <- "am_2005_liberia"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "attuquayefio_study_2003"] <- "da_2003_ghana"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "meinig_notes_1999"] <- "hm_1999_mali"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "mccullough_rapid_2005"] <- "jd_2005_ghana"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "kwaku_rapid_2006"] <- "pb_2006_ghana"
clean_citations$unique_id[clean_citations$BIBTEXKEY == "wright_rapid_2006"] <- "rn_2006_guinea"

final_citations <- clean_citations %>%
  select(unique_id, authors, full_citation) %>%
  drop_na(unique_id) %>%
  group_by(unique_id) %>%
  arrange(-str_length(full_citation)) %>%
  slice(1) %>%
  ungroup()


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
  mutate(across(.cols = everything(), .fns = str_to_sentence))
species <- read_rds(here("data_clean", "species_gbif.rds")) %>%
  drop_na(gbif_id)

# Add GBIF ID for species and genus identification to the raw data
rodent_classifications <- rodent_data %>%
  mutate(genus = str_to_sentence(genus)) %>%
  full_join(., genus %>%
              rename("genus_gbif" = gbif_id), by = "genus") %>%
  full_join(., genus_hierarchy, by = "genus") %>%
  full_join(., species %>%
              rename("species_gbif" = gbif_id), by = "classification") %>%
  mutate(gbif_id = ifelse(is.na(species_gbif), genus_gbif, species_gbif),
         iso3c = countrycode(as.character(country), "country.name", "iso3c")) %>%
  distinct() %>%
  drop_na(number) %>%
  group_by(unique_id) %>%
  mutate(path_link = paste0(unique_id, "_", row_number())) %>%
  ungroup()


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


final_columns <- c("occurrenceID", "associatedOccurrences", "associatedTaxa", "pathway", "basisOfRecord", "scientificName", "eventDate",
                   "countryCode", "taxonRank", "kingdom", "phylum", "class", "order", "family", "genus", "specificEpithet", "decimalLatitude", "decimalLongitude",
                   "geodeticDatum", "coordinateUncertaintyInMeters", "occurrenceStatus", "individualCount", "organismQuantity", "organismQuantityType",
                   "occurrenceRemarks", "dataGeneralizations", "country", "locality", "verbatimLocality", "identificationRemarks", "identifiedBy",
                   "datasetName", "eventID", "parentEventID", "bibliographicCitation", "rightsHolder", "license", "accessRights", "eventRemarks")

dataGeneralisations <- c("region" = "Centre of smallest associated administrative district used",
                         "regional" = "Centre of smallest associated administrative district used",
                         "Niakhar" = "Centre of smallest associated administrative district used",
                         "country" = "Centroid of the country used as no finer scale geographic information available",
                         "site" = "Coordinates reflect the location of the study site where the sample was collected")

prep_df <- rodent_classifications %>%
  bind_cols(as_tibble(matrix(nrow = nrow(rodent_classifications), ncol = length(gbif_col_names)), .name_repair = ~ gbif_col_names)) %>%
  relocate(colnames(gbif_template)) %>%
  mutate(occurrenceID = paste0(unique_id, "_", record_id),  # The occurrenceID submitted to GBIF will be a combination of the study identifier and record identifier
         basisOfRecord = "Human observation", # All rodents were trapped in either rodent traps or brought to the researchers by local hunters
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
         eventDate = case_when(is.na(year_trapping) ~ as.character(NA),
                               year_start == year_end & is.na(month_start) ~ year_start,
                               year_start != year_end & is.na(month_start) ~ paste0(year_start, "/", year_end),
                               year_start == year_end & is.na(month_end) ~ paste0(year_start, "-", month_start),
                               year_start == year_end & month_start == month_end ~ paste0(year_start, "-", month_start),
                               year_start == year_end & month_start != month_end ~ paste0(year_start, "-", month_start, "/", year_end, "-", month_end),
                               year_start != year_end ~ paste0(year_start, "-", month_start, "/", year_end, "-", month_end))) %>%
  select(-endDayOfYear, -year, -month, -day, -verbatimEventDate) %>%
  mutate(eventRemarks = as.character(NA)) %>%
  left_join(., species %>%
              rename(gbif_species = classification) %>%
              mutate(gbif_id = as.character(gbif_id)), by = "gbif_id") %>% # Merge GBIF IDs with accepted scientific name
  mutate(scientificName = str_to_sentence(case_when(!is.na(gbif_species) ~ gbif_species,
                                                    !is.na(genus) ~ genus,
                                                    TRUE ~ as.character(NA))),  # Use scientific name if available, otherwise genus association
         specificEpithet = case_when(nchar(str_split(scientificName, " ", simplify = TRUE)[, 2]) > 1 ~ str_split(scientificName, " ", simplify = TRUE)[, 2],
                                    TRUE ~ as.character(NA)),
         countryCode = countrycode(country, "country.name", "iso2c")) %>%
  select(-higherClassification, -infraspecificEpithet) %>%
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
         individualCount = number,
         occurrenceStatus = case_when(individualCount >= 1 ~ "Present",
                                      TRUE ~ "Absent")) %>%
  left_join(., study_titles, by = "unique_id") %>%
  left_join(., final_citations, by ="unique_id") %>%
  mutate(identifiedBy = authors,
         datasetName = title,
         parentEventID = unique_id,
         eventID = paste0(parentEventID, "_", if_else(is.na(eventDate), "", paste0(eventDate, "_")), locality),
         bibliographicCitation = full_citation,
         order = case_when(scientificName == "rodentia" ~ "rodentia",
                           TRUE ~ order),
         phylum = "chordata",
         kingdom = "animalia") %>%
  rename(license = license.y) %>%
  mutate(accessRights = case_when(str_detect(license, "all|distribution") ~ license,
                                  TRUE ~ NA),
         license = case_when(str_detect(license, "CC") ~ license,
                             TRUE ~ NA),
         rightsHolder = case_when(str_detect(rightsHolder, "authors") ~ authors,
                                  TRUE ~ rightsHolder))

# Add pathogen data -------------------------------

# Cleaning species names
pathogen_data <- read_rds(here("data_raw", "pathogen.rds")) %>%
  group_by(unique_id) %>%
  mutate(country = as_factor(country),
         path_link =  paste0(unique_id, "_", row_number()))

pathogen_tested <- c("path_1", "path_2", "path_3", "path_4", "path_5", "path_6")

path_all <- pathogen_data %>%
  ungroup() %>%
  dplyr::select(path_link, all_of(pathogen_tested), path_1_tested, path_2_tested, path_3_tested, path_4_tested, path_5_tested, path_6_tested,
                pcr_path_1_positive, pcr_path_2_positive, pcr_path_3_positive, pcr_path_4_positive, pcr_path_5_positive, pcr_path_6_positive,
                ab_ag_path_1_positive, ab_ag_path_2_positive, ab_ag_path_3_positive, ab_ag_path_4_positive, ab_ag_path_5_positive,
                culture_path_1_positive, culture_path_2_positive, culture_path_3_positive, histo_path_1_positive, histo_path_2_positive,
                histo_path_3_positive, histo_path_4_positive, histo_path_5_positive, histo_path_6_positive)

long_pathogen <- path_all %>%
  tibble() %>%
  pivot_longer(cols = c(path_1_tested, path_2_tested, path_3_tested, path_4_tested, path_5_tested, path_6_tested,
                        pcr_path_1_positive, pcr_path_2_positive, pcr_path_3_positive, pcr_path_4_positive, pcr_path_5_positive, pcr_path_6_positive,
                        ab_ag_path_1_positive, ab_ag_path_2_positive, ab_ag_path_3_positive, ab_ag_path_4_positive, ab_ag_path_5_positive,
                        culture_path_1_positive, culture_path_2_positive, culture_path_3_positive, histo_path_1_positive, histo_path_2_positive,
                        histo_path_3_positive, histo_path_4_positive, histo_path_5_positive, histo_path_6_positive)) %>%
  drop_na(value) %>%
  rename("assay" = name,
         "number" = value) %>%
  mutate(pathogen_tested = case_when(str_detect(assay, "path_1") ~ path_1,
                                     str_detect(assay, "path_2") ~ path_2,
                                     str_detect(assay, "path_3") ~ path_3,
                                     str_detect(assay, "path_4") ~ path_4,
                                     str_detect(assay, "path_5") ~ path_5,
                                     str_detect(assay, "path_6") ~ path_6,
                                     TRUE ~ "error")) %>%
  dplyr::select(-all_of(pathogen_tested))

pathogen_dictionary <- as.list(c("leishmania major", "lassa mammarenavirus", "toxoplasma gondii", "usutu virus", "taenia", "escherichia coli",
                                 "klebsiella pneumoniae", "arenaviridae", "trichuridae", "borrelia", "bartonella", "trypanosomatidae", "arenaviridae",
                                 "schistosoma", "borrelia crocidurae", "leishmania", "arenaviridae", "arenaviridae", "leptospira", "lassa mammarenavirus",
                                 "rift valley fever phlebovirus", "phlebovirus", "flavivirus", "hantaviridae", "lassa mammarenavirus", "lassa mammarenavirus",
                                 "strongyloides", "babesia", "eimeria", "plasmodium", "trypanosoma lewisi", "coxiella burnetii", "bartonella", "borrelia",
                                 "ehrlichia", "mycoplasma", "orentia", "rickettsia", "mammarenavirus", "orthopoxvirus", "lassa mammarenavirus", "orthopoxvirus",
                                 "mycobacterium", "borrelia", "schistosoma mansoni", "plagiorchis", "anaplasma", "mycobacterium", "lassa mammarenavirus",
                                 "lassa mammarenavirus"))
names(pathogen_dictionary) <- c("leishmania_major", "lassa_mammeranvirus", "toxoplasma_gondii", "usutu_virus", "hydatigera_species", "e_coli_esbl",
                                "k_pneumoniae_esbl", "arenaviruses", "trichuris_spp", "borrellia_spp", "bartonella_spp", "trypanosoma_spp", "arenaviridae_spp",
                                "schistosoma_spp", "borrelia_crocidurae", "leishmania_spp", "arenavirus", "arenavirus_spp", "leptospirosis_spp", "lassa_mammarenavirus",
                                "rift_valley_fever", "phleboviruses", "flavivirus", "hantavirus", "lassa_mammarenavirus_antigen", "lassa_mammarenavirus_antibody",
                                "strongyloides_spp", "babesia_spp", "eimeria_spp", "plasmodium_spp", "trypanasoma_lewisi", "coxiella_burnetii", "bartonella", "borrelia",
                                "ehrlichia", "mycoplasma_spp", "orentia", "rickettsia", "mammarenavirus", "orthopoxvirus", "lassa", "orthopoxvirus_spp",
                                "mycobacterium_spp", "borrelia_spp", "schistosoma_mansoni", "plagiorchis_species", "anaplasma", "mycobacteria_spp", "lassa_mammarenavirus_IgG_Ab",
                                "lassa_mammarenavirus_Ag")

path_species <- unique(unlist(pathogen_dictionary))

# Pull the hierarchy for each pathogen genus from gbif
if(!file.exists(here("data_clean", "pathogen_genus_hierarchy.rds"))){
  path_species[path_species == "trypanosoma lewisi"] = "Trypanosoma"
  path_species[path_species == "babesia"] = "Aconoidasida"
  path_species[path_species == "orentia"] = "Orientia"
  
  path_genera <- classification(snakecase::to_sentence_case(path_species), db = "gbif")
  path_classification <- as_tibble(do.call(rbind,c(path_genera, make.row.names = T)), rownames = "path") %>%
    mutate(path = str_remove_all(path, ".[1-9]")) %>%
    pivot_wider(id_cols = path, names_from = rank, values_from = name)
  
  path_classification$path[path_classification$path == "Trypanosoma"] = "Trypanosoma lewisi"
  path_classification$path[path_classification$path == "Aconoidasida"] = "Babesia"
  
  write_rds(path_classification %>%
              mutate(species = case_when(nchar(str_split(path, pattern = " ", simplify = TRUE)[, 2]) > 1 ~ path,
                                         TRUE ~ species)), here("data_clean", "pathogen_genus_hierarchy.rds"))

  } else {
  
  path_classification <- read_rds(here("data_clean", "pathogen_genus_hierarchy.rds"))
  
}

pathogen_prep <- long_pathogen %>%
  mutate(pathogen_name = recode(pathogen_tested, !!!pathogen_dictionary),
         group = case_when(str_detect(assay, "_1_") ~ 1,
                           str_detect(assay, "_2_") ~ 2,
                           str_detect(assay, "_3_") ~ 3,
                           str_detect(assay, "_4_") ~ 4,
                           str_detect(assay, "_5_") ~ 5,
                           str_detect(assay, "_6_") ~ 6)) %>%
  group_by(group) %>%
  group_split() %>%
  lapply(., function(x) {
  x %>%
    mutate(method = case_when(str_detect(assay, "_positive") ~ str_split(assay, "_", simplify = TRUE)[,1]),
           method = case_when(str_detect(method, "pcr") ~ "PCR",
                              str_detect(method, "ab|ag") ~ "Serology",
                              str_detect(method, "culture") ~ "Culture",
                              str_detect(method, "hist") ~ "Histology"),
           assay = case_when(str_detect(assay, "test") ~ "tested",
                             str_detect(assay, "pos") ~ "positive")) %>%
    group_by(path_link) %>%
    fill(method, .direction = "updown") %>%
    pivot_wider(names_from = assay, values_from = number)
}) %>%
  bind_rows() %>%
  filter(tested != 0) %>%
  group_by(path_link) %>%
  mutate(occurrenceID = case_when(n() > 1 ~ paste0(path_link, "_path_"),
                                  TRUE ~ paste0(path_link, "_path")),
         path_number = row_number(),
         occurrenceID = case_when(str_detect(occurrenceID, "_path_") ~ paste0(occurrenceID, path_number),
                                  TRUE ~ occurrenceID)) %>%
  ungroup() %>%
  mutate(pathway = "parasites on animals",
         basisOfRecord = "Material sample",
         Organism = str_to_sentence(pathogen_name)) %>%
  left_join(path_classification %>%
              rename(scientificName = species), by = c("Organism" = "path")) %>%
  mutate(taxonRank = case_when(!is.na(scientificName) ~ "species",
                               !is.na(genus) ~ "genus",
                               !is.na(family) ~ "family",
                               !is.na(order) ~ "order",
                               !is.na(class) ~ "class",
                               !is.na(phylum) ~ "phylum",
                               TRUE ~ kingdom),
         specificEpithet = case_when(!is.na(scientificName) ~ scientificName,
                                     TRUE ~ as.character(NA)),
         occurrenceStatus = case_when(positive >= 1 ~ "Present",
                                       TRUE ~ "Absent"),
         organismQuantity = positive,
         organismQuantityType = "samples",
         occurrenceRemarks = paste0("Number of individuals tested = ", tested),
         identificationRemarks = paste0("Assayed using ", method))
  

pathogen_rodent_link <- pathogen_prep %>%
  select(path_occurrenceID = occurrenceID, path_link, pathogen_name, occurrenceStatus) %>%
  left_join(prep_df %>%
              select(occurrenceID, path_link, scientificName), by = c("path_link")) %>%
  mutate(associatedOccurrences = case_when(occurrenceStatus == "Present" ~ paste0("parasites collected from '", occurrenceID, "'"),
                                           TRUE ~ as.character(NA)),
         associatedTaxa = case_when(occurrenceStatus == "Present" ~ paste0("parasite of ", str_to_lower(scientificName)),
                                    TRUE ~ as.character(NA))) %>%
  distinct(occurrenceID = path_occurrenceID, associatedOccurrences, associatedTaxa)

pathogen_prep_df <- pathogen_prep %>%
  left_join(pathogen_rodent_link) %>%
  left_join(prep_df %>%
              select(!any_of(names(pathogen_prep)[-1])),
            by  = "path_link") %>%
  select(any_of(final_columns)) %>%
  mutate(individualCount = NA)

rodent_pathogen_link <- prep_df %>%
  select(occurrenceID, path_link) %>%
  right_join(pathogen_prep %>%
               select(path_occurrenceID = occurrenceID, path_link, pathogen_name, method, organismQuantity)) %>%
  filter(organismQuantity != 0) %>%
  group_by(occurrenceID) %>%
  mutate(associatedOccurrences = paste0("host of '", path_occurrenceID, collapse = "', ", "'"),
         associatedTaxa = paste0("host of ", pathogen_name, collapse = ", ")) %>%
  distinct(occurrenceID, associatedOccurrences, associatedTaxa)

rodent_prep_df <- prep_df %>%
  left_join(rodent_pathogen_link) %>%
  select(any_of(final_columns))

# Build the GBIF dataframe ------------------------------------------------

# Add pathogen_prep_df to rodent_prep_df
final_df <- bind_rows(rodent_prep_df, pathogen_prep_df) %>%
  select(any_of(final_columns)) %>%
  mutate(scientificName = case_when(nchar(str_split(scientificName, pattern = " ", simplify = TRUE)[, 2]) >= 1 ~ scientificName,
                                          TRUE ~ as.character(NA)),
         kingdom = str_to_sentence(kingdom),
         phylum = str_to_sentence(phylum))

dir.create(here("final_data"))
write_rds(final_df, here("data_clean", "occurrence_df.rds"))
write_tsv(final_df, here("final_data", "taxa.txt"))
