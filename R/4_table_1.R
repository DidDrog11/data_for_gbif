source(here::here("R", "libraries.R"))

studies <- readRDS(here("data_clean", "studies.rds")) %>%
  filter(!is.na(reference_uid))

bib <- bib2df(here("citations", "include_final.bib")) %>%
  select(TITLE, JOURNALTITLE, SERIES, BIBTEXKEY, DOI, DOI.1, PMID, PMCID, ISSN, ISBN, URL) %>%
  mutate(reference = as.character(coalesce(DOI, PMID, ISSN, ISBN)))

tbl_1 <- left_join(studies, bib,
                         by = c("reference_uid" = "reference")) %>%
  select(year_publication, first_author, title, journal_name, reference_uid) %>%
  mutate(year_publication = as.integer(year_publication)) %>%
  select(year_publication, first_author, title, journal_name, reference_uid) %>%
  rename("Year publication" = year_publication,
         "Author" = first_author,
         "Title" = title,
         "Journal/Publication" = journal_name,
         "DOI/PMCID/ISSN/ISBN" = reference_uid) %>%
  arrange(`Year publication`, `Author`)

print(xtable(tbl_1, type = "latex"), file = here("manuscript", "table_1.tex"))
