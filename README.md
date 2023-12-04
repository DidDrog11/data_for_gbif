# Data for GBIF
[![DOI](https://zenodo.org/badge/599957741.svg)](https://zenodo.org/badge/latestdoi/599957741)

This repository is to clean and organise the data from a paper published on rodent trapping across [West Africa](https://github.com/DidDrog11/scoping_review).

The field names are available from a GBIF template and descriptions of the fields here (https://www.gbif.org/data-quality-requirements-occurrences).

This repository includes the code to pull down the data from the Google sheets document that was used for data extraction, clean the underlying data and format it in the structure required for GBIF for both the pathogen and rodent species'.

# Data abstraction

The following are the variables included in the Google sheets document.

## Study

The study sheet contained the following abstracted variables on information contained in the included studies.

| Variable name | Description | Values |
| --- | -- | -- |
| link | Digital Object Identifier (DOI) or persistent weblink to manuscript or report | URL |
| year_publication | Year of publication of manuscript or report | numeric |
| title | Title of manuscript or report | text |
| journal_name | Name of the journal or resource where the data was published | text |
| aim_1 - aim_3 | Stated aims of the study within the manuscript or report | text |
| first_author | The name of the first author of the publication or report | text |
| reference_uid | A unique reference (within this dataset) for the manuscript or report (typically the DOI, ISSN or ISBN) | text |
| unique_id | A unique identifier (within this dataset) for each study | text |
| metric | The reported outcome of small-mammal sampling within the study, typically abundance (number of detected individuals) or presence and absence | text |
| trap_types | The reported design of traps used within the study | text |
| trapping_method | The design of the trapping survey used within the report | text |
| repeated_visit | Whether the study utilised repeated sampling of the same locations | yes/no |
| geolocation_level | The resolution of geographic coordinates reported in the study | text |
| specieation | The method reported to identify detected individuals to species level | text |
| aim | The categorised aim of the study, whether to investigate rodent and small-mammal ecology or to investigate the risk of zoonoses | text |
| aim_detail | Further sub-categorisation of the study aims within the domains of rodent/small-mammal ecology or zoonosis risk | text |
| species_accumulation | Whether the study reported a species-accumulation function to estimate completeness of sampling | yes/no |
| diversity_measurement | Whether the study reported a metric of species diversity | yes/no |
| trapping_effort | Whether the study reported trapping effort at the same level as geolocation_level | yes/no/incomplete |
| pathogen | Whether the study included data on zoonotic pathogens or rodent pathogens | yes/no/rodent pathogen |

## Rodent

The rodent sheet contained the following abstracted variables on information contained in the included studies.

| Variable name | Description | Values |
| --- | -- | -- |
| unique_id | This unique identifier links the rodent species data to the source of this information (in the study sheet) | text |
| year_trapping | The year trapping was conducted | numeric |
| month_trapping | The month trapping was conducted (if reported) | text |
| country | The country in which trapping was conducted | text |
| region | The region in which trapping was conducted | text |
| town_village | The town or village in which trapping was conducted | text |
| latitude_dms | The latitude of sampling coordinates in degrees, minutes and seconds, separated by "_". Positive numbers indicate North of the equator | text |
| longitude_dms | The longitude of sampling coordinates in degrees, minutes and seconds, separated by "_". Positive numbers indicate West of the Greenwich Meridian | text |
| latitude_D_N | The latitude of sampling coordinates in decimal degrees. Positive numbers indicate North of the equator | text |
| longitude_D_E | The longitude of sampling coordinates in decimal degrees. Positive numbers indicate East of the Greenwich Meridian | text |
| UTM_coordinates | The Universal Transverse Mercator (UTM) coordinates of sampling | text |
| habitat | The reported habitat type in which sampling occurred | text |
| intensity_use | The intensity of anthropogenic pressure in the sampling location, categorised as intense if in or around areas of human habitation | text |
| genus | The reported genus of detected small mammals | text |
| species | The reported species of detected small mammals | text |
| number | The reported number of detected small mammals | numeric |
| trap_nights | The reported number of trap nights at the resolution of reported sampling | numeric |
| capture_rate | The reported capture rate or trap success at the resolution of reported sampling | numeric |
| trap_night_unit | The level of resolution for trap_nights reported in the study | text |
| study_nights | The number of nights of sampling at that location if number of trap_nights cannot be inferred | text |

## Pathogen

The pathogen sheet contained the following abstracted variables on information contained in the included studies.

| Variable name | Description | Values |
| --- | -- | -- |
| unique_id | This unique identifier links the pathogen data to the source of this information (in the study sheet) | text |
 year_trapping | The year trapping was conducted | numeric |
| month_trapping | The month trapping was conducted (if reported) | text |
| country | The country in which trapping was conducted | text |
| region | The region in which trapping was conducted | text |
| town_village | The town or village in which trapping was conducted | text |
| habitat | The reported habitat type in which sampling occurred | text |
| genus | The reported genus of detected small mammals | text |
| species | The reported species of detected small mammals | text |
| path_1-path_6 | The names of the pathogens tested in the study | text |
| latitude_dms | The latitude of sampling coordinates in degrees, minutes and seconds, separated by "_". Positive numbers indicate North of the equator | text |
| longitude_dms | The longitude of sampling coordinates in degrees, minutes and seconds, separated by "_". Positive numbers indicate West of the Greenwich Meridian | text |
| latitude_D_N | The latitude of sampling coordinates in decimal degrees. Positive numbers indicate North of the equator | text |
| longitude_D_E | The longitude of sampling coordinates in decimal degrees. Positive numbers indicate East of the Greenwich Meridian | text |
| UTM_coordinates | The Universal Transverse Mercator (UTM) coordinates of sampling | text |
| path_1_tested-path_6_tested | The number of individual small-mammals tested for the pathogen | numeric |
| pcr_path_1_positive-pcr_path_6_positive | The number of individual small-mammals positive for this pathogen by PCR assay | numeric |
| ab_ag_path_1_positive-ab_ag_path_5_positive | The number of individual small-mammals positive for this pathogen by antibody or antigen based assays | numeric |
| culture_path_1_positive-culture_path_3_positive | The number of individual small-mammals positive for this pathogen by culture based assays | numeric |
| histo_path_1_positive-histo_path_6_positive | The number of individual small-mammals positive for this pathogen by histological or pathological based assays | numeric |
