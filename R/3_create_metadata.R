david <- list(individualName = list(givenName = "David",
                                    surName = "Simons"),
              organizationName = list("The Royal Veterinary College, London, United Kingdom"),
              electronicMailAddress = list("dsimons19@rvc.ac.uk"))

lauren <- list(individualName = list(givenName = "Lauren A",
                                     surName = "Attfield"),
               organizationName = list("Imperial College London, London, United Kingdom"))

kate <-  list(individualName = list(givenName = "Kate E",
                                    surName = "Jones"),
              organizationName = list("University College London, London, United Kingdom"))

deborah <- list(individualName = list(givenName = "Deborah",
                                      surName = "Watson-Jones"),
                organizationName = list("London School of Hygiene and Tropical Medicine, London, United Kingdom"))

richard <-  list(individualName = list(givenName = "Deborah",
                                       surName = "Watson-Jones"),
                 organizationName = list("The Royal Veterinary College, London, United Kingdom"))

intellectualRights = "This work is licensed under a CC0 1.0 Universal (CC0 1.0) Public Domain Dedication"

abstract <- "This dataset contains records of the small-mammals, primarily Rodentia, trapped across West Africa and included in published literature that was obtained through a comprehensive review of rodent trapping studies across the region. Data were extracted from studies to produce detection, non-detection histories at described trapping locations. Geolocation were perfored on all included locations. Individual species that were detected elsewhere in the contributing studies records but were not reported from a trapping location were interepreted as non-detections and so explicitly recorded as 0 counts."

keywordSet <- list(
  list(keywordThesaurus = "GBIF Dataset Type Vocabulary",
       keyword = list("Occurrence"))
)

my_eml <- list(system = "http://gbif.org",
               scope = "system",
               schemaLocation="eml://ecoinformatics.org/eml-2.1.1 http://rs.gbif.org/schema/eml-gbif-profile/1.1/eml.xsd",
               packageId = "22p01224-12z5-12d4-91r2-32cp1a81f14g/eml-v1.0.xml",
               dataset = list(
                 title = "Data from: Rodent trapping studies as an overlooked information source for understanding endemic and novel zoonotic spillover",
                 creator = list(david, lauren, kate, deborah, richard),
                 pubDate = "2023",
                 intellectualRights = intellectualRights,
                 abstract = abstract,
                 keywordSet = keywordSet,
                 contact = david,
                 coverage = list(geographicCoverage = list(geographicDescription = "This dataset covers at least 1,661 sampling sites, from 14 West African countries.",
                                                           boundingCoordinates = list(westBoundingCoordinate = "-23.67",
                                                                                      eastBoundingCoordinate = "21.84",
                                                                                      northBoundingCoordinate = "23.9",
                                                                                      southBoundingCoordinate = "21.83")),
                                 taxonomicCoverage = list(generalTaxonomicCoverage = "This dataset contains occurrence records of species from four orders of Chordata, namel Afrosoicida, Erinaceomorpha, Rodentia and Soricomorpha",
                                                          taxonomicClassification = list(taxonRankName = "phylum",
                                                                                         taxonRankValue = "Chordata"))),
                 project = list(
                   title = "Data from: Rodent trapping studies as an overlooked information source for understanding endemic and novel zoonotic spillover",
                   personnel = list(individualName = list(givenName = david$individualName$givenName,
                                                          surName = david$individualName$surName),
                                    role = list("PrincipalInvestigator")),
                   funding = list("Data collection and analysis were funded through a a PhD award from the UK Biotechnology and Biological Sciences Research Council [BB/M009513/1]."),
                   studyAreaDescription = list(descriptor = list(name = "generic",
                                                                 citableClassificationSystem="false",
                                                                 descriptorValue = "Data were extracted from studies reporting rodent trapping in the West African region. The United Nations defined 16 countries within this region, including,  Benin, Burkina Faso, Cape Verde, The Gambia, Ghana, Guinea, Guinea-Bissau, Ivory Coast, Liberia, Mali, Mauritania, Niger, Nigeria, Senegal, Sierra Leone, and Togo, as well as Saint Helena, Ascension and Tristan da Cunha (United Kingdom Overseas Territory).")))))

write_eml(my_eml, here("final_data", "eml.xml"))
eml_validate(here("final_data", "eml.xml"))
