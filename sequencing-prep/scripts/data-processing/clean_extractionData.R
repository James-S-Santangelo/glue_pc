# Script to assess GLUE DNA extractions for useable plants to pass on to library
# prep
#
# Authors: James S. Santangelo and Sophie Koch

##################################################
#### LOAD AND CONCATENATE RAW EXTRACTION DATA ####
##################################################

# Read in data with Qubit concentrations for plants assayed by Beata
# Some plants were re-extracted and assayed twice
Qubit <- read_csv("data/raw/Qubit_extractions.csv", na = c("", "na")) %>% 
  select(row, column, continent, name, pop, site, individual, "Qubit_ng.ul Ext 1", "Qubit_ng.ul Ext 2 (SK)", "original comments") %>% 
  rename("qubit_1" = "Qubit_ng.ul Ext 1",
         "qubit_2" = "Qubit_ng.ul Ext 2 (SK)",
         "comments_1" = "original comments") %>% 
  mutate(pop = as.character(pop),
         individual = as.character(individual),
         site = as.character(site)) %>% 
  filter(name != "Almada") %>% # Not using Almada for now
  arrange(., name)

# Read in plants assayed by Sophie. 
new_ext <- read_csv("data/raw/New_extractions.csv", na = c("", "na")) %>% 
  dplyr::select(row, column, continent, name, pop, site, individual, "Qubit_ng.ul", "new comments") %>% 
  rename("qubit_3" = "Qubit_ng.ul",
         "comments_2" = "new comments") %>% 
  mutate(pop = as.character(pop),
         individual = as.character(individual),
         site = fct_recode(site, "u" = "U"),
         site = as.character(site)) %>% 
  # name = as.character(fct_recode(name, "Almada" = "Almada-2"))) %>%
  filter(name != "Almada")  %>%
  distinct() %>% 
  arrange(., name)

# Combine Beata and Sophie's Qubit data
all_ext <- full_join(Qubit, new_ext, by = c("name", "pop", "individual", "continent")) %>% 
  mutate(qubit_1 = as.numeric(qubit_1),
         qubit_2 = as.numeric(qubit_2),
         qubit_3 = as.numeric(qubit_3),
         
         # Handle NAs in site info for two cities
         site.x = case_when(name == "Wollongong" & pop %in% c("1", "2", "3", "4", "5") ~ "u",
                            TRUE ~ site.x),
         site.y = case_when(name == "Punta Arenas" & pop == "22" ~ "r",
                            TRUE ~ site.y),
         
         # Make sure urban and rural site infor matches between datasets
         site = case_when(
           (!(is.na(site.x)) & is.na(site.y)) ~ site.x,
           (!(is.na(site.y)) & is.na(site.x)) ~ site.y,
           site.x == site.y ~ site.x,
           TRUE ~ "ERR")) %>% 
  rename("row_firstPlating" = "row.x",
         "column_firstPlating" = "column.x",
         "row_secondPlating" = "row.y",
         "column_secondPlating" = "column.y",
         "city" = "name") %>% 
  
  # Remove plants with no Qubit info
  filter(!((is.na(qubit_1) & is.na(qubit_2) & is.na(qubit_3))) | !(is.na(individual))) %>% 
  mutate(city = str_replace(city, "ö", "o")) %>% 
  separate(city, into = c("city", "extra"), sep = "[,|;]") %>% 
  mutate(city = str_replace(city, pattern = " ", replacement = "_")) %>% 
  dplyr::select(-site.x, -site.y, -extra) %>% 
  dplyr::select(continent, city, pop, site, individual, row_firstPlating, column_firstPlating, 
         row_secondPlating, column_secondPlating, qubit_1, qubit_2, qubit_3, comments_1, comments_2) %>% 
  arrange(., city) %>% 
  mutate(plantID = paste(city, pop, individual, sep = "_")) %>% 
  dplyr::select(-contains("comment")) %>% 
  dplyr::select(continent, city, pop, individual, site, plantID, contains("qubit"), everything()) 

#############################################################################
#### ASSESS BREAKDOWN OF USEABLE PLANTS BY CITY, HABITAT, AND POPULATION ####
#############################################################################

# Load cline summary 
clineSummary <- read_csv("data/raw/linearClinesSummary.csv") %>% 
  dplyr::select(pvalHCN, city) %>% 
  mutate(city = as.character(fct_recode(city, "Munster" = "Muenster")))

# Useable plants with minimum concentration of 10 ng/uL
Goodplants_10 <- all_ext %>% 
  group_by(continent, city) %>% 
  mutate(is_good = case_when(qubit_1 > 10 | qubit_2 > 10 | qubit_3 > 10 ~ 1,
                             TRUE ~ 0 )) %>% 
  filter(is_good == 1) %>% 
  arrange(continent, city) %>% 
  ungroup()

# Number of good plants by site and population
numGoodplants_10 <- Goodplants_10 %>% 
  group_by(continent, city, site) %>% 
  summarize(total_plants = sum(is_good)) %>% 
  spread(key = site, value = total_plants) %>% 
  ungroup() %>% 
  left_join(., clineSummary, by = "city") %>% 
  left_join(., Goodplants_10 %>% 
              select(city, site, pop) %>% 
              arrange(city, as.numeric(pop)) %>% 
              group_by(city, site) %>%
              distinct() %>% 
              mutate(pop_reLabeled = 1:n(),
                     sitePop = paste(site, pop_reLabeled, sep = "_")) %>% 
              left_join(., Goodplants_10, by = c("city", "site", "pop")) %>% 
              group_by(city, sitePop) %>% 
              summarize(total_plants = sum(is_good)) %>% 
              spread(key = sitePop, value = total_plants),
            by = "city") %>% 
  mutate(significant = ifelse(pvalHCN < 0.05, "Yes", "No")) %>% 
  dplyr::select(-pvalHCN) %>% 
  arrange(continent, city) %>% 
  left_join(Goodplants_10 %>% 
              select(city, site, pop) %>% 
              arrange(city, as.numeric(pop)) %>% 
              group_by(city, site) %>%
              distinct() %>% 
              mutate(pop_reLabeled = 1:n(),
                     sitePop = paste(site, pop_reLabeled, sep = "_"),
                     sitePop_map = paste(sitePop, pop, sep = "=")) %>% 
              ungroup() %>% 
              group_by(city) %>% 
              summarise(sitePop_mapped = paste(sitePop_map, collapse = "; ")),
            by = "city") %>%
  dplyr::select(continent, city, significant, everything())

outpath <- 'data/clean/extractions/'
print(sprintf('Cleaned extraction data saved to %s', outpath))
write_csv(numGoodplants_10, paste0(outpath, 'numGoodplants_10.csv'))
write_csv(all_ext, paste0(outpath, 'allExtractions.csv'))

