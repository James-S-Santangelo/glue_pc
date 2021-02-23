# Script to clean city center lat/longs for merging with pop means. 
# City-center Lat/longs from https://www.latlong.net/ but may have been
#   manually moved if we disagreed with its placement. 
# Script also adds continent and country to dataset
#
# Author: James. S. Santangelo

# Load city center dataframe, which was manually created in Excel
city_centers <- read_csv("data/raw/latLongs_cityCenters.csv")

city_center_modified <- city_centers %>%
  
  # If coordinate was manually moved, use new location
  mutate(latitude_city = ifelse(is.na(Latitude_moved), 
                                Latitude_web, Latitude_moved),
         longitude_city = ifelse(is.na(Longitude_moved), 
                                 Longitude_web, Longitude_moved),
         
         # Remove/replace special characters (e.g., in city names)
         city = gsub("[;|,|(].*$", "", City),
         city = str_trim(city, side = "right"),
         city = str_replace_all(city, c(" " = "_", "\\." = "",
                                        'ü' = 'u', 'ï' = 'i',
                                        'ë' = 'e', 'ä' = 'a',
                                        'ö' = 'o')),
         city = fct_recode(city, "Washington" = "Washington_DC",
                           "Frankfurt" = "Frankfurt_am_Main",
                           "Sioux_Falls" = "Sioux_Falls_East")) %>%
  filter(city != "Sioux_Falls_West") %>% 
  dplyr::select(city, latitude_city, longitude_city, Country) %>%
  na.omit()

# Add continent
city_center_modified$continent <- coords2continent(city_center_modified %>% 
                                                   dplyr::select(longitude_city, latitude_city))

# Refactor levels and manually add continents for cities that were missed by function
city_center_modified <- city_center_modified %>%
  mutate(continent = fct_recode(continent, "South America" = "South America and the Caribbean",
                                "Oceania" = "Australia"),
         continent = case_when(city == "Charlottetown" ~ "North America",
                               city == "Halifax" ~ "North America",
                               city == "Lisbon" ~ "Europe",
                               city == "Punta_Arenas" ~ "South America",
                               city == "Quebec_City" ~ "North America",
                               city == "Saint_John" ~ "North America",
                               city == "Stockholm" ~ "Europe",
                               Country == "Mexico" ~ "North America",
                               Country == "USA" ~ "North America",
                               Country == "Tasmania" ~ "Oceania",
                               Country == "Colombia" ~ "South America",
                               TRUE ~ as.character(continent)))

# Write lat long dataset to disk
write_csv(city_center_modified, path = "data/clean/latLong_cityCenters_clean.csv",
          col_names = TRUE)
