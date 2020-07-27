# clear workspace
rm(ls=list())

# Load libraries
library(dplyr)
library(stringr)
library(tidyr)

# Load richness data
diversity_all_richness <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Load abundance data 
diversity_all_abundance <- readRDS(file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")

# Load Simpsons data 
diversity_all_simpson <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Remove columns of the abundance table that don't match with the other tables
diversity_all_abundance <- diversity_all_abundance[ , -c(25,26)]

# Merge all tables
coordinates_sites <- base::rbind.data.frame(diversity_all_richness, 
                                      diversity_all_abundance, 
                                      diversity_all_simpson) 

# Get a list of sites for nePlants
nePlants <-coordinates_sites %>% 
  
  # subset nePlants
  subset(Kingdom == "nePlantae") %>% 
  
  # droplevels()
  droplevels() %>%
  
  # Remove duplicates
  distinct(SSBS) 

# Get a list of sites for nePlants
Plants <-coordinates_sites %>% 
  
  # subset fleshy-fruited plants
  subset(Kingdom == "Plantae") %>% 
  
  # droplevels 
  droplevels() %>%
  
  # remove duplicates
  distinct(SSBS) 


# Get the sites where both types of plants were assesed
Both_plants <- Plants %>% 
  
  # subset studies that appear in bot nePlantae and Plantae 
  base::subset(SSBS %in% nePlants$SSBS) %>% 
  
  # create a vector
  pull() %>% 
  
  # transform into a character vector
  as.character()


coordinates_sites <- coordinates_sites %>% 
  
  
  # Create a new column to modify the names of the kingdoms
  mutate(
  newKingdom = base::paste(Kingdom),
  newKingdom = stringr::str_replace_all(newKingdom, c("Animalia" = "Frugivores",
                                                      "nePlantae" = "Abiotically-dispersed plants",
                                                 "Plantae" = "Fleshy-fruited plants"
                                                ))) %>%
  
  # Remove columns we're not interested in
  select(SSBS, Longitude, Latitude, Kingdom, newKingdom)


# Assign a new kindom name to the sites where both types of plants were assessed
for(i in 1:nrow(coordinates_sites)){
  if(coordinates_sites[i, 1] %in% Both_plants){
    coordinates_sites[i, 5] <- "Both types of plants"
  }
}


# Get unique coordinates for each site
coordinates_sites <- coordinates_sites %>%
  
  # Get one row for each study
  dplyr::distinct(SSBS, .keep_all = TRUE) 

# Export table
write.csv(coordinates_sites, "./output/map_files/location_sites.csv")
