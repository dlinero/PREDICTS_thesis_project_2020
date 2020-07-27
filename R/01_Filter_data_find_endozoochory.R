# 2. Function that updates the vector of species names that have endozoochorus dispersal syndrome by 
# taking a specific dataset and searching for the TRUE condition for endozoochory

Filter_endozoochory_terms <- function(DatasetName, FilterValue ) {
  c(Endoo_PlantSpecies, as.character(Dispersal_syndrome %>% # Update vector
                                       
                                       # Select specific dataset
                                       subset(Dataset == DatasetName) %>%
                                       
                                       # Filter out all of the species that have already been identified 
                                       # as having endozoochory
                                       subset(Updatedname %nin% Endoo_PlantSpecies) %>% 
                                       
                                       # Filter values with condition TRUE for endozoochory
                                       filter(str_detect(OrigValueStr, FilterValue)) %>%
                                       
                                       # Get unique species names
                                       distinct(Updatedname) %>%
                                       
                                       # Make a vector
                                       droplevels() %>% pull(Updatedname)))
}
