# 1. Function that returns the first record for each SSBS (SSBS = Concatenation of Source_ID, Study_number, Block and Site_number)

Number_sites_landuses <- function(dataframe){
  # Get a vector with the SSBS names 
  Sites <- dataframe %>% pull(SSBS) %>% as.character()
  # Find the index of the first match for each SSBS and eliminate the duplicates
  Indeces <- match(Sites, dataframe$SSBS) %>% unique()
  # Get a data frame with only the first matches of each SSBS
  justFirsts <- dataframe[Indeces, ] 
  return(justFirsts)
}


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
  






