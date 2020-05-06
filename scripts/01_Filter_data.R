rm(list = ls())
library(data.table)
library(dplyr)
library(ngram)
library(stringr)
library(taxize)
library(Hmisc)

# Import PREDICTS dataset 

diversity <- readRDS("./data/PREDICTS2020/PREDICTS2020.rds")

glimpse(diversity) # Take a look at the dimensions and structure of the data

# Filter tropics and biomes 

TropicsDiversity <- diversity %>%
  
  # Select the coordinates that fall within the tropics
  filter(between(Latitude, -23.43, 23.43)) %>% 
  
  # Select records located inside forest ecosystems
  subset(Biome == "Tropical & Subtropical Moist Broadleaf Forests" | Biome == "Tropical & Subtropical Dry Broadleaf Forests") %>%
  
  # drop to remove levels with zero records
  droplevels() 

# Export table
saveRDS(TropicsDiversity, file = "./output/intermediate_files/01_Filter_data_Tropical_records.rds")

# Import table
TropicsDiversity <- readRDS("./output/intermediate_files/01_Filter_data_Tropical_records.rds")

# Number of records per region
table(TropicsDiversity$Region)

# Number of records per biome
table(TropicsDiversity$Biome)

# Total number of records
length(TropicsDiversity$Source_ID)

# Total number of species
length(unique(TropicsDiversity$Best_guess_binomial))

# Record per land-use and intensity level
table(TropicsDiversity$Predominant_habitat, TropicsDiversity$Use_intensity)


# Filter frugivore animals in PREDICTS----------------------------------------------------------------
# Import mammal dataset:

# Eltontraits: https://figshare.com/articles/Data_Paper_Data_Paper/3559887
# http://www.esapubs.org/archive/ecol/E095/178/metadata.php#_ENREF_3

Mammal_Diet <- read.csv("./data/Eltontraits/MamFuncDat.csv", sep = ";")

Frugivore_mammals <- Mammal_Diet %>%
  
  # Filter mammals that have a diet composed of at least 50% fruit
  filter(Diet.Fruit >= 50) %>% 
  
  # Select the columns of scientific name, % of diet composed of fruits, and diet certainty level
  select(Scientific, Diet.Fruit, Diet.Certainty)  

# Number of frugivore mammal species
length(unique(Frugivore_mammals$Scientific))

# Import bird dataset

# EltonTraits: https://figshare.com/articles/Data_Paper_Data_Paper/3559887
# http://www.esapubs.org/archive/ecol/E095/178/metadata.php

Birds_Diet <- read.csv("./data/Eltontraits/BirdFuncDat.csv", sep = ";")

Frugivores <- Birds_Diet %>% 
  
  # Filter birds that have a diet composed of at least 50% fruits
  filter(Diet.Fruit >= 50) %>%
  
  # Select the columns of scientific name, % of diet composed of fruits, and diet certainty level
  select(Scientific, Diet.Fruit, Diet.Certainty)  %>%
  
  # Bind bird and mammal data
  rbind(Frugivore_mammals)  %>% 
  
  # Rename the scientific name column to match the PREDICTS dataset 
  rename(Best_guess_binomial = Scientific)

# Check frugivores nomenclature

# Get a character vector with the species names
Frugivore_Names <- as.character(Frugivores$Best_guess_binomial)

# Get the data source id for Catalogue of life
sources <- gnr_datasources() 

# Get the submitted species names and the corresponding name used in the Catalogue of Life (id = 1)
Frugivore_Names<- gnr_resolve(names = Frugivore_Names, data_source_ids = 1) %>% 
  
  # Identify records that have a different name in the Catalogue of Life (if the name is the same
  # the score will be 0.988)
  subset(score < 0.98) 

write.csv(Frugivore_Names, "./output/intermediate_files/01_Filter_data_Frugivores_nomenclature_check.csv")
# Get the records that belong to frugivore birds and mammals in the PREDICTS dataset 

# Merge both datasets and retain only the species of PREDICTS that match the species in the Frugivores table
PREDICTS_frugivores <- merge(x=TropicsDiversity, y=Frugivores, by = "Best_guess_binomial", all.x = FALSE) %>%
  droplevels()

# Export table
write.csv(PREDICTS_frugivores, "./output/cleaned_data/01_Filter_data_PREDICTS_Frugivores_records.csv")

# Import table
PREDICTS_frugivores <- read.csv("./output/cleaned_data/01_Filter_data_PREDICTS_Frugivores_records.csv")

# Get the number of species in the new dataset
length(unique(PREDICTS_frugivores$Best_guess_binomial))

# Number of records belonging to Birds and mammals
table(PREDICTS_frugivores$Class)

# Number of records belonging to different levels of uncertainty
table(PREDICTS_frugivores$Diet.Certainty, PREDICTS_frugivores$Class.y)

# Number of records per order
table(PREDICTS_frugivores$Order, PREDICTS_frugivores$Class.y)

# Number of records per order per region
table(PREDICTS_frugivores$Order, PREDICTSfrugivores$Region)

# Records per land-use and intensity level
table(PREDICTS_frugivores$Predominant_habitat, PREDICTS_frugivores$Use_intensity)

# Total number of sourceIDs
length(unique(PREDICTS_frugivores$Source_ID))

# Total number of studies
length(unique(PREDICTS_frugivores$SS))

# Total number of sites
length(unique(PREDICTS_frugivores$SSBS))

# Get the first record for each SSBS (SSBS = Concatenation of Source_ID, Study_number, Block and Site_number)
justFirstsFrugivores <- Number_sites_landuses(PREDICTS_frugivores)

# Number of sites per land-use type and land-use intensity
addmargins(table(justFirstsFrugivores$Predominant_habitat,justFirstsFrugivores$Use_intensity), 2)


# Export table of first matches
write.csv(justFirstsFrugivores, "./output/intermediate_files/01_Filter_data_First_matches_frugivores.csv")

# Filter plants with endozoochory in PREDICTS: ------------------------------------------------------

# Import TRY dataset
TRYdata <- fread("./data/TRY/9212.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T) 

# Check taxonomy 
# Determine the records that are identified at least at the species level

Plantspecies <- TRYdata %>% 
  
  # get only the records that have information for the trait
  subset(TraitID != "NA") %>% 
  
  # select only the columns with important information
  select(Dataset, SpeciesName, AccSpeciesName, TraitID, TraitName,
         DataName, OriglName, OrigValueStr, Reference, Comment) %>% 
  
  # drop levels
  droplevels() %>% 
  
  # create a column to store the information of level of classification
  mutate(Nwords = NA)


# Count the words of the species names to determine the level of classification
for (i in 1:nrow(Plantspecies)){
  
  # print the index to see the progress of the loop
  print(i)
  
  # count the words in the consolidated species name
  Plantspecies$Nwords[i] <- wordcount(Plantspecies$AccSpeciesName[i]) 
}

# Select the species that are identified at the species or subspecies level

Plantspecies <- Plantspecies  %>% 
  
  # Filter the rows that are classified at species, subspecies or variety level
  filter(Nwords == 2 | Nwords == 4) %>%
  
  #drop levels
  droplevels()

# Export table
write.csv(Plantspecies, "./output/intermediate_files/01_Filter_data_Plant_species_and_subspecies.csv")

# Import table
Plantspecies <- read.csv("./output/intermediate_files/01_Filter_data_Plant_species_and_subspecies.csv")

# Check taxonomy using the Catalogue of life

# Make a vector with the species names
PlantNames <- Plantspecies %>%
  
  # Get unique names
  distinct(AccSpeciesName) %>%
  
  # Make a character vector with the names of the species
  pull()


# Divide the vector because the gnr_resolve function does not work with large character vectors
a <- PlantNames[1:5500] 
b <- PlantNames[5501:11000]
c <- PlantNames[11001:15983]

# Get the data source id for Catalogue of life
sources <- gnr_datasources() 

# Get the matched name in the Catalogue of life 
aUpdated <- gnr_resolve(names= a, data_source_ids=1)   
bUpdated <- gnr_resolve(names= b, data_source_ids=1) 
cUpdated <- gnr_resolve(names= c, data_source_ids=1) 

# Export results
Updated_plant_names <- rbind(aUpdated, bUpdated, cUpdated)
write.csv(Updated_plant_names, "./output/intermediate_files/01_Filter_data_Plant_submitted_and_updated_name.csv")

# Import results
Updated_plant_names <- read.csv("./output/intermediate_files/01_Filter_data_Plant_submitted_and_updated_name.csv")


# Raise subspecies and varieties to species level
Updated_plant_names <- Updated_plant_names %>%  
  
  # Get the updated names and raise the subspecies and varieties to species by selecting the 
  # first two words of the string (this will include synonyms)
  
  mutate(Updatedname = word(matched_name, 1, 2)) %>% 
  
  # Select the original name and the updated name in the Catalogue of life
  select(user_supplied_name, Updatedname) %>%
  
  # Rename the column to match the original TRY dataset
  rename("AccSpeciesName" = user_supplied_name)


# Find out what species names (and what synonyms) are present in the PREDICTS data.

Updated_plant_names <- Updated_plant_names %>% 
  
  # Subset the species that are present in both datasets
  subset(Updatedname %in% TropicsDiversity$Best_guess_binomial) %>%
  
  # drop levels
  droplevels()

# How many species are present in both datasets
length(unique(Updated_plant_names$Updatedname))

# Create a table with the TRY records that belong to species, subspecies and varieties
# that are present (as species) in the PREDICTS database
Plantspecies <- merge(x=Plantspecies, y=Updated_plant_names, by = "AccSpeciesName", all.x = FALSE)


# Search for plants with endozoochory 

# Divide the datasets in the traits 28= "Dispersal syndrome", 99= "Fruit type"

# Get the table for dispersal syndrome
Dispersal_syndrome <- Plantspecies %>% 
  subset(TraitID == 28) %>% 
  droplevels()

# Export results
write.csv(Dispersal_syndrome, "output/intermediate_files/01_Filter_data_All_dispersal_syndrome_records.csv")

# Get the table for fruit type
Fruit_type <- Plantspecies %>% subset(TraitID == 99) %>% droplevels()

# Export results 
write.csv(Fruit_type, "output/intermediate_files/01_Filter_data_All_fruit_type_records.csv")

# See the factors of dispersal syndrome type
DS_levels <- as.data.frame(levels(Dispersal_syndrome$OrigValueStr))

# Remove all species that are explicitly classified as non-dispersed by endozoochory 

Not_endozoo <- c("a-wind", "Anemo", "anemochory", "Anemochory", "ant", "ant-dispersed",
                 "Ants", "autochor", "Autochory", "ballochor", "boleochor", "Cichlids", "Pacu",
                 "wind+water", "commerce", "wind", "dysochor", "epizoochor", "ethelochor", 
                 "exozoochory", "explosive", "Explosive", "Fish", "Floating", "generative",
                 "germinule", "machinery", "hemerochor", "man", "meteorochor", 
                 "nautochor", "ombrochor", "seed contamination", "ground", "	speirochor", "snail", 
                 "turtles", "Wind", "water", "Water", "agochor", "Anchorage", "terrapin",
                 "Pheidole", "car", "carried", "clothes", "environmental", "erosion", "Geochelone",
                 "rainwash", "reptiles", "roe", "speirochor", "unassisted", "Unassisted", 
                 "unspecialised")

Dispersal_syndrome <- Dispersal_syndrome %>% 
  
  # Filter out all of the records with values in the vector Not_endozoo
  filter(!str_detect(OrigValueStr, paste(Not_endozoo, collapse = "|"))) %>%
  
  # drop levels
  droplevels()

# Filter species that are explicitly classified as dispersed by endozoochory:

# Since the TRY dataset is composed by multiple datasets, I analyze each one separately
# because the terms used are different. 

# The Miombo dataset uses a 1 to specify species that have endozoochory
Endoo_PlantSpecies <- Dispersal_syndrome %>%
  
  # Select specific dataset
  subset(Dataset == "Miombo tree species - SLA, leaf and seed size") %>%
  
  # Filter values with condition TRUE for endozoochory
  filter(OriglName == "dispersal mode: endozoochory" & OrigValueStr == 1) %>%
  
  # Make a character vector
  droplevels() %>% pull(Updatedname) %>% as.character()

# Use the function Filter_endozoochory_terms to filter the rest of databases (to see the function, open document
# ./R/01_Filter_data_functions.R)

Endoo_PlantSpecies <- Filter_endozoochory_terms(DatasetName = "The LEDA Traitbase", FilterValue = "endozoochor")
Endoo_PlantSpecies <- Filter_endozoochory_terms(DatasetName = "BASECO: a floristic and ecological database of Mediterranean French flora",
                             FilterValue = "Endo-zoochory")
Endoo_PlantSpecies <- Filter_endozoochory_terms(DatasetName = "Global Seed Mass, Plant Height Database", FilterValue = "eaten")
Endoo_PlantSpecies <- Filter_endozoochory_terms(DatasetName = "KEW Seed Information Database (SID)", FilterValue = "eaten")


# Determine what species have an animal as the dispersion vector, but the 
# dataset does not specify whether it is endozoochory or epizoochory
Revise_species <-   Dispersal_syndrome %>% 
  
  # Subset all the species that are not present in the Endozoochory list
  subset(Updatedname %nin% Endoo_PlantSpecies) %>% 
  
  # Get unique species names
  distinct(Updatedname, .keep_all = TRUE) %>%
  
  # drop levels
  droplevels()

# Identify the species in the Revise list that have fleshy fruits
# and therefore may be dispersed by endozoochory. 

# See the factors of fruit type
FT_levels <- as.data.frame(levels(Fruit_type$OrigValueStr))

# Create a character vector with the factors of fleshy fruits
Fleshy <- c("berry", "Berry", "drupe", "Drupe")

# Identify the species that have fleshy fruits
Fruit_type <- Fruit_type %>% 
  
  # Filter all of the records with fleshy fruits
  filter(str_detect(OrigValueStr, paste(Fleshy, collapse = "|"))) %>%
  
  # Select only the name of the species and the fruit type
  select(Updatedname, OrigValueStr) %>%
  
  # Get the unique species names
  distinct(Updatedname, .keep_all = TRUE) %>% 
  
  # drop levels
  droplevels()

# Add to the list of Endozoochory plants the species that were already identified as being dispersed by animals
# and also have fleshy fruits. 
# I'm also going to add those species without a dispersal syndrome value, but with fleshy fruits 
Endoo_PlantSpecies <- c(Endoo_PlantSpecies, as.character(merge(x = Revise_species, y = Fruit_type, by= "Updatedname", all = TRUE) %>% # Merge the
                                                           # datasets of plants dispersed by animals, and plants with fleshy fruits 
                                                           
                                                           # Filter out all of the species that don't have information of fruit type or don't 
                                                           # have fleshy fruits
                                                           subset(OrigValueStr.y != "NA") %>% 
                                                           
                                                           # drop levels
                                                           droplevels() %>% 
                                                           
                                                           # Get unique species names
                                                           distinct(Updatedname) %>%
                                                           
                                                           # Make a vector
                                                           droplevels() %>% pull(Updatedname)))

# Eliminate repeated records between the Fruit type and Dispersal syndrome datasets
Endoo_PlantSpecies <- unique(Endoo_PlantSpecies)

# Update list to revise
Revise_species <-   Revise_species %>% 
  
  # Subset all the species that are not present in the Endozoochory list
  subset(Updatedname %nin% Endoo_PlantSpecies) %>% 
  
  # drop levels
  droplevels()

# Export table
write.csv(Revise_species, "./output/intermediate_files/01_Filter_data_Plant_species_to_revise.csv")

# Load species checked manually
Revised <- read.csv("./output/intermediate_files/01_Filter_data_Plant_species_revised.csv", sep = ";")

# Get the final list of endozoochorous plant species 
Endoo_PlantSpecies <- c(Endoo_PlantSpecies, as.character(Revised  %>% 
                                                           pull(Updatedname)))

# Export list of species
write.csv(Endoo_PlantSpecies, "./output/intermediate_files/01_Filter_data_List_endozoochory_species.csv")

# Import list of species
Endoo_PlantSpecies <- read.csv("./output/intermediate_files/01_Filter_data_List_endozoochory_species.csv", header = TRUE) %>%
  
  # rename column to match the PREDICTS dataset
  rename(Best_guess_binomial = x) %>%
  
  # get unique species names
  distinct(Best_guess_binomial)

# Get the PREDICTS records that belong to endozoochoric plants
PREDICTSendooPlants <- merge(x=TropicsDiversity, y= Endoo_PlantSpecies, by= "Best_guess_binomial", all.x = FALSE)

# Export table
write.csv(PREDICTSendooPlants, "./output/cleaned_data/01_Filter_data_PREDICTS_EndozooPlants_records.csv")

# Import table
PREDICTSendooPlants <- read.csv("./output/cleaned_data/01_Filter_data_PREDICTS_EndozooPlants_records.csv")

# Get the number of species in the new dataset
length(unique(PREDICTSendooPlants$Best_guess_binomial))

# Number of records belonging to different orders
table(PREDICTSendooPlants$Order)

# Number of records per region
table(PREDICTSendooPlants$Region)

# Total number of sourceIDs
length(unique(PREDICTSendooPlants$Source_ID))

# Total number of studies
length(unique(PREDICTSendooPlants$SS))

# Total number of sites
length(unique(PREDICTSendooPlants$SSBS))

# Get the first record for each SSBS (SSBS = Concatenation of Source_ID, Study_number, Block and Site_number)
justFirstsEndoPlants <- Number_sites_landuses(PREDICTSendooPlants)

# Number of sites per land-use type and land-use intensity
addmargins(table(justFirstsEndoPlants$Predominant_habitat,justFirstsEndoPlants$Use_intensity), 2)

# Export table of first matches
write.csv(justFirstsEndoPlants, "./output/intermediate_files/01_Filter_data_First_matches_endoPlants.csv")

# Merge PREDICTS frugivores data frame with the PREDICTS endozoochorous plants  -----------------------------------------

# Eliminate columns that don't match between the two data frames
PREDICTSendooPlants <- PREDICTSendooPlants[,c(-1,-3)]
PREDICTS_frugivores <- PREDICTS_frugivores[,c(-1 ,-92, -91)]

# Bind the rows of frugivores to the PREDICTS endozoochorous plants data frame
PREDICTS_frugivores_and_plants <- rbind(PREDICTSendooPlants, PREDICTS_frugivores)

# # Add class
# There are some rows that do not have information on the Kingdom or class, since class is going to
# be used as an interaction term in the model, I am going to complete that information. 

Noclass_kingdom  <- which(PREDICTS_frugivores_and_plants$Class == "" | PREDICTS_frugivores_and_plants$Kingdom == "")
List_sp <- as.data.frame(unique(PREDICTS_frugivores_and_plants[Noclass_kingdom, 1]))
PREDICTS_frugivores_and_plants[c(15717), c(1, 75, 76, 77, 78)]


# This loop searches for the index of species that do not have Kingdom or Class details, and 
# adds the information in the corresponding rows
for (i in Noclass_kingdom) {
  
  # Manually add the species that lack information and their kingdom and class
  if (PREDICTS_frugivores_and_plants[i, 1] %in% List_sp[,1]) { 
    PREDICTS_frugivores_and_plants[i, 75] <- "Plantae" 
    PREDICTS_frugivores_and_plants[i, 77] <- "Magnoliopsida"
  }
}


# Export table 
saveRDS(PREDICTS_frugivores_and_plants , "./output/cleaned_data/01_Filter_data_PREDICTS_Frugivores_and_Endoplants.rds")

# Import table
PREDICTS_frugivores_and_plants <- readRDS("./output/cleaned_data/01_Filter_data_PREDICTS_Frugivores_and_Endoplants.rds")


# Total number of sourceIDs
length(unique(PREDICTS_frugivores_and_plants$Source_ID))

# Total number of studies
length(unique(PREDICTS_frugivores_and_plants$SS))

# Total number of sites
length(unique(PREDICTS_frugivores_and_plants$SSBS))

# Recalculate number of sites

# Get the first record for each SSBS (SSBS = Concatenation of Source_ID, Study_number, Block and Site_number)
justFirsts_Frugivores_and_plants <- Number_sites_landuses(PREDICTS_frugivores_and_plants)

# Number of sites per land-use type and land-use intensity
addmargins(table(justFirsts_Frugivores_and_plants$Predominant_habitat,justFirsts_Frugivores_and_plants$Use_intensity), 2)

# Export table of first matches
write.csv(justFirsts_Frugivores_and_plants, "./output/intermediate_files/01_Filter_data_First_matches_Frugivores_and_endoPlants.csv")


