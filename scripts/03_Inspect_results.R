# Clear workspace
rm(list = ls())

# Load libraries

library(dplyr) # for handy data manipulation functions
library(tidyr) # ditto
library(lme4) # for mixed effects models
library(ggplot2) # for making graphs
library(influence.ME) # for detecting influential cases

### 1. Load data -------------------------------------------------------------------------------------------

# 1.1 Load site metrics data -------------------------------------------------------
# species richness 
diversity_all_richness <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")

# abundance 
diversity_all_abundance <- readRDS(file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")

# Transform rescaled abundance
diversity_all_abundance <- mutate(diversity_all_abundance, 
                                  logAbundance = log(RescaledAbundance + 1)
)


# Simpsons 
diversity_all_simpson <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")


# 1.2. Load data that contains sp names  -------------------------------------------------------
# I am going to upload the tables that have been:
# * Corrected for sampling effort
# * Filtered for more than 1 sp per study (in the case of sp. richness and simpsons diversity)
# * Removed sites with NaN values (only for simpsons and richness because the othe abundance
# one I saved it before doing that step)
# * Merge sites
# * Rename predominant habitat

# species richness 
sp_richness <- readRDS(file = "./output/intermediate_files/02_Statistical_Analysis_richness_species_all.rds")

# Abundance
sp_abundance <- readRDS(file = "./output/intermediate_files/02_Statistical_Analysis_Abundance_species_all.rds")

# Simpsons diversity
sp_simpson <- readRDS(file = "./output/intermediate_files/02_Statistical_Analysis_Simpson_species_all.rds")


### 2. Get the SSBS that we want to inspect -------------------------------------------------------------------

# 2.1. Drop sites that don't have the selected LUI (Land-use-intensity) --------------------------------------------------------------

# Drop the sites that we don't analyze with our selected land-use type/intensities 
# combination

# Call the function that merges lan-uses and intensities
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R")

# Create the vectors that hold the land-uses that we want to keep with 
# different use intensities 
land_uses_separate <- c("Primary","Cropland", "ISV", "Plantation forest")

# Create a vector with the land-uses where we want to merge the light 
# and intense use intensities
land_uses_light_intense <- c("Primary", "Cropland", "ISV", "Plantation forest")

# Merge LUI for species richness
diversity_all_richness <- Merge_landUses_and_intensities(dataset = diversity_all_richness, 
                                                     index = 0, 
                                                     land_uses_separate_intensities = land_uses_separate,
                                                     land_uses_merge_light_intense = land_uses_light_intense,
                                                     "Primary Minimal use")
# Drop NA values
diversity_richness <- drop_na(diversity_all_richness, 
                              Species_richness, LandUse.0)  %>%
  droplevels()

# Merge LUI for abundance
diversity_all_abundance <- Merge_landUses_and_intensities(dataset = diversity_all_abundance, 
                                                      index = 0, 
                                                      land_uses_separate_intensities = land_uses_separate,
                                                      land_uses_merge_light_intense = land_uses_light_intense,
                                                      "Primary Minimal use")

# Drop NA values
diversity_abundance <- drop_na(diversity_all_abundance, 
                               RescaledAbundance, LandUse.0)  %>%
  droplevels()

# Merge LUI for Simpson diversity
diversity_all_simpson  <- Merge_landUses_and_intensities(dataset = diversity_all_simpson,
                                                    index = 0, 
                                                    land_uses_separate_intensities = land_uses_separate,
                                                    land_uses_merge_light_intense = land_uses_light_intense,
                                                    "Primary Minimal use")
# Drop NA values
diversity_simpson <- drop_na(diversity_all_simpson, 
                               Simpson_diversity, LandUse.0)  %>%
  droplevels()


# 2.2. Check pasture SSBS-----------------------------------------------------------

# I am going to check why plants have greater abundance and no decrease in species richness
# and simpsons diversity for pasture compared to primary minimal use

# 2.2.1. Species richness -----------------------------------------------------------
# Get the SSBS of sites belonging to pasture and primary minimal use for Plants 

richness <- diversity_richness %>% 
  
  # Filter LUI
  dplyr::filter(LandUse.0 %in% c("Primary Minimal use", "Pasture All")) %>%
  
  # Filter plants
  dplyr::filter(Kingdom %in% c("Plantae")) %>%
  
  # drop levels
  droplevels() 

# Make a boxplot comparing LUI (remember glmm compares sites within studies)
ggplot(richness, aes(x=LandUse.0, y=Species_richness)) + 
  geom_boxplot() + ylim(0,20)

# Get some statistics of species richness for both land uses
richness %>%
  
  # Group by LUI
  group_by(LandUse.0) %>% 
  
  # Calculate statistics
  summarise(mean_richness = mean(Species_richness), 
                                               sd = sd(Species_richness),
                                               qt_1 = quantile(Species_richness, prob=c(0.25)),
                                               qt_2 = quantile(Species_richness, prob=c(0.5)),
                                               qt_3 = quantile(Species_richness, prob=c(0.75)),
                                               qt_4 = quantile(Species_richness, prob=c(1)))

# I am going to check the sites that have a species richness greater than the mean/median 
# for pasture

pasture_richness <- richness %>% 
  
  # Filter pasture
  filter(LandUse.0 == "Pasture All") %>%
  
  # Filter sites with species richness greater than the mean/median
  filter(Species_richness > 2) %>%
  
  # droplevels
  droplevels() %>%
  
  pull(SSBS) %>% as.character()

# 2.2.2. Abundance -----------------------------------------------------------
# Get the SSBS of sites belonging to pasture and primary minimal use for Plants 

abundance <- diversity_abundance %>% 
  
  # Filter LUI
  dplyr::filter(LandUse.0 %in% c("Primary Minimal use", "Pasture All")) %>%
  
  # Filter plants
  dplyr::filter(Kingdom %in% c("Plantae")) %>%
  
  # drop levels
  droplevels() 

# Make a boxplot comparing LUI (remember glmm compares sites within studies)
ggplot(abundance, aes(x=LandUse.0, y=RescaledAbundance)) + 
  geom_boxplot()

# Get some statistics of abundance for both land uses
abundance %>%
  
  # Group by LUI
  group_by(LandUse.0) %>% 
  
  # Calculate statistics
  summarise(mean = mean(RescaledAbundance), 
            sd = sd(RescaledAbundance),
            qt_1 = quantile(RescaledAbundance, prob=c(0.25)),
            qt_2 = quantile(RescaledAbundance, prob=c(0.5)),
            qt_3 = quantile(RescaledAbundance, prob=c(0.75)),
            qt_4 = quantile(RescaledAbundance, prob=c(1)))

# I am going to check the sites that have a rescaled abundance greater than the median 
# for pasture

pasture_abundance <- abundance %>% 
  
  # Filter pasture
  filter(LandUse.0 == "Pasture All") %>%
  
  # Filter sites with species richness greater than the median
  filter(RescaledAbundance > 0.07) %>%
  
  # droplevels
  droplevels() %>%
  
  pull(SSBS) %>% as.character()

# 2.2.3. Simpsons -----------------------------------------------------------
# Get the SSBS of sites belonging to pasture and primary minimal use for Plants 

simpson <- diversity_simpson %>% 
  
  # Filter LUI
  dplyr::filter(LandUse.0 %in% c("Primary Minimal use", "Pasture All")) %>%
  
  # Filter plants
  dplyr::filter(Kingdom %in% c("Plantae")) %>%
  
  # drop levels
  droplevels() 

# Make a boxplot comparing LUI (remember glmm compares sites within studies)
ggplot(simpson, aes(x=LandUse.0, y=Simpson_diversity)) + 
  geom_boxplot() + ylim(0, 5)

# Get some statistics of simpsons diversity for both land uses
simpson %>%
  
  # Group by LUI
  group_by(LandUse.0) %>% 
  
  # Calculate statistics
  summarise(mean = mean(Simpson_diversity), 
            sd = sd(Simpson_diversity),
            qt_1 = quantile(Simpson_diversity, prob=c(0.25)),
            qt_2 = quantile(Simpson_diversity, prob=c(0.5)),
            qt_3 = quantile(Simpson_diversity, prob=c(0.75)),
            qt_4 = quantile(Simpson_diversity, prob=c(1)))

# I am going to check the sites that have a simpson's index greater than the median 
# for pasture

pasture_simpson <- simpson %>% 
  
  # Filter pasture
  filter(LandUse.0 == "Pasture All") %>%
  
  # Filter sites with species richness greater than the median
  filter(Simpson_diversity > 1.35) %>%
  
  # droplevels
  droplevels() %>%
  
  pull(SSBS) %>% as.character()


# See how many sites are shared among the metrics 
length(which(pasture_richness %in% pasture_abundance))
length(which(pasture_simpson %in% pasture_abundance))
length(which(pasture_simpson %in% pasture_richness))


### 3. Inspect the sites for pasture ---------------------------------------------

# 3.1. Richness ---------------------------------------------------------------------
pasture_richness_inspect <- sp_richness %>% 
  
  # subset the sites we want to inspect
  base::subset(SSBS %in% pasture_richness) %>%
  
  # Subset only plants dispersed by endozoochory
  base::subset(Kingdom == "Plantae") %>%
  
  # drop levels
  droplevels() %>%
  
  # select columns we are interested in
  dplyr::select(Country, Source_ID, SS, Diversity_metric_type, Diversity_metric_unit,
        Habitat_as_described, Predominant_land_use, Use_intensity, SSBS, 
         Best_guess_binomial, Taxon_name_entered, Measurement, Kingdom) %>%
  
  # Select species present(1) or with abundance greater than 0
  dplyr::filter(Measurement != 0)

# Resulting number of species
length(unique(pasture_richness_inspect$Best_guess_binomial))

# Resulting number of studies
length(unique(pasture_richness_inspect$SS))

# Export results
write.csv(pasture_richness_inspect, "./output/intermediate_files/03_Inspect_results_Pasture_Species_richness.csv")

# According to what I checked, I've got some species that I doubt are dispersed 
# by endozoochory, like:
#* Juncus effusus (LEDA trait base)
#* Juncus inflexus (LEDA trait base)
#* Plantago major (LEDA trait base)
#* Oxalis corniculata (BASECO database)

# So now, I'm going to directly load these databases (that means I'm not going to use 
# the TRY dataset, but download these databases directly from their respective sources)

# LEDA traitbase
LEDA_traitbase <- read.csv("./data/LEDA traitbase/dispersal_type.csv", header = TRUE, sep = ";")

LEDA_traitbase <- LEDA_traitbase %>% 
  
  # Subset the species of interest
  subset(SBS.name %in% c("Juncus effusus", "Juncus inflexus", "Plantago major"))


#* Juncus effusus (LEDA trait base) = can be abiotically-dispersed (e.g. hay making
# machinery), but also has zoochory (both epi and endoo, with dispersal agents like:
# cattle, deer, hare, rabbits, sheeps and domestic animals).

#* Juncus inflexus (LEDA trait base) = can be abiotically-dispersed, but also has zoochory 
# (both epi and endoo, with dispersal agents like: rabbits)

#* Plantago major (LEDA trait base) = can be abiotically-dispersed, but also has zoochory 
# (both epi and endoo, with dispersal agents like: sheep, domestic animals, roe, horse, goat, deer, cattle)

# BASECO traitbase (I did not find the original table)
# Upload disppersal syndrome traits from the revised TRY dataset
Dispersal_syndrome <- read.csv("output/intermediate_files/01_Filter_data_All_dispersal_syndrome_records.csv")

Baseco_traitbase <- Dispersal_syndrome %>% 
  
  # Subset records belonging to that database
  subset(Dataset == "BASECO: a floristic and ecological database of Mediterranean French flora") %>% 
  
  # Subset the species of interest
  subset(Updatedname %in% c("Oxalis corniculata")) %>% droplevels()

#* Oxalis corniculata (BASECO database) = endozoochory. But according to "Shapes of ballistic seed dispersal distributions:
# a comparison of Oxalis corniculata with a theoretical model" it can be abiotically dispersed
# According to "Seed dispersal by chacma baboons and syntopic ungulates in southern African savannas" it can 
# be consumed and dispersed by baboons, elands, and impalas in savannas. 


# 3.2. Abundance ---------------------------------------------------------------------

pasture_abundance_inspect <- sp_abundance %>% 
  
  # subset the sites we want to inspect
  base::subset(SSBS %in% pasture_abundance) %>%
  
  # Subset only plants dispersed by endozoochory
  base::subset(Kingdom == "Plantae") %>%
  
  # drop levels
  droplevels() %>%
  
  # select columns we are interested in
  dplyr::select(Country, Source_ID, SS, Diversity_metric_type, Diversity_metric_unit,
                Habitat_as_described, Predominant_land_use, Use_intensity, SSBS, 
                Best_guess_binomial, Taxon_name_entered, Measurement, Kingdom) %>%
  
  # Select species present(1) or with abundance greater than 0
  dplyr::filter(Measurement != 0)

# Resulting number of species
length(unique(pasture_abundance_inspect$Best_guess_binomial))

# Resulting number of studies
length(unique(pasture_abundance_inspect$SS))

# Export results
write.csv(pasture_abundance_inspect, "./output/intermediate_files/03_Inspect_results_Pasture_Abundance.csv")

# 3.3. Simpson ---------------------------------------------------------------------

pasture_simpson_inspect <- sp_simpson %>% 
  
  # subset the sites we want to inspect
  base::subset(SSBS %in% pasture_simpson) %>%
  
  # Subset only plants dispersed by endozoochory
  base::subset(Kingdom == "Plantae") %>%
  
  # drop levels
  droplevels() %>%
  
  # select columns we are interested in
  dplyr::select(Country, Source_ID, SS, Diversity_metric_type, Diversity_metric_unit,
                Habitat_as_described, Predominant_land_use, Use_intensity, SSBS, 
                Best_guess_binomial, Taxon_name_entered, Measurement, Kingdom) %>%
  
  # Select species present(1) or with abundance greater than 0
  dplyr::filter(Measurement != 0)

# Resulting number of species
length(unique(pasture_simpson_inspect$Best_guess_binomial))

# Resulting number of studies
length(unique(pasture_simpson_inspect$SS))

# Export results
write.csv(pasture_simpson_inspect, "./output/intermediate_files/03_Inspect_results_Pasture_Simpson.csv")

