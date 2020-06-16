# Clear workspace
rm(list = ls())

# Load libraries
library(yarg) # useful functions for dealing with PREDICTS data
library(roquefort) # useful PREDICTS functions, particularly for plotting
library(raster) # for dealing with spatial data
library(dplyr) # for handy data manipulation functions
library(tidyr) # ditto
library(lme4) # for mixed effects models
library(car) # for getting anova tables with significance values
library(DHARMa) # for model criticism plots
library(MuMIn) # for checking explanatory power of mixed effects models
library(stringr) # to replace text
library(lmerTest) # to get p-values for estimates in LMMs

# --- Description -----------------------------------------------------------------

# Model Within-sample species richness as the number of total species sampled in eahc site. 


# ---- 1. Load data ------------------------------------------------------------------

diversitS <- readRDS("./output/cleaned_data/01_Filter_data_PREDICTS_Filtered_table.rds")

# ---- 2. Select studies with more than one species ---------------------------------

# According to Phillips et al., (2017) "2.	Studies where the sampling focused on a 
# single species or a predetermined list of species (rather than recording any species
# within the focal taxonomic or ecological group that was sampled) were removed to 
# avoid biasing species-richness estimates". 

# So, for now, I'm just going to remove studies that focused on a single species

# Filter the data needed to run the model
diversity1S <- diversity  %>% 

  # By grouping by SS we are comparing sites with the same diversity metric type
  # either abundance or occurence 
  group_by(SS) %>%
  
  # Create a new column to calculate the number of species sampled per study 
  # n_distinct = length(unique())
   mutate(N_species_sampled = n_distinct(Taxon_name_entered)) %>%
  
  # ungroup data frame  
  ungroup() %>%
  
  # Filter the studies that sampled more than 1 species
  filter(N_species_sampled > 1) %>%
  
  droplevels()

# Minimum number of species sampled
min(diversity1$N_species_sampled)

# Maximum number of species sampled
max(diversity1$N_species_sampled)

# Remove N_species_sampled_column
diversity1S <- diversity1S %>% select(-N_species_sampled)

# ---- 3. Merge Sites -------------------------------------------------------------

# We merge sites within studies that have identical coordinates, start and end dates(and land-use 
# type and intensity). We do this because sometimes authors record different points in a transect as
# different sites, which might not be meaningful if they share land-use, intensity and coordinates (the
# coordinates have an error) OR maybe it is, but we are having a conservative approach

diversity2S <- yarg::MergeSites(diversity1S,
                               silent = TRUE,
                               merge.extra = "Wilderness_area")

# ----4. Rename Predominant habitat --------------------------------------------------

# Rename the column predominant habitat, as the dataset is actually refering to land use
diversity2S <- rename(diversity2S,
                     Predominant_land_use = Predominant_habitat)

# ----5. Remove sites that will produce NaN values --------------------------------------------

diversity2S <- diversity2S %>%  subset(SSBS %nin% c("MJ1_2009__Lehouck 2 Fururu 10", 
                                                    "MJ1_2009__Lehouck 2 Fururu 11", 
                                                    "MJ1_2009__Lehouck 2 Macha 12", 
                                                    "MJ1_2009__Lehouck 2 Macha 13")) %>%
  droplevels()

# ----5.  Calculate diversity metrics -----------------------------------------------

diversity3S <- diversity2S %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Kingdom")) 


# Export table
saveRDS(diversity3S, file = "./output/cleaned_data/02_Statistical_Analysis_Site_metrics_animals_endoo_richness.rds")

# Import table
diversity4 <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Site_metrics.rds")






