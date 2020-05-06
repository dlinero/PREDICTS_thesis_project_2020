library(dplyr)

# Total abundance model ---------------------------------------------------------------------------

# Import table
PREDICTS_frugivores_and_plants <- readRDS("./output/cleaned_data/01_Filter_data_PREDICTS_Frugivores_and_Endoplants.rds")

# Filter the data needed to run the model
Abundances <- PREDICTS_frugivores_and_plants %>% 
  
  # Filter the studies that report abundance measures 
  filter(Diversity_metric_type == "Abundance") %>%
  
  # Create a new column to calculate the number of species sampled per study 
  # n_distinct = length(unique())
  group_by(SS) %>% mutate(N_species_sampled = n_distinct(Taxon_name_entered)) %>%
  
  # ungroup data frame  
  ungroup() %>%
  
  # Filter the studies that sampled more than 1 species
  filter(N_species_sampled > 1) %>%
  
  droplevels()

# Explore the number of resulting sites
Sites_abundances <- Number_sites_landuses(Abundances)
length(Sites_abundances$SSBS)

# Total number of sourceIDs
length(unique(Sites_abundances$Source_ID))

# Total number of studies
length(unique(Sites_abundances$SS))

# Number of sites per land-use type and land-use intensity
addmargins(table(Sites_abundances$Predominant_habitat,Sites_abundances$Use_intensity), 2)

# Number of sites per land-use type and land-use intensity
addmargins(table(Sites_abundances$Predominant_habitat,Sites_abundances$Use_intensity), 2)

# Add Class: apparently 
# Number of sites per land-use type and land-use intensity for animals


# Number of sites per land-use type and land-use intensity for plants


