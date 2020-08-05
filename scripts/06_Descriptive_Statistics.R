

diversity <- readRDS(file = "./output/cleaned_data/01_Filter_data_frugi_endooPlants_notEndooPlants_records.rds")

# ---- 2. Correct abundance measures using Sampling effort ---------------------------------

# The CorrectSamplingEffort function groups the dataset by SS and :
# 1. Finds the maximum value of Sampling effort
# 2. Divides every sampling effort value by the maximum value within each study. 
# 3. Then it takes the abundance measure and divides it by the rescaled sampling effort. 

diversity1 <- yarg::CorrectSamplingEffort(diversity) 


# Descriptive statistics: 
# Centrality: medians are better for distributions with outliers
# Spread: variance or range

# Number of sources
length(unique(diversity1$Source_ID))
# Number of studies
length(unique(diversity1$SS))
# Number of sites
length(unique(diversity1$SSBS))
# Number of countries 
length(unique(diversity1$Country))


# Number of species
frugi  <- diversity1 %>% subset(Kingdom == "Animalia") %>%
  select(Best_guess_binomial, Kingdom, Order) %>% 
  distinct(Best_guess_binomial, .keep_all = TRUE) %>% 
  droplevels()

table(frugi$Order)

endoplants <- diversity1 %>% subset(Kingdom == "Plantae") %>%
  select(Best_guess_binomial, Kingdom, Order) %>% 
  distinct(Best_guess_binomial, .keep_all = TRUE) %>% 
  droplevels()

abioplants <- diversity1 %>% subset(Kingdom == "nePlantae") %>%
  select(Best_guess_binomial, Kingdom, Order) %>% 
  distinct(Best_guess_binomial, .keep_all = TRUE) %>% 
  droplevels()
  
# Year
# Tranform the date into a character
diversity1$Sample_end_latest<- as.character(diversity1$Sample_end_latest)
# Select only the year
diversity1$year <- trimws(sapply(strsplit(diversity1$Sample_end_latest, "-"), `[[`, 1))
# Calculate median year
median(diversity1$year)
range(diversity1$year)

# How many studies and sites in average in each paper

# Mean number of species assessed in each paper

# Abundance
# Load abundance data 
diversity_all_abundance <- readRDS(file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")
median(diversity_all_abundance$Total_abundance, na.rm= TRUE)
range(diversity_all_abundance$Total_abundance, na.rm= TRUE)
length(which(diversity_all_abundance$Kingdom == "Animalia"))
length(which(diversity_all_abundance$Kingdom == "Plantae"))
length(which(diversity_all_abundance$Kingdom == "nePlantae"))

# Species richness
diversity_all_richness <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")
hist(diversity_all_richness$Species_richness)
median(diversity_all_richness$Species_richness, na.rm= TRUE)
range(diversity_all_richness$Species_richness, na.rm= TRUE)

# Simpsons diversity
diversity_all_simpson <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")
hist(diversity_all_simpson$Simpson_diversity)
median(diversity_all_simpson$Simpson_diversity, na.rm= TRUE)
range(diversity_all_simpson$Simpson_diversity, na.rm= TRUE)

