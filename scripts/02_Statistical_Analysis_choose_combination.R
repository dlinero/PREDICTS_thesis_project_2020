# Load libraries
library(dplyr)
library(lme4)

# Try all the combinations of possible land-use types and land-use categories

# In order to know which levels to use in all of the models:
# I will the run the abundance, species richness and Simpson´s diversity models with 
# the same land-use/intensities levels, add their AICs and compare which option is better.

# Call the function that merges lan-uses and intensities
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R")

# Load richness data
diversity_all_richness <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Load abundance data 
diversity_all_abundance <- readRDS(file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")
# Transform rescaled abundance
diversity_all_abundance <- mutate(diversity_all_abundance, 
                                  logAbundance = log(RescaledAbundance + 1)
)

# Load Simpsons data 
diversity_all_simpson <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")
# Transform Simpson's diversity index
diversity_all_simpson <- diversity_all_simpson %>%
  
  # Create a new column to transform the response variable with log
  mutate(log_one_over_D = log(Simpson_diversity))


# --- First option ---------------------------------------------------------------------------------------------
# For the first attempt I will:

# Primary can be divided into the three level intensities in all cases
#	Cropland can be divided in minimal use and light/intense use 
# ISV can be divided in minimal use and light/intense use 
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense 
#	YSV has to be merged 

# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate <- c("Primary", "Cropland", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense <- c("Cropland", "ISV", "Plantation forest")


diversity_all_richness <- Merge_landUses_and_intensities(dataset = diversity_all_richness, index = 1, 
                                                land_uses_separate_intensities = land_uses_separate,
                                                land_uses_merge_light_intense = land_uses_light_intense,
                                                "Primary Minimal use")

diversity_all_abundance <- Merge_landUses_and_intensities(dataset = diversity_all_abundance, index = 1, 
                                                          land_uses_separate_intensities = land_uses_separate,
                                                          land_uses_merge_light_intense = land_uses_light_intense,
                                                          "Primary Minimal use")

diversity_all_simpson <- Merge_landUses_and_intensities(dataset = diversity_all_simpson, index = 1, 
                                                          land_uses_separate_intensities = land_uses_separate,
                                                          land_uses_merge_light_intense = land_uses_light_intense,
                                                          "Primary Minimal use")

# 1. Richness

# Get complete cases, in order to compare models that use the same rows
diversity_richness <- drop_na(diversity_all_richness, Species_richness, LandUse.1)  %>% droplevels()

#Check levels
levels(diversity_richness$LandUse.1)


# Model richness with those land-use classes
m1_richness <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_richness, family = poisson, 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))



# Get results of the model
summary(m1_richness)
# AIC
AIC(m1_richness)
#Residual plots
simulationOutput_m1_richness <- simulateResiduals(fittedModel = m1_richness,
                                                  plot = TRUE)


# 2. Abundance
# Get complete cases, in order to compare models that use the same rows
diversity_abundance <- drop_na(diversity_all_abundance, RescaledAbundance, LandUse.1)

#Check levels
levels(diversity_abundance$LandUse.1)

# Model abundance with those land-use classes
m1_abundance <- lmer(logAbundance ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                               (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_abundance)

# Get results of the model
summary(m1_abundance)
# AIC
AIC(m1_abundance)
#Residual plots
simulationOutput_m1_abundance <- simulateResiduals(fittedModel = m1_abundance,
                                                   plot = TRUE)

# 3. Simpson Diversity
diversity_simpson <- drop_na(diversity_all_simpson, Simpson_diversity, LandUse.1)

#Check levels
levels(diversity_simpson$LandUse.1)

# Model simpson's index with those land-use classes
m1_simpson <- lmer(log_one_over_D ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                       (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_simpson)

# Get results of the model
summary(m1_simpson)
# AIC
AIC(m1_simpson)
#Residual plots
simulationOutput_m1_simpson <- simulateResiduals(fittedModel = m1_simpson,
                                                 plot = TRUE)


# --- Second option ---------------------------------------------------------------------------------------------

# For the second attempt I will:

# Primary can be divided into the three level intensities in all cases
#	Cropland has to be merged
# ISV can be divided in minimal use and light/intense use 
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense 
#	YSV has to be merged 


# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate_2 <- c("Primary", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense_2 <- c("ISV", "Plantation forest")


diversity_richness <- Merge_landUses_and_intensities(dataset = diversity_richness, index = 2, 
                                                                 land_uses_separate_intensities = land_uses_separate_2,
                                                                 land_uses_merge_light_intense = land_uses_light_intense_2,
                                                                 "Primary Minimal use")

diversity_abundance <- Merge_landUses_and_intensities(dataset = diversity_abundance, index = 2, 
                                                          land_uses_separate_intensities = land_uses_separate_2,
                                                          land_uses_merge_light_intense = land_uses_light_intense_2,
                                                          "Primary Minimal use")

diversity_simpson <- Merge_landUses_and_intensities(dataset = diversity_simpson, index = 2, 
                                                     land_uses_separate_intensities = land_uses_separate_2,
                                                     land_uses_merge_light_intense = land_uses_light_intense_2,
                                                     "Primary Minimal use")

# 1. Richness

#Check levels
levels(diversity_richness$LandUse.2)


# Model richness with those land-use classes
m2_richness <- glmer(Species_richness ~ LandUse.2 + Kingdom + LandUse.2:Kingdom +
                       (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_richness, family = poisson, 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))



# Get results of the model
summary(m2_richness)
# AIC
AIC(m2_richness)
#Residual plots
simulationOutput_m2_richness <- simulateResiduals(fittedModel = m2_richness, plot = TRUE)


# 2. Abundance

#Check levels
levels(diversity_abundance$LandUse.2)

# Model abundance with those land-use classes
m2_abundance <- lmer(logAbundance ~ LandUse.2 + Kingdom + LandUse.2:Kingdom +
                       (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_abundance)

# Get results of the model
summary(m2_abundance)
# AIC
AIC(m2_abundance)
#Residual plots
simulationOutput_m2_abundance <- simulateResiduals(fittedModel = m2_abundance, plot = TRUE)

# 3. Simpson Diversity

#Check levels
levels(diversity_simpson$LandUse.2)

# Model simpson's index with those land-use classes
m2_simpson <- lmer(log_one_over_D ~ LandUse.2 + Kingdom + LandUse.2:Kingdom +
                     (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_simpson)

# Get results of the model
summary(m2_simpson)
# AIC
AIC(m2_simpson)
#Residual plots
simulationOutput_m2_simpson <- simulateResiduals(fittedModel = m2_simpson, plot = TRUE)

# --- Third option ---------------------------------------------------------------------------------------------

# For the third attempt I will:

# Primary can be divided into the two level intensities in all cases
#	Cropland has to be merged
# ISV can be divided in minimal use and light/intense use 
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense 
#	YSV has to be merged 


# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate_3 <- c("Primary", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense_3 <- c("Primary", "ISV", "Plantation forest")


diversity_richness <- Merge_landUses_and_intensities(dataset = diversity_richness, index = 3, 
                                                     land_uses_separate_intensities = land_uses_separate_3,
                                                     land_uses_merge_light_intense = land_uses_light_intense_3,
                                                     "Primary Minimal use")

diversity_abundance <- Merge_landUses_and_intensities(dataset = diversity_abundance, index = 3, 
                                                      land_uses_separate_intensities = land_uses_separate_3,
                                                      land_uses_merge_light_intense = land_uses_light_intense_3,
                                                      "Primary Minimal use")

diversity_simpson <- Merge_landUses_and_intensities(dataset = diversity_simpson, index = 3, 
                                                    land_uses_separate_intensities = land_uses_separate_3,
                                                    land_uses_merge_light_intense = land_uses_light_intense_3,
                                                    "Primary Minimal use")

# 1. Richness

#Check levels
levels(diversity_richness$LandUse.3)


# Model richness with those land-use classes
m3_richness <- glmer(Species_richness ~ LandUse.3 + Kingdom + LandUse.3:Kingdom +
                       (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_richness, family = poisson, 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))



# Get results of the model
summary(m3_richness)
# AIC
AIC(m3_richness)
#Residual plots
simulationOutput_m3_richness <- simulateResiduals(fittedModel = m3_richness, plot = TRUE)


# 2. Abundance

#Check levels
levels(diversity_abundance$LandUse.3)

# Model abundance with those land-use classes
m3_abundance <- lmer(logAbundance ~ LandUse.3 + Kingdom + LandUse.3:Kingdom +
                       (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_abundance)

# Get results of the model
summary(m3_abundance)
# AIC
AIC(m3_abundance)
#Residual plots
simulationOutput_m3_abundance <- simulateResiduals(fittedModel = m3_abundance, plot = TRUE)

# 3. Simpson Diversity

#Check levels
levels(diversity_simpson$LandUse.3)

# Model richness with those land-use classes
m3_simpson <- lmer(log_one_over_D ~ LandUse.3 + Kingdom + LandUse.3:Kingdom +
                     (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_simpson)

# Get results of the model
summary(m3_simpson)
# AIC
AIC(m3_simpson)
#Residual plots
simulationOutput_m3_simpson <- simulateResiduals(fittedModel = m3_simpson, plot = TRUE)

# --- Fourth option ---------------------------------------------------------------------------------------------

# For the fourth attempt I will:

# Primary can be divided into the two level intensities in all cases
#	Cropland can be divided into the two level intensities in all cases
# ISV can be divided in minimal use and light/intense use 
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense 
#	YSV has to be merged 


# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate_4 <- c("Primary","Cropland", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense_4 <- c("Primary", "Cropland", "ISV", "Plantation forest")


diversity_richness <- Merge_landUses_and_intensities(dataset = diversity_richness, index = 4, 
                                                     land_uses_separate_intensities = land_uses_separate_4,
                                                     land_uses_merge_light_intense = land_uses_light_intense_4,
                                                     "Primary Minimal use")

diversity_abundance <- Merge_landUses_and_intensities(dataset = diversity_abundance, index = 4, 
                                                      land_uses_separate_intensities = land_uses_separate_4,
                                                      land_uses_merge_light_intense = land_uses_light_intense_4,
                                                      "Primary Minimal use")

diversity_simpson <- Merge_landUses_and_intensities(dataset = diversity_simpson, index = 4, 
                                                    land_uses_separate_intensities = land_uses_separate_4,
                                                    land_uses_merge_light_intense = land_uses_light_intense_4,
                                                    "Primary Minimal use")

# 1. Richness

#Check levels
levels(diversity_richness$LandUse.4)


# Model richness with those land-use classes
m4_richness <- glmer(Species_richness ~ LandUse.4 + Kingdom + LandUse.4:Kingdom +
                       (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_richness, family = poisson, 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))



# Get results of the model
summary(m4_richness)
# AIC
AIC(m4_richness)
#Residual plots
simulationOutput_m4_richness <- simulateResiduals(fittedModel = m4_richness, plot = TRUE)


# 2. Abundance

#Check levels
levels(diversity_abundance$LandUse.4)

# Model abundance with those land-use classes
m4_abundance <- lmer(logAbundance ~ LandUse.4 + Kingdom + LandUse.4:Kingdom +
                       (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_abundance)

# Get results of the model
summary(m4_abundance)
# AIC
AIC(m4_abundance)
#Residual plots
simulationOutput_m4_abundance <- simulateResiduals(fittedModel = m4_abundance, plot = TRUE)

# 3. Simpson Diversity

#Check levels
levels(diversity_simpson$LandUse.4)

# Model simpson's with those land-use classes
m4_simpson <- lmer(log_one_over_D ~ LandUse.4 + Kingdom + LandUse.4:Kingdom +
                     (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_simpson)

# Get results of the model
summary(m4_simpson)
# AIC
AIC(m4_simpson)
#Residual plots
simulationOutput_m4_simpson <- simulateResiduals(fittedModel = m_simpson, plot = TRUE)

# --- Choose bets combination ---------------------------------------------------------------------------------------------

# First option
sum(AIC(m1_richness) + AIC(m1_abundance), AIC(m1_simpson))

# Second option
sum(AIC(m2_richness) + AIC(m2_abundance), AIC(m2_simpson))

# Third option
sum(AIC(m3_richness) + AIC(m3_abundance), AIC(m3_simpson))

# Fourth option
sum(AIC(m4_richness) + AIC(m4_abundance), AIC(m4_simpson))
