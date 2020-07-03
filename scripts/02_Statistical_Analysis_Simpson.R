# Clear workspace
rm(list = ls())

# Load libraries
library(yarg) # useful functions for dealing with PREDICTS data
library(roquefort) # useful PREDICTS functions, particularly for plotting
library(dplyr) # for handy data manipulation functions
library(tidyr) # ditto
library(lme4) # for mixed effects models
library(car) # for getting anova tables with significance values
library(DHARMa) # for model criticism plots
library(MuMIn) # for checking explanatory power of mixed effects models
library(stringr) # to replace text
library(lmerTest) # to get p-values for estimates in LMMs


# --- Description -----------------------------------------------------------------

# Model Within-sample Simpson diversity.
# According to the SiteMetrics function the Simpson's Diversity Index is calculated as
# as 1/D, where D = sum((a/A)^2), where a is the abundance of a single taxon, and A is the 
# total abundance of all taxa at a site.

# ---- 1. Load data ------------------------------------------------------------------

# Load the table that contains the abudance and occurrence measurements for frugivores
#  along with the plants dispersed and not dispersed by animals. 

# Import table
diversityS <- readRDS(file = "./output/cleaned_data/01_Filter_data_frugi_endooPlants_notEndooPlants_records.rds")

# ---- 2. Correct abundance measures using Sampling effort ---------------------------------

# The CorrectSamplingEffort function groups the dataset by SS and :
# 1. Finds the maximum value of Sampling effort
# 2. Divides every sampling effort value by the maximum value within each study. 
# 3. Then it takes the abundance measure and divides it by the rescaled sampling effort. 

diversityS <- yarg::CorrectSamplingEffort(diversityS) 

# ---- 3. Select studies with more than one species ---------------------------------

# According to Phillips et al., (2017) "Studies where the sampling focused on a 
# single species or a predetermined list of species (rather than recording any species
# within the focal taxonomic or ecological group that was sampled) were removed to 
# avoid biasing species-richness estimates". 

# So, for now, I'm just going to remove studies that focused on a single species

# Create a dataset with the SS that assessed more than 1 species
list <- diversityS  %>% 
  
  # By grouping by SS we are comparing sites with the same diversity metric type
  # either abundance or occurence 
  dplyr::group_by(SS) %>%
  
  # Create a new column to calculate the number of species sampled per study 
  # n_distinct = length(unique())
  dplyr::mutate(N_species_sampled = dplyr::n_distinct(Taxon_name_entered)) %>%
  
  # ungroup data frame  
  dplyr::ungroup() %>%
  
  # Filter the studies that sampled more than 1 species
  dplyr::filter(N_species_sampled > 1) %>%
  
  base::droplevels()

# Minimum number of species sampled
min(list$N_species_sampled)

# Maximum number of species sampled
max(list$N_species_sampled)

# Create a character vector of unique SS that assessed more than 1 species
list <- as.character(unique(list$SS))

# Filter the studies present in the list
diversityS <- diversityS %>% base::subset(SS %in% list) %>% droplevels()

# ----4. Remove sites that will produce NaN values --------------------------------------------
# Doing the abundance models, I identified some sites where (due to the filtering of specific species)
# recorded only zero abundances. I am going to remove those sites as they will produce 
# Non infinite values. 

diversityS <- diversityS %>%  base::subset(SSBS %nin% c("MJ1_2009__Lehouck 2 Fururu 10", 
                                                        "MJ1_2009__Lehouck 2 Fururu 11", 
                                                        "MJ1_2009__Lehouck 2 Macha 12", 
                                                        "MJ1_2009__Lehouck 2 Macha 13")) %>%
  base::droplevels()

# ---- 5. Split datasets --------------------------------------------------------

# I am going to separate the records that belong to plants not dispersed by animals
diversityS_notEndo <- diversityS %>% base::subset(Kingdom == "nePlantae") %>% base::droplevels()

# Get the table without plants dispersed by animals
diversityS_frugi_endo <- diversityS %>% base::subset(Kingdom != "nePlantae") %>% base::droplevels()

# ---- 6. Merge Sites for animals and plants-------------------------------------------------------------

# We merge sites within studies that have identical coordinates, start and end dates(and land-use 
# type and intensity). We do this because sometimes authors record different points in a transect as
# different sites, which might not be meaningful if they share land-use, intensity and coordinates (the
# coordinates have an error) OR maybe it is, but we are having a conservative approach

diversityS_frugi_endo <- yarg::MergeSites(diversityS_frugi_endo, silent = TRUE, 
                                          merge.extra = "Wilderness_area")

diversityS_notEndo <- yarg::MergeSites(diversityS_notEndo, silent = TRUE, 
                                       merge.extra = "Wilderness_area")

# ----7. Rename Predominant habitat --------------------------------------------------

# Rename the column predominant habitat, as the dataset is actually refering to land use
diversityS_frugi_endo <-  dplyr::rename(diversityS_frugi_endo,
                                        Predominant_land_use = Predominant_habitat)

diversityS_notEndo <-  dplyr::rename(diversityS_notEndo,
                                     Predominant_land_use = Predominant_habitat)


# ----8.  Calculate diversity metrics -----------------------------------------------

# Calculate diversity metrics for animals and endoozoocoric plants
diversity1S_frugi_endo <- diversityS_frugi_endo %>%
  
  # add Diversity_metric_is_valid column
  dplyr::mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Kingdom")) 

# Calculate diversity metrics for plants not dispersed by animals
diversity1S_notendo <- diversityS_notEndo %>%
  
  # add Diversity_metric_is_valid column
  dplyr::mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Kingdom")) 

# Merge the site metrics for all organisms
diversity_all <- base::rbind.data.frame(diversity1S_frugi_endo, diversity1S_notendo)

# ---9. Check the results--------------------------------------------------------------------

# Check results

diversityS_frugi_endo   %>% 
  
  # subset one site to check
  base::subset(SSS == "CM1_2009__Letcher 1 9") %>% base::droplevels() %>%
  
  # Calculate each species abundance over total abundance at the site
  dplyr::mutate(Proportion = Measurement/sum(Measurement)) %>%
  
  # Sum the proportions and squared them
  dplyr::mutate(Simpson_D = sum(Proportion^2)) %>%
  
  # Calculate 1/D
  dplyr::mutate(one_over_D = 1/Simpson_D) %>% 
  
  # Select only columns we're interested in
  dplyr::select(SSS, Diversity_metric_type, Taxon_name_entered, Best_guess_binomial, Measurement, Proportion, Simpson_D, one_over_D)

# Compare with the SiteMetrics result
diversity_all %>% base::subset(SSS == "CM1_2009__Letcher 1 9") %>% dplyr::select(SSS, Simpson_diversity)

# Export table
saveRDS(diversity_all, file = "./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Import table
diversity_all <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")


# ---- 10. Check the type of measure for the abundance records --------------------------------------------------------

# Merge tables of frugivores, endoplants and not endoplants
diversityS_combined <- rbind.data.frame(diversityS_frugi_endo, diversityS_notEndo)

# Get the unique diversity metric type and unit for each study
Diversity_metric_unit <- diversityS_combined %>% 
  
  dplyr::distinct(SS, .keep_all = TRUE) %>%
  
  dplyr::select(Source_ID, SS, Diversity_metric_type, Diversity_metric_unit, Kingdom, 
                Study_common_taxon)

# Filter the measurement table in order to get more information on single papers. 
check <- diversityS_combined %>% subset(Source_ID == "YY1_2016__Mohandass") %>%
  select(Source_ID, Predominant_habitat, Use_intensity, SS, SSS, Diversity_metric_type,
         Best_guess_binomial, Diversity_metric, Diversity_metric_unit,
         Diversity_metric_is_effort_sensitive, Measurement, Habitat_as_described, Sampling_target, Sampling_effort, Sampling_effort_unit,
         Country, Location)

# After checking the type of abundance measures I will continue with the model

############################## MODEL TESTING #################################

# ---- 11. Count number of sites--------------------------------------------------------

# Count number of sites for the SiteMetrics table in order to know the sample size
# for each land-use type and intensity.

# Remove sites that don't have a simpson diversity value or land-use type
diversity_simpson <- drop_na(diversity_all, 
                             Simpson_diversity, Predominant_land_use) %>% droplevels()

# Number of sites
addmargins(table(diversity_simpson$Predominant_land_use, diversity_simpson$Use_intensity, diversity_simpson$Kingdom), 2)


# ---- 12. Try an initial combination--------------------------------------------------------

# I am going to try the combination that allows at least 30 sites in each combination
# except MSV for nePlants that has in total 20 sites 

# Primary can be divided into the two level intensities in all cases
#	Cropland has to be merged
# ISV can be divided in minimal use and light/intense use 
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal (27 animals) use and light/intense 
#	YSV has to be merged 

# Call the function that merges lan-uses and intensities
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R")

# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate_1 <- c("Primary", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense_1 <- c("Primary", "ISV", "Plantation forest")

diversity_simpson <- Merge_landUses_and_intensities(dataset = diversity_simpson, index = 1, 
                                                    land_uses_separate_intensities = land_uses_separate_1,
                                                    land_uses_merge_light_intense = land_uses_light_intense_1,
                                                    "Primary Minimal use")


# --- 12.1 Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(diversity_simpson[ , c("LandUse.1", "Kingdom")])

# ---12.2 Complete cases --------------------------------------------------------------------------------

# Drop sites that don't have index measures or land-use data
diversity_simpson1 <- drop_na(diversity_simpson, 
                             Simpson_diversity, LandUse.1) %>% droplevels()


# Check number of sites
addmargins(table(diversity_simpson1$LandUse.1, diversity_simpson1$Kingdom), 2)

# ---12.3 Choose random effects structure  -----------------------------------------------------------------

# To select the random-effects structure, we use the method recommended by (Zuur et al., 2009) 
# of taking the most complex fixed-effects structure, including all interactions, that will be 
# tested in the second stage of modelling, and use it to compare different random-effects 
# structures.

m1 <- lmer(Simpson_diversity ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
             (1|SS) + (1|SSB), data = diversity_simpson1) #Is Singular

summary(m1)
simulation_simpson1 <- simulateResiduals(fittedModel = m1, plot = TRUE)
plot(m1)
hist(diversity_simpson1$Simpson_diversity)
# The simplest model does not converge

# ---- 13. Try a second combination--------------------------------------------------------

# I am going to merge all categories

# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate_2 <- "NA"
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense_2 <- "NA"

diversity_simpson <- Merge_landUses_and_intensities(dataset = diversity_simpson, index = 2, 
                                                    land_uses_separate_intensities = land_uses_separate_2,
                                                    land_uses_merge_light_intense = land_uses_light_intense_2,
                                                    "Primary All")

# --- 13.1 Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(diversity_simpson[ , c("LandUse.2", "Kingdom")])

# ---13.2 Complete cases --------------------------------------------------------------------------------

# Drop sites that don't have index measures or land-use data
diversity_simpson2 <- drop_na(diversity_simpson, 
                              Simpson_diversity, LandUse.2) %>% droplevels()


# Check number of sites
addmargins(table(diversity_simpson2$LandUse.2, diversity_simpson2$Kingdom), 2)

# ---13.3 Choose random effects structure  -----------------------------------------------------------------

m2 <- lmer(Simpson_diversity ~ LandUse.2 + Kingdom + LandUse.2:Kingdom +
             (1|SS) + (1|SSB), data = diversity_simpson2) #Is Singular

simulation_simpson2 <- simulateResiduals(fittedModel = m2, plot = TRUE)
plot(m2)
hist(diversity_simpson2$Simpson_diversity)
# The simplest model does not converge



# ---- 14. Try with log(1/D)--------------------------------------------------------

# ---14.1 Choose random effects structure  -----------------------------------------------------------------


diversity_simpson <- diversity_simpson %>%
  
  # Create a new column to transform the response variable with log
  mutate(log_one_over_D = log(Simpson_diversity))


m1_log <- lmer(log_one_over_D ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
             (1|SS) + (1|SSB), data = diversity_simpson) 

m2_log <- lmer(log_one_over_D ~ LandUse.1 + Kingdom + LandUse.1:Kingdom + 
             (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_simpson)

m3_log <- lmer(log_one_over_D ~ LandUse.1 + Kingdom + LandUse.1:Kingdom + 
             (1+LandUse.1|SS) + (1|SSB), data = diversity_simpson)
# fit is singular

m4_log <- lmer(log_one_over_D ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +  
             (1+Predominant_land_use|SS) + (1|SSB), data = diversity_simpson)
# fit is singular

# Compare the models that converged
AIC(m1_log, m2_log)

# ---14.2. Choose fixed effects structure  -----------------------------------------------------------------

# See the significance of the terms in the model 
Anova(m2_log)

# ---14.3. Plot residuals: ----------------------------------------------------------------

simulationOutput2 <- simulateResiduals(fittedModel = m2_log, plot = TRUE)
plot(m2_log)

# ---14.4. Estimates: ----------------------------------------------------------------
summary(m2_log)

# ---14.4. Plot results: ----------------------------------------------------------------


# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of the plot labels
PlotErrBar_interactions(model = m2_log, resp = "Abundance", Effect1 = "LandUse.1", Effect2 = "Kingdom",
                        ylims = c(-2.5,2.5), pointtype = c(16,17,18),blackwhite = FALSE)



# Plot the differences between estimates 
PlotErrBar_interactions_modi(model =  m2_log, resp = "Abundance", Effect1 = "LandUse.1", Effect2 = "Kingdom",
                             ylims = c(-1,1), pointtype = c(16,17, 18),blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10), y = -0.9, labels = c("Primary Min",
                                         "Cropland All",
                                         "ISV LI",
                                         "ISV Min",
                                         "MSV All",
                                         "Pasture All",
                                         "Plantation LI",
                                         "Plantation Min",
                                         "Primary LI",
                                         "YSV A"), srt = 18, cex= 0.7)





# ---------------------RESULTS WITH THE BEST COMBINATION -------------------------------------------

# 15. Load table -------------------------------------------------------------------------------------

# Load Simpsons data that has already been filtered for studies that assessed more than 1 species
# and has corrected abundance measures

diversity_all_simpson <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Transform Simpson's diversity index
diversity_all_simpson <- diversity_all_simpson %>%
  
  # Create a new column to transform the response variable with log
  mutate(log_one_over_D = log(Simpson_diversity))

# 15.1. Merge land-use intensities ------------------------------------------------------------------------

# Call the function that merges lan-uses and intensities
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R")

# Create the vectors that hold the land-uses that we want to keep with 
# different use intensities 
land_uses_separate_final <- c("Primary","Cropland", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light 
# and intense use intensities
land_uses_light_intense_final <- c("Primary", "Cropland", "ISV", "Plantation forest")


# Create the levels of the best combination for the tables that have all the records 
diversity_all_simpson <- Merge_landUses_and_intensities(dataset = diversity_all_simpson,
                                                          index = 0, 
                                                          land_uses_separate_intensities = land_uses_separate_final,
                                                          land_uses_merge_light_intense = land_uses_light_intense_final,
                                                          "Primary Minimal use")

# ---15.2. Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(diversity_all_simpson[ , c("LandUse.0", "Kingdom")])

# ---15.3. Complete cases --------------------------------------------------------------------------------

# Create a table with complete cases
diversity_all_simpson <- drop_na(diversity_all_simpson, 
                                   Simpson_diversity, LandUse.0) %>% droplevels()

# Check number of sites
table(diversity_all_simpson$LandUse.0, diversity_all_simpson$Kingdom)

# ---15.4 Model with log: -------------------------------------------------------------- 

# I will log 1/D because otherwise the model does not converge
m1_final <- lmer(log_one_over_D ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                   (1|Source_ID) + (1|SS) + (1|SSB),
                 data = diversity_all_simpson)

# test significance of fixed effects 
Anova(m1_final)

# Export estimates
require(broom)
out <- tidy(m1_final)

# ---15.5 Plot residuals: ----------------------------------------------------------------

plot(m1_final)

#----15.6 Compare with model that has all the land-use intensities merge -----------------------------

# Merge all the intensity levels for all land-use categories 

# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate_null1 <- "NA"
# Create a vector with the lad-uses that we want to merge
land_uses_light_intense_null1 <- "NA"


# This function creates new columns with the land-uses and use intensities that we want, then it then merges both columns 
# into a LanUse column
diversity_all_simpson <- Merge_landUses_and_intensities(diversity_all_simpson, 
                                                          1,
                                                          land_uses_separate_null1, 
                                                          land_uses_light_intense_null1,
                                                          "Primary All")
# Check levels
levels(diversity_all_simpson$LandUse.1)


# Model richness with those land-use classes
m2_final <- lmer(log_one_over_D ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                   (1|Source_ID) + (1|SS) + (1|SSB),
                 data = diversity_all_simpson)

# Check they have the same sample size with the complex model
summary(m2_final)

# Compare AICs
AIC(m1_final, m2_final)

# Next we will check if we've lost a significant amount of explanatory power 
# by removing this interaction. If we have, we want to keep the more complex
# model. If we haven't lost a significant amount of explanatory power, then we
# can keep the simpler model. 

anova(m1_final, m2_final)

#--15.7. Compare with null model of plants -----------------------------------------------

# For the null model I am going to replace the nePlantae for Plantae 
diversity_all_simpson <- diversity_all_simpson %>%
  
  mutate(
    
    # Create a new kingdom column as the copy of the kingdom column we used
    # in the first model
    Kingdom.1 = paste(Kingdom),
    
    # Replace nePlantae for Plantae
    Kingdom.1 = recode_factor(Kingdom.1, "nePlantae" = "Plantae"), 
    
    # set reference level
    Kingdom.1 = factor(Kingdom.1),
    Kingdom.1 = relevel(Kingdom.1, ref = "Animalia"))

# Check the levels
levels(diversity_all_simpson$Kingdom.1)

m3_final <- lmer(log_one_over_D ~ LandUse.0 + Kingdom.1 + LandUse.0:Kingdom.1 +
                   (1|Source_ID) + (1|SS) + (1|SSB),
                 data = diversity_all_simpson)

# Check they have the same sample size with the complex model
summary(m3_final)

# check if we've lost a significant amount of explanatory power 
anova(m1_final, m3_final, test = "F")


# ---15.8. Plot results: ----------------------------------------------------------------

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of the plot labels
PlotErrBar_interactions(model = m1_final, resp = "Simpson's diversity index", Effect1 = "LandUse.0", Effect2 = "Kingdom",
                        ylims = c(-3,2.5), pointtype = c(16,17,18),blackwhite = FALSE)


# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m1_final,
                             resp = "Simpson's diversity index",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-1, 1),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -0.9, labels = c("Primary Minimal",
                          "Cropland Light-Intense",
                          "Cropland Minimal",
                          "ISV Light-Intense",
                          "ISV Minimal",
                          "MSV All",
                          "Pasture All",
                          "Plantation Light-Intense",
                          "Plantation Minimal",
                          "Primary Light-Intense",
                          "YSV All"),
     srt = 18, cex= 0.7)
