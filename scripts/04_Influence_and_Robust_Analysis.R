# Load libraries
library(influence.ME)
library(dplyr)
library(data.table)
library(broom)
library(lme4) 
library(robustlmm)

# Load the modified influence function
source("./R/04_Analise_Influence_modified_influenceME_function.R") 


# 1. Try to understand influence.ME -------------------------------------------------------------------------

# Load data to try solutions
sleepstudy<- lme4::sleepstudy

# Run a simple model
m1 <- lmer(Reaction ~ Days + (Days | Subject),
           sleepstudy)

# Get the coefficients for the model
summary(m1)

# Firts I'm going to explore what the influence function of the influence.ME does
INFME <- influence(m1, "Subject", count = TRUE)
INFME$alt.fixed

# Now, I'm going to see what the cooks.distance output is like
influence.ME::cooks.distance.estex(INFME)

# Now I'm going to explore what happens if I use the argument select
INFME2 <- influence(m1, group = "Subject", select = c("308"), count = TRUE) # It only exludes the subject I'm indicating
influence.ME::cooks.distance.estex(INFME2) # Since it only excludes one subject, it only gives one cook distance

# Now I'm going to see if I can select multiple levels to be excluded
INFME3 <- influence(m1, group = "Subject", select = c("308","309"), count = TRUE) 
influence.ME::cooks.distance.estex(INFME3) # The function excludes the levels but all at the same time, so I only get one value of cooks

# Now, I've modified the influence function of the package influence.ME,
# so that now I provide a vector with the levels that I don't want to exclude 
# in any of iterations
exclude <- c("308", "309", "310")

# I'm going to try the modified function
INFME4 <- influence_modified(m1, "Subject", count = TRUE) # It returns altered fixed parameters for 15 levels (out of a total of 18 levels)

# Get the fixed estimates of the original model (based on the full data)
INFME4$or.fixed # They are the same compared to the summary of the initial model m1

# Get the Standard Error of the estimates of the original model
INFME4$or.se # They are the same compared to the summary of the initial model m1

# Matrix of the fixed parameters estimate and SE, after iteratively subsets of data are removed
INFME4$alt.fixed
INFME4$alt.se

# Im going to compare if the results of the modified function are the sames if I do the iterations
# manually

# Update the original model, removing one subject each time
m1.1 <- update(m1, subset = Subject != "330")
m1.2 <- update(m1, subset = Subject != "331")
m1.3 <- update(m1, subset = Subject != "332")
m1.4 <- update(m1, subset = Subject != "333")
m1.5 <- update(m1, subset = Subject != "334")
m1.6 <- update(m1, subset = Subject != "335")
m1.7 <- update(m1, subset = Subject != "337")
m1.8 <- update(m1, subset = Subject != "349")

# Compare the coefficients of all of the models
compareCoefs(m1, m1.1, m1.2, m1.3, m1.4, m1.5, m1.6, m1.7, m1.8) # It works


# 2. Abundance data ------------------------------------------------------------------------------------

# Now I'm going to search for those source ID's that should be excluded beacuse without them 
# the models are rank deficient

# Load abundance data
abundance <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")

# Create the logAbundance column
abundance <- mutate(abundance, logAbundance = log(RescaledAbundance + 1))

# Add the Land use data
# Load function to merge LUI
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R") 

# Vector of the land uses where intensity levels will remain separate
land_uses_separate <- c("Primary", "Cropland", "ISV", "Plantation forest")

# Vector of land uses where light and intense uses will be merged
land_uses_light_intense <- c("Primary", "Cropland", "ISV", "Plantation forest")

# Merge LUI
abundance <- Merge_landUses_and_intensities(dataset = abundance, 
                                       index = 0, 
                                       land_uses_separate_intensities = land_uses_separate, 
                                       land_uses_merge_light_intense = land_uses_light_intense, 
                                       reference = "Primary Minimal use")


# Get complete cases
abundance <- drop_na(abundance, RescaledAbundance, LandUse.0) %>% droplevels()

# Organize levels of LUI
abundance$LandUse.0 <- factor(abundance$LandUse.0, 
                                          levels = c("Primary Minimal use",
                                                     "Primary Light-intense use",
                                                     "MSV All",
                                                     "ISV Minimal use",
                                                     "ISV Light-intense use",
                                                     "YSV All",
                                                     "Plantation forest Minimal use",
                                                     "Plantation forest Light-intense use",
                                                     "Pasture All",
                                                     "Cropland Minimal use",
                                                     "Cropland Light-intense use"))

# Run model
full_model_abundance <- lmer(logAbundance ~ LandUse.0 + Kingdom + LandUse.0*Kingdom +
                     (1|Source_ID) + (1|SS) + (1|SSB), abundance)

summary(full_model_abundance)

# Get the levels of Source ID
levels_abundance <- as.data.frame(influence.ME::grouping.levels(full_model_abundance, "Source_ID"))

# Apply original influence function to locate source IDs that we can't remove
influence.ME::influence(full_model_abundance, group = "Source_ID", count = TRUE)

# Remove Source ID 18: it  has 900 sites:
length(which(data$Source_ID == "DI1_2013__deLima"))

# It has more than half sites for cropland light-intense use and ISV light-intense use
deLima <-  data %>% base::subset(Source_ID == "DI1_2013__deLima") %>% droplevels() 
table(deLima$LandUse.0, deLima$Kingdom)

# Try to fit a model without that source ID
full_model.1_abundance <- update(full_model_abundance, subset = Source_ID != "DI1_2013__deLima")
summary(full_model.1)

# Use the modified function to search for more Source ID's that shouldn't be removed
exclude <- "DI1_2013__deLima"
INFM1_abundance <- influence_modified(full_model_abundance, group = "Source_ID", count = TRUE)

# Export list
saveRDS(INFM1_abundance, file = "./output/intermediate_files/04_Analysis_influence_Abundance.Rds")
INFM1_abundance <- readRDS(file = "./output/intermediate_files/04_Analysis_influence_Abundance.Rds")


# Now we're going to see if there are any influential cases: 

# The basic rationale behind measuring influential cases is that when iteratively single units 
# are omitted from the data, models based on these data should not produce substantially 
# different estimates. To standardize the assessment of how influential a (single group of)
# observation(s) is, several measures of influence are common practice. First, DFBETAS is a
# standardized measure of the absolute difference between the estimate with a particular case 
# included and the estimate without that particular case. Second, Cook's distance provides 
# an overall measurement of the change in all parameter estimates, or a selection thereof.

# First, I am going to calculate the Cooks distance

# Make a dataframe with the results of cooks distance
cooks_abundance <- as.data.frame(influence.ME::cooks.distance.estex(INFM1_abundance, sort = TRUE))
# Turn Source ID into a column
cooks_abundance <- setDT(cooks_abundance, keep.rownames = "Source_ID")
# Count number of levels 
length(influence.ME::grouping.levels(full_model_abundance, "Source_ID"))
# Calculate cut off for cook's distance, excluding one source ID that we're not evaluating
4/85
# Select the values of all columns that are greater for the cut off
cook_sources_abundance <- cooks_abundance %>% filter(V1 > 0.05) %>% droplevels
cook_Sources_abundance <- cooks_abundance %>% filter(V1 > 0.05) %>% droplevels %>% pull(Source_ID) %>% as.character()
# Calculate number of sites that belong to those source IDs
length(which(abundance$Source_ID %in% cook_Sources_abundance))


# Then the DF betas:
# DFBETAS (standardized difference of the beta) is a measure that standardizes the absolute 
# difference in parameter estimates between a (mixed effects) regression model based on a full 
# set of data, and a model from which a (potentially influential) subset of data is removed.
# A value for DFBETAS is calculated for each parameter in the model separately.

# Make a dataframe with the results of dfbetas
dfbetas_abundance <- as.data.frame(influence.ME::dfbetas.estex(INFM1_abundance))
# Turn Source ID into a column
dfbetas_abundance <- setDT(dfbetas_abundance, keep.rownames = "Source_ID")
# Count number of levels of the grouping variables
length(influence.ME::grouping.levels(full_model_abundance, "Source_ID"))
# Calcultae cut off
2/sqrt(85)
# Select the values of all columns that are greater for the cut off
dfbetas_sources_abundance <-dfbetas_abundance %>%  
  filter_if(is.numeric, any_vars(. > 0.22))
# Get the names of the Sources IDs that have a DFBETAS greater than the cutoff for any parameters
dfbetas_Sources_abundance <- dfbetas_sources_abundance%>% pull(Source_ID) %>% as.character
# Calculate number of Source IDs that are shared between both measures of change
length(which(cook_Sources_abundance %in% dfbetas_Sources_abundance))
# Try a different cutoff 
dfbetas_sources_abundance.2 <-dfbetas_abundance %>%  
  filter_if(is.numeric, any_vars(. > 1))


# Plot the results

# COOK'S distance
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for cooksdistance
# is 4/n, with n: number of groups in the grouping factor under evaluation
4/85
influence.ME::plot.estex(x=INFM1_abundance, which = "cook", cutoff = 0.047)


# DFBETAS
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for DFBETAS
# is 2/sqrt(n), with n: number of groups in the grouping factor under evaluation
2/(sqrt(85))
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 1:4)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 5:8)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 9:12)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 13:16)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 17:20)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 21:24)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 25:28)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 0.22, parameters = 29:30)

influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 1:4)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 5:8)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 9:12)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 13:16)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 17:20)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 21:24)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 25:28)
influence.ME::plot.estex(x=INFM1_abundance, which = "dfbetas", cutoff = 1, parameters = 29:30)

# 3. Species richness data ------------------------------------------------------------------------------------

# Now I'm going to search for those source ID's that should be excluded beacuse without them 
# the models are rank deficient

# Load species richness data
sp_richness <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Vector of the land uses where intensity levels will remain separate
land_uses_separate <- c("Primary","Cropland", "ISV", "Plantation forest")

# Vector of land uses where light and intense uses will be merged
land_uses_light_intense <- c("Primary", "Cropland", "ISV", "Plantation forest")

# Merge LUI
sp_richness <- Merge_landUses_and_intensities(dataset = sp_richness, 
                                       index = 0, 
                                       land_uses_separate_intensities = land_uses_separate, 
                                       land_uses_merge_light_intense = land_uses_light_intense, 
                                       reference = "Primary Minimal use")


# Get complete cases
sp_richness <- drop_na(sp_richness, 
                         Species_richness, LandUse.0) %>% droplevels()
# Organize levels of LUI
sp_richness$LandUse.0 <- factor(sp_richness$LandUse.0, 
                                levels = c("Primary Minimal use",
                                           "Primary Light-intense use",
                                           "MSV All",
                                           "ISV Minimal use",
                                           "ISV Light-intense use",
                                           "YSV All",
                                           "Plantation forest Minimal use",
                                           "Plantation forest Light-intense use",
                                           "Pasture All",
                                           "Cropland Minimal use",
                                           "Cropland Light-intense use"
                                ))

# Run model
full_model_richness <- glmer(Species_richness ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                    (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = sp_richness, family = poisson, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))


summary(full_model_richness)

# Get the levels of Source ID
levels_richness <- as.data.frame(influence.ME::grouping.levels(full_model_richness, "Source_ID"))

# Apply original influence function to locate source IDs that we can't remove
influence.ME::influence(full_model_richness, group = "Source_ID", count = TRUE)


# Use the modified function to search for more Source ID's that shouldn't be removed
exclude <- "DI1_2013__deLima"
INFM1_richness <- influence_modified(full_model_richness, group = "Source_ID", count = TRUE)

# Export list
saveRDS(INFM1_richness, file = "./output/intermediate_files/04_Analysis_influence_Species_richness.Rds")
INFM1_richness <- readRDS(file = "./output/intermediate_files/04_Analysis_influence_Species_richness.Rds")


# Make a dataframe with the results of cooks distance
cook_richness <- as.data.frame(influence.ME::cooks.distance.estex(INFM1_richness, sort = TRUE))
# Turn Source ID into a column
cook_richness <- setDT(cook_richness, keep.rownames = "Source_ID")
# Count number of grouping levels
length(influence.ME::grouping.levels(full_model_richness, "Source_ID"))
# Calculate cut off for cook's distance
4/84
# Select the values of all columns that are greater for the cut off
cook_sources_richness <- cook_richness %>% filter(V1 > 0.05) %>% droplevels
cook_Sources_richness <- cook_richness %>% filter(V1 > 0.05) %>% droplevels %>% pull(Source_ID) %>% as.character()
# Calculate number of sites that belong to those source IDs
length(which(sp_richness$Source_ID %in% cook_Sources_richness))


# Then the DF betas:
# Make a dataframe with the results of dfbetas
dfbetas_richness <- as.data.frame(influence.ME::dfbetas.estex(INFM1_richness))
# Turn Source ID into a column
dfbetas_richness <- setDT(dfbetas_richness, keep.rownames = "Source_ID")
# Calculate cut off for dfbetas
2/(sqrt(84))
# Select the values of all columns that are greater for the cut off
dfbetas_sources_richness <-dfbetas_richness %>%  
  filter_if(is.numeric, any_vars(. > 0.22))
# Get the names of the Sources IDs that have a DFBETAS greater than the cutoff for any parameters
dfbetas_Sources_richness <- dfbetas_sources_richness %>% pull(Source_ID) %>% as.character
# Count number of sources repeated in cook's and dfbetas
length(which(cook_Sources_richness %in% dfbetas_Sources_richness))

# Plot the results

# COOK'S distance
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for cooksdistance
# is 4/n, with n: number of groups in the grouping factor under evaluation
4/84
influence.ME::plot.estex(x=INFM1_richness, which = "cook", cutoff = 0.05)


# DFBETAS
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for DFBETAS
# is 2/sqrt(n), with n: number of groups in the grouping factor under evaluation
2/(sqrt(84))
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 1:4)
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 5:8)
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 9:12)
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 13:16)
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 17:20)
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 21:24)
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 25:28)
influence.ME::plot.estex(x=INFM1_richness, which = "dfbetas", cutoff = 0.22, parameters = 29:33)

# 4. Simpson's diversity ------------------------------------------------------------------------------------

# Now I'm going to search for those source ID's that should be excluded beacuse without them 
# the models are rank deficient

# Load species richness data
simpson <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Simpson_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Vector of the land uses where intensity levels will remain separate
land_uses_separate <- c("Primary","Cropland", "ISV", "Plantation forest")

# Vector of land uses where light and intense uses will be merged
land_uses_light_intense <- c("Primary", "Cropland", "ISV", "Plantation forest")

# Merge LUI
simpson <- Merge_landUses_and_intensities(dataset = simpson, 
                                              index = 0, 
                                              land_uses_separate_intensities = land_uses_separate, 
                                              land_uses_merge_light_intense = land_uses_light_intense, 
                                              reference = "Primary Minimal use")


# Get complete cases
simpson <- drop_na(simpson, 
                       Simpson_diversity, LandUse.0) %>% droplevels()

# Transform Simpson's diversity index
simpson  <- simpson  %>%
  
  # Create a new column to transform the response variable with log
  mutate(log_one_over_D = log(Simpson_diversity))

# Organize levels of LUI
simpson$LandUse.0 <- factor(simpson$LandUse.0, 
                            levels = c("Primary Minimal use",
                                       "Primary Light-intense use",
                                       "MSV All",
                                       "ISV Minimal use",
                                       "ISV Light-intense use",
                                       "YSV All",
                                       "Plantation forest Minimal use",
                                       "Plantation forest Light-intense use",
                                       "Pasture All",
                                       "Cropland Minimal use",
                                       "Cropland Light-intense use"
                            ))

# Run model
full_model_simpson <- lmer(log_one_over_D ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                             (1|Source_ID) + (1|SS) + (1|SSB),
                           data = simpson)

summary(full_model_simpson)

# Get the levels of Source ID
levels_simpson <- as.data.frame(influence.ME::grouping.levels(full_model_simpson, "Source_ID"))

# Apply original influence function to locate source IDs that we can't remove
influence.ME::influence(full_model_simpson, group = "Source_ID", count = TRUE)


# Use the modified function to search for more Source ID's that shouldn't be removed
exclude <- "DI1_2013__deLima"
INFM1_simpson <- influence_modified(full_model_simpson, group = "Source_ID", count = TRUE)

# Export list
saveRDS(INFM1_simpson, file = "./output/intermediate_files/04_Analysis_influence_Simpson.Rds")
INFM1_simpson <- readRDS(file = "./output/intermediate_files/04_Analysis_influence_Simpson.Rds")


# Make a dataframe with the results of cooks distance
cook_simpson <- as.data.frame(influence.ME::cooks.distance.estex(INFM1_simpson, sort = TRUE))
# Turn Source ID into a column
cook_simpson <- setDT(cook_simpson, keep.rownames = "Source_ID")
# Count number of grouping levels
length(influence.ME::grouping.levels(full_model_simpson, "Source_ID"))
# Calculate cut off for cook's distance
4/73
# Select the values of all columns that are greater for the cut off
cook_sources_simpson <- cook_simpson %>% filter(V1 > 0.06) %>% droplevels
cook_Sources_simpson <- cook_simpson %>% filter(V1 > 0.06) %>% droplevels %>% pull(Source_ID) %>% as.character()
# Calculate number of sites that belong to those source IDs
length(which(simpson$Source_ID %in% cook_Sources_simpson))


# Then the DF betas:
# Make a dataframe with the results of dfbetas
dfbetas_simpson <- as.data.frame(influence.ME::dfbetas.estex(INFM1_simpson))
# Turn Source ID into a column
dfbetas_simpson <- setDT(dfbetas_simpson, keep.rownames = "Source_ID")
# Calculate cut off for dfbetas
2/(sqrt(73))
# Select the values of all columns that are greater for the cut off
dfbetas_sources_simpson <-dfbetas_simpson %>%  
  filter_if(is.numeric, any_vars(. > 0.23))
# Get the names of the Sources IDs that have a DFBETAS greater than the cutoff for any parameters
dfbetas_Sources_simpson <- dfbetas_sources_simpson %>% pull(Source_ID) %>% as.character
# Count number of sources repeated in cook's and dfbetas
length(which(cook_Sources_simpson %in% dfbetas_Sources_simpson))

# Plot the results

# COOK'S distance
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for cooksdistance
# is 4/n, with n: number of groups in the grouping factor under evaluation
4/73
influence.ME::plot.estex(x=INFM1_simpson, which = "cook", cutoff = 0.06)


# DFBETAS
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for DFBETAS
# is 2/sqrt(n), with n: number of groups in the grouping factor under evaluation
2/(sqrt(73))
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 1:4)
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 5:8)
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 9:12)
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 13:16)
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 17:20)
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 21:24)
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 25:28)
influence.ME::plot.estex(x=INFM1_simpson, which = "dfbetas", cutoff = 0.23, parameters = 29:33)


# 5. Removing source IDs ------------------------------------------------------------------------------------

# 5.1 Abundance --------------------------------------------------------------------------------------

# I'm going to start by deleting the studies present in both cook's analysis and dfbetas
remo_abundance <- merge(x= cook_sources_abundance, y = dfbetas_sources_abundance, by = "Source_ID")
print(remo_abundance$Source_ID)

# Try to fit a model without those sources
m_remo_abundance <- update(full_model_abundance, subset = Source_ID %nin% 
                          c("SH1_2013__CIFORcameroon", "YY1_2009__Selwyn", "SH1_2002__Sheil",
                            "MH1_2011__Phalan", "CM1_2012__Katovai", "YY1_2018__Guillemot",
                            "TN1_2007__ODea", "DG1_2012__Ding", "HP1_2008__Hanson",
                            "DL1_2010__Luskin", "YY1_2016__Sunil", "YP1_2011__Wang",
                            "CM1_2009__Letcher"))

# Sources that create rank defficiency: "SH1_2013__CIFORphilippines"
summary(m_remo_abundance)

# Calculate estex object with full and modified model
INFM1_abundance_remove <- influence.ME::influence(model = full_model_abundance, group = "Source_ID",
                                                  select =  c("SH1_2013__CIFORcameroon", "YY1_2009__Selwyn", "SH1_2002__Sheil",
                                                              "MH1_2011__Phalan", "CM1_2012__Katovai", "YY1_2018__Guillemot",
                                                              "TN1_2007__ODea", "DG1_2012__Ding", "HP1_2008__Hanson",
                                                              "DL1_2010__Luskin", "YY1_2016__Sunil", "YP1_2011__Wang",
                                                              "CM1_2009__Letcher"))
INFM1_abundance_remove$or.fixed
dfbetas.estex(INFM1_abundance_remove)
cooks.distance.estex(INFM1_abundance_remove)
#  See what estimates are not longer significant
summary(full_model_abundance)
summary(m_remo_abundance)

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of the plot labels
PlotErrBar_interactions(model = full_model_abundance, resp = "Abundance", Effect1 = "LandUse.0", Effect2 = "Kingdom",
                        ylims = c(-3,2.5), pointtype = c(16,17,18),blackwhite = FALSE)


# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m_remo_abundance,
                             resp = "Rescaled total abundance",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-0.2, 0.45),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -0.18, labels = c("Primary Minimal",
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

# Try to run a robust model
library(robustlmm)

robust_abundance <- rlmer(logAbundance ~ LandUse.0 + Kingdom + LandUse.0*Kingdom +
                            (1|Source_ID) + (1|SS) + (1|SSB), abundance)

# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = robust_abundance,
                             resp = "Rescaled total abundance",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-0.2, 0.45),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -0.18, labels = c("Primary Minimal",
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



# 5.2 Species richness --------------------------------------------------------------------------------------


# I'm going to start by deleting the studies present in both cook's analysis and dfbetas
remo_richness <- merge(x= cook_sources_richness, y = dfbetas_sources_richness, by = "Source_ID", all.x = TRUE)
print(remo_richness$Source_ID)

# Try to fit a model without those sources
m_remo_richness <- update(full_model_richness, subset = Source_ID %nin% 
                             c("CM2_2012__Cicuzza", "SH1_2002__Sheil", "TN1_2007__ODea", 
                               "MH1_2011__Phalan", "CM1_2009__Letcher", "DG1_2012__Ding",
                               "SH1_2013__CIFORcameroon", "DI1_2013__Azhar", 
                               "CM2_2012__Sambuichi", "HP1_2009__Kessler", "SC1_2006__Mayfield",
                               "YY1_2015__Mandle", "CM2_2007__Sambuichi", "SC1_2011__Bakayoko",
                               "YY1_2016__Sunil", "YY1_2018__Guillemot",
                               "CM1_2012__Katovai"))

summary(m_remo_richness)

# Calculate estex object with full and modified model
INFM1_richness_remove <- influence.ME::influence(model = full_model_richness, group = "Source_ID",
                                                  select =  c("CM2_2012__Cicuzza", "SH1_2002__Sheil", "TN1_2007__ODea", 
                                                              "MH1_2011__Phalan", "CM1_2009__Letcher", "DG1_2012__Ding",
                                                              "SH1_2013__CIFORcameroon", "DI1_2013__Azhar", 
                                                              "CM2_2012__Sambuichi", "HP1_2009__Kessler", "SC1_2006__Mayfield",
                                                              "YY1_2015__Mandle", "CM2_2007__Sambuichi", "SC1_2011__Bakayoko",
                                                              "YY1_2016__Sunil", "YY1_2018__Guillemot",
                                                              "CM1_2012__Katovai"))
INFM1_richness_remove$or.fixed
dfbetas.estex(INFM1_richness_remove)
cooks.distance.estex(INFM1_richness_remove)
#  See what estimates are not longer significant
summary(full_model_richness)
summary(m_remo_richness)

# stargazer(full_model_richness, m_remo_richness, title="Results", align=TRUE)
# influence.ME::sigtest(INFM1_richness_remove)

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m_remo_richness,
                             resp = "Species richness",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-3,4),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -3, labels = c("Primary Minimal",
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

# 5.3 Simpson diversity  --------------------------------------------------------------------------------------


# I'm going to start by deleting the studies present in both cook's analysis and dfbetas
remo_simpson <- merge(x= cook_sources_simpson, y = dfbetas_sources_simpson, by = "Source_ID", all.x = TRUE)
print(remo_simpson$Source_ID)

# Try to fit a model without those sources
m_remo_simpson <- update(full_model_simpson, subset = Source_ID %nin% 
                            c("SH1_2002__Sheil", "TN1_2007__ODea", "MH1_2011__Phalan",
                              "CM1_2009__Letcher", "DG1_2012__Ding", "DI1_2013__Azhar",
                              "SC1_2011__Bakayoko", "YY1_2009__Selwyn", "HP1_2008__Presley",
                              "SC1_2006__Mayfield", "SH1_2013__CIFORcameroon"))

summary(m_remo_simpson)

# Calculate estex object with full and modified model
INFM1_simpson_remove <- influence.ME::influence(model = full_model_simpson, group = "Source_ID",
                                                 select =  c("SH1_2002__Sheil", "TN1_2007__ODea", "MH1_2011__Phalan",
                                                             "CM1_2009__Letcher", "DG1_2012__Ding", "DI1_2013__Azhar",
                                                             "SC1_2011__Bakayoko", "YY1_2009__Selwyn", "HP1_2008__Presley",
                                                             "SC1_2006__Mayfield", "SH1_2013__CIFORcameroon"))
INFM1_simpson_remove$or.fixed
dfbetas.estex(INFM1_simpson_remove)
cooks.distance.estex(INFM1_simpson_remove)
#  See what estimates are not longer significant
summary(full_model_simpson)
summary(m_remo_simpson)

# stargazer(full_model_richness, m_remo_richness, title="Results", align=TRUE)
# influence.ME::sigtest(INFM1_richness_remove)

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m_remo_simpson,
                             resp = "Simpson",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-1, 1),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -1, labels = c("Primary Minimal",
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


## 6. Use studies instead of Sources ---------------------------------------------------------------------------------

# 6.1. Abundance -----------------------------------------------------------------------
full_model_abundance <- lmer(logAbundance ~ LandUse.0 + Kingdom + LandUse.0*Kingdom +
                               (1|Source_ID) + (1|SS) + (1|SSB), abundance)

# Get the levels of Source ID
levels_abundance_SS <- as.data.frame(influence.ME::grouping.levels(full_model_abundance, "SS"))


# Apply original influence function to locate source IDs that we can't remove
INFM2_abundance <- influence.ME::influence(full_model_abundance, group = "SS", count = TRUE)


# Remove SS ID 20: it  has 600 sites:
length(which(data$SS == "DI1_2013__deLima 2"))

# Use the modified function to search for more Source ID's that shouldn't be removed
exclude <- "DI1_2013__deLima 2"
INFM2_abundance <- influence_modified(full_model_abundance, group = "SS", count = TRUE)

# Export list
saveRDS(INFM2_abundance, file = "./output/intermediate_files/04_Analysis_influence_Abundance_SS.Rds")
INFM2_abundance <- readRDS(file = "./output/intermediate_files/04_Analysis_influence_Abundance_SS.Rds")

# Make a dataframe with the results of cooks distance
cooks_abundance_SS <- as.data.frame(influence.ME::cooks.distance.estex(INFM2_abundance, sort = TRUE))
# Turn SS into a column
cooks_abundance_SS <- setDT(cooks_abundance_SS, keep.rownames = "SS")
# Count number of levels 
length(influence.ME::grouping.levels(full_model_abundance, "SS"))
# Calculate cut off for cook's distance, excluding one study that we're not evaluating
4/102
# Select the values of all columns that are greater for the cut off
cook_ss_abundance <- cooks_abundance_SS %>% filter(V1 > 0.04) %>% droplevels
cook_SS_abundance <- cooks_abundance_SS %>% filter(V1 > 0.04) %>% droplevels %>% pull(SS) %>% as.character()
# Calculate number of sites that belong to those studies
length(which(abundance$SS %in% cook_SS_abundance))


# Then the DF betas:
# Make a dataframe with the results of dfbetas
dfbetas_abundance_SS <- as.data.frame(influence.ME::dfbetas.estex(INFM2_abundance))
# Turn SS into a column
dfbetas_abundance_SS <- setDT(dfbetas_abundance_SS, keep.rownames = "SS")
# Calculate cut off
2/sqrt(102)
# Select the values of all columns that are greater for the cut off
dfbetas_ss_abundance <-dfbetas_abundance_SS %>%  
  filter_if(is.numeric, any_vars(. > 0.2))
# Get the names of the studies that have a DFBETAS greater than the cutoff for any parameters
dfbetas_SS_abundance <- dfbetas_ss_abundance%>% pull(SS) %>% as.character
# Calculate number of studies that are shared between both measures of change
length(which(cook_SS_abundance %in% dfbetas_SS_abundance))
# Try a different cutoff 
dfbetas_ss_abundance.2 <-dfbetas_abundance_SS %>%  
  filter_if(is.numeric, any_vars(. > 1))


# Plot the results

# COOK'S distance
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for cooksdistance
# is 4/n, with n: number of groups in the grouping factor under evaluation
4/102
influence.ME::plot.estex(x=INFM2_abundance, which = "cook", cutoff = 0.04)


# DFBETAS
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for DFBETAS
# is 2/sqrt(n), with n: number of groups in the grouping factor under evaluation
2/(sqrt(102))
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 1:4)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 5:8)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 9:12)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 13:16)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 17:20)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 21:24)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 25:28)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 0.2, parameters = 29:33)

influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 1:4)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 5:8)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 9:12)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 13:16)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 17:20)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 21:24)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 25:28)
influence.ME::plot.estex(x=INFM2_abundance, which = "dfbetas", cutoff = 1, parameters = 29:30)

# 6.2. Species richness  -----------------------------------------------------------------------

# Get the levels of Source ID
levels_richness_SS <- as.data.frame(influence.ME::grouping.levels(full_model_richness, "SS"))


# Apply original influence function to locate source IDs that we can't remove
INFM2_richness <- influence.ME::influence(full_model_richness, group = "SS", count = TRUE)

# Use the modified function to search for more Source ID's that shouldn't be removed
exclude <- "DI1_2013__deLima 2"
INFM2_richness <- influence_modified(full_model_richness, group = "SS", count = TRUE)

# Export list
saveRDS(INFM2_richness, file = "./output/intermediate_files/04_Analysis_influence_Richness_SS.Rds")
INFM2_richness <- readRDS(file = "./output/intermediate_files/04_Analysis_influence_Richness_SS.Rds")

# Make a dataframe with the results of cooks distance
cooks_richness_SS <- as.data.frame(influence.ME::cooks.distance.estex(INFM2_richness, sort = TRUE))
# Turn SS into a column
cooks_richness_SS <- setDT(cooks_richness_SS, keep.rownames = "SS")
# Count number of levels 
length(influence.ME::grouping.levels(full_model_richness, "SS"))
# Calculate cut off for cook's distance, excluding one study that we're not evaluating
4/108
# Select the values of all columns that are greater for the cut off
cook_ss_richness <- cooks_richness_SS %>% filter(V1 > 0.04) %>% droplevels
cook_SS_richness <- cooks_richness_SS %>% filter(V1 > 0.04) %>% droplevels %>% pull(SS) %>% as.character()
# Calculate number of sites that belong to those studies
length(which(sp_richness$SS %in% cook_SS_richness))


# Then the DF betas:
# Make a dataframe with the results of dfbetas
dfbetas_richness_SS <- as.data.frame(influence.ME::dfbetas.estex(INFM2_richness))
# Turn SS into a column
dfbetas_richness_SS <- setDT(dfbetas_richness_SS, keep.rownames = "SS")
# Calculate cut off
2/sqrt(102)
# Select the values of all columns that are greater for the cut off
dfbetas_ss_richness <-dfbetas_richness_SS %>%  
  filter_if(is.numeric, any_vars(. > 0.2))
# Get the names of the studies that have a DFBETAS greater than the cutoff for any parameters
dfbetas_SS_richness <- dfbetas_ss_richness %>% pull(SS) %>% as.character
# Calculate number of studies that are shared between both measures of change
length(which(cook_SS_richness %in% dfbetas_SS_richness))
# Try a different cutoff 
dfbetas_ss_richness.2 <-dfbetas_richness_SS %>%  
  filter_if(is.numeric, any_vars(. > 1))


# Plot the results

# COOK'S distance
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for cooksdistance
# is 4/n, with n: number of groups in the grouping factor under evaluation
4/102
influence.ME::plot.estex(x=INFM2_richness, which = "cook", cutoff = 0.04)


# DFBETAS
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for DFBETAS
# is 2/sqrt(n), with n: number of groups in the grouping factor under evaluation
2/(sqrt(90))
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 1:4)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 5:8)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 9:12)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 13:16)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 17:20)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 21:24)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 25:28)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 29:33)



# 6.3. Simpson -----------------------------------------------------------------------

# Get the levels of Source ID
levels_simpson_SS <- as.data.frame(influence.ME::grouping.levels(full_model_simpson, "SS"))


# Apply original influence function to locate source IDs that we can't remove
INFM2_simpson <- influence.ME::influence(full_model_simpson, group = "SS", count = TRUE)

# Use the modified function to search for more Source ID's that shouldn't be removed
exclude <- "DI1_2013__deLima 2"
INFM2_simpson <- influence_modified(full_model_simpson, group = "SS", count = TRUE)

# Export list
saveRDS(INFM2_simpson, file = "./output/intermediate_files/04_Analysis_influence_Simpson_SS.Rds")
INFM2_abundance <- readRDS(file = "./output/intermediate_files/04_Analysis_influence_Simpson_SS.Rds")

# Make a dataframe with the results of cooks distance
cooks_simpson_SS <- as.data.frame(influence.ME::cooks.distance.estex(INFM2_simpson, sort = TRUE))
# Turn SS into a column
cooks_simpson_SS <- setDT(cooks_simpson_SS, keep.rownames = "SS")
# Count number of levels 
length(influence.ME::grouping.levels(full_model_simpson, "SS"))
# Calculate cut off for cook's distance, excluding one study that we're not evaluating
4/90
# Select the values of all columns that are greater for the cut off
cook_ss_simpson <- cooks_simpson_SS %>% filter(V1 > 0.044) %>% droplevels
cook_SS_simpson <- cooks_simpson_SS %>% filter(V1 > 0.044) %>% droplevels %>% pull(SS) %>% as.character()
# Calculate number of sites that belong to those studies
length(which(simpson$SS %in% cook_SS_simpson))


# Then the DF betas:
# Make a dataframe with the results of dfbetas
dfbetas_simpson_SS <- as.data.frame(influence.ME::dfbetas.estex(INFM2_simpson))
# Turn SS into a column
dfbetas_simpson_SS <- setDT(dfbetas_simpson_SS, keep.rownames = "SS")
# Calculate cut off
2/sqrt(90)
# Select the values of all columns that are greater for the cut off
dfbetas_ss_simpson <-dfbetas_simpson_SS %>%  
  filter_if(is.numeric, any_vars(. > 0.21))
# Get the names of the studies that have a DFBETAS greater than the cutoff for any parameters
dfbetas_SS_simpson <- dfbetas_ss_simpson %>% pull(SS) %>% as.character
# Calculate number of studies that are shared between both measures of change
length(which(cook_SS_simpson %in% dfbetas_SS_simpson))
# Try a different cutoff 
dfbetas_ss_simpson.2 <-dfbetas_simpson_SS %>%  
  filter_if(is.numeric, any_vars(. > 1))


# Plot the results

# COOK'S distance
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for cooksdistance
# is 4/n, with n: number of groups in the grouping factor under evaluation
4/90
influence.ME::plot.estex(x=INFM2_simpson, which = "cook", cutoff = 0.044)


# DFBETAS
# According to Nieuwenhuis et al., (2012), as a rule of thumb a cutoff value for DFBETAS
# is 2/sqrt(n), with n: number of groups in the grouping factor under evaluation
2/(sqrt(90))
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 1:4)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 5:8)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 9:12)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 13:16)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 17:20)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 21:24)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 25:28)
influence.ME::plot.estex(x=INFM2_simpson, which = "dfbetas", cutoff = 0.21, parameters = 29:33)

## 7. Robust analysis  ---------------------------------------------------------------------------------

# Robust linear mixed effects models:

# ------------------------------------- Abundance -------------------------------------------------

# Run the original model
full_model_abundance <- lmer(logAbundance ~ LandUse.0 + Kingdom + LandUse.0*Kingdom +
                               (1|Source_ID) + (1|SS) + (1|SSB), data= abundance)

# Run the robust analysis
full_model_abundance_robust <- rlmer(logAbundance ~ LandUse.0 + Kingdom + LandUse.0*Kingdom +
                                       (1|Source_ID) + (1|SS) + (1|SSB), data = abundance)

# Plot the residuals
plot(full_model_abundance_robust)

# Get the robustness weights for the sites
we_abundance <- getME(full_model_abundance_robust, "w_e")
# histogram of the robustness weights of sites
histogram(we_abundance)
# Get the robustness weights for the groups
wb_abundance <- getME(full_model_abundance_robust, "w_b")
# histogram of the robustness weights of source IDs
wb_abundance[["Source_ID"]] %>% pull() %>% histogram()
# histogram of the robustness weights of studies
wb_abundance[["SS"]] %>% pull() %>% histogram()
# histogram of the robustness weights of SSB
wb_abundance[["SSB"]] %>% pull() %>% histogram()

# Get the estimates
summary(full_model_abundance)
summary(full_model_abundance_robust)

# Plot the results
# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")


# Export pdf with the graph
pdf(file = "./output/figures/04_Influence_and_Robust_Analysis_Abundance_results.pdf", width = 10)

# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = full_model_abundance_robust,
                             resp = "log(Total abundance + 1)",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-0.25, 0.3),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -0.18, labels = c("Primary Minimal use",
                           "Primary Light-intense use",
                           "MSV All",
                           "ISV Minimal use",
                           "ISV Light-intense use",
                           "YSV All",
                           "Plantation forest Minimal use",
                           "Plantation forest Light-intense use",
                           "Pasture All",
                           "Cropland Minimal use",
                           "Cropland Light-intense use"),
     srt = 18, cex= 0.7)


# Clear 
dev.off()


# ------------------------------------- Sp. richness -------------------------------------------------

# Run the original model
full_model_richness <- glmer(Species_richness ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                               (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = sp_richness, family = poisson, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))


# Run the robust analysis:
# Since the package only works for gaussian family, I'm going to model sp. richness as log(sp. richness + 1) 
# Create the logAbundance column
sp_richness <- mutate(sp_richness, logRichness = log(Species_richness + 1))

full_model_richness_robust <- rlmer(logRichness ~ LandUse.0 + Kingdom + LandUse.0*Kingdom +
                                      (1|Source_ID) + (1|SS) + (1|SSB), data = sp_richness)
# Run the lme4 model but with log
full_model_richness_log <- lmer(logRichness ~ LandUse.0 + Kingdom + LandUse.0*Kingdom +
                                      (1|Source_ID) + (1|SS) + (1|SSB), data = sp_richness)


# Plot the residuals
plot(full_model_richness_robust)

# Get the robustness weights for the sites
we_richness <- getME(full_model_richness_robust, "w_e")
# histogram of the robustness weights of sites
histogram(we_richness)
# Get the robustness weights for the groups
wb_richness <- getME(full_model_richness_robust, "w_b")
# histogram of the robustness weights of source IDs
wb_richness[["Source_ID"]] %>% pull() %>% histogram()
# histogram of the robustness weights of studies
wb_richness[["SS"]] %>% pull() %>% histogram()
# histogram of the robustness weights of SSB
wb_richness[["SSB"]] %>% pull() %>% histogram()

# Get the estimates
summary(full_model_richness)
summary(full_model_richness_robust)

# Plot the results
# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Export pdf with the graph
pdf(file = "./output/figures/04_Influence_and_Robust_Analysis_Species_Richness_results.pdf", width = 10)

# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = full_model_richness_robust,
                             resp = "log(Species Richness + 1)",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-3, 3),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -3, labels = c("Primary Minimal use",
                        "Primary Light-intense use",
                        "MSV All",
                        "ISV Minimal use",
                        "ISV Light-intense use",
                        "YSV All",
                        "Plantation forest Minimal use",
                        "Plantation forest Light-intense use",
                        "Pasture All",
                        "Cropland Minimal use",
                        "Cropland Light-intense use"),
     srt = 18, cex= 0.7)

# Clear
dev.off()

# Export pdf with the log graph
pdf(file = "./output/figures/04_Influence_and_Robust_Analysis_Species_Richness_results_log.pdf", width = 10)

# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = full_model_richness_log,
                             resp = "log(Species Richness + 1)",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-3, 3),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -3, labels = c("Primary Minimal use",
                        "Primary Light-intense use",
                        "MSV All",
                        "ISV Minimal use",
                        "ISV Light-intense use",
                        "YSV All",
                        "Plantation forest Minimal use",
                        "Plantation forest Light-intense use",
                        "Pasture All",
                        "Cropland Minimal use",
                        "Cropland Light-intense use"),
     srt = 18, cex= 0.7)

# Clear
dev.off()


# ------------------------------------- Simpson diversity -------------------------------------------------

# Run the original model
full_model_simpson <- lmer(log_one_over_D ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                             (1|Source_ID) + (1|SS) + (1|SSB),
                           data = simpson)

# Run the robust analysis:
full_model_simpson_robust <- rlmer(log_one_over_D ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                                     (1|Source_ID) + (1|SS) + (1|SSB),
                                   data = simpson)

# Plot the residuals
plot(full_model_simpson_robust)

# Get the robustness weights for the sites
we_simpson <- getME(full_model_simpson_robust, "w_e")
# histogram of the robustness weights of sites
histogram(we_simpson)
# Get the robustness weights for the groups
wb_simpson <- getME(full_model_simpson_robust, "w_b")
# histogram of the robustness weights of source IDs
wb_simpson[["Source_ID"]] %>% pull() %>% histogram()
# histogram of the robustness weights of studies
wb_simpson[["SS"]] %>% pull() %>% histogram()
# histogram of the robustness weights of SSB
wb_simpson[["SSB"]] %>% pull() %>% histogram()

# Get the estimates
summary(full_model_simpson)
summary(full_model_simpson_robust)

# Plot the results
# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Export pdf with the graph
pdf(file = "./output/figures/04_Influence_and_Robust_Analysis_Simpson_results.pdf", width = 10)


# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = full_model_simpson_robust,
                             resp = "log (Simpson's diversity)",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-1, 1),
                             pointtype = c(16,17, 18),
                             blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), 
     y = -1, labels = c("Primary Minimal use",
                        "Primary Light-intense use",
                        "MSV All",
                        "ISV Minimal use",
                        "ISV Light-intense use",
                        "YSV All",
                        "Plantation forest Minimal use",
                        "Plantation forest Light-intense use",
                        "Pasture All",
                        "Cropland Minimal use",
                        "Cropland Light-intense use"),
     srt = 18, cex= 0.7)

# Clear
dev.off()
