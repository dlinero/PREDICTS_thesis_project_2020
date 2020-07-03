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
library(aods3) # to check for overdispersion



# --- Description -----------------------------------------------------------------

# Model Within-sample species richness as the number of total species sampled in each site. 


# ---- 1. Load data ------------------------------------------------------------------

# Load the table that contains the frugivore information along with the plants dispersed
# and not dispersed by animals 

diversityS <- readRDS(file = "./output/cleaned_data/01_Filter_data_frugi_endooPlants_notEndooPlants_records.rds")

# ---- 2. Select studies with more than one species ---------------------------------

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

# Filter the dataset
diversityS <- diversityS %>% base::subset(SS %in% list) %>% droplevels()
# ----3. Remove sites that will produce NaN values --------------------------------------------
# Doing the abundance models, I identified some sites where (due to the filtering of specific species)
# recorded only zero abundances. I am going to remove those sites as they will produce 
# Non infinite values. 

diversityS <- diversityS %>%  base::subset(SSBS %nin% c("MJ1_2009__Lehouck 2 Fururu 10", 
                                                        "MJ1_2009__Lehouck 2 Fururu 11", 
                                                        "MJ1_2009__Lehouck 2 Macha 12", 
                                                        "MJ1_2009__Lehouck 2 Macha 13")) %>%
  base::droplevels()

# ---- 4. Split datasets --------------------------------------------------------

# I am going to separate the records that belong to plants not dispersed by animals
diversityS_notEndo <- diversityS %>% base::subset(Kingdom == "nePlantae") %>% base::droplevels()

# Get the table without plants dispersed by animals
diversityS_frugi_endo <- diversityS %>% base::subset(Kingdom != "nePlantae") %>% base::droplevels()

# ---- 5. Merge Sites for animals and Endooplants-------------------------------------------------------------

# We merge sites within studies that have identical coordinates, start and end dates(and land-use 
# type and intensity). We do this because sometimes authors record different points in a transect as
# different sites, which might not be meaningful if they share land-use, intensity and coordinates (the
# coordinates have an error) OR maybe it is, but we are having a conservative approach

diversityS_frugi_endo <- yarg::MergeSites(diversityS_frugi_endo, silent = TRUE, 
                                          merge.extra = "Wilderness_area")

diversityS_notEndo <- yarg::MergeSites(diversityS_notEndo, silent = TRUE, 
                                       merge.extra = "Wilderness_area")

# ----6. Rename Predominant habitat --------------------------------------------------

# Rename the column predominant habitat, as the dataset is actually refering to land use
diversityS_frugi_endo <-  dplyr::rename(diversityS_frugi_endo,
                                        Predominant_land_use = Predominant_habitat)

diversityS_notEndo <-  dplyr::rename(diversityS_notEndo,
                                     Predominant_land_use = Predominant_habitat)


# ----7.  Calculate diversity metrics -----------------------------------------------

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


# ---8. Check the results--------------------------------------------------------------------

# Check the results for studies that assessed abundance
diversityS_frugi_endo %>% 
  
  # subset the table with all of the species records
  base::subset(SSS == "BS1_2010__Page 1 1") %>%
  
  # Select species that have an abundance greater than zero
  dplyr::filter(Measurement > 0) %>%
  
  # Create a column with the number of species present
  dplyr::mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  dplyr::select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_frugi_endo %>% base::subset(SSS == "BS1_2010__Page 1 1") %>% dplyr::select(SSS, Species_richness)

# Check the results for studies that assessed abundance
diversityS_notEndo %>% 
  
  # subset the table with all of the species records
  base::subset(SSS == "CM1_2007__MarinSpiotta 1 1") %>%
  
  # Select species that have an abundance greater than zero
  dplyr::filter(Measurement > 0) %>%
  
  # Create a column with the number of species present
  dplyr::mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  dplyr::select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_notendo %>% base::subset(SSS == "CM1_2007__MarinSpiotta 1 1") %>% dplyr::select(SSS, Species_richness)

# Check the results with studies that occurrence
diversityS_frugi_endo %>% 
  
  # subset the table with all of the species records
  base::subset(SSS == "DL1_2011__Latta 2 1") %>%
  
  # Select species that are present in the site 
  dplyr::filter(Measurement == 1) %>% 
  
  # Create a column with the number of species present
  dplyr::mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  dplyr::select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_frugi_endo %>% base::subset(SSS == "DL1_2011__Latta 2 1") %>% dplyr::select(SSS, Species_richness)

# Check the results with studies that occurrence
diversityS_notEndo %>% 
  
  # subset the table with all of the species records
  base::subset(SSS == "CM2_2012__Cicuzza 1 1") %>%
  
  # Select species that are present in the site 
  dplyr::filter(Measurement == 1) %>% 
  
  # Create a column with the number of species present
  dplyr::mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  dplyr::select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_notendo %>% base::subset(SSS == "CM2_2012__Cicuzza 1 1") %>% dplyr::select(SSS, Species_richness)

# ---9. Merge site metrics -------------------------------------------------------------------
# Merge the site metrics for all organisms

diversity_all <- base::rbind.data.frame(diversity1S_frugi_endo, diversity1S_notendo)

# Export table
saveRDS(diversity_all, file = "./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")

# Import table
diversity_all <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")


############################## MODEL TESTING #################################

# ---10. Simple model ----------------------------------------------------------------------

# What I'm going to do is a Plant contrast, that means I'm firts going to model the complex model
# (the one that assumme the reponse between endo and not-endo plants is significantly different), in this model 
# I will  put the Kingdom interaction that has 3 levels: animals and endo plants and not endo plants. Then I am going
# to model the null model that assumes there is no difference between endo and not-endo plants. In this model
# the Kingdom will have 2 levels: animals and plants. Then I will compare both models using an ANOVA

# 10.1. Checking number of sites -------------------------------------------------------------

# I will fisrt explore the number of sites that have a known land use category

# Drop sites that don't have richness measures or land-use data
diversity_all <- drop_na(diversity_all, 
                         Species_richness, Predominant_land_use) %>% droplevels()

# Check number of sites
addmargins(table(diversity_all$Predominant_land_use, diversity_all$Use_intensity, diversity_all$Kingdom), 2)


# 10.2. Merging land-use intensities  -------------------------------------------------------------
# According to the number of sites, we can split the land-use types and intensities in:

# Primary can be divided into the three level intensities in all cases
#	Cropland can be divided in minimal use and light/intense use in animal, endo plant and non endo plant
# ISV can be divided in minimal use and light/intense use in animal, endo plant and non endo plant
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense use in animal, endo plant and non endo plant
#	YSV has to be merged 


# Call the function that merges lan-uses and intensities
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R")
# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate <- c("Primary", "Cropland", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense <- c("Cropland", "ISV", "Plantation forest")


diversity_all <- Merge_landUses_and_intensities(dataset = diversity_all, index = 1, 
                                                land_uses_separate_intensities = land_uses_separate,
                                                land_uses_merge_light_intense = land_uses_light_intense,
                                                "Primary Minimal use")

# Check levels
levels(diversity_all$LandUse.1)

# ---10.3. Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(diversity_all[ , c("LandUse.1", "Kingdom")])

# ---10.4. Complete cases --------------------------------------------------------------------------------

# Check number of sites
table(diversity_all$LandUse.1, diversity_all$Kingdom)

# ---10.5. Choose between GLMM or LMM -----------------------------------------------------------------

# Species richness was modelled with Poisson error distribution and log-link function;
# there was evidence of significant overdispersion in these models so an observation-level
# random effect was included to account for this (i.e., a Poisson-lognormal model)(De Palma et al., 2016).

hist(diversity_all$Species_richness)

# ---10.6 Choose random effects structure  -----------------------------------------------------------------

# To select the random-effects structure, we use the method recommended by (Zuur et al., 2009) 
# of taking the most complex fixed-effects structure, including all interactions, that will be 
# tested in the second stage of modelling, and use it to compare different random-effects 
# structures.

m1 <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
              (1|SS) + (1|SSB), data = diversity_all, family = poisson)

# Warning message:
#In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.0399178 (tol = 0.002, component 1)

m1 <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
              (1|SS) + (1|SSB), data = diversity_all, family = poisson, 
            
            # I'm increasing the number of iterations as the model initially 
            # didn't converge
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))


m2 <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
              (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_all, family = poisson, 
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

m3 <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
              (1+LandUse.1|SS) + (1|SSB), data = diversity_all, family = poisson, 
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
# Is singular

m4 <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
              (1+Predominant_land_use.1|SS) + (1|SSB), data = diversity_all, family = poisson, 
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
# Is singular

# Compare the models that converged
AIC(m1, m2)


# ---10.7. Test for overdispersion: ----------------------------------------------------------------

# I am going to use the Ben Bolker's function available at:
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(m2)

# Or using the aods3 package
gof(m2)


# ---10.8. Model with quasipoisson: -------------------------------------------------------------- 
# I am going to add an observation-level random effects to account for the overdispersion of the data

m2_q <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_all, family = poisson, 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))


# ---10.11 Plot residuals: ----------------------------------------------------------------

simulationOutputq <- simulateResiduals(fittedModel = m2_q)
# Acces the qq plot
plotQQunif(simulationOutputq)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutputq)

# Test for overdispersion
overdisp_fun(m2_q)


#--------------------11. See what combination of land-use and land-use intensities is better -----------------------------
# In order to know if it's better to use the same levels used in the abundance model I will fit the 
# abundance and richness with the same land-use/intensities levels, and add their AIC and compare which option is better.


# --- First option ---------------------------------------------------------------------------------------------
# For the first attempt I will:

# Primary can be divided into the three level intensities in all cases
#	Cropland can be divided in minimal use and light/intense use (for the abundance endo plants have 19 sites in cropland )
# ISV can be divided in minimal use and light/intense use 
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense 
#	YSV has to be merged 

# Richness

# Call the function that merges lan-uses and intensities
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R")
# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate <- c("Primary", "Cropland", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light and intense use intensities
land_uses_light_intense <- c("Cropland", "ISV", "Plantation forest")


diversity_all <- Merge_landUses_and_intensities(dataset = diversity_all, index = 1, 
                                                land_uses_separate_intensities = land_uses_separate,
                                                land_uses_merge_light_intense = land_uses_light_intense,
                                                "Primary Minimal use")

# Get complete cases, in order to compare models that use the same rows
diversity_richness <- drop_na(diversity_all, Species_richness, LandUse.1)  %>% droplevels()

#Check levels
levels(diversity_richness$LandUse.1)


# Model richness with those land-use classes
m2_q <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_richness, family = poisson, 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))



summary(m2_q)
AIC(m2_q)

simulationOutputq <- simulateResiduals(fittedModel = m2_q)
# Acces the qq plot
plotQQunif(simulationOutputq)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutputq)

# Test for overdispersion
overdisp_fun(m2_q)

# Abundance

# Import table that has the site metrics calculated for the abudance model (that means including animals and both types of plants)
diversity_all_abundance <- readRDS( file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")


# This function creates new columns with the land-uses and use intensities that we want, then it then merges both columns 
# into a LanUse column
diversity_all_abundance <- Merge_landUses_and_intensities(diversity_all_abundance, 
                                                          1,
                                                          land_uses_separate, 
                                                          land_uses_light_intense,
                                                          "Primary Minimal use")

# Get complete cases, in order to compare models that use the same rows
diversity_abundance <- drop_na(diversity_all_abundance, RescaledAbundance, LandUse.1)  %>% droplevels()

# Check levels
levels(diversity_abundance$LandUse.1)

# Transform rescaled abundance
diversity_abundance <- mutate(diversity_abundance, 
                              logAbundance = log(RescaledAbundance + 1)
)


# Model abundance with those land-use classes
m2_a <- lmer(logAbundance ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_abundance)
summary(m2_a)
AIC(m2_a)

simulationOutputa <- simulateResiduals(fittedModel = m2_a)
# Acces the qq plot
plotQQunif(simulationOutputa)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutputa)

# --- Second option ---------------------------------------------------------------------------------------------
# For the second attempt I will:

# Primary can be divided into the three level intensities in all cases
#	Cropland has to be merged
# ISV can be divided in minimal use and light/intense use 
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense 
#	YSV has to be merged 

# Richness

# The only difference between both models is that we are going to merge the intensities for cropland
# in the second option

# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate2 <- c("Primary","ISV", "Plantation forest")
# Create a vector with the lad-uses that we want to merge
land_uses_light_intense2 <- c("ISV", "Plantation forest")


# This function creates new columns with the land-uses and use intensities that we want, then it then merges both columns 
# into a LanUse column
diversity_richness <- Merge_landUses_and_intensities(diversity_richness, 
                                                     2,
                                                     land_uses_separate2, 
                                                     land_uses_light_intense2,
                                                     "Primary Minimal use")
# Check levels
levels(diversity_richness$LandUse.2)

# Model richness with those land-use classes
m2_q2 <- glmer(Species_richness ~ LandUse.2 + Kingdom + LandUse.2:Kingdom +
                 (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_richness, family = poisson, 
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

summary(m2_q2)
AIC(m2_q2)

simulationOutputq_2 <- simulateResiduals(fittedModel = m2_q2)
# Acces the qq plot
plotQQunif(simulationOutputq_2)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutputq_2)

# Test for overdispersion
overdisp_fun(m2_q2)


# Abundance

# This function creates new columns with the land-uses and use intensities that we want, then it merges both columns 
# into a LandUse column
diversity_abundance <- Merge_landUses_and_intensities(diversity_abundance, 
                                                      2,
                                                      land_uses_separate2, 
                                                      land_uses_light_intense2, 
                                                      "Primary Minimal use")
# Check levels
levels(diversity_abundance$LandUse.2)

# Model abundance with those land-use classes
m2_a_2 <- lmer(logAbundance ~ LandUse.2 + Kingdom + LandUse.2:Kingdom +
                 (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_abundance)
summary(m2_a_2)
AIC(m2_a_2)

simulationOutputa_2 <- simulateResiduals(fittedModel = m2_a_2)
# Acces the qq plot
plotQQunif(simulationOutputa_2)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutputa_2)

# Compare AICs

# For the first combination
AIC(m2_q)
AIC(m2_a)
AIC(m2_q) + AIC(m2_a)
# For the second combination
AIC(m2_q2)
AIC(m2_a_2)
AIC(m2_q2) + AIC(m2_a_2)

# According to these results I should keep the land-use levels that the species richness models
# allowed me to separate
# That means my more complex model will be 
m2_q <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_all, family = poisson, 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))


levels(diversity_all$LandUse.1)

#--------------------12. Compare with model that has all the land-use intensities merge -----------------------------

# Complete cases: get the rows that have information for the levels that we choose in the more complex model
diversity_all <- drop_na(diversity_all, 
                         Species_richness, LandUse.1) %>% droplevels()


# Merge all the intensity levels for all land-use categories 

# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate3 <- "NA"
# Create a vector with the lad-uses that we want to merge
land_uses_light_intense3 <- "NA"


# This function creates new columns with the land-uses and use intensities that we want, then it then merges both columns 
# into a LanUse column
diversity_all <- Merge_landUses_and_intensities(diversity_all, 
                                                3,
                                                land_uses_separate3, 
                                                land_uses_light_intense3,
                                                "Primary All")
# Check levels
levels(diversity_all$LandUse.3)


# Model richness with those land-use classes
m2_q_3 <- glmer(Species_richness ~ LandUse.3 + Kingdom + LandUse.3:Kingdom +
                  (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_all, family = poisson, 
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

AIC(m2_q, m2_q_3)

# Next we will check if we've lost a significant amount of explanatory power 
# by removing this interaction. If we have, we want to keep the more complex
# model. If we haven't lost a significant amount of explanatory power, then we
# can keep the simpler model. 

anova(m2_q, m2_q_3)

#--------------------13. Compare with null model of plants -----------------------------

# ---13.1. Create null model: -----------------------------------------------

# For the null model I am going to replace the nePlantae for Plantae 
diversity_all <- diversity_all %>%
  
  
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
levels(diversity_all$Kingdom.1)

# ---13.2. Compare complex and null models: -----------------------------------------------

m2_q_4 <- glmer(Species_richness ~ LandUse.1 + Kingdom.1 + LandUse.1:Kingdom.1 +
                  (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_all, family = poisson, 
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
summary(m2_q_4)

anova(m2_q_4, m2_q, test = "F")


# ---13.4 Plot results: ----------------------------------------------------------------

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of the plot labels
PlotErrBar_interactions(model = m2_q, resp = "Abundance", Effect1 = "LandUse.1", Effect2 = "Kingdom",
                        ylims = c(-2.5,2.5), pointtype = c(16,17,18),blackwhite = FALSE)



# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m2_q, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                             ylims = c(-3,2.5), pointtype = c(16,17, 18),blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), y = -3, labels = c("Primary M",
                                         "Cropland LI",
                                         "Cropland M",
                                         "ISV LI",
                                         "ISV M",
                                         "MSV A",
                                         "Pasture A",
                                         "Plantation LI",
                                         "Plantation M",
                                         "Primary I",
                                         "Primary L",
                                         "YSV A"), srt = 18, cex= 0.7)


# ---------------------RESULTS WITH THE BEST COMBINATION -------------------------------------------

# 14.1. Load table -------------------------------------------------------------------------------------
# Import table that has already been filteres to include studies that assessed more than 1 species
diversity_all <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Richness_Site_metrics_animals_endooPlants_notendooPlants.rds")

# 14.2. Merge land-use intensities ------------------------------------------------------------------------

# Call the function that merges lan-uses and intensities
source("./R/02_Statistical_Analysis_merge_LandUses_Intensities.R")

# Create the vectors that hold the land-uses that we want to keep with 
# different use intensities 
land_uses_separate_final <- c("Primary","Cropland", "ISV", "Plantation forest")
# Create a vector with the land-uses where we want to merge the light 
# and intense use intensities
land_uses_light_intense_final <- c("Primary", "Cropland", "ISV", "Plantation forest")


# Create the levels of the best combination for the tables that have all the records 
diversity_all <- Merge_landUses_and_intensities(dataset = diversity_all,
                                                         index = 0, 
                                                         land_uses_separate_intensities = land_uses_separate_final,
                                                         land_uses_merge_light_intense = land_uses_light_intense_final,
                                                         "Primary Minimal use")

# ---14.3 Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(diversity_all[ , c("LandUse.0", "Kingdom")])

# ---14.4. Complete cases --------------------------------------------------------------------------------

# Create a table with complete cases
diversity_all <- drop_na(diversity_all, 
                              Species_richness, LandUse.0) %>% droplevels()

# Check number of sites
table(diversity_all$LandUse.0, diversity_all$Kingdom)

# ---14.5. Choose between GLMM or LMM -----------------------------------------------------------------
 
# Species richness was modelled with Poisson error distribution and log-link function;
# there was evidence of significant overdispersion in these models so an observation-level
# random effect was included to account for this (i.e., a Poisson-lognormal model)(De Palma et al., 2016).

m1_final <- glmer(Species_richness ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                    (1|Source_ID) + (1|SS) + (1|SSB) , data = diversity_all, family = poisson, 
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

# ---14.6. Test for overdispersion: ----------------------------------------------------------------

# I am going to use the Ben Bolker's function available at:
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(m1_final) # It's overdispersed


# ---14.7. Model with quasipoisson: -------------------------------------------------------------- 
# I am going to add an observation-level random effects to account for 
# the overdispersion of the data

m2_final <- glmer(Species_richness ~ LandUse.0 + Kingdom + LandUse.0:Kingdom +
                (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_all, family = poisson, 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

# test overdispersion
overdisp_fun(m2_final) # Not overdispersed

# test significance of fixed effects 
Anova(m2_final)

# Export estimates
require(broom)
out <- tidy(m2_final)

# ---14.8 Plot residuals: ----------------------------------------------------------------

simulationOutputq_final <- simulateResiduals(fittedModel = m2_final)
# Acces the qq plot
plotQQunif(simulationOutputq_final)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutputq_final)

#----14.9 Compare with model that has all the land-use intensities merge -----------------------------

# Merge all the intensity levels for all land-use categories 

# Create the vectors that hold the land-uses that we want to keep with different use intensities 
land_uses_separate_null1 <- "NA"
# Create a vector with the lad-uses that we want to merge
land_uses_light_intense_null1 <- "NA"


# This function creates new columns with the land-uses and use intensities that we want, then it then merges both columns 
# into a LanUse column
diversity_all <- Merge_landUses_and_intensities(diversity_all, 
                                                1,
                                                land_uses_separate_null1, 
                                                land_uses_light_intense_null1,
                                                "Primary All")
# Check levels
levels(diversity_all$LandUse.1)


# Model richness with those land-use classes
m3_final <- glmer(Species_richness ~ LandUse.1 + Kingdom + LandUse.1:Kingdom +
                  (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_all, family = poisson, 
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

# Check they have the same sample size with the complex model
summary(m3_final)

# Compare AICs
AIC(m2_final, m3_final)

# Next we will check if we've lost a significant amount of explanatory power 
# by removing this interaction. If we have, we want to keep the more complex
# model. If we haven't lost a significant amount of explanatory power, then we
# can keep the simpler model. 

anova(m2_final, m3_final)

#--15. Compare with null model of plants -----------------------------------------------

# For the null model I am going to replace the nePlantae for Plantae 
diversity_all <- diversity_all %>%
  
  
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
levels(diversity_all$Kingdom.1)

m4_final <- glmer(Species_richness ~ LandUse.0 + Kingdom.1 + LandUse.0:Kingdom.1 +
                  (1|Source_ID) + (1|SS) + (1|SSB) + (1|SSBS), data = diversity_all, family = poisson, 
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

# Check they have the same sample size with the complex model
summary(m4_final)

# check if we've lost a significant amount of explanatory power 
anova(m2_final, m4_final, test = "F")


# ---16. Plot results: ----------------------------------------------------------------

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of the plot labels
PlotErrBar_interactions(model = m2_final, resp = "Species richness", Effect1 = "LandUse.0", Effect2 = "Kingdom",
                        ylims = c(-3,2.5), pointtype = c(16,17,18),blackwhite = FALSE)


# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m2_final,
                             resp = "Species Richness",
                             Effect1 = "LandUse.0", 
                             Effect2 = "Kingdom",
                             ylims = c(-3,2.5),
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



