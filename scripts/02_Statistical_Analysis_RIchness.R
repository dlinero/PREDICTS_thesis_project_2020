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
library(aods3) # to check for overdispersion
library(lmerTest) # to get p-values for estimates in LMMs


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

# Filter the data needed to run the model
diversityS <- diversityS  %>% 

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
min(diversityS$N_species_sampled)

# Maximum number of species sampled
max(diversityS$N_species_sampled)

# Remove N_species_sampled_column
diversityS <- diversityS %>% select(-N_species_sampled)

# ----3. Remove sites that will produce NaN values --------------------------------------------
# Doing the abundance models, I identified some sites where (due to the filtering of specific species)
# recorded only zero abundances. I am going to remove those sites as they will produce 
# Non infinite values. 

diversityS <- diversityS %>%  subset(SSBS %nin% c("MJ1_2009__Lehouck 2 Fururu 10", 
                                                                          "MJ1_2009__Lehouck 2 Fururu 11", 
                                                                          "MJ1_2009__Lehouck 2 Macha 12", 
                                                                          "MJ1_2009__Lehouck 2 Macha 13")) %>%
  droplevels()

# ---- 4. Split datasets --------------------------------------------------------

# I am going to separate the records that belong to plants not dispersed by animals
diversityS_notEndo <- diversityS %>% subset(Kingdom == "nePlantae") %>% droplevels()

# Get the table without plants dispersed by animals
diversityS_frugi_endo <- diversityS %>% subset(Kingdom != "nePlantae") %>% droplevels()

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
diversityS_frugi_endo <-  rename(diversityS_frugi_endo,
                      Predominant_land_use = Predominant_habitat)

diversityS_notEndo <-  rename(diversityS_notEndo,
                                 Predominant_land_use = Predominant_habitat)


# ----7.  Calculate diversity metrics -----------------------------------------------

# Calculate diversity metrics for animals and endoozoocoric plants
diversity1S_frugi_endo <- diversityS_frugi_endo %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Kingdom")) 

# Calculate diversity metrics for plants not dispersed by animals
diversity1S_notendo <- diversityS_notEndo %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Kingdom")) 


# ---8. Check the results--------------------------------------------------------------------

# Check the results for studies that assessed abundance
diversityS_frugi_endo %>% 
  
  # subset the table with all of the species records
  subset(SSS == "BS1_2010__Page 1 1") %>%
  
  # Select species that have an abundance greater than zero
  filter(Measurement > 0) %>%
  
  # Create a column with the number of species present
  mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_frugi_endo %>% subset(SSS == "BS1_2010__Page 1 1") %>% select(SSS, Species_richness)

# Check the results for studies that assessed abundance
diversityS_notEndo %>% 
  
  # subset the table with all of the species records
  subset(SSS == "CM1_2007__MarinSpiotta 1 1") %>%
  
  # Select species that have an abundance greater than zero
  filter(Measurement > 0) %>%
  
  # Create a column with the number of species present
  mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_notendo %>% subset(SSS == "CM1_2007__MarinSpiotta 1 1") %>% select(SSS, Species_richness)

# Check the results with studies that occurrence
diversityS_frugi_endo %>% 
  
  # subset the table with all of the species records
  subset(SSS == "DL1_2011__Latta 2 1") %>%
  
  # Select species that are present in the site 
  filter(Measurement == 1) %>% 
  
  # Create a column with the number of species present
  mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_frugi_endo %>% subset(SSS == "DL1_2011__Latta 2 1") %>% select(SSS, Species_richness)

# Check the results with studies that occurrence
diversityS_notEndo %>% 
  
  # subset the table with all of the species records
  subset(SSS == "CM2_2012__Cicuzza 1 1") %>%
  
  # Select species that are present in the site 
  filter(Measurement == 1) %>% 
  
  # Create a column with the number of species present
  mutate(Species_richness = n_distinct(Taxon_name_entered)) %>%
  
  select(SSS, Taxon_name_entered, Measurement, Species_richness) 

# Compare with the SiteMetrics result
diversity1S_notendo %>% subset(SSS == "CM2_2012__Cicuzza 1 1") %>% select(SSS, Species_richness)

# ---9. Merge site metrics -------------------------------------------------------------------
# Merge the site metrics for all organisms

diversity_all <- rbind.data.frame(diversity1S_frugi_endo, diversity1S_notendo)

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

diversity_all <- diversity_all %>%
  # make a level of Primary minimal
  mutate(
    
    # Create a new land-use column to store the transformations we are going to make
    Predominant_land_use_1 = paste(Predominant_land_use), 
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use_1 = recode_factor(Predominant_land_use_1, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide are transformed into NA, urban too 
    Predominant_land_use_1 = na_if(Predominant_land_use_1, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use_1 = na_if(Predominant_land_use_1, "Cannot decide"),
    Predominant_land_use_1 = na_if(Predominant_land_use_1, "Urban"),
    
    
    # Give a shorter name to some land use categories 
    Predominant_land_use_1 = str_replace_all(Predominant_land_use_1, pattern = c("Young secondary vegetation" = "YSV",
                                                                             "Intermediate secondary vegetation" = "ISV", 
                                                                             "Mature secondary vegetation" = "MSV")),
    
    Predominant_land_use_1 = factor(Predominant_land_use_1),
    Predominant_land_use_1 = relevel(Predominant_land_use_1, ref = "Primary")
  )


# Drop sites that don't have richness measures or land-use data
diversity_all <- drop_na(diversity_all, 
                                  Species_richness, Predominant_land_use_1) %>% droplevels()

# Check number of sites
addmargins(table(diversity_all$Predominant_land_use_1, diversity_all$Use_intensity, diversity_all$Kingdom), 2)


# 10.2. Merging land-use intensities  -------------------------------------------------------------
# According to the number of sites, we can split the land-use types and intensities in:

# Primary can be divided into the three level intensities in all cases
#	Cropland can be divided in minimal use and light/intense use in animal, endo plant and non endo plant
# ISV can be divided in minimal use and light/intense use in animal, endo plant and non endo plant
# MSV has to be merged
# Pasture has to be merged 
# Plantation forest can be divided in minimal use and light/intense use in animal, endo plant and non endo plant
#	YSV has to be merged 

# Create a vector with the land uses categories that have enough sites to separate them 
# into different intensities
LandUse_divide_r <- c("Primary", "Cropland", "ISV", "Plantation forest")


diversity_all <- diversity_all %>%
  
  mutate(
    
    # Create a new column for Use intensity to not overwrite the original
    Use_intensity_1 = paste(Use_intensity),
    
    # Drop the Cannot decide intensity levels for the land-use categories that have 
    # enough sites for the minimal, and light/intense
    Use_intensity_1 = ifelse(Predominant_land_use_1 %in% LandUse_divide_r & Use_intensity_1 == "Cannot decide",
                           NA,
                           paste(Use_intensity_1)), 
    
    # Join the intensity levels of light and intense for cropland, plantation forest and ISV 
    Use_intensity_1 = ifelse((Predominant_land_use_1 == "Plantation forest" |
                              Predominant_land_use_1 == "ISV" |
                              Predominant_land_use_1 == "Cropland") & 
                             (Use_intensity_1 == "Intense use" | 
                                Use_intensity_1 == "Light use"),
                           str_replace_all(Use_intensity_1,
                                           pattern = c("Intense use" = "Light-intense use", 
                                                       "Light use" = "Light-intense use")), 
                           paste(Use_intensity_1)),
    
    # Merge all the intensity levels for those land-use categories that don't have enough sites in each land-use type/intensity combination
    Use_intensity_1 = ifelse(Predominant_land_use_1 %nin% LandUse_divide_r,
                           str_replace_all(Use_intensity_1, pattern = c("Intense use" = "All", 
                                                                      "Light use" = "All", 
                                                                      "Minimal use" = "All", 
                                                                      "Cannot decide" = "All")), 
                           paste(Use_intensity_1)),
    
    
    # Paste the land-use classes and intensity levels
    LandUse = ifelse(Predominant_land_use_1 != "NA" & Use_intensity_1 != "NA",
                     paste(Predominant_land_use_1, Use_intensity_1),
                     NA),
    
    # set reference level
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary Minimal use")
  )

# ---10.3. Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(diversity_all[ , c("LandUse", "Kingdom")])

# ---10.4. Complete cases --------------------------------------------------------------------------------

# Create a table with complete cases
diversity_all <- drop_na(diversity_all, 
                              Species_richness, LandUse) %>% droplevels()

# Check number of sites
table(diversity_all$LandUse, diversity_all$Kingdom)

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

m1 <- glmer(Species_richness ~ LandUse + Kingdom + LandUse:Kingdom +
              (1|SS) + (1|SSB), data = diversity_all, family = poisson)

# Warning message:
#In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.0399178 (tol = 0.002, component 1)

m1 <- glmer(Species_richness ~ LandUse + Kingdom + LandUse:Kingdom +
              (1|SS) + (1|SSB), data = diversity_all, family = poisson, 
            
            # I'm increasing the number of iterations as the model initially 
            # didn't converge
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))


m2 <- glmer(Species_richness ~ LandUse + Kingdom + LandUse:Kingdom +
              (1|SS) + (1|SSB) + (1|Source_ID), data = diversity_all, family = poisson, 
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

m3 <- glmer(Species_richness ~ LandUse + Kingdom + LandUse:Kingdom +
              (1+LandUse|SS) + (1|SSB), data = diversity_all, family = poisson, 
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
# Is singular

m4 <- glmer(Species_richness ~ LandUse + Kingdom + LandUse:Kingdom +
              (1+Predominant_land_use_1|SS) + (1|SSB), data = diversity_all, family = poisson, 
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
# Is singular

# Compare the models that converged
AIC(m1, m2)

# ---10.6 Choose fixed effects structure  -----------------------------------------------------------------

# See the significance of the terms in the model.
Anova(m2)

# The Land_use:kingdom interaction is significant, so I won’t remove this term 
# from the model. 

# ---10.7 Model estimates  -----------------------------------------------------------------
summary(m2)

# Export summary
require(broom)
out <- tidy(m4_2)

# ---10.8 Plot residuals: ----------------------------------------------------------------

simulationOutput <- simulateResiduals(fittedModel = m2)
# Acces the qq plot
plotQQunif(simulationOutput)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput)

# ---10.9. Test for overdispersion: ----------------------------------------------------------------

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


# ---10.10. Model with quasipoisson: -------------------------------------------------------------- 
# I am going to add an observation-level random effects

m2_q <- glmer(Species_richness ~ LandUse + Kingdom + LandUse:Kingdom +
                (1|SS) + (1|SSB) + (1|Source_ID) + (1|SSBS), data = diversity_all, family = poisson, 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

# ---10.11 Plot residuals: ----------------------------------------------------------------

simulationOutputq <- simulateResiduals(fittedModel = m2_q)
# Acces the qq plot
plotQQunif(simulationOutputq)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutputq)

overdisp_fun(m2_q)

# ---12.10. Plot results: ----------------------------------------------------------------

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of the plot labels
PlotErrBar_interactions(model = m2_q, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
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
