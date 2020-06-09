
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

# ---- 1. Load data ------------------------------------------------------------------

diversity <- readRDS("./output/cleaned_data/01_Filter_data_PREDICTS_Filtered_table.rds")

# ---- 2. Correct abundance measures using Sampling effort ---------------------------------

# The CorrectSamplingEffort function groups the dataset by SS and :
# 1. Finds the maximum value of Sampling effort
# 2. Divides every sampling effort value by the maximum value within each study. 
# 3. Then it takes the measurement and divides it by the rescaled sampling effort. 

diversity1 <- yarg::CorrectSamplingEffort(diversity) 

# ---- 3. Merge Sites -------------------------------------------------------------

# We merge sites within studies that have identical coordinates, start and end dates(and land-use 
# type and intensity). We do this because sometimes authors record different points in a transect as
# different sites, which might not be meaningful if they share land-use, intensity and coordinates (the
# coordinates have an error) OR maybe it is, but we are having a conservative approach

diversity2 <- yarg::MergeSites(diversity1,
                               silent = TRUE,
                               merge.extra = "Wilderness_area")

# ----4. Rename Predominant habitat --------------------------------------------------

# Rename the column predominant habitat, as the dataset is actually refering to land use
diversity2 <- rename(diversity2,
                     Predominant_land_use = Predominant_habitat)


# See abundance measures
animals <- diversity2 %>% subset(Diversity_metric_type == "Abundance" & Kingdom == "Animalia") %>% droplevels()
levels(animals$Diversity_metric_unit)

plants <- diversity2 %>% subset(Diversity_metric_type == "Abundance" & Kingdom == "Plantae") %>% droplevels()
levels(plants$Diversity_metric_unit)

# ----5.  Calculate diversity metrics -----------------------------------------------

diversity3 <- diversity2 %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Kingdom")) %>% 
  
  # calculate the maximum abundance within each study
  group_by(SS) %>%
  mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance",
                               max(Total_abundance),
                               NA)) %>%
  ungroup() %>%
  
  # now calculate the rescaled abundance (abundance divided by the maximum within each study)
  mutate(RescaledAbundance = ifelse(Diversity_metric_type == "Abundance",
                                    Total_abundance/MaxAbundance,
                                    NA))

# ----6.  Check non- infinite values -----------------------------------------------

# I am going to check for non-infinite values, as the function is giving a warning

check <- diversity3[which(is.nan(diversity3$RescaledAbundance)), ]
# I am going to exclude those values from further analyzes
diversity4 <- diversity3[-which(is.nan(diversity3$RescaledAbundance)),]

# Export table
saveRDS(diversity4, file = "./output/cleaned_data/02_Statistical_Analysis_Site_metrics.rds")

# Import table
diversity4 <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Site_metrics.rds")

############################## MODEL TESTING #################################

# ---7. Simple model ----------------------------------------------------------------------

# So, in this first attempt I'm just going to model sqrt(abundance), without any interactions 
# with kingdom. I'm also going to join all the light and intense uses. 

# Check the number of sites per land use and use intensity combination
table(diversity4$Predominant_land_use, diversity4$Use_intensity)

first_model_data <- diversity4 %>%
  # make a level of Primary minimal
  mutate(
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide are transformed into NA, urban too because it has only 40 sites. The use
    # intensity cannot decide will also equal to NA
    Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"),
    Predominant_land_use = na_if(Predominant_land_use, "Urban"),
    Use_intensity = na_if(Use_intensity, "Cannot decide"),
    
    # Join the intensity levels of light and intense. 
    Use_intensity = str_replace_all(Use_intensity, pattern = c("Intense use" = "Light-intense use",
                                                               "Light use" = "Light-intense use")),
    
    # Give a shorter name to some land use categories 
    Predominant_land_use = str_replace_all(Predominant_land_use, pattern = c("Young secondary vegetation" = "YSV",
                                            "Intermediate secondary vegetation" = "ISV", 
                                           "Mature secondary vegetation" = "MSV")),
    
    # Paste the land-use classes and intensity levels
    LandUse = ifelse(Predominant_land_use != "NA" & Use_intensity != "NA",
                     paste(Predominant_land_use, Use_intensity),
                     NA),
    
    
    # set reference level
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary Minimal use"),
  )

# Check the resulting levels of land use
levels(first_model_data$LandUse)


# ---7.1. Test for collinearity  ----------------------------------------------------------------------

# Since in this case I only have one explanatory variable, I don't need to test for collinearity.

# ---7.2. Complete cases --------------------------------------------------------------------------------

# I am going to exclude those sites that do not have a value for the response and explanatory variables.
# In this way, I will also exclude those sites that do not have a total abundance value because the diversity
# metric type was occurrence instead of abudance. 

first_model_data1 <- drop_na(first_model_data, 
                       Total_abundance, LandUse) %>% droplevels()

# Check the number of sites per land use and use intensity combination
table(first_model_data1$LandUse)

# I am going to drop the Pasture minimal use level, because the number of sites is < 30.
first_model_data1 <- first_model_data1 %>% subset(LandUse != "Pasture Minimal use") %>% droplevels()

# ---7.3. Choose between GLMM or LMM -----------------------------------------------------------------

# Abundance data usually has a nonnormal error distribution because it has a positive mean-variance
# relationship and zero-inflation (Purvis et al., 2018). Given that some abudance measures are not integers
# (some are relative abudance or densities), I am not going to model the abudance with a Poisson distribution,
# but I’m going to transform it in order to meet the assumptions of linear mixed models.


# Distribution of the original variable
hist(first_model_data1$RescaledAbundance) 

# Transform RescaledAbundance.
first_model_data1 <- mutate(first_model_data1, 
                      logAbundance = log(RescaledAbundance + 1),
                      sqrtAbundance = sqrt(RescaledAbundance)
)

# Plot histograms of the rescaled abundance
ggplot(first_model_data1, aes(x=LandUse, y= RescaledAbundance)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 20))

# ---7.4. Choose random effects structure  -----------------------------------------------------------------

# # To select the random-effects structure, we use the method recommended by (Zuur et al., 2009) 
# of taking the most complex fixed-effects structure, including all interactions, that will be 
# tested in the second stage of modelling, and use it to compare different random-effects 
# structures.

# Study identity is always included as a random intercept because of the differences in the 
# diversity metrics that will be caused by the fundamental differences in methods, sampling effort
# etc. among different studies. We also tend to include block as a random intercept, to account 
# for the spatial configuration of sites. We can also include random slopes within study, to 
# allow the effects of explanatory variables to vary among studies. 

m1 <- lmer(sqrtAbundance ~ LandUse +
             (1|SS) + (1|SSB), data = first_model_data1)

m2 <- lmer(sqrtAbundance ~ LandUse + 
             (1|SS) + (1|SSB) + (1|Source_ID), data = first_model_data1)

m3 <- lmer(sqrtAbundance ~ LandUse + 
             (1+LandUse|SS) + (1|SSB), data = first_model_data1)
# fit is singular

m4 <- lmer(sqrtAbundance ~ LandUse + 
             (1+LandUse|SS) + (1|SSB) + (1|Source_ID), data = first_model_data1)
# fit is singular 

m5 <- lmer(sqrtAbundance ~ LandUse + 
             (1+Predominant_land_use|SS) + (1|SSB), data = first_model_data1)
# fit is singular

m6 <- lmer(sqrtAbundance ~ LandUse + 
           (1+Predominant_land_use|SS) + (1|SSB) + (1|Source_ID), data = first_model_data1)

m7 <- lmer(sqrtAbundance ~ LandUse + 
             (1+Use_intensity|SS) + (1|SSB), data = first_model_data1)

m8 <- lmer(sqrtAbundance ~ LandUse + 
             (1+Use_intensity|SS) + (1|SSB) + (1|Source_ID), data = first_model_data1)


# compare the models that converged using Akaike's Information Criterion (AIC)
AIC(m1, m2, m7, m8)

# ---7.5. Choose fixed effects structure  -----------------------------------------------------------------

# See the significance of the terms in the model (that in this case is just one)
Anova(m8)
Anova

# ---7.6. Model estimates  -----------------------------------------------------------------
summary(m8)

# ---7.7. Plot residuals: ----------------------------------------------------------------

simulationOutput1 <- simulateResiduals(fittedModel = m8)
# Acces the qq plot
plotQQunif(simulationOutput1)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput1)

# ---7.8. Plot results: ----------------------------------------------------------------

roquefort::PlotErrBar(model = m8,
                      data = first_model_data1,
                      responseVar = "sqrt(Rescaled Abundance)",
                      seMultiplier = 1.96,
                      secdAge = FALSE,
                      logLink = "n",
                      catEffects = c("LandUse"),
                      forPaper = TRUE,
                      plotLabels = FALSE,
                      plotLandUses = FALSE,
                      xtext.srt = 25, 
                      ylim = c(-65, 60), 
                      order = c(1, 11, 7, 6, 5, 4, 13, 12, 11, 10, 8, 3, 2))

text(x = c(1:14), y = -63, labels = c("Primary minimal", 
                                      "Primary light-intense", 
                                      "MSV minimal",
                                      "MSV light-intense",
                                      "ISV minimal",
                                      "ISV light-intense",
                                      "YSV minimal",
                                      "YSV light-intense", 
                                      "Plantation forest minimal",
                                      "Plantation forest light-intense",
                                      "Pasture light-intense",
                                      "Cropland minimal",
                                      "Cropland light-intense"), srt = 20, cex= 0.7)

# 8. --- Model with kingdom interaction ---------------------------------------------------------------------------

# I will fisrt explore the number of sites that have a known land use category

second_model_data <- diversity4 %>%
  # make a level of Primary minimal
  mutate(
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide are transformed into NA, urban too because it has only 40 sites. The use
    # intensity cannot decide will also equal to NA
    Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"),
    Predominant_land_use = na_if(Predominant_land_use, "Urban"),
   
    
     # Give a shorter name to some land use categories 
    Predominant_land_use = str_replace_all(Predominant_land_use, pattern = c("Young secondary vegetation" = "YSV",
                                                                             "Intermediate secondary vegetation" = "ISV", 
                                                                             "Mature secondary vegetation" = "MSV")),
    
    Predominant_land_use = factor(Predominant_land_use),
    Predominant_land_use = relevel(Predominant_land_use, ref = "Primary"),
  )


# Drop sites that don't have abundance measures or land-use data
second_model_data1 <- drop_na(second_model_data, 
                             Total_abundance, Predominant_land_use) %>% droplevels()


# Check the number of sites for animals and plants
addmargins(table(second_model_data1$Predominant_land_use, second_model_data1$Use_intensity, second_model_data1$Kingdom), 2)

# Create a vector with the land uses categories that have enough sites to separate them 
# into different intensities
LandUse_divide <- c("Primary", "ISV", "Plantation forest")


# Merge all categories that don't have enough sites on their own 
second_model_data2 <- second_model_data1 %>%
  
  mutate(
    
    # Drop the Cannot decide intensity levels for the land-use categories that have 
    # enough sites for the minimal, and light/intense
    Use_intensity = ifelse(Predominant_land_use %in% LandUse_divide & Use_intensity == "Cannot decide",
                           NA,
                           paste(Use_intensity)), 
    
    # Join the intensity levels of light and intense for Plantation forest and ISV 
    Use_intensity = ifelse((Predominant_land_use == "Plantation forest" | Predominant_land_use == "ISV") & (Use_intensity == "Intense use" | Use_intensity == "Light use"),
      str_replace_all(Use_intensity, pattern = c("Intense use" = "Light-intense use", "Light use" = "Light-intense use")), paste(Use_intensity)),
    
    # Merge all the intensity levels for those land-use categories that don't have enough sites in each land-use type/intensity combination
    Use_intensity = ifelse(Predominant_land_use %nin% LandUse_divide,
                           str_replace_all(Use_intensity, pattern = c("Intense use" = "All", "Light use" = "All", "Minimal use" = "All", "Cannot decide" = "All")), 
                           paste(Use_intensity)),
    

    # Paste the land-use classes and intensity levels
    LandUse = ifelse(Predominant_land_use != "NA" & Use_intensity != "NA",
                     paste(Predominant_land_use, Use_intensity),
                     NA),
    
    # set reference level
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary Minimal use")
  )


# ---8.1. Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(second_model_data2[ , c("LandUse", "Kingdom")])

# ---8.2. Complete cases --------------------------------------------------------------------------------

# Create a table with complete cases
second_model_data2 <- drop_na(second_model_data2, 
                              Total_abundance, LandUse) %>% droplevels()

# Check number of sites
table(second_model_data2$LandUse, second_model_data2$Kingdom)

# Check number of studies
SS_second_model2 <- second_model_data2 %>% group_by(Kingdom, LandUse) %>% 
  summarise(Number_studies = length(unique(SS))) %>%  
  ungroup() 

sum(SS_second_model2$Number_studies)

# ---8.3. Choose between GLMM or LMM -----------------------------------------------------------------

# Distribution of the original variable
hist(second_model_data2$RescaledAbundance) 

# Transform RescaledAbundance.
second_model_data2 <- mutate(second_model_data2, 
                            logAbundance = log(RescaledAbundance + 1),
                            sqrtAbundance = sqrt(RescaledAbundance)
)

# Plot histograms of the rescaled abundance
ggplot(second_model_data2, aes(x=LandUse, y= RescaledAbundance, color= Kingdom)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 20))

# ---8.4. Choose random effects structure  -----------------------------------------------------------------

# To select the random-effects structure, we use the method recommended by (Zuur et al., 2009) 
# of taking the most complex fixed-effects structure, including all interactions, that will be 
# tested in the second stage of modelling, and use it to compare different random-effects 
# structures.

m2_1 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1|SS) + (1|SSB), data = second_model_data2)

m2_2 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1|SS) + (1|SSB) + (1|Source_ID), data = second_model_data2)

m2_3 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1+LandUse|SS) + (1|SSB), data = second_model_data2)
# IsSingular

m2_4 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1+LandUse|SS) + (1|SSB) + (1|Source_ID), data = second_model_data2)

m2_5 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1+Predominant_land_use|SS) + (1|SSB), data = second_model_data2)
# IsSingular

m2_6 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1+Predominant_land_use|SS) + (1|SSB) + (1|Source_ID), data = second_model_data2)

m2_7 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1+Use_intensity|SS) + (1|SSB), data = second_model_data2)
# IsSingular

m2_8 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
             (1+Use_intensity|SS) + (1|SSB) + (1|Source_ID), data = second_model_data2)


# compare the models that converged using Akaike's Information Criterion (AIC)
AIC(m2_1, m2_2)

# anova-like table of random effects via likelihood ratio tests
ranova(m2_2)

# ---8.5. Choose fixed effects structure  -----------------------------------------------------------------

# See the significance of the terms in the model.
Anova(m2_2)

# The Land_use:kingdom interaction is significant, so I won’t remove this term 
# from the model. 

# ---8.6. Model estimates  -----------------------------------------------------------------
summary(m2_2)


# Export summary
require(broom)
out <- tidy(m2_2)
# ---8.7. Plot residuals: ----------------------------------------------------------------

simulationOutput2 <- simulateResiduals(fittedModel = m2_2)
# Acces the qq plot
plotQQunif(simulationOutput2)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput2)

# ---8.8. Plot results: ----------------------------------------------------------------

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m2_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                        ylims = c(-0.4,0.4), pointtype = c(16,17),blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), y = -0.38, labels = c("Primary minimal use", 
                                      "Cropland all uses", 
                                      "ISV light-intense use",
                                      "ISV minimal use",
                                      "MSV all uses",
                                      "Pasture all uses",
                                      "Plantation forest light-intense use",
                                      "Plantation forest minimal use", 
                                      "Primary intense use",
                                      "Primary light use",
                                      "YSV all uses"), srt = 15, cex= 0.9)

# ---8.9. Run the models with log: ----------------------------------------------------------------

m2l_1 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB), data = second_model_data2)

m2l_2 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = second_model_data2)

AIC(m2l_1, m2l_2)

simulationOutput2l <- simulateResiduals(fittedModel = m2l_2)
# Acces the qq plot
plotQQunif(simulationOutput2l)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput2l)

# The residual plots look better with sqrt.

