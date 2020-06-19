
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

# ---- Description --------------------------------------------------------------------

# Model rescaled total abundance against land use type and intensity

# ---- 1. Load data ------------------------------------------------------------------

diversity <- readRDS("./output/cleaned_data/01_Filter_data_PREDICTS_Filtered_table.rds")

# ---- 2. Correct abundance measures using Sampling effort ---------------------------------

# The CorrectSamplingEffort function groups the dataset by SS and :
# 1. Finds the maximum value of Sampling effort
# 2. Divides every sampling effort value by the maximum value within each study. 
# 3. Then it takes the abundance measure and divides it by the rescaled sampling effort. 

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
# Since I filtered specific species, there are sites where those species had an abundance
# of zero and the maximum abundance within a site was also zero.
# I am going to exclude those values from further analyzes
diversity4 <- diversity3[-which(is.nan(diversity3$RescaledAbundance)),]

# Export table
saveRDS(diversity4, file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_frugi_endooPlants.rds")

# Import table
diversity4 <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_frugi_endooPlants.rds")

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
    
    # indeterminate secondary veg and cannot decide are transformed into NA, urban too because it has only 40 sites. 
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

# Check number of studies /// I think this code is not working properly
SS_second_model2 <- second_model_data2 %>% group_by(Kingdom, LandUse) %>% 
  summarise(Number_studies = length(unique(SS))) %>%  
  ungroup() 

# Total number of studies
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


# compare the models that converged using Akaike's Information Criterion (AIC)
AIC(m2_1, m2_2)

# anova-like table of random effects via likelihood ratio tests
ranova(m2_2)
ranova(m2_3)

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

# 9. --- Model with kingdom interaction simplification ---------------------------------------------------------------------------

Third_model_data <- diversity4 %>%
  # make a level of Primary minimal
  mutate(
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide are transformed into NA, urban too because it has only 40 sites. 
    Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"),
    Predominant_land_use = na_if(Predominant_land_use, "Urban"),
    
    
    # Give a shorter name to some land use categories 
    Predominant_land_use = str_replace_all(Predominant_land_use, pattern = c("Young secondary vegetation" = "YSV",
                                                                             "Intermediate secondary vegetation" = "ISV", 
                                                                             "Mature secondary vegetation" = "MSV")),
    
   
     # Merge all the intensity levels for all land-use categories 
    Use_intensity = ifelse(Use_intensity != "NA", "All", "NA"),
    
    # Paste the land-use classes and intensity levels
    LandUse = ifelse(Predominant_land_use != "NA" & Use_intensity != "NA",
                     paste(Predominant_land_use, Use_intensity),
                     NA),
    
    # set reference level
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary All")
  )


# ---9.1. Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(Third_model_data[ , c("LandUse", "Kingdom")])

# ---9.2. Complete cases --------------------------------------------------------------------------------

# Create a table with complete cases
Third_model_data2 <- drop_na(Third_model_data, 
                              Total_abundance, LandUse) %>% droplevels()

# Check number of sites
table(Third_model_data2$LandUse, Third_model_data2$Kingdom)

# Check number of studies /// I think this code is not working properly 
SS_third_model<- Third_model_data2 %>% group_by(Kingdom, LandUse) %>% 
  summarise(Number_studies = length(unique(SS))) %>%  
  ungroup() 

# Total number of studies
sum(SS_third_model$Number_studies)

# ---9.3. Choose between GLMM or LMM -----------------------------------------------------------------

# Distribution of the original variable
hist(Third_model_data2$RescaledAbundance) 

# Transform RescaledAbundance.
Third_model_data2 <- mutate(Third_model_data2, 
                             logAbundance = log(RescaledAbundance + 1),
                             sqrtAbundance = sqrt(RescaledAbundance)
)

# Plot histograms of the rescaled abundance
ggplot(Third_model_data2, aes(x=LandUse, y= RescaledAbundance, color= Kingdom)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 20))

# ---9.4. Choose random effects structure  -----------------------------------------------------------------

# To select the random-effects structure, we use the method recommended by (Zuur et al., 2009) 
# of taking the most complex fixed-effects structure, including all interactions, that will be 
# tested in the second stage of modelling, and use it to compare different random-effects 
# structures.

m3_1 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB), data = Third_model_data2)

m3_2 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = Third_model_data2)

m3_3 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+LandUse|SS) + (1|SSB), data = Third_model_data2)
# IsSingular

m3_4 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+LandUse|SS) + (1|SSB) + (1|Source_ID), data = Third_model_data2)
# IsSingular

m3_5 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+Predominant_land_use|SS) + (1|SSB), data = Third_model_data2)
# IsSingular

m3_6 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+Predominant_land_use|SS) + (1|SSB) + (1|Source_ID), data = Third_model_data2)
# IsSingular

# Compare the models that converged
AIC(m3_1, m3_2)

# anova-like table of random effects via likelihood ratio tests
ranova(m3_2)

# ---9.5. Choose fixed effects structure  -----------------------------------------------------------------

# See the significance of the terms in the model.
Anova(m3_2)

# The Land_use:kingdom interaction is significant, so I won’t remove this term 
# from the model. 

# ---9.6. Model estimates  -----------------------------------------------------------------
summary(m3_2)


# Export summary
require(broom)
out <- tidy(m3_2)

# ---9.7. Plot residuals: ----------------------------------------------------------------

simulationOutput3 <- simulateResiduals(fittedModel = m3_2)
# Acces the qq plot
plotQQunif(simulationOutput3)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput3)

# ---9.8. Plot results: ----------------------------------------------------------------

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

par(mar = c(10, 2, 2, 2))
PlotErrBar_interactions(model = m3_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                             ylims = c(-0.4,0.4), pointtype = c(16,17),blackwhite = FALSE)



# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m3_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                             ylims = c(-0.4,0.4), pointtype = c(16,17),blackwhite = FALSE)

# Plot the x label
text(x = c(0.8, 1.8, 2.8, 3.8, 4.7, 5.7, 6.8), y = -0.38, labels = c("Primary", 
                                            "Cropland", 
                                            "ISV",                                          
                                            "MSV",
                                            "Pasture",
                                            "Plantation forest",
                                            "YSV all"), srt = 15, cex= 0.9)

# --9.9. Run the models with log: ----------------------------------------------------------------

m3l_1 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
                (1|SS) + (1|SSB), data = Third_model_data2)

m3l_2 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
                (1|SS) + (1|SSB) + (1|Source_ID), data = Third_model_data2)

AIC(m3l_1, m3l_2)

simulationOutput3l <- simulateResiduals(fittedModel = m3l_2)
# Acces the qq plot
plotQQunif(simulationOutput3l)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput3l)

# The residual plots look better with sqrt.

# Plot results
PlotErrBar_interactions_modi(model = m3l_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                             ylims = c(-0.2,0.2), pointtype = c(16,17),blackwhite = FALSE)

# Plot the x label
text(x = c(0.8, 1.8, 2.8, 3.8, 4.7, 5.7, 6.8), y = -0.2, labels = c("Primary", 
                                                                     "Cropland", 
                                                                     "ISV",                                          
                                                                     "MSV",
                                                                     "Pasture",
                                                                     "Plantation forest",
                                                                     "YSV all"), srt = 15, cex= 0.9)

# 10. Choose between mergin all use -intensities -------------------------------------------------------
# for all land-uses or only for those that do not have enough sites for each combination

# I have to make the models with only the rows of data that are the same between both tables. 

# make a copy of the database: 
WH_merged <- second_model_data2

# Model with separated land-use intensities
m2_2 <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = WH_merged)
summary(m2_2)
# Same model but with log 
m2_2l <- lmer(logAbundance~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = WH_merged)


# Check overdispersion
simulateResiduals_m2_2 <- simulateResiduals(fittedModel = m2_2) 
testDispersion(simulateResiduals_m2_2)

# Merge land-use intensities for the second_model_data2
WH_merged <- WH_merged %>% mutate(
  
  # Merge all the intensity levels for all land-use categories 
  Use_intensity = ifelse(Use_intensity != "NA", "All", "NA"),

  # Paste the land-use classes and intensity levels
  LandUse = ifelse(Predominant_land_use != "NA" & Use_intensity != "NA",
                 paste(Predominant_land_use, Use_intensity),
                 NA),

# set reference level
LandUse = factor(LandUse),
LandUse = relevel(LandUse, ref = "Primary All")
)

# Check number of sites
table(second_model_data2_merged$LandUse, second_model_data2_merged$Kingdom)

# Model with merged land-uses
m2_2_m <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = WH_merged)
summary(m2_2_m)

# same model but with log
m2_2_ml <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
                 (1|SS) + (1|SSB) + (1|Source_ID), data = WH_merged)

# Check overdispersion
simulateResiduals_m2_2_m <- simulateResiduals(fittedModel = m2_2_m) 
testDispersion(simulateResiduals_m2_2_m)


# Test what model is better
# Akaike Information Criteria
AIC(m2_2, m2_2_m)

# Next, we will check if we’ve lost a significant amount of explanatory 
# power by removing this interaction. If we have, we want to keep the more
# complex model. If we haven’t lost a significant amount of explanatory power,
# then we can keep this simpler model. The test we’ll use is a likelihood ratio
# test (LRT).

anova(m2_2_m, m2_2, test = "F")

# The LRT is not significant, so by removing the interaction from the model, we didn't
# lost a significant amount of explanatory power. This means that we should keep the more 
# simpler model. 

AIC(m2_2l, m2_2_ml)
anova(m2_2l, m2_2_ml, test = "F")

# 11. Models including plants that are not endozoochorus -------------------------------------------------- 

# Load PREDICTS data that has been filtered to include:
# * Only tropical records
# * Only plants that do not depend on animals for their dispersal

# Import
PREDICTS_not_endoo_plants <- readRDS(file = "./output/cleaned_data/01_Filter_data_Plants_not_dispersed_by_animals.rds")

# ---- 11.2. Correct abundance measures using sampling effort ---------------------------------

# The CorrectSamplingEffort function groups the dataset by SS and :
# 1. Finds the maximum value of Sampling effort
# 2. Divides every sampling effort value by the maximum value within each study. 
# 3. Then it takes the measurement and divides it by the rescaled sampling effort. 

PREDICTS_not_endoo_plants1 <- yarg::CorrectSamplingEffort(PREDICTS_not_endoo_plants) 

# ---- 11.3. Merge Sites -------------------------------------------------------------

# We merge sites within studies that have identical coordinates, start and end dates(and land-use 
# type and intensity). We do this because sometimes authors record different points in a transect as
# different sites, which might not be meaningful if they share land-use, intensity and coordinates (the
# coordinates have an error) OR maybe it is, but we are having a conservative approach

PREDICTS_not_endoo_plants2 <- yarg::MergeSites(PREDICTS_not_endoo_plants1,
                               silent = TRUE,
                               merge.extra = "Wilderness_area")

# ----11.4. Rename Predominant habitat --------------------------------------------------

# Rename the column predominant habitat, as the dataset is actually refering to land use
PREDICTS_not_endoo_plants2 <- rename(PREDICTS_not_endoo_plants2,
                     Predominant_land_use = Predominant_habitat)


plants <- PREDICTS_not_endoo_plants2 %>% subset(Diversity_metric_type == "Abundance") %>% droplevels()
levels(plants$Diversity_metric_unit)

# ----11.5.  Calculate diversity metrics -----------------------------------------------

PREDICTS_not_endoo_plants3 <- PREDICTS_not_endoo_plants2 %>%
  
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

# Export table
saveRDS(PREDICTS_not_endoo_plants3, file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_not_endoo_plants.rds")

# Import table
PREDICTS_not_endoo_plants3 <- readRDS("./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_not_endoo_plants.rds")

# Check number of sites
addmargins(table(PREDICTS_not_endoo_plants3$Predominant_land_use, PREDICTS_not_endoo_plants3$Use_intensity), 2)


############################## MODEL TESTING INCLUDING MORE PLANTS #################################

# ---12. Simple model ----------------------------------------------------------------------

# What I'm going to do is a Plant contrast, that means I'm firts going to model the complex model
# (the one that assumme the reponse between endo and not-endo plants is significantly different), in this model 
# I will  put the Kingdom interaction that has 3 levels: animals and endo plants and not endo plants. Then I am going
# to model the null model that assumes there is no difference between endo and not-endo plants. In this model
# the Kingdom will have 2 levels: animals, plants. Then I will
# compare both models using an ANOVA

# 12.1. Checking number of sites -------------------------------------------------------------

# I will fisrt explore the number of sites that have a known land use category for the plants not dispersed
# by animals

NotEndoo_plants_sites <- PREDICTS_not_endoo_plants3 %>%
  # make a level of Primary minimal
  mutate(
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide are transformed into NA, urban too 
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
NotEndoo_plants_sites1 <- drop_na(NotEndoo_plants_sites, 
                              Total_abundance, Predominant_land_use) %>% droplevels()


# Check the number of sites 
addmargins(table(NotEndoo_plants_sites1$Predominant_land_use, NotEndoo_plants_sites1$Use_intensity), 2)


# 12.2. Merge the metrics of plants not dispersed by animals, plants dispersed by endozoochory and
# animals ------------------------------------------------------------------------------------------

# Add an identifier to differentiate between endo and not-endo plants
PREDICTS_not_endoo_plants3[, "Kingdom"] <- "nePlantae"
# Merge the site metrics tables
diversity5 <- rbind(diversity4, PREDICTS_not_endoo_plants3)

# Export table 
saveRDS(diversity5, file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")

# Import table 
diversity5 <- readRDS( file = "./output/cleaned_data/02_Statistical_Analysis_Abundance_Site_metrics_combined_animals_endooPlants_notEndooPlants.rds")

# Since the not_endoo plants have enough sites to separate primary into light, minimal and intense, and 
# also YSV and plantation forest into minimal and light/intense. I will try to run in with
# those classes and then compare if it's better to merge everything.

# Create a vector with the land uses categories that have enough sites to separate them 
# into different intensities
LandUse_divide <- c("Primary", "ISV", "Plantation forest")

fourth_model_data <- diversity5 %>%
  # make a level of Primary minimal
  mutate(
    
    # Create a new column of predominant habitat, to not overwrite the original
    Predominant_land_use_copy = paste(Predominant_land_use),
    
    # Create a new column of use intensity, to not overwrite the original
    Use_intensity_copy = paste(Use_intensity),
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use_copy = recode_factor(Predominant_land_use_copy, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide are transformed into NA, urban too because it has only 40 sites. 
    Predominant_land_use_copy = na_if(Predominant_land_use_copy, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use_copy = na_if(Predominant_land_use_copy, "Cannot decide"),
    Predominant_land_use_copy = na_if(Predominant_land_use_copy, "Urban"),
    
    
    # Give a shorter name to some land use categories 
    Predominant_land_use_copy = str_replace_all(Predominant_land_use_copy, pattern = c("Young secondary vegetation" = "YSV",
                                                                             "Intermediate secondary vegetation" = "ISV", 
                                                                             "Mature secondary vegetation" = "MSV")),
    
    # Drop the Cannot decide intensity levels for the land-use categories that have 
    # enough sites for the minimal, and light/intense
    Use_intensity_copy = ifelse(Predominant_land_use_copy %in% LandUse_divide & Use_intensity_copy == "Cannot decide",
                           NA,
                           paste(Use_intensity_copy)), 
    
    # Join the intensity levels of light and intense for Plantation forest and ISV 
    Use_intensity_copy = ifelse((Predominant_land_use_copy == "Plantation forest" | Predominant_land_use_copy == "ISV") & (Use_intensity_copy == "Intense use" | Use_intensity_copy == "Light use"),
                           str_replace_all(Use_intensity_copy, pattern = c("Intense use" = "Light-intense use", "Light use" = "Light-intense use")), paste(Use_intensity_copy)),
    
    # Merge all the intensity levels for those land-use categories that don't have enough sites in each land-use type/intensity combination
    Use_intensity_copy = ifelse(Predominant_land_use_copy %nin% LandUse_divide,
                           str_replace_all(Use_intensity_copy, pattern = c("Intense use" = "All", "Light use" = "All", "Minimal use" = "All", "Cannot decide" = "All")), 
                           paste(Use_intensity_copy)),
    
    
    # Paste the land-use classes and intensity levels
    LandUse = ifelse(Predominant_land_use_copy != "NA" & Use_intensity_copy != "NA",
                     paste(Predominant_land_use_copy, Use_intensity_copy),
                     NA),
    
    # set reference level
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary Minimal use")
    
  )

# Check landuse levels
levels(fourth_model_data$LandUse)

# ---12.3. Test for collinearity  ----------------------------------------------------------------------

# Since I'm going to explore the collinearity between categorical variables, I'm going to use 
# the Generalized variance Inflation Factors function provided by Zuur et al., (2009)

# Get the function
source("https://highstat.com/Books/Book2/HighstatLibV10.R")

# Calculate the VIF
corvif(fourth_model_data[ , c("LandUse", "Kingdom")])

# --- 12.4 Complete cases --------------------------------------------------------------------------------

# Drop sites that don't have abundance measures or land-use data
fourth_model_data1 <- drop_na(fourth_model_data, 
                              Total_abundance, LandUse) %>% droplevels()


# Check the number of sites for animals and plants
addmargins(table(fourth_model_data1$LandUse, fourth_model_data1$Kingdom), 2)

# ---12.5 Choose between GLMM or LMM -----------------------------------------------------------------

# Distribution of the original variable
hist(fourth_model_data1$RescaledAbundance) 

# Transform RescaledAbundance.
fourth_model_data2 <- mutate(fourth_model_data1, 
                            logAbundance = log(RescaledAbundance + 1),
                            sqrtAbundance = sqrt(RescaledAbundance)
)

# Plot histograms of the rescaled abundance
ggplot(fourth_model_data2, aes(x=LandUse, y= RescaledAbundance, color= Kingdom)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 20))

# ---12.6 Choose random effects structure  -----------------------------------------------------------------

# To select the random-effects structure, we use the method recommended by (Zuur et al., 2009) 
# of taking the most complex fixed-effects structure, including all interactions, that will be 
# tested in the second stage of modelling, and use it to compare different random-effects 
# structures.

m4_1 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB), data = fourth_model_data2)

m4_2 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = fourth_model_data2)

m4_3 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+LandUse|SS) + (1|SSB), data = fourth_model_data2)
# IsSingular

m4_4 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+LandUse|SS) + (1|SSB) + (1|Source_ID), data = fourth_model_data2)
# IsSingular

m4_5 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+Predominant_land_use|SS) + (1|SSB), data = fourth_model_data2)
# IsSingular

m4_6 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1+Predominant_land_use|SS) + (1|SSB) + (1|Source_ID), data = fourth_model_data2)
# IsSingular

# Compare the models that converged
AIC(m4_1, m4_2)

# anova-like table of random effects via likelihood ratio tests
ranova(m4_2)

# ---12.7. Choose fixed effects structure  -----------------------------------------------------------------

# See the significance of the terms in the model.
Anova(m4_2)

# The Land_use:kingdom interaction is significant, so I won’t remove this term 
# from the model. 

# ---12.8. Model estimates  -----------------------------------------------------------------
summary(m4_2)

# Export summary
require(broom)
out <- tidy(m4_2)

# ---12.9. Plot residuals: ----------------------------------------------------------------

simulationOutput4 <- simulateResiduals(fittedModel = m4_2)
# Acces the qq plot
plotQQunif(simulationOutput4)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput4)

# ---12.10. Plot results: ----------------------------------------------------------------

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of the plot labels
PlotErrBar_interactions(model = m4_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                        ylims = c(-0.4,0.4), pointtype = c(16,17,18),blackwhite = FALSE)



# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m4_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                             ylims = c(-0.25,0.4), pointtype = c(16,17, 18),blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), y = -0.18, labels = c("Primary minimal use",
                                            "Cropland all uses",
                                            "ISV light-intense use",
                                            "ISV minimal use",
                                            "MSV all uses",
                                            "Pasture all uses",
                                            "Plantation forest light-intense use",
                                            "Plantation forest minimal use",
                                            "Primary intense use",
                                            "Primary light use",
                                            "YSV all uses"), srt = 18, cex= 0.7)

# ---12.11. Merge all land use intensities: ----------------------------------------------------------------

# Check overdispersion of the first model
testDispersion(simulationOutput4)

# Merge land-use intensities for the fourth_model_data2
fourth_model_data2 <- fourth_model_data2 %>% mutate(
  
  # Create a new column as the copy of the use intensities we used in the first
  # selected model
  Use_intensity_copy_2 = paste(Use_intensity_copy),
  
  # Merge all the intensity levels for all land-use categories 
  Use_intensity_copy_2 = ifelse(Use_intensity_copy_2 != "NA", "All", "NA"),

  # Paste the land-use classes and intensity levels in a new land use column
  LandUse_copy = ifelse(Predominant_land_use_copy != "NA" & Use_intensity_copy_2 != "NA",
                   paste(Predominant_land_use_copy, Use_intensity_copy_2),
                   NA),
  
  # set reference level
  LandUse_copy = factor(LandUse_copy),
  LandUse_copy = relevel(LandUse_copy, ref = "Primary All")
)

# Model with merged land uses
m4_2_ml <- lmer(logAbundance ~ LandUse_copy + Kingdom + LandUse_copy:Kingdom +
                  (1|SS) + (1|SSB) + (1|Source_ID), data = fourth_model_data2)

# Check overdispersion
simulateResiduals_m4_2_m <- simulateResiduals(fittedModel = m4_2_ml) 
testDispersion(simulateResiduals_m4_2_m)

# Next, we will check if we’ve lost a significant amount of explanatory 
# power by removing this interaction. If we have, we want to keep the more
# complex model. If we haven’t lost a significant amount of explanatory power,
# then we can keep this simpler model. The test we’ll use is a likelihood ratio
# test (LRT).

anova(m4_2_ml, m4_2, test = "F")

# The LRT is significant, so by removing the interaction from the model, we did
# lose a significant amount of explanatory power. This means that we should keep the more 
# complex model. 

# ---12.12. Compare with null model: ----------------------------------------------------------------

# ---12.12.1. Create null model: -----------------------------------------------

# For the null model I am going to replace the nePlantae for Plantae 
fourth_model_data2 <- fourth_model_data2 %>%
  
  
  mutate(
    
    # Create a new kingdom column as the copy of the kingdom column we used
    # in the first model
    Kingdom_copy = paste(Kingdom),
      
      # Replace nePlantae for Plantae
    Kingdom_copy = recode_factor(Kingdom, "nePlantae" = "Plantae"), 

# set reference level
LandUse = factor(LandUse),
LandUse = relevel(LandUse, ref = "Primary Minimal use"), 
Kingdom_copy = factor(Kingdom_copy),
Kingdom_copy = relevel(Kingdom_copy, ref = "Animalia"))

# Check the levels
levels(fourth_model_data2$Kingdom_copy)

# ---12.12.1. Compare complex and null models: -----------------------------------------------

m4_2_mk <- lmer(logAbundance ~ LandUse + Kingdom_copy + LandUse:Kingdom_copy +
               (1|SS) + (1|SSB) + (1|Source_ID), data = fourth_model_data2)
summary(m4_2_mk)

anova(m4_2_mk, m4_2, test = "F")

# ---12.13. Plot the results of the selected abundance model ----------------------------------------

# Selected model
m4_2 <- lmer(logAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
               (1|SS) + (1|SSB) + (1|Source_ID), data = fourth_model_data2)

# Selected model with sqrt
m4_2_sqrt <- lmer(sqrtAbundance ~ LandUse + Kingdom + LandUse:Kingdom +
              (1|SS) + (1|SSB) + (1|Source_ID), data = fourth_model_data2)

# Read r code from a file which contains the function to make the plot
source("./R/PlotErrBar_interactions.R")
source("./R/PlotErrBar_interactions_modified.R")

# Check the order of land use factors 
PlotErrBar_interactions(model = m4_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                        ylims = c(-0.4,0.4), pointtype = c(16,17,18),blackwhite = FALSE)



# Plot the differences between estimates 
PlotErrBar_interactions_modi(model = m4_2, resp = "Abundance", Effect1 = "LandUse", Effect2 = "Kingdom",
                             ylims = c(-0.25,0.25), pointtype = c(16,17,18),blackwhite = FALSE)

# Plot the x label
text(x = c(0.8:10.8), y = -0.18, labels = c("Primary minimal use",
                                            "Cropland all uses",
                                            "ISV light-intense use",
                                            "ISV minimal use",
                                            "MSV all uses",
                                            "Pasture all uses",
                                            "Plantation forest light-intense use",
                                            "Plantation forest minimal use",
                                            "Primary intense use",
                                            "Primary light use",
                                            "YSV all uses"), srt = 18, cex= 0.9)

