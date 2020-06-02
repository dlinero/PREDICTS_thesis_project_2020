rm(list = ls())
# Install some functions specifically developed for handling PREDICTS data.

# install the PREDICTS yarg package
install.packages("yarg_0.1-14.tar.gz", repos = NULL, type = "source")
install.packages("glmmADMB", repos = "http://R-Forge.R-project.org")
install.packages("roquefort_0.1-2.tar.gz", repos = NULL, type = "source")

# Load some packages for manipulating and modelling the data.

library(yarg) # useful functions for dealing with PREDICTS data
library(roquefort) # useful PREDICTS functions, particularly for plotting
library(raster) # for dealing with spatial data
library(dplyr) # for handy data manipulation functions
library(tidyr) # ditto
library(lme4) # for mixed effects models
library(car) # for getting anova tables with significance values
library(DHARMa) # for model criticism plots
library(MuMIn) # for checking explanatory power of mixed effects models
library(stringr)

# If you’re working with a raw extract from the PREDICTS database, it will have diversity and the date it was extracted in the name.

# ---- 1. Load data ------------------------------------------------------------------

diversity <- readRDS("./output/cleaned_data/01_Filter_data_PREDICTS_Filtered_table.rds")

# ---- 2. Correct abundance measures using Sampling effort 

# The next thing is correct for differences in sampling effort across sites within a single study.
# This rescales sampling effort for each study to have a maximum value of 1, and then divides any 
# diversity measurements that are sensitive to sampling effort by this rescaled measure of relative effort.

# What the CorrectSamplingEffort function does is grouping the dataset by SS, findind the maximum value
# of Sampling effort, dividing every sampling effort value by the maximum value within each study. This gives the
# rescaled sampling effort. Then it takes the measurement and divides it by the rescaled sampling effort. 

diversity1 <- yarg::CorrectSamplingEffort(diversity) 

# Export_table
saveRDS(diversity1, "output/intermediate_files/02_Analysis_walkthrough_Effort_corrected_measures")

# ---- 3. Merge Sites -------------------------------------------------------------

# Next we’ll merge any sites that are within the same land-use type and that have identical
# coordinates, start and end dates.

# We merge sites within studies that have identical coordinates, start and end dates(and land-use 
# type and intensity). We do this because sometimes authors record different points in a transect as
# different sites, which might not be meaningful if they share land-use, intensity and coordinates (the
# coordinates have an error) OR maybe it is, but we are having a conservative approach

diversity2 <- yarg::MergeSites(diversity1,
                              silent = TRUE,
                              merge.extra = "Wilderness_area")

# Export table
saveRDS(diversity2, "output/intermediate_files/02_Analysis_walkthrough_Merged_sites")

# ----4. Rename Predominant habitat --------------------------------------------------

# Rename the Predominant_habitat column, since it’s not really “habitat” that we’re looking at.
# This will also keep things a little more consistent with the public version of the PREDICTS database.

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

# Check results
# Check the results of the function SiteMetrics for a specific study
diversity3 %>% filter(SSS == "CM1_2012__Katovai 1 43") %>% select(SSS, Total_abundance, Species_richness, 
                                                                  MaxAbundance, RescaledAbundance)

# Calculate total abundance and species richness for one study and see if it matches
# the results of the SiteMetrics function
diversity2 %>%
  filter(SSS == "CM1_2012__Katovai 1 43") %>%
  select(Taxon_name_entered, Measurement, SSS) %>%
  mutate(Total_Abundance = sum(Measurement), Species_richness = length(unique(Taxon_name_entered)))

# Export table
saveRDS(diversity3, file = "./output/intermediate_files/02_Analysis_walkthrough_Site_metrics.rds")

# Import table

diversity3 <- readRDS("./output/intermediate_files/02_Analysis_walkthrough_Site_metrics.rds")


# take a look at the LandUse/Use intensity split
table(diversity3$Predominant_land_use, diversity3$Use_intensity, diversity3$Kingdom)

# ----6. Tidy up land-use column ---------------------------------------------------------------

# The Predominant_land_use column holds the land use categories. They often need some tidying up.

diversity4 <- diversity3 %>%
  
  mutate(
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide get NA, urban too beacuse there are only 40 sites
    Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"),
    Predominant_land_use = na_if(Predominant_land_use, "Urban"),

    
    # set reference levels
    Predominant_land_use = factor(Predominant_land_use),
    Predominant_land_use = relevel(Predominant_land_use, ref = "Primary"),
    Use_intensity = factor(Use_intensity),
    Use_intensity = relevel(Use_intensity, ref = "Minimal use")
  )


# take a look at the LandUse/Use intensity split
table(diversity4$Predominant_land_use, diversity4$Use_intensity, diversity4$Kingdom)

########################################### FIRST ATTEMPT #################################

# So, in this first attempt I'm just going to model sqrt(abundance), without any interactions 
# with kingdom. I'm also going to join all the light and intense uses. 


diversity5 <- diversity4 %>%  mutate(Use_intensity = str_replace_all(Use_intensity, pattern = c("Intense use" = "Intense light use",
                                                                    "Light use" = "Intense light use")), 
         Use_intensity = na_if(Use_intensity, "Cannot decide"),
         Use_intensity = factor(Use_intensity),
         Use_intensity = relevel(Use_intensity, ref = "Minimal use")) 

table(diversity5$Predominant_land_use, diversity5$Use_intensity)

# Collinearity: ---------------------------------------------------------------

# Ok, so let’s check whether any of the candidate explanatory variables are strongly correlated 
# with one another. You can look at the linear correlations between variables (using cor.test)
# or correlations between continuous and categorical variables (using aov). However, we tend to
# have multiple categorical varialbes. A good alternative is to look at Generalized Variance 
# Inflation Factors. You can use a function provided by the Zuur book.

source("https://highstat.com/Books/Book2/HighstatLibV10.R")

corvif(diversity5[ , c("Predominant_land_use", "Use_intensity")])

# There are different ‘rules of thumb’ about how high is too high when it comes to GVIFs. 
# Under 3 is great and under 5 is ok. Some people also say under 10 is acceptable, but I 
# think that’s a bit too high personally.

# If collinearity is a problem, then you’ll need to consider which variables are most important
# to include and what ones can be dropped.

# Complete cases: ------------------------------------------------------------

# By dropping the sites that have a total abundance of NA, we are excluding the sites with a 
# diversity metric type equal to occurrence. 

model_data5 <- drop_na(diversity5, 
                      Total_abundance, Predominant_land_use, Use_intensity) %>% droplevels()

table(model_data5$Predominant_land_use, model_data5$Use_intensity)

# Starting maximal model: ----------------------------------------------------------

# We’ll start by modelling Total_abundance. The maximal model is the most complex model you want 
# to test. There’s a limit to how complicated this can or should be. Think carefully about what 
# variables you want to test and what variables you might need to control for, and think carefully
# about what interactions are likely to be important.

# If you’re analysing count data, you might also find that there are far too many zeros in your data
# than you’d expect given a poisson distribution. If this is the case, you will need to consider 
# using a zero-inflated poisson distribution or a hurdle model. Neither of these are available
# in lme4; if you want to use lme4, you can use a ‘modified’ hurdle model by first modelling 
# the probability of species occurrence and then modelling the abundance of species that are 
# present.

# Let's do some data explorations
hist(model_data5$RescaledAbundance) # Bound to zero

# First, let’s transform RescaledAbundance.

model_data5 <- mutate(model_data5, 
                     logAbundance = log(RescaledAbundance + 1),
                     sqrtAbundance = sqrt(RescaledAbundance)
)

hist(model_data5$logAbundance) 
hist(model_data5$sqrtAbundance)


ggplot(model_data5, aes(x=Predominant_land_use, y= RescaledAbundance, color= Use_intensity)) + 
  geom_boxplot()

# Now try some models

m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity +
             Predominant_land_use:Use_intensity + 
             (1|SS) + (1|SSB), data = model_data5)

# Choose the random effects: -------------------------------------------------------------

# To select the random-effects structure, we use the method recommended by (Zuur et al., 2009) 
# of taking the most complex fixed-effects structure, including all interactions, that will be 
# tested in the second stage of modelling, and use it to compare different random-effects 
# structures.

# Study identity is always included as a random intercept because of the differences in the 
# diversity metrics that will be caused by the fundamental differences in methods, sampling effort
# etc. among different studies. We also tend to include block as a random intercept, to account 
# for the spatial configuration of sites. We can also include random slopes within study, to 
# allow the effects of explanatory variables to vary among studies. You may also wish you 
# include Source as a random intercept.

m2 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity +
             Predominant_land_use:Use_intensity + 
             (1+Predominant_land_use|SS) + (1|SSB), data = model_data5)
# boundary (singular) fit: see ?isSingular
# isSingular(m2) == TRUE

m3 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + 
             Predominant_land_use:Use_intensity + 
             (1+Use_intensity|SS) + (1|SSB), data = model_data5)

# compare the models that converged using Akaike's Information Criterion (AIC)
AIC(m1, m2, m3)

# The lowest AIC value is the best model of the bunch. In this case, m2. However m2 has a warning
# so I'm going to choose m3. This has Study (SS) and Block (SSB) as random intercepts and a
# random slope of Use_intensity.
# This means that we’re allowing the effect of use_intensity on abundance to vary among studies.
# We’re not going to go into detail about this here.

# Choose the best fixed effect structure: -------------------------------------------------------------
# Let’s assess the fixed effects in this model.

Anova(m3)

# The Predominant_land_use:Use_intensity interaction is significant, so we can’t remove this term 
# from the model. This means that even if the main effects (LandUse and Use_intensity) weren’t
# significant, we wouldn’t want to remove them from the model.

# The same is true for the Predominant_land_use:loghpd interaction. But let’s imagine it wasn’t 
# significant. Let’s go through the model simplification process. We start with the most
# complicated, least significant effect and try to remove it from the model. In this case, it is
# Predominant_land_use:loghpd.

# m3.1 <- update(m3,~.- Predominant_land_use:loghpd)

# Next, we will check if we’ve lost a significant amount of explanatory power by removing this 
# interaction. If we have, we want to keep the more complex model. If we haven’t lost a 
# significant amount of explanatory power, then we can keep this simpler model. 
# The test we’ll use is a likelihood ratio test (LRT).

# anova(m3.1, m3)

# Note that the Likelihood Ratio Test refits the models with Maximimum Likelihood. 
# This is best when testing different fixed effects, but when testing different random effects,
# the models should be fit using Restricted Maximum Likelihood (REML, the default option).

# The LRT is significant, so by removing the interaction from the model, we have lost a significant 
# amount of explanatory power. This means that we should keep the more complex model. 
# If it wasn’t significant, we would go with the simpler model and look at the Anova table
# again and remove the next most complicated, least significant term and repeat the process
# until everything left in the model is significant (this is backwards stepwise model simplification).
# The remaining model is our minimum adequate model.

# Interpretation: -------------------------------------------------------------

# Now let’s look at the model estimates of our mimumum adequate model 
# (which in our case is also our maximal model).

summary(m3)

# Plot residuals: ----------------------------------------------------------------

plot(m3)
#  for a correctly specified model we would expect assymptotically:
# *a uniform (flat) distribution of the scaled residuals
# * uniformity in y direction if we plot against any predictor.

# Simulate the residuals and plot them
simulationOutput <- simulateResiduals(fittedModel = m3, plot = T)

# Acces the qq plot
plotQQunif(simulationOutput)
# Plot the residuals against the predicted value 
plotResiduals(simulationOutput)

# Plot the results: -------------------------------------------------------------


roquefort::PlotErrBar(model = m3,
                      data = model_data5,
                      responseVar = "sqrt(Rescaled Abundance)",
                      seMultiplier = 1.96,
                      secdAge = TRUE,
                      logLink = "n",
                      catEffects = c("Predominant_land_use"),
                      forPaper = TRUE,
                      plotLabels = FALSE)


# Note that we usually collapse Land Use and Use intensity into a single factor - this makes plotting
# the results a little simpler, but it can make simplifying the model a little more complicated.

# You can also use PlotContEffects in the roquefort package to plot PREDICTS models to see how human
# population density responds to land-use change.

###########################################################################################################

 

# with log
m12 <- lmer(logAbundance ~ Predominant_land_use + Use_intensity +
             Predominant_land_use:Use_intensity + 
             (1|SS) + (1|SSB), data = model_data5)
m22 <- lmer(logAbundance ~ Predominant_land_use + Use_intensity +
             Predominant_land_use:Use_intensity + 
             (1+Predominant_land_use|SS) + (1|SSB), data = model_data5)
# boundary (singular) fit: see ?isSingular
# isSingular(m2) == TRUE

m32 <- lmer(logAbundance ~ Predominant_land_use + Use_intensity + 
             Predominant_land_use:Use_intensity + 
             (1+Use_intensity|SS) + (1|SSB), data = model_data5)

# compare the models that converged using Akaike's Information Criterion (AIC)
AIC(m12, m22, m32)

Anova(m32)

summary(m32)
plot(m32)

# Simulate the residuals and plot them
simulationOutput <- simulateResiduals(fittedModel = m32, plot = T)

# Acces the qq plot
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

roquefort::PlotErrBar(model = m32,
                      data = model_data5,
                      responseVar = "log Abundance",
                      seMultiplier = 1.96,
                      secdAge = TRUE,
                      logLink = "n",
                      catEffects = c("Predominant_land_use", "Use_intensity"),
                      forPaper = TRUE,
                      plotLabels = FALSE)
