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

# So, in this first attempt I'm just going to model (log(Abundance) + 1), without any interactions 
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
# diversity metric type of occurrence. 

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

hist(model_data5$RescaledAbundance) # Bound to zero

# First, let’s transform RescaledAbundance.

model_data5 <- mutate(model_data5, 
                     logAbundance = log(RescaledAbundance + 1),
                     sqrtAbundance = sqrt(RescaledAbundance)
)

hist(model_data5$logAbundance) 
hist(model_data5$sqrtAbundance)

m1 <- lmer(logAbundance ~ Predominant_land_use + Use_intensity + 
             Predominant_land_use*Use_intensity + 
             (1|SS) + (1|SSB), data = model_data5)

summary(m1)










###########################################################################################################

 
# Coarsen factor levels
model_data1 <- model_data %>% filter(Predominant_land_use != "Urban") %>% droplevels() %>%
  mutate(Use_intensity = str_replace_all(Use_intensity, pattern = c("Intense use" = "Intense light use",
                                                "Light use" = "Intense light use")), 
         Use_intensity = factor(Use_intensity),
         Use_intensity = relevel(Use_intensity, ref = "Minimal use")) %>% filter(sqrtAbundance != 0) %>% droplevels()

table(model_data1$Predominant_land_use, model_data1$Use_intensity)

m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + 
             Predominant_land_use:Use_intensity + 
             (1|SS) + (1|SSB), data = model_data1)





