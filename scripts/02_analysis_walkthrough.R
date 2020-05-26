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

# Im going to upload my filtered table

diversity <- readRDS("./output/cleaned_data/01_Filter_data_PREDICTS_Frugivores_and_Endoplants.rds")
diversity$Sample_start_earliest <- as.Date(diversity$Sample_start_earliest)
diversity$Sample_end_latest <- as.Date(diversity$Sample_end_latest)
diversity$Site_area <- as.numeric(diversity$Site_area)
diversity$Site_area_unit <- as.factor(diversity$Site_area_unit)

diversity_all <- readRDS("./data/PREDICTS2020/PREDICTS2020.rds")
diversity_all <- diversity_all[1:10000, ]

# The next thing is correct for differences in sampling effort across sites within a single study.
# This rescales sampling effort for each study to have a maximum value of 1, and then divides any 
# diversity measurements that are sensitive to sampling effort by this rescaled measure of relative effort.


# What the CorrectSamplingEffort function does is grouping the dataset by SS, findind the maximum value
# of Sampling effort, dividing every sampling effort value by the maximum value within each study. This gives the
# rescaled sampling effort. The it takes the measurement and divides it by the rescaled sampling effor. 

diversity <- yarg::CorrectSamplingEffort(diversity) 
diversity_all <- yarg::CorrectSamplingEffort(diversity_all)

# Next we’ll merge any sites that are within the same land-use type and that have identical
# coordinates, start and end dates.

# -----------------------------------------------------------------------
# I'm getting an error: 
# Error in yarg::MergeSites(diversity, silent = TRUE, merge.extra = "Wilderness_area") : 
# length(overwrite) == length(measurements) is not TRUE

# synthesize.cols and setdiff(colnames(diversity), c(match.cols, site.level.merge.cols,  .... are not equal:
# Lengths (8, 37) differ (string compare on first 8)
# 8 string mismatches


setdiff(diversity, diversity_all)
diversity <- yarg::MergeSites(diversity,
                              silent = TRUE,
                              merge.extra = "Wilderness_area")

diversity_all <- yarg::MergeSites(diversity_all,
                              silent = TRUE,
                              merge.extra = "Wilderness_area")

# --------------------------------------------------------------------------

# Rename the Predominant_habitat column, since it’s not really “habitat” that we’re looking at.
# This will also keep things a little more consistent with the public version of the PREDICTS database.

diversity <- rename(diversity,
                    Predominant_land_use = Predominant_habitat)

### 1. Calculate diversity metrics -----------------------------------------------

sites <- diversity %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
 
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use")) %>% 
  
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

### --------------------------- Total abundance try ---------------------------------------------
# I'm getting this warning message when running the SiteMetrics Warning message:
# In min(x[x > 0]) : no non-missing arguments to min; returning Inf


# I'm not getting the same results if I do the operations only for one site

TotaAbundance <- diversity %>% filter(SSS == "BS1_2010__Page 1 1") %>% select(Best_guess_binomial, Diversity_metric, 
                                                                              Diversity_metric_unit, Measurement, SS, SSS) %>%
  mutate(Total_Abundance = sum(Measurement), Species_richness = length(unique(Best_guess_binomial))) 

TotaAbundance2 <- diversity %>% filter(SSS == "BS1_2010__Page 1 2") %>% select(Best_guess_binomial, Diversity_metric, 
                                                                       Diversity_metric_unit, Measurement, SS, SSS) %>%
  mutate(Total_Abundance = sum(Measurement), Species_richness = length(unique(Best_guess_binomial))) 

### ---------------------------------------------------------------------------------------------------

# The Predominant_land_use column holds the land use categories. They often need some tidying up.

sites <- sites %>%
  
  mutate(
    
    # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary", 
                                         "Primary non-forest" = "Primary"),
    
    # indeterminate secondary veg and cannot decide get NA
    Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"),
    Use_intensity =  as.character(Use_intensity), 
    Use_intensity = na_if(Use_intensity, "Cannot decide"),

    
    # set reference levels
    Predominant_land_use = factor(Predominant_land_use),
    Predominant_land_use = relevel(Predominant_land_use, ref = "Primary"),
    Use_intensity = factor(Use_intensity),
    Use_intensity = relevel(Use_intensity, ref = "Minimal use")
  )


# take a look at the LandUse/Use intensity split
table(sites$Predominant_land_use, sites$Use_intensity)

# Number of studies where the abundance is not an integer
f <- which(sites$Total_abundance < 1 & sites$Total_abundance > 0)


### 2. Model Site level diversity --------------------------------------------------------------

# Collinearity: ---------------------------------------------------------------

# Ok, so let’s check whether any of the candidate explanatory variables are strongly correlated 
# with one another. You can look at the linear correlations between variables (using cor.test)
# or correlations between continuous and categorical variables (using aov). However, we tend to
# have multiple categorical varialbes. A good alternative is to look at Generalized Variance 
# Inflation Factors. You can use a function provided by the Zuur book.

source("https://highstat.com/Books/Book2/HighstatLibV10.R")

corvif(sites[ , c("Predominant_land_use", "Use_intensity")])

# There are different ‘rules of thumb’ about how high is too high when it comes to GVIFs. 
# Under 3 is great and under 5 is ok. Some people also say under 10 is acceptable, but I 
# think that’s a bit too high personally.

# If collinearity is a problem, then you’ll need to consider which variables are most important
# to include and what ones can be dropped.

# Complete cases: ------------------------------------------------------------

model_data <- drop_na(sites, 
                      Total_abundance, Predominant_land_use, Use_intensity)

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

hist(model_data$RescaledAbundance) # Bound to zero

# First, let’s transform RescaledAbundance.

model_data <- mutate(model_data, 
                     logAbundance = log(RescaledAbundance + 1),
                     sqrtAbundance = sqrt(RescaledAbundance)
)

hist(model_data$logAbundance) 
hist(model_data$sqrtAbundance)

m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + 
             Predominant_land_use:Use_intensity + 
             (1|SS) + (1|SSB), data = model_data)

# fixed-effect model matrix is rank deficient so dropping 12 columns / coefficients
# boundary (singular) fit: see ?isSingular

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

