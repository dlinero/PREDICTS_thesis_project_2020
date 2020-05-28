####----------------- Filter tropics -------------------
TropicsDiversity1 <- yarg::CorrectSamplingEffort(TropicsDiversity) 

TropicsDiversity2 <- yarg::MergeSites(TropicsDiversity1,
                                         silent = TRUE,   merge.extra = "Wilderness_area")


TropicsDiversity3 <-  TropicsDiversity2  %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_habitat")) %>% 
  
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

TropicsDiversity_correct <- TropicsDiversity2 %>%
  
  filter(SSS == "BS1_2010__Page 1 1") %>%
  
  select(Taxon_name_entered, Measurement, SSS) %>%
  
  mutate(Total_Abundance = sum(Measurement), 
         Species_richness = length(unique(Taxon_name_entered)))

###------------------- Filter frugivores -----------------
frugivores <- read.csv("output/intermediate_files/01_Filter_data_Frugivores_names.csv")
frugivores <- as.character(frugivores$x)


PREDICTS_frugivores <- TropicsDiversity %>% 
  
  subset(Best_guess_binomial %in% frugivores) %>% 
  
  droplevels()

PREDICTS_frugivores1 <- yarg::CorrectSamplingEffort(PREDICTS_frugivores) 

PREDICTS_frugivores2 <- yarg::MergeSites(PREDICTS_frugivores1,
                                      silent = TRUE,   merge.extra = "Wilderness_area")

PREDICTS_frugivores3 <-  PREDICTS_frugivores2  %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_habitat")) %>% 
  
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

PREDICTS_frugivores_correct <- PREDICTS_frugivores2 %>%
  
  filter(SSS == "DI1_2010__SaldanaVazquez 1 1") %>%
  
  select(Taxon_name_entered, Measurement, SSS) %>%
  
  mutate(Total_Abundance = sum(Measurement), 
         Species_richness = length(unique(Taxon_name_entered)))

###------------------- Filter plants -----------------

plants <- read.csv("output/intermediate_files/01_Filter_data_List_endozoochory_species.csv")
plants <- as.character(plants$x)

PREDICTS_plants <- TropicsDiversity %>% 
  
  subset(Best_guess_binomial %in% plants) %>% 
  
  droplevels()

PREDICTS_plants1 <- yarg::CorrectSamplingEffort(PREDICTS_plants) 

PREDICTS_plants2 <- yarg::MergeSites(PREDICTS_plants1,
                                         silent = TRUE,   merge.extra = "Wilderness_area")

PREDICTS_plants3 <-  PREDICTS_plants2  %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_habitat")) %>% 
  
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

PREDICTS_plants_correct <- PREDICTS_plants2 %>%
  
  filter(SSS == "BS1_2010__Page 1 1") %>%
  
  select(Taxon_name_entered, Measurement, SSS) %>%
  
  mutate(Total_Abundance = sum(Measurement), 
         Species_richness = length(unique(Taxon_name_entered)))

### --------------------------------- Filtered table -----------------------

frugivores_and_plants <- c(frugivores, plants)

PREDICTS_plants_and_frugivores <- TropicsDiversity %>% 
  
  subset(Best_guess_binomial %in% frugivores_and_plants) %>% 
  
  droplevels()

Noclass_kingdom  <- which(PREDICTS_plants_and_frugivores$Class == "" | PREDICTS_plants_and_frugivores$Kingdom == "")
List_sp <- as.data.frame(unique(PREDICTS_plants_and_frugivores[Noclass_kingdom, "Best_guess_binomial"]))

# This loop searches for the index of species that do not have Kingdom or Class details, and 
# adds the information in the corresponding rows
for (i in Noclass_kingdom) {
  
  # add  information of kingdom and class 
  if (PREDICTS_plants_and_frugivores[i, "Best_guess_binomial"] %in% List_sp[,1]) { 
    PREDICTS_plants_and_frugivores[i, "Kingdom"] <- "Plantae" 
    PREDICTS_plants_and_frugivores[i, "Class"] <- "Magnoliopsida"
  }
}

PREDICTS_plants_and_frugivores1 <- yarg::CorrectSamplingEffort(PREDICTS_plants_and_frugivores) 

PREDICTS_plants_and_frugivores2 <- yarg::MergeSites(PREDICTS_plants_and_frugivores1,
                                     silent = TRUE,   merge.extra = "Wilderness_area")

PREDICTS_plants_and_frugivores3 <-  PREDICTS_plants_and_frugivores2  %>%
  
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  
  # The extra.cols parameter is used for columns that we want to 
  # transferred to the final site-level data frame and that the function 
  # does not add  automatically
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_habitat")) %>% 
  
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


PREDICTS_plants_and_frugivores_correct <- PREDICTS_plants_and_frugivores2 %>%
  
  filter(SSS == "BS1_2010__Page 1 1") %>%
  
  select(Taxon_name_entered, Measurement, SSS) %>%
  
  mutate(Total_Abundance = sum(Measurement), 
         Species_richness = length(unique(Taxon_name_entered)))

# Export table 
saveRDS(PREDICTS_plants_and_frugivores , "./output/cleaned_data/01_Filter_data_PREDICTS_Filtered_table.rds")

# Total number of sourceIDs
length(unique(PREDICTS_plants_and_frugivores$Source_ID))

# Total number of studies
length(unique(PREDICTS_plants_and_frugivores$SS))

# Total number of sites
length(unique(PREDICTS_plants_and_frugivores$SSBS))
