
Merge_landUses_and_intensities <- function(dataset, 
                                           index,
                                           land_uses_separate_intensities,
                                           land_uses_merge_light_intense, 
                                           reference) {
  
  # Create the column names that will be the current column names plus an index
  name_predominant_habitat_column <- base::paste("Predominant_land_use", index, sep = ".")
  name_use_intensity_column <- base::paste("Use_intensity", index , sep = ".")
  name_LandUse <- base::paste("LandUse", index, sep = ".")
  
  dataset %>%
    # make a level of Primary minimal
    dplyr::mutate(
      
      # Create a new land-use column to store the transformations we are going to make, the name
      # comes from the previous step where we paste the original name and the index
      !!name_predominant_habitat_column := base::paste(Predominant_land_use),
      
      # collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
      !!name_predominant_habitat_column := dplyr::recode_factor(!!as.name(name_predominant_habitat_column), 
                                                                "Primary forest" = "Primary", 
                                                                "Primary non-forest" = "Primary"),
      
      # indeterminate secondary veg and cannot decide are transformed into NA, urban too 
      !!name_predominant_habitat_column := dplyr::na_if(!!as.name(name_predominant_habitat_column), "Secondary vegetation (indeterminate age)"),
      !!name_predominant_habitat_column := dplyr::na_if(!!as.name(name_predominant_habitat_column), "Cannot decide"),
      !!name_predominant_habitat_column := dplyr::na_if(!!as.name(name_predominant_habitat_column), "Urban"),
      
      
      # Give a shorter name to some land use categories 
      !!name_predominant_habitat_column := stringr::str_replace_all(!!as.name(name_predominant_habitat_column), pattern = c("Young secondary vegetation" = "YSV",
                                                                                                                            "Intermediate secondary vegetation" = "ISV", 
                                                                                                                            "Mature secondary vegetation" = "MSV")),
      
      # Create a new use_intensity column to store the transformations we are going to make, the name
      # comes from the previous step where we paste the original name and the index
      !!name_use_intensity_column := base::paste(Use_intensity),
      
      # Drop the Cannot decide intensity levels for the land-use categories that have 
      # enough sites for the minimal, and light/intense
      !!name_use_intensity_column := ifelse(!!as.name(name_predominant_habitat_column) %in% land_uses_separate_intensities & 
                                              !!as.name(name_use_intensity_column) == "Cannot decide",
                                            NA,
                                            base::paste(!!as.name(name_use_intensity_column))), 
      
      # Join the intensity levels of light and intense for those land-uses where it is necessary 
      !!name_use_intensity_column := ifelse((!!as.name(name_predominant_habitat_column) %in% land_uses_merge_light_intense) & 
                                              (!!as.name(name_use_intensity_column) == "Intense use" | 
                                                 !!as.name(name_use_intensity_column) == "Light use"),
                                            stringr::str_replace_all(!!as.name(name_use_intensity_column),
                                                                     pattern = c("Intense use" = "Light-intense use", 
                                                                                 "Light use" = "Light-intense use")), 
                                            paste(!!as.name(name_use_intensity_column))),
      
      # Merge all the intensity levels for those land-use categories that don't have enough sites in each land-use type/intensity combination
      !!name_use_intensity_column := ifelse(!!as.name(name_predominant_habitat_column) %nin% land_uses_separate_intensities,
                                            stringr::str_replace_all(!!as.name(name_use_intensity_column), pattern = c("Intense use" = "All", 
                                                                                                                       "Light use" = "All", 
                                                                                                                       "Minimal use" = "All", 
                                                                                                                       "Cannot decide" = "All")), 
                                            base::paste(!!as.name(name_use_intensity_column))),
      
      
      # Paste the land-use classes and intensity levels
      !!name_LandUse := ifelse(!!as.name(name_predominant_habitat_column) != "NA" & !!as.name(name_use_intensity_column) != "NA",
                               base::paste(!!as.name(name_predominant_habitat_column), !!as.name(name_use_intensity_column)),
                               NA),
      
      # set reference level
      !!name_LandUse := base::factor(!!as.name(name_LandUse), ordered = FALSE),
      !!name_LandUse := relevel(!!as.name(name_LandUse), ref = reference)
    )
}
