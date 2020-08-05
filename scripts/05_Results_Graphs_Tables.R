# Load libraries
library(ggplot2)


# Create a vector with the total number of sources, influential sources
# according to Cook´s distance and DFBETAS for each of the response variables
number_sources <- c(85, 14, 54, 84, 17, 40, 73, 11, 41)

# Create a vector with the name of the variable
statistic <- c("Total number", "Cook´s distance", "DFBETAS")

# Create a vector with the name of the model
model <- rep(c("Total abundance", "Species richness", "Simpson´s diversity"), each = 3)

# Create ID
ID <- seq(from = 0, to = 4, by =.5)

# Combine into a dataframe
number_influential <- data.frame(model, statistic, number_sources, ID)

# Export the graph 
pdf(file = "./output/graphs_tables/05_Results_Graphs_Tables_Influence_Sources.pdf", width = 10)

# Plot number of source IDs
ggplot(number_influential, aes(number_sources, ID, color = statistic)) +
  
   # Add points
   geom_point() + 
  
  # Add lines to the points
  geom_segment(aes(x = 0, y = ID, xend = number_sources, yend = ID)) + 
  
  # Flip the graph
  coord_flip() +
  
  # Assing the name of the labels
  labs(x = "Number of Source IDs", y = "Model") +
  
  # Change background and font size
  theme(text = element_text(size = 12)) +
  theme_bw() 
  
# Clear plot
dev.off()
