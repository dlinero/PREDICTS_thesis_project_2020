# Load libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)


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

# Orders -------------------------------------------------------------------------

# Import table
frugi_endooPlants_notEndooPlants <- readRDS(file = "./output/cleaned_data/01_Filter_data_frugi_endooPlants_notEndooPlants_records.rds")

# Calculate number of records per order for frugivores
frugi <- frugi_endooPlants_notEndooPlants %>% subset(Kingdom == "Animalia") %>%
  droplevels() 
table(frugi$Order)

# Calculate number of records per order for endozoochorous plants
endo <- frugi_endooPlants_notEndooPlants %>% subset(Kingdom == "Plantae") %>%
  droplevels() 
table(endo$Order)

# Calculate number of records per order for abiotically-dispersed plants
abi <- frugi_endooPlants_notEndooPlants %>% subset(Kingdom == "nePlantae") %>%
  droplevels() 
table(abi$Order)

# Calculate number of species per order for frugivores
frugi <- frugi %>% distinct(Best_guess_binomial, .keep_all = TRUE)
table(frugi$Order)

# Calculate number of species per order for endo
endo <- endo %>% distinct(Best_guess_binomial, .keep_all = TRUE)
table(endo$Order)

# Calculate number of species per order for abi
abi <- abi %>% distinct(Best_guess_binomial, .keep_all = TRUE)
table(abi$Order)

# Import summary table 
orders <- read.csv("./output/graphs_tables/orders.csv", sep= ";")

# Ordering the functional groups
orders$group <- factor(orders$group, 
                       levels = c("Aves", "Mammalia", 
                                  "Endozoochorous plants", 
                                  "Abiotically-dispersed plants"))
# Eliminate spaces between words
orders$order <-trimws(orders$order)

# Export graph
pdf(file = "./output/graphs_tables/05_Results_Graphs_Tables_Orders.pdf", width = 14)


ggplot(orders) + 
  # barplot
  geom_bar(
    # Organize the bars in descending order and give different colors to the orders
    aes(x = reorder(order, -species), y = species, fill = order, group = order), 
    stat='identity', position = 'dodge'
  ) +
  # Add a number with the number of records for every order
  geom_text(
    aes(x = order, y = species, label = records, group = order),
    position = position_dodge(width = 1),
    vjust = -0.3, size = 3, fontface = "bold"
  ) + 
  # Remove legend and give axis labels
  theme_bw() + theme(legend.position = "none") + labs(y = "Number of species",
                                                      x = "Order") +
  # Faceted plot of groups, allow the order names to be different in every facet
  facet_wrap( ~ group, scales = "free_y", strip.position = 'top') +
  # Rotate graph
  theme(strip.placement = "outside", axis.text=element_text(size=7)) + coord_flip()

# Clear
dev.off()
