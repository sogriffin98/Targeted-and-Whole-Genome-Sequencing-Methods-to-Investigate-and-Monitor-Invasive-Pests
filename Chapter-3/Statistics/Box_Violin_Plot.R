# Set working directory
setwd("~/PhD/Thesis/chapter 1b distance between nests/Distance Stats Test 2nd Dec")

# Load libraries
library(ggplot2)

# Read the CSV
data <- read.csv("final_dataset.csv", stringsAsFactors = FALSE)

# Make sure MeanLR is numeric
data$MeanLR <- as.numeric(data$MeanLR)

# Check for any NA values created during conversion
sum(is.na(data$MeanLR))  # Should be 0 ideally

# Set the order of Kin_relationship
data$Kin_relationship <- factor(data$Kin_relationship, 
                                levels = c("Full Sibling", "Cousin", "Aunt-Niece", "Great-Aunt-Niece"))

# Distance_KM violin in grayscale
ggplot(data, aes(x = Kin_relationship, y = Distance_KM, fill = Kin_relationship)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_grey(start = 0.8, end = 0.2) +  # light to dark grey
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Distance_KM by Kin Relationship",
       x = "Kin Relationship",
       y = "Distance (KM)") +
  theme(legend.position = "none") +
  ylim(0, max(data$Distance_KM, na.rm = TRUE))

# MeanLR violin in grayscale
ggplot(data, aes(x = Kin_relationship, y = MeanLR, fill = Kin_relationship)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of MeanLR by Kin Relationship",
       x = "Kin Relationship",
       y = "MeanLR") +
  theme(legend.position = "none") +
  ylim(0, max(data$MeanLR, na.rm = TRUE))

# Simple black & white boxplot
ggplot(pair_data, aes(x = Kin_relationship, y = Distance_KM)) +
  geom_boxplot(fill = "white", colour = "black", outlier.shape = 21, outlier.size = 2) +
  labs(
    title = "Distribution of Nest-to-Nest Distances by Kinship Category",
    x = "Kin relationship",
    y = "Distance between nests (km)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )


