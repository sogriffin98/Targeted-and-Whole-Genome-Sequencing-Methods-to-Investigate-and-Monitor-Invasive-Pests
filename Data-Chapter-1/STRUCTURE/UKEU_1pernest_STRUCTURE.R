# Set working directory
setwd("~/PhD/Study Species/Vespa velutina/STRUCTURE/ylh_uk_eu_3rdJuly")

# Load libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Read the CSV
structure_data <- read.csv("structure_summary_uk_eu.csv")

# Ensure Delta_K is numeric
structure_data$Delta_K <- as.numeric(structure_data$Delta_K)

# Common theme for both plots
clean_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),       # remove grid lines
    axis.line = element_line(color = "black"),  # add axis lines
    axis.ticks = element_line(color = "black"), # add tick marks
    plot.title = element_blank()        # no title
  )

# Plot 1: Mean LnP(D) with error bars
plot_ln <- ggplot(structure_data, aes(x = K, y = Mean_LnP)) +
  geom_errorbar(aes(ymin = Mean_LnP - Stdev_LnP, ymax = Mean_LnP + Stdev_LnP),
                width = 0.2, color = "grey40") +
  geom_point(size = 2) +
  geom_line(color = "steelblue", size = 1) +
  labs(x = "K", y = "Mean LnP(D)") +
  clean_theme

# Plot 2: Delta K (remove NA rows)
plot_delta <- ggplot(filter(structure_data, !is.na(Delta_K)), aes(x = K, y = Delta_K)) +
  geom_point(size = 2) +
  geom_line(color = "darkred", size = 1) +
  labs(x = "K", y = expression(Delta~K)) +
  clean_theme

# Show the plots separately
print(plot_ln)
print(plot_delta)
