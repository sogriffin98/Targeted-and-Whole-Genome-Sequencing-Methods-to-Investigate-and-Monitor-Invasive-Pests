# Set working directory
setwd("~/PhD/Thesis/chapter 1b distance between nests/Distance Stats Test 2nd Dec")

# Load Packages
library(ggplot2)
library(dplyr)
library(lme4)
library(FSA)

# Read in data
pair_data <- read.csv("final_dataset.csv", stringsAsFactors = FALSE)
pair_data$Kin_relationship <- as.factor(pair_data$Kin_relationship)

# Create a dataset without full siblings
pair_no_sibs <- pair_data %>% 
  filter(Kin_relationship != "Full Sibling")

# Scatterplot with all samples
plot_all <- ggplot(pair_data, aes(x = Distance_KM, y = MeanLR,
                                  colour = Kin_relationship)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 1) +
  scale_color_manual(values = c(
    "Full Sibling" = "#4477AA",      
    "Cousin" = "#EE6677",            
    "Aunt-Niece" = "#228833",        
    "Great Aunt-Niece" = "#CCBB44"   
  )) +
  labs(
    title = "All Samples",
    x = "Distance between nests (km)",
    y = "Mean Lynch & Ritland relatedness",
    colour = "Kin relationship"
  ) +
  theme_minimal(base_size = 14)

plot_all

# Scatterplot without full siblings
plot_no_sibs <- ggplot(pair_no_sibs, aes(x = Distance_KM, y = MeanLR,
                                         colour = Kin_relationship)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 1) +
  scale_color_manual(values = c(
    "Cousin" = "#EE6677",            
    "Aunt-Niece" = "#228833",        
    "Great Aunt-Niece" = "#CCBB44"   
  )) +
  labs(
    title = "Without Full Siblings",
    x = "Distance between nests (km)",
    y = "Mean Lynch & Ritland relatedness",
    colour = "Kin relationship"
  ) +
  theme_minimal(base_size = 14)

plot_no_sibs

# Combined plots
pair_data$Dataset <- "All samples"
pair_no_sibs$Dataset <- "Without full siblings"

combined <- rbind(pair_data, pair_no_sibs)
combined$Dataset <- factor(combined$Dataset,
                           levels = c("All samples", "Without full siblings"))

kin_colors <- c(
  "Full Sibling" = "#4477AA",      # blue
  "Cousin" = "#EE6677",            # vermillion
  "Aunt-Niece" = "#228833",        # teal
  "Great Aunt-Niece" = "#CCBB44"   # yellowish purple
)

# For the panel without full siblings, "Full_sibling" will simply not appear in the legend
plot_combined <- ggplot(combined, aes(x = Distance_KM, y = MeanLR,
                                      colour = Kin_relationship)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 1) +
  scale_color_manual(values = kin_colors) +
  facet_wrap(~ Dataset) +
  labs(
    title = "Relatedness vs Distance (Comparison)",
    x = "Distance between nests (km)",
    y = "Mean Lynch & Ritland relatedness",
    colour = "Kin relationship"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

plot_combined

# Run statistical tests on the datasets
run_all_stats <- function(df, label = "Dataset") {
  
  cat("\n=====================================================\n")
  cat(" STATISTICS FOR:", label, "\n")
  cat("=====================================================\n\n")
  
# Correlations
  cat("Pearson correlation:\n")
  print(cor.test(df$MeanLR, df$Distance_KM, method = "pearson"))
  cat("\nSpearman correlation:\n")
  print(cor.test(df$MeanLR, df$Distance_KM, method = "spearman"))
  
# MLPE mixed model
  cat("\nMLPE mixed model:\n")
  m <- lmer(MeanLR ~ Distance_KM +
              (1 | NestID1) + (1 | NestID2),
            data = df)
  print(summary(m))
  
# ANOVA + Tukey
  cat("\nANOVA (if assumptions met):\n")
  a <- aov(MeanLR ~ Kin_relationship, data = df)
  print(summary(a))
  
  cat("\nTukey HSD:\n")
  print(TukeyHSD(a))
  
# Kruskal-Wallis + Dunn
  cat("\nKruskal-Wallis:\n")
  print(kruskal.test(MeanLR ~ Kin_relationship, data = df))
  
  cat("\nDunn Test (Bonferroni):\n")
  print(dunnTest(MeanLR ~ Kin_relationship,
                 data = df, method = "bonferroni"))
  
# Distance ~ relationship
  cat("\nDistance difference by relationship:\n")
  print(kruskal.test(Distance_KM ~ Kin_relationship, data = df))
  
  cat("\nDunn Test (distance):\n")
  print(dunnTest(Distance_KM ~ Kin_relationship,
                 data = df, method = "bonferroni"))
  
  cat("\n=====================================================\n\n")
}

# Run the stats tests for both datasets
run_all_stats(pair_data, label = "All Samples")
run_all_stats(pair_no_sibs, label = "Without Full Siblings")



###############################################################
### 4️⃣ COMBINED 2-PANEL SCATTERPLOT WITH CUSTOM SHAPES
###############################################################

# Add dataset labels
pair_data$Dataset <- "A) All samples"
pair_no_sibs$Dataset <- "B) Without full siblings"

combined <- rbind(pair_data, pair_no_sibs)
combined$Dataset <- factor(combined$Dataset,
                           levels = c("All samples", "Without full siblings"))

# Colour-blind-safe palette (Tol Bright)
kin_colors <- c(
  "Full Sibling" = "#CCBB44",      # blue
  "Cousin" = "#EE6677",            # vermillion
  "Aunt-Niece" = "#228833",        # teal
  "Great Aunt-Niece" = "#4477AA"   # yellowish purple
)

# Custom shape mapping
kin_shapes <- c(
  "Full Sibling" = 17,     # triangle
  "Cousin" = 16,           # circle
  "Aunt-Niece" = 15,       # square
  "Great Aunt-Niece" = 4   # cross
)

plot_combined <- ggplot(combined, aes(x = Distance_KM, y = MeanLR)) +
  # Points only get colour and shape
  geom_point(aes(colour = Kin_relationship, shape = Kin_relationship),
             alpha = 0.7, size = 3) +
  # Trend line for all points in each facet (black, original)
  geom_smooth(data = subset(combined, Dataset == "All samples"),
              method = "lm", se = TRUE, colour = "black", size = 1) +
  geom_smooth(data = subset(combined, Dataset == "Without full siblings"),
              method = "lm", se = TRUE, colour = "black", size = 1) +
  scale_color_manual(values = kin_colors) +
  scale_shape_manual(values = kin_shapes) +
  facet_wrap(~ Dataset) +
  labs(
    title = "Relatedness vs Distance (Comparison)",
    x = "Distance between nests (km)",
    y = "Mean Lynch & Ritland relatedness",
    colour = "Kin relationship",
    shape = "Kin relationship"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  )

plot_combined




# Read data
pair_data <- read.csv("final_dataset.csv", stringsAsFactors = FALSE)
pair_data$Kin_relationship <- as.factor(pair_data$Kin_relationship)

# Black & white / grayscale palette
bw_palette <- c(
  "Full Sibling" = "black",
  "Cousin" = "grey20",
  "Aunt-Niece" = "grey50",
  "Great Aunt-Niece" = "grey80"
)

# Boxplot (black & white)
ggplot(pair_data, aes(x = Kin_relationship, y = Distance_KM, fill = Kin_relationship)) +
  geom_boxplot(alpha = 1, outlier.shape = 21, outlier.size = 2, colour = "black") +
  scale_fill_manual(values = bw_palette) +
  labs(
    title = "Distribution of Nest-to-Nest Distances by Kinship Category",
    x = "Kin relationship",
    y = "Distance between nests (km)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )



