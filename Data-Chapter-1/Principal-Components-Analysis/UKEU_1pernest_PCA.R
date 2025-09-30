# Set working directory
setwd("B:/combined_run_data/PCA/Original_WORKED/UK_EU_1pernest/vcf")

# Install packages
install.packages("ggplot2")
install.packages("readr")

# load packages
library(readr)
library(ggplot2)

# important plink output files 
PCA_file_covMat <- read.table("UKEU_pca_results.eigenvec", header=FALSE) 
# import the .cov file
eigenvals <- read.table("UKEU_pca_results.eigenval", header=FALSE)
# remove double ID column
pca <- PCA_file_covMat[,-1]
#pca <- PCA_file_covMat
head(pca)
# set header names
names(pca)[1] <- "ind"
ncolpca <- ncol(pca)
PCs <- ncolpca -1
names(pca)[2:ncol(pca)] <- paste0("PC", 1:PCs)
head (pca)

# load in the sample ID assignments 
meta<- read.csv("sample_ID_UKEU_1pernest.csv", header=FALSE)
head(meta)
pop_assings <- meta[,2]
levs <- as.factor(pop_assings)
levs2 <- levels(levs)
# make dataframe for plot
PC1 <- pca$PC1
PC2 <- pca$PC2
PC3 <- pca$PC3
PC4 <- pca$PC4

factor_colors <- as.character(as.factor(meta$V3))
col_levs <- unique(factor_colors)
factor_shapes <- as.numeric(meta$V4)

#factor_colors <- colors[as.numeric((as.factor(pop_assings)))]
#factor_shapes <- shapes[(as.factor(pop_assings))]

plot <- cbind(pop_assings, PC1,PC2,PC3,PC4,factor_colors,factor_shapes)
head(plot)
plot2 <- as.data.frame(plot)
head(plot2)

PC1_axis <- round(eigenvals[1,], digits = 2)
PC2_axis <- round(eigenvals[2,], digits = 2)
PC3_axis <- round(eigenvals[3,], digits = 2)
PC4_axis <- round(eigenvals[4,], digits = 2)

# Create the scatter plot for PC1 vs PC2

21.1006 + 29.79838

par(mar = c(5, 4, 4, 8), xpd = TRUE)

plot(plot2$PC1, plot2$PC2, 
     col = plot2$factor_colors,          
     pch = as.numeric(plot2$factor_shapes),
     xlab = paste("PC 1 (", PC1_axis, "%)", sep = ""),  # Combining "PC 1" with the value of PC1
     ylab = paste("PC 2 (", PC2_axis, "%)", sep = ""))

legend <- read.csv("legend.csv", header=FALSE)
colnames(legend) <- c("Pop", "Color", "Shape")

# Add legend to the right outside the plot
legend("topright",                                  
       inset = c(-0.1, 0),           #Moves the legend outside
       legend = legend$Pop, 
       col = legend$Color, 
       pch = as.numeric(legend$Shape),
       cex = 0.75,
       bty = "n")  

# Create the scatter plot for PC1 vs PC3

21.1006 + 29.79838

par(mar = c(5, 4, 4, 8), xpd = TRUE)

plot(plot2$PC1, plot2$PC3, 
     col = plot2$factor_colors,          
     pch = as.numeric(plot2$factor_shapes),
     xlab = paste("PC 1 (", PC1_axis, "%)", sep = ""),  # Combining "PC 1" with the value of PC1
     ylab = paste("PC 3 (", PC2_axis, "%)", sep = ""))

legend <- read.csv("legend.csv", header=FALSE)
colnames(legend) <- c("Pop", "Color", "Shape")

# Add legend to the right outside the plot
legend("topright",                                  
       inset = c(-0.1, 0),           #Moves the legend outside
       legend = legend$Pop, 
       col = legend$Color, 
       pch = as.numeric(legend$Shape),
       cex = 0.75,
       bty = "n")  

# Create the scatter plot for PC1 vs PC4

21.1006 + 29.79838

par(mar = c(5, 4, 4, 8), xpd = TRUE)

plot(plot2$PC1, plot2$PC4, 
     col = plot2$factor_colors,          
     pch = as.numeric(plot2$factor_shapes),
     xlab = paste("PC 1 (", PC1_axis, "%)", sep = ""),  # Combining "PC 1" with the value of PC1
     ylab = paste("PC 4 (", PC2_axis, "%)", sep = ""))

legend <- read.csv("legend.csv", header=FALSE)
colnames(legend) <- c("Pop", "Color", "Shape")

# Add legend to the right outside the plot
legend("topright",                                  
       inset = c(-0.1, 0),           #Moves the legend outside
       legend = legend$Pop, 
       col = legend$Color, 
       pch = as.numeric(legend$Shape),
       cex = 0.75,
       bty = "n")  

# Create the scatter plot for PC2 vs PC3

21.1006 + 29.79838

par(mar = c(5, 4, 4, 8), xpd = TRUE)

plot(plot2$PC2, plot2$PC3, 
     col = plot2$factor_colors,          
     pch = as.numeric(plot2$factor_shapes),
     xlab = paste("PC 2 (", PC1_axis, "%)", sep = ""),  # Combining "PC 1" with the value of PC1
     ylab = paste("PC 3 (", PC2_axis, "%)", sep = ""))

legend <- read.csv("legend.csv", header=FALSE)
colnames(legend) <- c("Pop", "Color", "Shape")

# Add legend to the right outside the plot
legend("topright",                                  
       inset = c(-0.1, 0),           #Moves the legend outside
       legend = legend$Pop, 
       col = legend$Color, 
       pch = as.numeric(legend$Shape),
       cex = 0.75,
       bty = "n")  

# Create the scatter plot for PC2 vs PC4

21.1006 + 29.79838

par(mar = c(5, 4, 4, 8), xpd = TRUE)

plot(plot2$PC2, plot2$PC4, 
     col = plot2$factor_colors,          
     pch = as.numeric(plot2$factor_shapes),
     xlab = paste("PC 2 (", PC1_axis, "%)", sep = ""),  # Combining "PC 1" with the value of PC1
     ylab = paste("PC 4 (", PC2_axis, "%)", sep = ""))

legend <- read.csv("legend.csv", header=FALSE)
colnames(legend) <- c("Pop", "Color", "Shape")

# Add legend to the right outside the plot
legend("topright",                                  
       inset = c(-0.1, 0),           #Moves the legend outside
       legend = legend$Pop, 
       col = legend$Color, 
       pch = as.numeric(legend$Shape),
       cex = 0.75,
       bty = "n")  

# Create the scatter plot for PC3 vs PC4

21.1006 + 29.79838

par(mar = c(5, 4, 4, 8), xpd = TRUE)

plot(plot2$PC3, plot2$PC4, 
     col = plot2$factor_colors,          
     pch = as.numeric(plot2$factor_shapes),
     xlab = paste("PC 3 (", PC1_axis, "%)", sep = ""),  # Combining "PC 1" with the value of PC1
     ylab = paste("PC 4 (", PC2_axis, "%)", sep = ""))

legend <- read.csv("legend.csv", header=FALSE)
colnames(legend) <- c("Pop", "Color", "Shape")

# Add legend to the right outside the plot
legend("topright",                                  
       inset = c(-0.1, 0),           #Moves the legend outside
       legend = legend$Pop, 
       col = legend$Color, 
       pch = as.numeric(legend$Shape),
       cex = 0.75,
       bty = "n")  
