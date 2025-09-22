# Set working directory
setwd("~/Desktop/microsats_graph")

# Load packages
library(ggplot2)
library(adegenet)
library(related)
library(viridis)

# Load genotype data
genotype_data <- read.table("raw_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Rename columns: replace _A -> .1 and _B -> .2
colnames(genotype_data)[-1] <- gsub("_A$", ".1", colnames(genotype_data)[-1])
colnames(genotype_data)[-1] <- gsub("_B$", ".2", colnames(genotype_data)[-1])

# Extract allele data (all except ID)
allele_data <- genotype_data[, -1]

# Get unique locus names by removing the ".1" or ".2" suffixes
loci <- unique(gsub("\\.(1|2)$", "", colnames(allele_data)))

# Create a new data frame with one column per locus,
# combining the two alleles with a separator "/"
combined_alleles <- data.frame(matrix(nrow = nrow(allele_data), ncol = length(loci)))
colnames(combined_alleles) <- loci

for (loc in loci) {
  col1 <- paste0(loc, ".1")
  col2 <- paste0(loc, ".2")
  combined_alleles[[loc]] <- paste(allele_data[[col1]], allele_data[[col2]], sep = "/")
}

# Convert combined allele columns to factors
combined_alleles[] <- lapply(combined_alleles, factor)

# Now convert to genind object,
# specifying sep = "/" (the allele separator),
# and ncode = 3 (number of characters per allele, adjust if needed)
file_genind <- df2genind(combined_alleles, ploidy = 2,
                         ind.names = genotype_data$ID,
                         sep = "/", ncode = 3, type = "codom")

# Convert genind to data frame and save
data_df <- genind2df(file_genind, sep = "\t")
write.table(data_df, file = "genind_dataframe_100iterations.txt", sep = "\t", row.names = TRUE, col.names = FALSE)

# Prepare data for related package

# Format related package input
# related package expects data in a matrix or data frame where
# first column = individual IDs
# followed by 2 columns per locus (alleles)

# Combine IDs with allele data
related_input <- cbind(genotype_data$ID, allele_data)
colnames(related_input)[1] <- "ID"

# Convert to data.frame, ensure correct types
related_input <- as.data.frame(related_input, stringsAsFactors = FALSE)

# related package expects all genotype data to be numeric (alleles)
# Convert allele columns to numeric
for(i in 2:ncol(related_input)){
  related_input[, i] <- as.numeric(as.character(related_input[, i]))
}

# Select only necessary columns (assuming 15 loci)
# 15 loci × 2 alleles + 1 ID column = 31 columns
related_input_trim <- related_input[, 1:(2 * 15 + 1)]

# Convert to matrix (without ID) for readgenotypedata()
genotypes_matrix <- as.matrix(related_input_trim[, -1])

# Use IDs separately
individual_ids <- related_input_trim[, 1]

# Combine into one data frame:
related_input <- cbind(individual_ids, genotypes_matrix)

# Write to a file, e.g. tab-delimited with no row names or column names:
write.table(related_input, "related_input_100iterations.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Now read into related package:
related_data <- readgenotypedata("related_input_100iterations.txt")

# Compare estimators
sink("estimator_output_microsatellites_100iterations.txt")
compareestimators(related_data, 100)
sink()

# Calculate relatedness values
genotypes <- read.table("related_input_100iterations.txt", header = FALSE)
relatedness_output <- coancestry(genotypes, error.rates = 0.05, lynchrd = 1)
write.table(relatedness_output$relatedness, "relatedness_output_100iterations.txt", sep = ",", row.names = TRUE, col.names = TRUE)

# Import data
related_results <- read.table("relatedness_output_100iterations.txt", header = TRUE, sep = ",", row.names = 1)

# Custom ggplot theme
coolbeans <- theme(
  legend.position = "right",
  legend.title = element_text(size = 15),
  legend.text = element_text(size = 12),
  axis.text.x = element_text(size = 12, colour = "black", face = "bold"),
  axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
  axis.title.x = element_text(size = 15, colour = "black", vjust = +0.5, face = "bold"),
  axis.title.y = element_text(size = 15, colour = "black", vjust = +0.5, face = "bold"),
  panel.background = element_rect(fill = "white", colour = NULL),
  axis.line = element_line(colour = "black", size = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# Plot histogram of relatedness
ggplot(relatedness_output$relatedness, aes(x = lynchrd)) +
  geom_histogram(binwidth = 0.025, colour = "black", fill = "steelblue") +
  scale_x_continuous(expand = c(0, 0), name = "Lynch & Ritland Relatedness") +
  scale_y_continuous(expand = c(0, 0), name = "Density") +
  coolbeans

# Save relatedness results to file
write.table(relatedness_output$relatedness, "relatedness_output_1000iterations.txt", sep = ",", row.names = TRUE, col.names = TRUE)

pdf("relatedness_histogram.pdf", width = 8, height = 6)  # open PDF device with filename and size
ggplot(relatedness_output$relatedness, aes(x = lynchrd)) +
  geom_histogram(binwidth = 0.025, colour = "black", fill = "steelblue") +
  scale_x_continuous(expand = c(0, 0), name = "Lynch & Ritland Relatedness") +
  scale_y_continuous(expand = c(0, 0), name = "Density") +
  coolbeans
dev.off()  # close PDF device and write file

# You can write the relatedness values to a text document:
write.table(relatedness_output$relatedness,"relatedness_100iterations.txt",sep=",",row.names = TRUE, col.names = TRUE)

# Simulate related individuals based on allele frequencies from our data:
sim <- familysim(related_data$freqs, 100)
sim_output  <- coancestry(sim , quellergt = 1)
simrel  <- cleanuprvals(sim_output$relatedness , 100)
relvalues  <- simrel[, 10]
label1  <- rep("PO", 100)
label2  <- rep("Full", 100)
label3  <- rep("Half", 100)
label4  <- rep("Unrelated", 100)
Relationship<- c(label1 , label2 , label3 , label4)
newdata  <- as.data.frame(cbind(Relationship , relvalues))
newdata$relvalues  <- as.numeric(as.character(newdata$relvalues))

# Plot simulated data:
ggplot(newdata, aes(x= relvalues, colour = Relationship, fill = Relationship))+ 
  geom_density(alpha = 0.3, linewidth = 0.75)+
  scale_x_continuous(name = "Lynch & Ritland relatedness")+
  scale_y_continuous(name = "density")+
  scale_colour_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE, option = "D")+
  coolbeans
