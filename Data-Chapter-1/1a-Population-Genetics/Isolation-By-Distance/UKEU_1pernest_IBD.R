# Set working directory
setwd("B:/combined_run_data/IBD/uk_eu_1_per_nest")

# Load required libraries
required_packages <- c("adegenet", "geosphere", "vegan", "ggplot2")
installed_packages <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed_packages)
if(length(to_install)) install.packages(to_install, dependencies = TRUE)
lapply(required_packages, library, character.only = TRUE)

# Load genetic data (PLINK raw format)
gl <- read.PLINK("UKEU_1pernest_data_plink_output.raw")  # adjust filename if needed

# Load coordinates
coords <- read.csv("UKEU_1pernest_coords.csv")
coords <- as.matrix(coords[, c("longitude", "latitude")])

# Check that individuals match between genotype and coords
stopifnot(nrow(coords) == nInd(gl))

# Calculate genetic distance using Euclidean distance on genotype matrix
geno_matrix <- as.matrix(gl)
gen_dist <- dist(geno_matrix)  # returns a 'dist' object

# Calculate geographic distances (km)
geo_dist <- distm(coords, fun = distHaversine) / 1000
geo_dist <- as.dist(geo_dist)

# Perform Mantel test
mantel_result <- mantel(gen_dist, geo_dist, method = "pearson", permutations = 9999)
print(mantel_result)

# Prepare data for plotting
df <- data.frame(
  Geographic_Distance_km = as.vector(geo_dist),
  Genetic_Distance = as.vector(gen_dist)
)

# Plot IBD
plot(as.vector(geo_dist), as.vector(gen_dist),
     xlab = "Geographic Distance (km)",
     ylab = "Genetic Distance",
     main = "Isolation by Distance")
abline(lm(as.vector(gen_dist) ~ as.vector(geo_dist)), col = "red")
