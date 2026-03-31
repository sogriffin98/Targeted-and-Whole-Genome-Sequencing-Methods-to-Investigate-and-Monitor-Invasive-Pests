# Install and load packages
install.packages(c("vcfR","adegenet","ape","poppr","RColorBrewer", "promises"))

library(promises)
library(poppr)
library(vcfR)
library(adegenet)
library(ape)
library(RColorBrewer)

# Set working directory
setwd("~/PhD/Thesis/Corrections/trees")

# Function to read VCF, clean sample names, calculate distance, and plot NJ tree
plot_nj_tree <- function(
    vcf_file,
    title = "",              # <-- default: no title
    metadata_file = NULL,
    add_legend = FALSE,
    name_map_file = NULL
){
  
  # Read VCF
  vcf <- read.vcfR(vcf_file)
  
  # Convert VCF to genind object
  gen <- vcfR2genind(vcf)
  
  # Clean sample names
  ind_names_clean <- gsub("^\\.\\/bam\\/", "", indNames(gen))
  ind_names_clean <- gsub("\\.mq30\\.maxinsert\\.primaryalignment\\.bam$", "", ind_names_clean)
  indNames(gen) <- ind_names_clean
  
  # Calculate distance matrix
  dist_matrix <- bitwise.dist(gen)
  
  # Build Neighbor Joining tree
  nj_tree <- nj(dist_matrix)
  
  # --- OPTIONAL: rename tip labels using CSV ---
  if(!is.null(name_map_file)){
    name_map <- read.csv(name_map_file, stringsAsFactors = FALSE)
    lookup <- setNames(name_map$New, name_map$Current)
    
    nj_tree$tip.label <- ifelse(
      nj_tree$tip.label %in% names(lookup),
      lookup[nj_tree$tip.label],
      nj_tree$tip.label
    )
  }
  # ------------------------------------------------
  
  # Default: black tips
  tip_colors <- rep("black", nInd(gen))
  
  # If metadata provided, assign colours
  if(!is.null(metadata_file)){
    meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
    meta$Sample <- trimws(meta$Sample)
    
    match_idx <- match(indNames(gen), meta$Sample)
    countries <- meta$Country[match_idx]
    countries[is.na(countries)] <- "Unknown"
    
    unique_countries <- unique(countries)
    other_countries <- setdiff(unique_countries, "UK")
    
    n_other <- length(other_countries)
    if(n_other <= 8){
      other_colors <- brewer.pal(n_other, "Dark2")
    } else {
      other_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_other)
    }
    
    country_colors <- c("UK" = "black", setNames(other_colors, other_countries))
    
    if("Unknown" %in% names(country_colors) == FALSE & any(countries == "Unknown")){
      country_colors["Unknown"] <- "grey"
    }
    
    tip_colors <- country_colors[countries]
  }
  
  # Scale tip labels
  tip_cex <- 1.5 / log10(nInd(gen) + 1)
  
  # --- Plotting section ---
  
  # Collapse top margin if no title
  top_margin <- ifelse(title == "", 1, 4)
  
  par(mar = c(5, 4, top_margin, 10))
  
  plot(nj_tree,
       tip.color = tip_colors,
       cex = tip_cex,
       main = title,   # empty string = no title
       x.lim = c(0, max(node.depth.edgelength(nj_tree)) * 1.3))
  
  add.scale.bar(
    x = 0,
    y = -0.5,
    length = signif(max(nj_tree$edge.length) / 5, 3)
  )
  
  if(!is.null(metadata_file) && add_legend){
    par(xpd = TRUE)
    legend("right",
           inset = c(-0.25, 0),
           legend = names(country_colors),
           col = country_colors,
           pch = 19,
           cex = 0.7,
           bty = "n")
  }
  
  return(nj_tree)
}

##### UK/EU Data (with colours + legend, NO TITLE)
nj_tree_uk_eu <- plot_nj_tree(
  vcf_file = "AH_104inds_minDP6_mac0.05_biallelic.recode.vcf",
  title = "",
  metadata_file = "metadata.csv",
  add_legend = TRUE,
  name_map_file = "new_names.csv"
)

##### UK Data (NO TITLE)
nj_tree_uk <- plot_nj_tree(
  vcf_file = "AH_UK_89inds_minDP6_mac0.05_biallelic.recode.vcf",
  title = "",
  name_map_file = "new_names.csv"
)
