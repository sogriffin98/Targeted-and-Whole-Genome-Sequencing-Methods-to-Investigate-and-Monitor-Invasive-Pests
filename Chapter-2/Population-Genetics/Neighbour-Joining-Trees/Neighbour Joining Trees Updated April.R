# Set working directory
setwd("~/PhD/Thesis/Corrections/trees")

# Load packages
library(vcfR)
library(adegenet)
library(ape)
library(poppr)
library(RColorBrewer)

# Bootstrap (500 replicates)
bootstrap_nj <- function(gen, nboot = 500){
  
  X <- gen$tab
  loci_names <- locNames(gen)
  ploidy_val <- ploidy(gen)
  
  boot_fun <- function(xx){
    sampled <- sample(ncol(xx), replace = TRUE)
    X_boot <- xx[, sampled, drop = FALSE]
    
    # Convert back to genind
    gen_boot <- new("genind",
                    tab = X_boot,
                    ploidy = ploidy_val,
                    type = "codom",
                    pop = pop(gen),
                    strata = gen@strata,
                    loc.names = loci_names[sampled])
    
    dist_boot <- bitwise.dist(gen_boot)
    nj(dist_boot)
  }
  
  boot.phylo(
    phy = nj(bitwise.dist(gen)),
    x = X,
    FUN = boot_fun,
    B = nboot
  )
}

# Plotting the nj tree
plot_nj_tree <- function(
    vcf_file,
    title = "",
    metadata_file = NULL,
    add_legend = FALSE,
    name_map_file = NULL,
    colour_by = "country"
){
  
  # Read VCF
  vcf <- read.vcfR(vcf_file)
  gen <- vcfR2genind(vcf)
  
  # Clean sample names
  ind_names_clean <- gsub("^\\.\\/bam\\/", "", indNames(gen))
  ind_names_clean <- gsub("\\.mq30\\.maxinsert\\.primaryalignment\\.bam$", "", ind_names_clean)
  indNames(gen) <- ind_names_clean
  
  # Distance matrix and NJ tree
  dist_matrix <- bitwise.dist(gen)
  nj_tree <- nj(dist_matrix)
  
  # Rename sample names from VCF to uniform names using the new_names.csv
  if(!is.null(name_map_file)){
    name_map <- read.csv(name_map_file, stringsAsFactors = FALSE)
    
    # Only rename samples that actually exist in the tree
    valid <- name_map$Current %in% nj_tree$tip.label
    lookup <- setNames(name_map$New[valid], name_map$Current[valid])
    
    nj_tree$tip.label <- ifelse(
      nj_tree$tip.label %in% names(lookup),
      lookup[nj_tree$tip.label],
      nj_tree$tip.label
    )
  }
  
  # Colour the sample names by country or year
  tip_colors <- rep("black", nInd(gen))
  legend_labels <- NULL
  legend_colors <- NULL
  
  if(!is.null(metadata_file)){
    meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
    meta$Sample <- trimws(meta$Sample)
    match_idx <- match(indNames(gen), meta$Sample)
    
    if(colour_by == "country"){
      
      countries <- meta$Country[match_idx]
      countries[is.na(countries)] <- "Unknown"
      
      unique_countries <- unique(countries)
      other_countries <- setdiff(unique_countries, "UK")
      
      n_other <- length(other_countries)
      other_colors <- if(n_other <= 8) brewer.pal(n_other, "Dark2") else colorRampPalette(brewer.pal(8, "Dark2"))(n_other)
      
      country_colors <- c("UK" = "black", setNames(other_colors, other_countries))
      
      if(!"Unknown" %in% names(country_colors) && any(countries == "Unknown")){
        country_colors["Unknown"] <- "grey"
      }
      
      tip_colors <- country_colors[countries]
      legend_labels <- names(country_colors)
      legend_colors <- country_colors
      
    } else if(colour_by == "year"){
      
      years <- meta$Year[match_idx]
      
      # Sort years numerically
      unique_years <- sort(unique(years))
      
      # Generate a palette in chronological order
      year_palette <- colorRampPalette(brewer.pal(8, "Set1"))
      year_colors <- setNames(year_palette(length(unique_years)), unique_years)
      
      # Assign colours
      tip_colors <- year_colors[as.character(years)]
      
      # Legend in chronological order
      legend_labels <- unique_years
      legend_colors <- year_colors
      
    }
  }
  
  # Bootstrap support (500 replicates)
  
  boot_values <- bootstrap_nj(gen, nboot = 500)
  nj_tree$node.label <- boot_values
  
  # Plot the NJ trees
  tip_cex <- 1.5 / log10(nInd(gen) + 1)
  top_margin <- ifelse(title == "", 1, 4)
  
  par(mar = c(5, 4, top_margin, 10))
  
  plot(nj_tree,
       tip.color = tip_colors,
       cex = tip_cex,
       main = title,
       x.lim = c(0, max(node.depth.edgelength(nj_tree)) * 1.3))
  
  add.scale.bar(x = 0, y = -0.5,
                length = signif(max(nj_tree$edge.length) / 5, 3))
  
  nodelabels(text = nj_tree$node.label,
             cex = 0.6,
             frame = "none",
             adj = c(1.2, -0.2))
  
  if(add_legend && !is.null(legend_labels)){
    par(xpd = TRUE)
    legend("right",
           inset = c(-0.4, 0),
           legend = legend_labels,
           col = legend_colors,
           pch = 19,
           cex = 0.7,
           bty = "n")
  }
  
  return(nj_tree)
}

# Plot both trees
# UK + EU tree (colour by country)
nj_tree_uk_eu <- plot_nj_tree(
  vcf_file = "AH_104inds_minDP6_mac0.05_biallelic.recode.vcf",
  metadata_file = "metadata.csv",
  add_legend = TRUE,
  name_map_file = "new_names.csv",
  colour_by = "country"
)

# UK-only tree (colour by year)
nj_tree_uk <- plot_nj_tree(
  vcf_file = "AH_UK_89inds_minDP6_mac0.05_biallelic.recode.vcf",
  metadata_file = "metadata.csv",
  add_legend = TRUE,
  name_map_file = "new_names.csv",
  colour_by = "year"
)

# Save trees as pdfs
pdf("UK_EU_tree.pdf", width = 16, height = 10)   # wide canvas

par(mar = c(5, 4, 4, 20), xpd = TRUE)   # ← big right margin + allow drawing outside

plot_nj_tree(
  vcf_file = "AH_104inds_minDP6_mac0.05_biallelic.recode.vcf",
  metadata_file = "metadata.csv",
  add_legend = TRUE,
  name_map_file = "new_names.csv",
  colour_by = "country"
)
dev.off()

pdf("UK_tree.pdf", width = 16, height = 10)

par(mar = c(5, 4, 4, 20), xpd = TRUE)

plot_nj_tree(
  vcf_file = "AH_UK_89inds_minDP6_mac0.05_biallelic.recode.vcf",
  metadata_file = "metadata.csv",
  add_legend = TRUE,
  name_map_file = "new_names.csv",
  colour_by = "year"
)

dev.off()


