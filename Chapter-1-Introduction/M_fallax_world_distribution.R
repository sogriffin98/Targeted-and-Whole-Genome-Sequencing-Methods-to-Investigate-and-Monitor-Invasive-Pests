# Set working directory
setwd("~/PhD/Thesis/Figures/Meloidogyne fallax world distribution")

# Load libraries (install once if needed)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Read species presence CSV
species_countries <- read.csv("species_presence.csv",
                              stringsAsFactors = FALSE)

# Fix country naming issues
species_countries <- species_countries %>%
  mutate(country = case_when(
    country == "England" ~ "United Kingdom",
    TRUE ~ country
  )) %>%
  distinct(country) %>%          # remove duplicates
  mutate(present = "Present")

# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Join species presence to world map
world_species <- world %>%
  left_join(species_countries, by = c("name_long" = "country")) %>%
  mutate(present = ifelse(is.na(present), "Absent", present))

# Check for any remaining mismatches
unmatched <- setdiff(species_countries$country, world$name_long)
if (length(unmatched) > 0) {
  message("Unmatched country names:")
  print(unmatched)
}

# Plot map
species_map <- ggplot(world_species) +
  geom_sf(aes(fill = present), color = "gray40", size = 0.2) +
  scale_fill_manual(
    values = c(
      "Present" = "darkgreen",
      "Absent" = "lightgray"
    )
  ) +
  theme_minimal() +
  labs(
    title = "Global distribution of Meloidogyne fallax",
    fill = "Presence"
  )

# Display and Save Map
print(species_map)
ggsave(
  "Meloidogyne_fallax_world_distribution.png",
  species_map,
  width = 10,
  height = 6,
  dpi = 300
)
