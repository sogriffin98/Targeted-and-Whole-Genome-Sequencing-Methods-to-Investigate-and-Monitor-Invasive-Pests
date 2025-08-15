# Set working directory
setwd("~/PhD/Study Species/Vespa velutina/Mapping of outbreaks/R code 23rd June 2025")

# Load required libraries
install.packages(c("tidyverse", "sf", "readxl", "janitor", "tidygeocoder", "rnaturalearthdata", "ggspatial", "cowplot"))

library(tidyverse)
library(sf)
library(readxl)
library(janitor)
library(tidygeocoder)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(cowplot)

# Load and clean Excel data
data <- read_excel("mapping_of_hornets.xlsx") %>%
  clean_names()

# Geocode unique locations (OpenStreetMap)
unique_locations <- data %>%
  select(location) %>%
  distinct() %>%
  geocode(address = location, method = "osm", lat = latitude, long = longitude)

# Join geocoded coordinates back to main data
data_geo <- left_join(data, unique_locations, by = "location")

# Save coordinates to avoid re-geocoding later
write_csv(unique_locations, "geocoded_locations.csv")

# Convert to sf object (spatial)
data_sf <- st_as_sf(data_geo, coords = c("longitude", "latitude"), crs = 4326)
data_sf$year <- as.factor(data_sf$year)
data_sf$type <- as.factor(data_sf$type)

# Load base UK map
uk_map <- ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf")

# Define South East bounding box for inset
se_bbox <- st_bbox(c(xmin = 0, xmax = 1.5, ymin = 50.8, ymax = 52), crs = st_crs(4326))

# Crop for inset
se_data_sf <- st_crop(data_sf, se_bbox)
se_uk_map <- st_crop(uk_map, se_bbox)

# Define a fixed color scale for years to use in both maps
color_scale <- scale_color_manual(
  values = c(
    "2016" = "#F8766D", "2017" = "#D89000", "2018" = "#A3A500",
    "2019" = "#39B600", "2020" = "#00BF7D", "2021" = "#00BFC4",
    "2022" = "#00B0F6", "2023" = "#9590FF", "2024" = "#E76BF3",
    "2025" = "#FF62BC"
  )
)

# Create main map WITHOUT title
main_map <- ggplot() +
  geom_sf(data = uk_map, fill = "grey96", color = "gray60") +
  geom_sf(data = data_sf, aes(color = year, shape = type), size = 2, alpha = 0.8) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-8, 2), ylim = c(49.5, 59), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(color = "Year", shape = "Type") +
  color_scale

inset_map <- ggplot() +
  geom_sf(data = se_uk_map, fill = "grey96", color = "gray60") +
  geom_sf(data = se_data_sf, aes(color = year, shape = type), size = 2, alpha = 0.8) +
  coord_sf(xlim = c(0, 1.5), ylim = c(50.8, 52)) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(title = "South East England") +
  color_scale

# Display the two maps separately
print(main_map)
print(inset_map)

# Save maps
ggsave("hornet_main_map.png", plot = main_map, width = 10, height = 7, dpi = 300)
ggsave("hornet_inset_map.png", plot = inset_map, width = 7, height = 7, dpi = 300)