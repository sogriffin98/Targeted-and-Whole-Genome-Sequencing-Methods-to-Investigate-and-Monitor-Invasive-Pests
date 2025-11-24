# Install required packages
install.packages(c("tidygeocoder", "geosphere", "dplyr"))

# Load libraries
library(tidygeocoder)
library(geosphere)
library(dplyr)

# 1. Read in CSV file
df <- read.csv("locations_all.csv", stringsAsFactors = FALSE)

# 2. Geocode both origin and destination using OpenStreetMap (free, no API key)
df <- df %>%
  geocode(origin_town, method = "osm", lat = origin_lat, long = origin_lon) %>%
  geocode(destination_town, method = "osm", lat = dest_lat, long = dest_lon)

# 3. Compute great-circle ("as-the-crow-flies") distance
df$distance_m <- distHaversine(
  p1 = df[, c("origin_lon", "origin_lat")],
  p2 = df[, c("dest_lon", "dest_lat")]
)

# Convert to kilometers
df$distance_km <- df$distance_m / 1000

# 4. Save the results
write.csv(df, "locations_all_with_distances.csv", row.names = FALSE)

# Preview
head(df)
