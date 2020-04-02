library(tidyverse)
library(dbplyr)
library(RSQLite)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
library(stringi)

# Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
setwd(wdir)

# Read data
obs_df <- read_delim(paste0(wdir,'gbif_pinopsida.csv'), delim = '\t')
tree_db = paste0(wdir, 'tree_ring_data_V2.db')

# Connect to database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)

# Open tables
species_db = as.data.frame(tbl(conn, 'species'))
tree_db = as.data.frame(tbl(conn, 'trees'))
site_db = as.data.frame(tbl(conn, "sites"))
obs_df$geom <- st_sfc(obs_df$decimalLatitude, obs_df$decimalLongitude)

spp_lookup <- read.csv(paste0(wdir, "itrdb_species_list.csv"))
spp_lookup <- spp_lookup %>%
  mutate(spp = paste0(genus, " ", species),
         code = str_replace_all(code, "*", "")) %>%
  arrange(spp)

# Explore itrdb species frequency
itrdb_species <- species_db %>%
  arrange(species_id)

site_count <- tree_db %>%
  select(species_id, site_id) %>%
  distinct() %>%
  count(species_id) %>%
  arrange(desc(n))

# Define species
spp <- "Pseudotsuga menziesii"
spp_code <- "psme"

spp <- "Pinus sylvestris"
spp_code <- "pisy"

spp <- "Cedrus libani"
spp_code <- "cdli"

spp <- "Pinus lambertiana"
spp_code <- "pila"

# Pull relevant ITRDB sites
site_list <- tree_db %>%
  filter(species_id == spp_code) %>%
  pull(site_id) %>%
  unique()

sites <- site_db %>%
  filter(site_id %in% site_list) %>%
  select(site_id, latitude, longitude) %>%
  drop_na() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Pull relevant GBIF observations
obs_spp <- obs_df %>%
  filter(species == spp) %>%
  select(species, decimalLatitude, decimalLongitude) %>%
  drop_na() %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)

# Plot species ranges
world <- ne_coastline(scale = "medium", returnclass = "sf")
ggplot() +
  geom_sf(data = world) +
  geom_sf(data = obs_spp, alpha = .15) + 
  geom_sf(data = sites, color = 'red', fill = 'red', alpha = .5) +
  coord_sf(crs = st_crs(54019), expand = FALSE)
 # coord_sf(xlim = c(-125, -110), ylim = c(30, 45), expand = FALSE) ## Western US
 # coord_sf(xlim = c(-10, 60), ylim = c(30, 55), expand = FALSE) ## Europe
# geom_hex(data = obs_spp, aes(x = decimalLongitude, y = decimalLatitude)) +
