#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 4/29/19
# Purpose: Characterize climate niche for different species
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
library(stringr)
library(raster)
select <- dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. ITRDB data
# obs_df <- read_delim(paste0(wdir,'gbif_pinopsida.csv'), delim = '\t')
# obs_df$geom <- st_sfc(obs_df$decimalLatitude, obs_df$decimalLongitude)
tree_db <- paste0(wdir, 'tree_ring_data_V2.db')
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)
species_db = as.data.frame(tbl(conn, 'species'))
tree_db = as.data.frame(tbl(conn, 'trees'))
site_db = as.data.frame(tbl(conn, "sites"))

# 2. Historic climate raster
clim_file <- paste0(wdir, 'HistoricCWD_AETGrids.Rdat')
load(clim_file)
cwd_historic <- sum(cwd_historic)
aet_historic <- sum(aet_historic)
pet_historic <- aet_historic + cwd_historic

# 3. Species range maps
range_file <- paste0(wdir, 'range_maps//merged_ranges.shp')
range_sf <- st_read(range_file)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize ITRDB species frequencies ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spp_lookup <- read.csv(paste0(wdir, "itrdb_species_list.csv"))
spp_lookup <- spp_lookup %>%
  mutate(spp = paste0(genus, " ", species),
         code = tolower(str_replace_all(code, "[*]", "")),
         ncode = nchar(code)) %>%
  arrange(spp) %>%
  as_tibble()

itrdb_species <- species_db %>%
  arrange(species_id)

site_count <- tree_db %>%
  select(species_id, site_id) %>%
  distinct() %>%
  count(species_id) %>%
  arrange(desc(n)) %>%
  select(species_id, n_sites = n)

site_count <- site_count %>%
  left_join(spp_lookup, by = c("species_id" = "code"))

tree_count <- tree_db %>%
  select(species_id, site_id) %>%
  count(species_id) %>%
  arrange(desc(n)) %>%
  select(species_id, n_trees = n)

site_count <- site_count %>%
  left_join(tree_count, by = "species_id")

site_count <- site_count %>%
  select(species_id, genus, species, spp, common_name, n_sites, n_trees)
write.csv2(site_count, paste0(wdir, "species_summary.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Map overlap between ITRDB and range ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define species
spp_code <- 'pilo'

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

# Pull relevant range map
sp_range <- range_sf %>%
  filter(sp_code == spp_code)

# Plot species ranges
world <- ne_coastline(scale = "medium", returnclass = "sf")
ggplot() +
  geom_sf(data = world) +
  geom_sf(data = sp_range, fill = 'blue', alpha = .9) +
  geom_sf(data = sites, color = 'red', fill = 'red', alpha = .5) +
  # coord_sf(crs = st_crs(54019), expand = FALSE)
  coord_sf(xlim = c(-125, -110), ylim = c(30, 45), expand = FALSE) ## Western US
 # coord_sf(xlim = c(-10, 60), ylim = c(30, 55), expand = FALSE) ## Europe

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species niches -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(cwd_historic)
plot(pet_historic)
plot(aet_historic)

cwd_crop <- crop(cwd_historic, sp_range)

cwd_vals <- raster::extract(cwd_historic, sp_range) %>% 
  unlist()
mean_cwd <- mean(cwd_vals)
sd_cwd <- sd(cwd_vals)

aet_vals <_ raster::extract(aet_historic, sp_range) %>% 
  unlist()
mean_aet <- mean(aet_vales)
sd_aet <_ sd(aet_vals)

ggplot(data = cwd_vals %>% as.data.frame(), aes(x=.)) +
  geom_histogram()
