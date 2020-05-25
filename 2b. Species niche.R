#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Characterize climate niche for different species
#
# Input files:
# - HistoricCWD_AETGrids.Rdat: Rasters describing historic CWD and AET
#     Generated using historic_cwdraster.R
# - 
#
# Todo ideas:
# - Could use both spatial and annual variation in CWD / AET to characterize niche
#
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

# 1. Historic climate raster
clim_file <- paste0(wdir, 'HistoricCWD_AETGrids.Rdat')
load(clim_file)
cwd_historic <- sum(cwd_historic)
aet_historic <- sum(aet_historic)
pet_historic <- aet_historic + cwd_historic

# 2. Species range maps
range_file <- paste0(wdir, 'range_maps//merged_ranges.shp')
range_sf <- st_read(range_file)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species niches -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull and organize climate distribution for species
pull_clim <- function(spp_code){
  # Pull relevant range map
  sp_range <- range_sf %>%
    filter(sp_code == spp_code)
  
  # Pull cwd and aet values
  cwd_vals <- raster::extract(cwd_historic, sp_range) %>% 
    unlist()
  aet_vals <- raster::extract(aet_historic, sp_range) %>% 
    unlist()
  
  # Combine into tibble
  clim_vals <- data.frame(cwd_vals, aet_vals)
  names(clim_vals) <- c('cwd', 'aet')
  clim_vals <- clim_vals %>% 
    mutate(pet = aet+cwd) %>% 
    as_tibble()
  return(clim_vals)
}

clim_df <- range_sf %>%
  pull(sp_code) %>% 
  unique() %>% 
  enframe(name = NULL) %>% 
  rename(sp_code = value)

clim_df$clim_vals <- map(clim_df$sp_code, pull_clim)

clim_df <- clim_df %>% 
  unnest()

write.csv(clim_df, paste0(wdir, "clim_niche.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CODE BELOW THIS POINT SHOULD BE MOVED TO DIFFERENT SCRIPT ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize ITRDB species frequencies ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
spp_code <- 'psme'

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


xlims <- c(-132, -102)
ylims <- c(26, 56)

# Plot species ranges
world <- ne_coastline(scale = "medium", returnclass = "sf")
map <- ggplot() +
  geom_sf(data = world) +
  geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  geom_sf(data = sites, color = 'darkblue', fill = 'red', alpha = .8) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US
map


library(raster)
library(rgdal)
library(viridis)


cwd_historic_df <- cwd_historic %>% 
  # mask(sp_range) %>%
  as.data.frame(xy = TRUE) %>% 
  drop_na() %>% 
  filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
  mutate(CWD = layer,
         cwd.mean = mean(layer),
         cwd.sd = sd(layer),
         std_cwd = (layer - cwd.mean) / cwd.sd)
map <- ggplot() +
  geom_raster(data = cwd_historic_df, aes(x = x, y = y, fill = CWD)) +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  scale_fill_viridis_c() +
  geom_sf(data = sites, color = 'red', alpha = 1) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  geom_sf(data = world, color = 'white') +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US

map

flm_df <- read_csv(paste0(wdir, 'first_stage\\log_cwd_pet.csv')) %>%
  select(-X1)

sp_regs <- flm_df %>% 
  filter(species_id == spp_code) %>% 
  select(site_id, estimate_cwd.an, std.error_cwd.an, p.value_cwd.an)


site_regs <- sites %>% 
  left_join(sp_regs, by = c("site_id")) %>% 
  mutate(ln_cwd = log(estimate_cwd.an)) %>% 
  drop_na()
  


map <- ggplot() +
  # geom_raster(data = cwd_historic_df, aes(x = x, y = y, fill = std_cwd)) +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  geom_sf(data = site_regs, aes(color = ln_cwd), alpha = 1) +
  scale_fill_viridis_c() +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  geom_sf(data = world, color = 'white', size = 1) +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US

map



# coord_sf(crs = st_crs(54019), expand = FALSE)
# coord_sf(xlim = c(-10, 60), ylim = c(30, 55), expand = FALSE) ## Europe

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot climate overlap of observations and species range -----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(cwd_historic)
plot(pet_historic)
plot(aet_historic)


ggplot(data = aet_vals %>% as.data.frame(), aes(x=.)) +
  geom_histogram()
ggplot(data = cwd_vals %>% as.data.frame(), aes(x=.)) +
  geom_histogram()

hex <- clim_df %>% 
  ggplot(aes(x = cwd, y = aet)) +
  geom_hex()
hex