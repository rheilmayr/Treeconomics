#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Quantify missing data from multiple data sources
#
# Input files:
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import libraries ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(sf)
library(naniar)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# Load ITRDB database
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)
tree_db = as.data.frame(tbl(conn, 'trees'))
spp_db = as.data.frame(tbl(conn,'species'))
site_db = as.data.frame(tbl(conn, "sites"))
obs_db = tbl(conn, 'observations_new')

# Load climate data
cwd_csv = paste0(wdir, 'CRU//essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) 
cwd_sites <- cwd_df %>%
  select(site_id) %>% 
  distinct() %>% 
  mutate(valid_clim = 1)

# Load species niches
range_file <- paste0(wdir, 'range_maps//merged_ranges.shp')
range_sf <- st_read(range_file)
st_geometry(range_sf) = NULL
range_species <- range_sf %>%
  select(species_id = sp_code) %>% 
  distinct() %>% 
  mutate(valid_range = 1)

# Load list of valid dendrochronologies
dendro_sites <- read_csv(paste0(wdir, 'out\\rwi_data\\2_valid_sites.csv')) %>% 
  select(site_id) %>% 
  mutate(valid_dendro = 1)

# Load first stage results
first_stage <- read_csv(paste0(wdir, 'out\\first_stage\\', 'log_cwd_pet.csv'))
first_stage <- first_stage %>% 
  select(site_id, species_id) %>% 
  mutate(valid_fs = 1)
 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge different datasets ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create list of valid species / site possibilities
index <- tree_db %>%
  select(species_id, site_id) %>% 
  distinct()

# Merge datasets to index
data_inventory <- index %>% 
  merge(cwd_sites, by = "site_id", all = TRUE) %>% 
  merge(range_species, by = 'species_id', all = TRUE) %>% 
  merge(dendro_sites, by = "site_id", all = TRUE) %>% 
  merge(first_stage, by = c("site_id", "species_id"), all = TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore possible data problems ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim(index)
data_inventory %>% gg_miss_upset()



