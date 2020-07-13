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
library(stringr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# Load list of ITRDB sites
itrdb_sites <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv')) %>% 
  select(collection_id, sp_id) %>% 
  distinct() %>% 
  mutate(sp_id = str_to_lower(sp_id)) %>% 
  drop_na()

# Load climate data
cwd_csv = paste0(wdir, 'out//climate//essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("collection_id" = as.character(site))
cwd_sites <- cwd_df %>%
  select(collection_id) %>% 
  distinct() %>% 
  mutate(valid_clim = 1)

# Load species niches
range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
range_sf <- st_read(range_file)
st_geometry(range_sf) = NULL
range_species <- range_sf %>%
  select(species_id = sp_code) %>% 
  distinct() %>% 
  mutate(valid_range = 1)

# Identify valid dendrochronologies
dendro_df <- read_csv(paste0(wdir, 'out\\dendro\\rwi_long.csv')) %>% 
  select(collection_id) %>% 
  distinct() %>% 
  mutate(valid_dendro = 1)

# Load first stage results
first_stage <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv'))
first_stage <- first_stage %>% 
  select(collection_id) %>%
  distinct() %>% 
  mutate(valid_fs = 1)
 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge different datasets ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge datasets to index
data_inventory <- itrdb_sites %>% 
  merge(cwd_sites, by = "collection_id", all = TRUE) %>% 
  merge(range_species, by.x = "sp_id", by.y = 'species_id', all = TRUE) %>% 
  merge(dendro_df, by = "collection_id", all = TRUE) %>% 
  merge(first_stage, by = "collection_id", all = TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore possible data problems ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim(itrdb_sites)
data_inventory %>% gg_miss_upset()
miss_var_summary(data_inventory)



#### NOTES
# Missing dendro: Largely due to strange format of ITRDB files, leading to parsing errors
# Missing first stage: For ones with valid dendro, missing due to no variation in CWD / PET
# Missing range: A few species dominate combined with sites that don't specify beyond genera
# Climate: Should dig in to this again, think I am missing some data here

