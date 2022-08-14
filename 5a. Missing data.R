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

# Load species ranges
range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
range_sf <- st_read(range_file)
st_geometry(range_sf) = NULL
range_species <- range_sf %>%
  select(sp_id = sp_code) %>% 
  distinct() %>% 
  mutate(valid_range = 1) %>% 
  drop_na()

# Identify valid dendrochronologies
dendro_df <- read_csv(paste0(wdir, 'out\\dendro\\rwi_long.csv')) %>% 
  select(collection_id) %>% 
  distinct() %>% 
  mutate(valid_dendro = 1)

# Load cwd data
cwd_csv = paste0(wdir, 'out//climate//essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("collection_id" = as.character(site))
cwd_sites <- cwd_df %>%
  select(collection_id) %>% 
  distinct() %>% 
  mutate(valid_clim = 1)

# Load site annual past climate
an_clim <- read_rds(paste0(wdir, "out\\climate\\site_an_clim.gz"))
an_clim <- an_clim %>% 
  select(collection_id) %>%
  distinct() %>% 
  mutate(valid_an_clim = 1)

# Load site average past climate
ave_clim <- read_rds(paste0(wdir, "out\\climate\\site_an_clim.gz"))
ave_clim <- ave_clim %>% 
  select(collection_id) %>%
  distinct() %>% 
  mutate(valid_ave_clim = 1)

# Load first stage results
first_stage <- read_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv'))
first_stage <- first_stage %>% 
  select(collection_id) %>%
  distinct() %>% 
  mutate(valid_fs = 1)
# file.info(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv'))$mtime

# Load site future climate
fut_clim <- read_rds(paste0(wdir, "out/climate/sp_clim_predictions.", compress = "gz"))
fut_clim <- fut_clim %>% 
  select(sp_id = sp_code) %>%
  distinct() %>% 
  mutate(valid_fut_clim = 1)

 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge different datasets ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge datasets to index
data_inventory <- itrdb_sites %>% 
  merge(cwd_sites, by = "collection_id", all = TRUE) %>% 
  merge(range_species, by = "sp_id", all = TRUE) %>% 
  merge(dendro_df, by = "collection_id", all = TRUE) %>% 
  merge(first_stage, by = "collection_id", all = TRUE) %>% 
  merge(ave_clim, by = "collection_id", all = TRUE) %>% 
  merge(an_clim, by = "collection_id", all = TRUE) %>% 
  merge(fut_clim, by = "sp_id", all = TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore possible data problems ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim(itrdb_sites)
data_inventory %>% gg_miss_upset()
miss_var_summary(data_inventory)

## Total number of complete sites in analysis
data_inventory %>% complete.cases() %>% sum()

#### NOTES
# Missing dendro: Largely due to strange format of ITRDB files, leading to parsing errors
# Missing first stage: For ones with valid dendro, missing due to no variation in CWD / PET
# Missing range: A few species dominate combined with sites that don't specify beyond genera

