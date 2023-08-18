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
summarize <- dplyr::summarize


# Convert site summary valid values to na/1
replace_valids = function(value){
  value = ifelse(value == TRUE, 1, NA)
  value = as.double(value)
  return(value)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# Load list of ITRDB sites
itrdb_sites <- read_csv(paste0(wdir, '1_input_processed/dendro/site_summary.csv')) %>% 
  mutate(sp_id = str_to_lower(sp_id))

itrdb_sites <- itrdb_sites %>% 
  filter(datasource %in% c("ritrdb_7.13", "itrdb_7.22"))  # Narrow to set of sites included in ITRDB or rITRDB

# Load species ranges
range_file <- paste0(wdir, '1_input_processed/species_ranges/merged_ranges_dissolve.shp')
range_sf <- st_read(range_file)
st_geometry(range_sf) = NULL
range_species <- range_sf %>%
  select(sp_id = sp_code) %>% 
  distinct() %>% 
  mutate(valid_range = 1) %>% 
  drop_na()

# Identify valid dendrochronologies
dendro_df <- read_csv(paste0(wdir, '2_output/dendro/rwi_long.csv')) %>% 
  select(collection_id) %>% 
  distinct() %>% 
  mutate(valid_dendro = 1)

# Load cwd data
cwd_csv = paste0(wdir, '2_output/climate/essentialcwd_data.csv')
cwd_df <- read_csv(cwd_csv)
cwd_df <- cwd_df %>% 
  mutate("collection_id" = as.character(site))
cwd_sites <- cwd_df %>%
  select(collection_id) %>% 
  distinct() %>% 
  mutate(valid_clim = 1)

# Load site annual climate
an_clim <- read_rds(paste0(wdir, "2_output/climate/site_an_clim.gz"))
vary_clim <- an_clim %>% 
  group_by(collection_id) %>% 
  summarize(cwd_n = n_distinct(cwd.an.spstd),
            pet_n = n_distinct(pet.an.spstd)) %>% 
  mutate(valid_weather_vary = pet_n > 1 & cwd_n > 1,
         valid_weather_vary = replace_valids(valid_weather_vary)) %>% 
  select(collection_id, valid_weather_vary)

an_clim <- an_clim %>% 
  select(collection_id) %>%
  distinct() %>% 
  mutate(valid_an_clim = 1)

# Load site average past climate
ave_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz"))
ave_clim <- ave_clim %>% 
  select(collection_id) %>%
  distinct() %>% 
  mutate(valid_ave_clim = 1)

# Load first stage results
first_stage <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_std.csv'))
first_stage <- first_stage %>% 
  select(collection_id) %>%
  distinct() %>% 
  mutate(valid_fs = 1)

# Load site future climate
fut_clim <- read_rds(paste0(wdir, "2_output/climate/sp_clim_predictions.", compress = "gz"))
fut_clim <- fut_clim %>% 
  select(sp_id = sp_code) %>%
  distinct() %>% 
  mutate(valid_fut_clim = 1)

 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge different datasets ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itrdb_sites <- itrdb_sites %>% 
  mutate(valid_dates = map_dbl(valid_dates, replace_valids),
         valid_parse = map_dbl(valid_parse, replace_valids),
         valid_rwl = map_dbl(valid_rwl, replace_valids)) %>% 
  select(collection_id, datasource, valid_dates, valid_parse, valid_rwl, sp_id)

# Merge datasets to index
data_inventory <- itrdb_sites %>% 
  left_join(cwd_sites, by = "collection_id") %>% 
  left_join(range_species, by = "sp_id") %>% 
  left_join(dendro_df, by = "collection_id") %>% 
  left_join(first_stage, by = "collection_id") %>% 
  left_join(ave_clim, by = "collection_id") %>% 
  left_join(an_clim, by = "collection_id") %>% 
  left_join(fut_clim, by = "sp_id") %>% 
  left_join(vary_clim, by = "collection_id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore possible data problems ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Number of sites included in either rITRDB 7.13 or ITRDB 7.22 where metadata 
# includes location, indicates they span into the study period (>1901) and have species label: 4362
data_inventory %>% 
  pull(collection_id) %>% 
  length()

# Number of sites remaining after RWL parse and second check for data from the study period: 4312 
data_inventory %>% 
  filter(valid_dates == 1, valid_parse == 1, valid_rwl == 1) %>% 
  pull(collection_id) %>% 
  length()

# Confirm this matches the output dendro data: 4312
dendro_df %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

# Sites dropped due to missing range maps: 448
# Sites remaining for anlaysis: 3864
data_inventory %>% 
  filter(valid_dates == 1, valid_parse == 1, valid_rwl == 1) %>% 
  group_by(valid_range) %>% 
  tally()


# Sites dropped due to missing raw climate/topo data: 20
# Sites remaining for analysis: 3844
data_inventory %>% 
  filter(valid_dates == 1, valid_parse == 1, valid_rwl == 1, valid_range == 1) %>% 
  group_by(valid_clim) %>% 
  tally()


# Sites dropped due to missing species-standardized climate data: None
# Sites remaining for analysis: 3844
data_inventory %>% 
  filter(valid_dates == 1, valid_parse == 1, valid_rwl == 1, valid_range == 1, valid_clim==1) %>% 
  group_by(valid_ave_clim, valid_an_clim, valid_fut_clim) %>% 
  tally()


# Sites dropped due to errors in first stage regression: 
# - Error 1: No variation in PET or CWD: 61
# - Error 2: Too few observations to identify regression: 8
# Sites remaining for analysis: 3775
data_inventory %>% 
  filter(valid_dates == 1, valid_parse == 1, valid_rwl == 1, valid_range == 1, valid_clim==1) %>% 
  group_by(valid_fs, valid_weather_vary) %>% 
  tally()


# Number of species / trees / observations included
flm_df$collection_id %>% unique() %>%  length()
flm_df$species_id %>% unique() %>% length()
flm_df$ntrees %>% sum()
flm_df$nobs %>% sum()





# Sites dropped in primary specification as outliers in first stage regression estimates:


data_inventory %>% gg_miss_upset()
miss_var_summary(data_inventory)

## Total number of complete sites in analysis
data_inventory %>% complete.cases() %>% sum()


data_inventory %>% 
  group_by(valid_rwl, valid_dendro, valid_range, valid_fs) %>% 
  tally()

