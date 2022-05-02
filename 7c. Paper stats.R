#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 7/24/21
# Purpose: Calculate metrics presented in paper
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(raster)
library(sf)
library(tidyverse)
# library(dbplyr)
# library(RSQLite)
# library(ggplot2)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(rgeos)
# library(stringr)
# library(rgdal)
# library(viridis)
# library(patchwork)
# library(Hmisc)
# library(prediction)
# library(colorspace)
# library(ggnewscale)



select <- dplyr::select
summarize <- dplyr::summarize

options(scipen=999)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv")) %>%
  select(-X1)

# 2. Species range maps
range_file <- paste0(wdir, 'in/species_ranges/merged_ranges.shp')
range_sf <- st_read(range_file)

# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)
site_loc <- site_smry %>% 
  select(collection_id, latitude, longitude)
flm_df <- flm_df %>% 
  left_join(site_loc, by = "collection_id") %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))

trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

# 5. Prediction rasters
rwi_list <- list.files(paste0(wdir, "out/predictions/sp_rwi_pred/"), pattern = ".rds", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))
# sp_predictions <- readRDS(paste0(wdir, "out/predictions/sp_predictions.rds"))


# 6. Second stage model
mod_df <- trim_df
cwd_mod <- readRDS(paste0(wdir, "out/second_stage/cwd_mod.rds"))
cwd_vcov <- readRDS(paste0(wdir, "out/second_stage/cwd_mod_vcov.rds"))
pet_mod <- readRDS(paste0(wdir, "out/second_stage/pet_mod.rds"))
int_mod <- readRDS(paste0(wdir, "out/second_stage/int_mod.rds"))

# 7. Genus model predictions
genus_predictions <- readRDS(paste0(wdir, "out/second_stage/genus_mods.rds"))


# 8. Historic climate raster
clim_file <- paste0(wdir, 'in/CRUData/historic_raster/HistoricCWD_AETGrids.Rdat')
load(clim_file)
cwd_historic <- sum(cwd_historic)
names(cwd_historic) = "cwd"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Abstract and intro --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Summary stats of data included in analysis
(n_species <- flm_df %>% 
  pull(species_id) %>% 
  unique() %>% 
  length())

(n_sites <- flm_df %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length())

(n_trees <- flm_df %>% 
  pull(ntrees) %>% 
  sum())

(n_obs <- flm_df %>% 
  pull(nobs) %>% 
  sum(na.rm = TRUE))
### TODO: Some rows have nobs = NA but have estimated sensitivity; need to go back to original processing to resolve this 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Impacts of water and energy on tree growth  ----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## In the median plot, a one standard deviation increase in CWD was associated with an XX percent decline in annual growth"
(cwd_te_median <- flm_df %>%
  pull(estimate_cwd.an) %>% 
  median())

## The site-level, marginal effect of CWD was significantly (p<0.05) negative for XX percent of sites, while only XX 
## percent of sites exhibited a significantly positive relationship. 
cwd_te_count_neg <- flm_df %>%
    filter(p.value_cwd.an<0.05,
           estimate_cwd.an<0) %>% 
    pull(collection_id) %>% 
    unique() %>% 
    length()

(cwd_te_shr_neg <- cwd_te_count_neg / n_sites)

cwd_te_count_pos <- flm_df %>%
  filter(p.value_cwd.an<0.05,
         estimate_cwd.an>0) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

(cwd_te_shr_pos <- cwd_te_count_pos / n_sites)


##  In contrast, a one standard deviation increase in PET was, on average, associated with an XX percent increase in growth. 
(cwd_te_median <- flm_df %>%
    pull(estimate_pet.an) %>% 
    median())

## The relationship between PET and growth was significantly positive for XX percent of sites, and negative for XX percent of sites.
pet_te_count_neg <- flm_df %>%
  filter(p.value_pet.an<0.05,
         estimate_pet.an<0) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

(pet_te_shr_neg <- pet_te_count_neg / n_sites)

pet_te_count_pos <- flm_df %>%
  filter(p.value_pet.an<0.05,
         estimate_pet.an>0) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

(pet_te_shr_pos <- pet_te_count_pos / n_sites)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Relatively dry sites exhibit less sensitivity to drought  --------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


