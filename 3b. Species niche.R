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
# - Could use both spatial and annual variation in CWD / AET to characterize niche - Waiting on new data from Fran (5/25/20)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(ggplot2)
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
clim_file <- paste0(wdir, 'in//CRUData//historic_raster//HistoricCWD_AETGrids_Annual.Rdat')
load(clim_file)
pet_historic <- aet_historic + cwd_historic

# 2. Species range maps
range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
range_sf <- st_read(range_file)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species niches -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull and organize climate distribution for species
pull_clim <- function(spp_code){
  print(spp_code)
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
  rename(sp_code = value) %>% 
  drop_na()

clim_df <- clim_df %>% 
  mutate(clim_vals = map(clim_df$sp_code, pull_clim))

clim_df <- clim_df %>% 
  unnest(clim_vals)

clim_df <- clim_df %>% 
  drop_na() %>% 
  group_by(sp_code) %>% 
  summarize(pet_mean = mean(pet),
            pet_sd = sd(pet),
            cwd_mean = mean(cwd),
            cwd_sd = sd(cwd))

write.csv(clim_df, paste0(wdir, "out//climate//clim_niche.csv"))




