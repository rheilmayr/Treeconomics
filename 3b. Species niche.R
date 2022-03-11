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
cwd_historic <- mean(cwd_historic)
aet_historic <- mean(aet_historic)
pet_historic <- aet_historic + cwd_historic
names(cwd_historic) = "cwd"
names(pet_historic) = "pet"
clim_historic <- raster::brick(list(cwd_historic, pet_historic))

# 2. Species range maps
range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
range_sf <- st_read(range_file)

# 3. Climate projections from CMIP5
cmip <- load(paste0(wdir, 'in\\CMIP5 CWD\\cmip5_cwdaet_end.Rdat'))
pet_raster <- aet_raster + cwd_raster
pet_future <- pet_raster
cwd_future <- cwd_raster
rm(pet_raster)
rm(cwd_raster)
rm(aet_raster)


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
  mutate(clim_vals = map(sp_code, pull_clim))

clim_df <- clim_df %>% 
  unnest(clim_vals)

niche_df <- clim_df %>% 
  drop_na() %>% 
  group_by(sp_code) %>% 
  summarize(pet_mean = mean(pet),
            pet_sd = sd(pet),
            cwd_mean = mean(cwd),
            cwd_sd = sd(cwd))

write.csv(niche_df, paste0(wdir, "out//climate//clim_niche.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create species-standardized historic climate ------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rasterize_spstd <- function(spp_code, clim_raster){
  sp_niche <- niche_df %>%
    drop_na() %>% 
    filter(sp_code == spp_code) %>% 
    select(-sp_code) 
  
  sp_range <- range_sf %>% 
    filter(sp_code == spp_code)
  pet_sp <- clim_raster %>%
    subset("pet") %>% 
    mask(sp_range)
  pet_spstd <- (pet_sp - sp_niche$pet_mean) / sp_niche$pet_sd
  names(pet_spstd) = "pet.spstd"
  
  cwd_sp <- clim_raster %>% 
    subset("cwd") %>% 
    mask(sp_range)
  cwd_spstd <- (cwd_sp - sp_niche$cwd_mean) / sp_niche$cwd_sd
  names(cwd_spstd) = "cwd.spstd"
  
  clim.spstd <- raster::brick(list(cwd_spstd, pet_spstd))
  return(clim.spstd)
}

## Create tibble of species
species_list <- niche_df %>% 
  select(sp_code)

## Create tibble of historic climates
sp_historic <- species_list %>% 
  mutate(clim_historic_sp = map(.x = sp_code, clim_raster = clim_historic, .f = rasterize_spstd)) %>% 
  as_tibble()

## Export historic
saveRDS(sp_historic, paste0(wdir, "out//climate//sp_clim_historic.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull CMIP projections -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pull_cmip_model <- function(cmip_idx){
  cwd_lyr <- cwd_future[[cmip_idx]]
  names(cwd_lyr) = "cwd"
  pet_lyr <- pet_future[[cmip_idx]]
  names(pet_lyr) = "pet"
  future_clim <- brick(list(cwd_lyr, pet_lyr))
  return(future_clim)
}

n_cmip_mods <- dim(pet_future)[3]
cmip_projections <- tibble(cmip_idx = 1:n_cmip_mods)
cmip_projections <- cmip_projections %>% 
  mutate(cmip_rast = map(cmip_idx, pull_cmip_model))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize future climate for each species ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Cross species list with cmip models
sp_fut_clim <- species_list %>% 
  crossing(cmip_idx = seq(n_cmip_mods)) %>% 
  left_join(cmip_projections, by = "cmip_idx")


## Rescale each CMIP model based on each species' climate
sp_fut_clim <- sp_fut_clim %>% 
  mutate(clim_future_sp = pmap(list(spp_code = sp_code, 
                                    clim_raster = cmip_rast), 
                               .f = rasterize_spstd))

## Select down to base columns
sp_fut_clim <- sp_fut_clim %>% 
  select(sp_code, cmip_idx, clim_future_sp)

## Export predictions
saveRDS(sp_fut_clim, paste0(wdir, "out//climate//sp_clim_predictions.rds"))

