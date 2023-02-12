#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Characterize climate niche for different species
#
# Input files:
# - HistoricCWD_AETGrids.Rdat: Rasters describing historic CWD and AET
#     Generated using historic_cwdraster.R
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
library(readr)
library(tmap)
library(tictoc)
select <- dplyr::select


library(furrr)
n_cores <- 6
future::plan(multisession, workers = n_cores)


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

# 2. Site-specific historic climate data
site_clim_csv <- paste0(wdir, 'out\\climate\\essentialcwd_data.csv')
site_clim_df <- read_csv(site_clim_csv)
site_clim_df <- site_clim_df %>% 
  mutate("site_id" = as.character(site)) %>% 
  rename(collection_id = site_id)

# 3. Load species information for sites
site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id) %>% 
  mutate(sp_code = tolower(sp_id)) %>% 
  select(-sp_id)


# 4. Species range maps
range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
range_sf <- st_read(range_file)

# 5. Climate projections from CMIP5
cmip_end <- load(paste0(wdir, 'in\\CMIP5 CWD\\cmip5_cwdaet_end.Rdat'))
pet_cmip_end <- aet_raster + cwd_raster
cwd_cmip_end <- cwd_raster
names(cwd_cmip_end) <- NULL # Resetting this due to strange names in file from CMIP processing
rm(cwd_raster)
rm(aet_raster)

cmip_start <- load(paste0(wdir, 'in\\CMIP5 CWD\\cmip5_cwdaet_start.Rdat'))
pet_cmip_start <- aet_raster + cwd_raster
cwd_cmip_start <- cwd_raster
names(cwd_cmip_start) <- NULL # Resetting this due to strange names in file from CMIP processing
rm(cwd_raster)
rm(aet_raster)


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Visually inspect data -----------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# tmap_mode("view")
# tm_shape(cwd_cmip_end) +
#   tm_raster() +
#   tm_facets(as.layers = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species niches -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull and organize climate distribution for species
pull_clim <- function(spp_code){
  print(spp_code)
  # Pull relevant range map
  sp_range <- range_sf %>%
    filter(sp_code == spp_code) %>% 
    rasterize(cwd_historic, getCover=TRUE)
  sp_range[sp_range==0] <- NA

  # Pull cwd and aet values
  cwd_vals <- cwd_historic %>% 
    mask(sp_range) %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  
  pet_vals <- pet_historic %>% 
    mask(sp_range) %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  
  # Combine into tibble
  clim_vals <- cwd_vals %>% 
    left_join(pet_vals, by = c("x", "y"))

  return(clim_vals)
}


species_list <- range_sf %>%
  pull(sp_code) %>% 
  unique() %>% 
  enframe(name = NULL) %>% 
  rename(sp_code = value) %>% 
  arrange(sp_code) %>% 
  drop_na()

clim_df <- species_list %>% 
  mutate(clim_vals = future_map(sp_code, 
                                .f = pull_clim,
                                .options = furrr_options(packages = c( "dplyr", "raster", "sf")),
                                .progress = TRUE))


## Summarize mean and sd of each species' climate
niche_df <- clim_df %>% 
  unnest(clim_vals) %>% 
  group_by(sp_code) %>% 
  summarize(pet_mean = mean(pet),
            pet_sd = sd(pet),
            cwd_mean = mean(cwd),
            cwd_sd = sd(cwd))


## Export species niche description
write.csv(niche_df, paste0(wdir, "out//climate//clim_niche.csv"))
# niche_df <- read.csv(paste0(wdir, "out//climate//clim_niche.csv")) %>%
#   select(-X)


# ## Export historic climates
# write_rds(clim_df, paste0(wdir, "out//climate//sp_clim_historic.gz"), compress = "gz")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Standardize historic climate -------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_standardize <- function(val, sp_mean, sp_sd){
  std_val <- (val - sp_mean) / sp_sd
  return(std_val)
}

sp_std_historic_df <- function(hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd){
  hist_clim_vals <- hist_clim_vals %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd))
  return(hist_clim_vals)
}


sp_std_future_df <- function(cmip_df, hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd){
  valid_locations <- hist_clim_vals %>% select(x,y)
  cmip_df <- valid_locations %>% 
    left_join(cmip_df, by = c("x", "y"))
  cmip_df <- cmip_df %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd))
  return(cmip_df)
}


clim_df <- clim_df %>% 
  left_join(niche_df, by = "sp_code")

clim_df <- clim_df %>% 
  mutate(clim_historic_sp = future_pmap(list(hist_clim_vals = clim_vals,
                                             pet_mean = pet_mean,
                                             pet_sd = pet_sd,
                                             cwd_mean = cwd_mean,
                                             cwd_sd = cwd_sd),
                                        .f = sp_std_historic_df,
                                        .options = furrr_options(packages = c( "dplyr"))))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull CMIP projections -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pet_end_df <- pet_cmip_end %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "pet_cmip_end"), 
              starts_with('layer.'))

cwd_end_df <- cwd_cmip_end %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "cwd_cmip_end"), 
              starts_with('layer.'))


pet_start_df <- pet_cmip_start %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "pet_cmip_start"), 
              starts_with('layer.'))

cwd_start_df <- cwd_cmip_start %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "cwd_cmip_start"), 
              starts_with('layer.'))


# ## Illustrate forawrd/backward conversion between df and raster
# crs_template <- crs(cwd_future)
# raster_template <- cwd_df %>% select(x,y)
# cwd_df <- cwd_df %>% 
#   drop_na()
# cwd_df2 <- raster_template %>% 
#   left_join(cwd_df, by = c("x", "y"))
# cwd_rast2 <- rasterFromXYZ(cwd_df2, crs = crs_template)

## Combine PET and CWD projections
cmip_df <- cwd_end_df %>% 
  full_join(pet_end_df, by = c("x", "y")) %>% 
  full_join(cwd_start_df, by = c("x", "y")) %>% 
  full_join(pet_start_df, by = c("x", "y"))
  
## Nest CMIP data
cmip_df <- cmip_df %>%
  mutate(idx = 1) %>% 
  group_by(idx) %>% 
  nest() %>% 
  ungroup() %>% 
  select(data)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize cmip climate for each species ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Cross species list with nested cmip data
sp_cmip_clim <- clim_df %>% 
  mutate(cmip_df = cmip_df$data)



sp_cmip_clim <- sp_cmip_clim %>% 
  mutate(clim_cmip_sp = future_pmap(list(cmip_df = cmip_df,
                                    hist_clim_vals = clim_vals,
                                    pet_mean = pet_mean,
                                    pet_sd = pet_sd,
                                    cwd_mean = cwd_mean,
                                    cwd_sd = cwd_sd), 
                               .f = sp_std_future_df,
                               .options = furrr_options(packages = c( "dplyr")))) %>% 
  select(-cmip_df)


# ## Check final result as raster
# species = "acsh"
# test_clim <- (sp_cmip_clim %>% filter(sp_code == species) %>% pull(clim_cmip_sp))[[1]]
# crs_template <- crs(cwd_cmip_end)
# raster_template <- cwd_cmip_end %>% as.data.frame(xy = TRUE) %>% select(x,y)
# test_clim <- raster_template %>%
#   left_join(test_clim, by = c("x", "y"))
# test_clim <- rasterFromXYZ(test_clim, crs = crs_template)
# range <- range_sf %>% filter(sp_code == species)
# tmap_mode("view")
# 
# tm_shape(test_clim$cwd_cmip_end1) +
#   tm_raster(palette = "-RdYlGn") +
#   tm_facets(as.layers = TRUE) +
#   tm_shape(range) + 
#   tm_fill(col = "lightblue")
# 
# tm_shape(test_clim$cwd_cmip_end1) +
#   tm_raster(palette = "-RdYlGn") +
#   tm_facets(as.layers = TRUE) +
#   tm_shape(test_clim$cwd_cmip_start1) +
#   tm_raster(palette = "-RdYlGn") +
#   tm_facets(as.layers = TRUE) +
#   tm_shape(range) + 
#   tm_fill(col = "lightblue")


## Export predictions
write_rds(sp_cmip_clim, paste0(wdir, "out/climate/sp_clim_predictions.", compress = "gz"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Apply species standardization to site-level data -----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# niche_df <- niche_df %>% 
#   select(sp_code = sp_code, sp_pet_mean = pet_mean, sp_pet_sd = pet_sd, sp_cwd_mean = cwd_mean, sp_cwd_sd = cwd_sd)


# Calculate site-level annual climate
site_clim_df = site_clim_df %>%
  group_by(collection_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum((aet+cwd)),
            .groups = "drop")


### Calculate site-level, average, historic, relative climate (for second stage)
ave_site_clim_df <- site_clim_df %>% 
  filter(year < 1980) %>% 
  group_by(collection_id) %>% 
  summarise(cwd.ave = mean(cwd.an, na.rm = TRUE), ## TODO: Shouldn't have to include na.rm here. Why do some sites have incomplete climate data?
            pet.ave = mean(pet.an, na.rm = TRUE)) %>% 
  drop_na() %>% 
  ungroup()

ave_site_clim_df <- ave_site_clim_df %>% 
  left_join(site_smry, by = "collection_id") %>% 
  group_by(sp_code) %>% 
  nest() %>% 
  left_join(niche_df, by = ("sp_code")) %>% 
  drop_na()

ave_site_clim_df <- ave_site_clim_df %>% 
  mutate(site_clim = future_pmap(list(hist_clim_vals = data,
                                      pet_mean = pet_mean,
                                      pet_sd = pet_sd,
                                      cwd_mean = cwd_mean,
                                      cwd_sd = cwd_sd),
                                 .f = sp_std_historic_df,
                                 .options = furrr_options(packages = c( "dplyr"))))

ave_site_clim_df <- ave_site_clim_df %>% 
  unnest(site_clim) %>% 
  rename(cwd.spstd = cwd.ave, pet.spstd = pet.ave) %>% 
  ungroup() %>% 
  select(-pet_mean, -pet_sd, -cwd_mean, -cwd_sd, -data, -sp_code)

write_rds(ave_site_clim_df, 
          paste0(wdir, "out/climate/site_ave_clim.", compress = "gz"))



### Calculate site-level, annual, historic, relative climate (for first stage) 
an_site_clim_df <- site_clim_df %>% 
  left_join(site_smry, by = "collection_id") %>% 
  group_by(sp_code) %>% 
  nest() %>% 
  left_join(niche_df, by = "sp_code") %>% 
  drop_na()

an_site_clim_df <- an_site_clim_df %>% 
  mutate(site_clim = future_pmap(list(hist_clim_vals = data,
                                             pet_mean = pet_mean,
                                             pet_sd = pet_sd,
                                             cwd_mean = cwd_mean,
                                             cwd_sd = cwd_sd),
                                        .f = sp_std_historic_df,
                                        .options = furrr_options(packages = c( "dplyr"))))

an_site_clim_df <- an_site_clim_df %>% 
  unnest(site_clim) %>% 
  rename(cwd.an.spstd = cwd.an, pet.an.spstd = pet.an) %>% 
  ungroup() %>% 
  select(-aet.an, -pet_mean, -pet_sd, -cwd_mean, -cwd_sd, -data, -sp_code)
 
write_rds(an_site_clim_df, 
          paste0(wdir, "out/climate/site_an_clim.", compress = "gz"))


# ## Exploring source of dropped sites - seems to be entirely driven by sites for species with no range maps
# an_site_clim_df %>% pull(collection_id) %>% unique() %>% length()
# site_clim_df %>% pull(collection_id) %>% unique() %>% length()
# # clim_sites <- clim_df %>% pull(collection_id) %>% unique()
# test_sites <- test %>% pull(collection_id) %>% unique()
# an_site_clim_df %>% unnest(cols = c(data)) %>% pull(collection_id) %>% unique() %>% length()
# an_site_clim_df %>% unnest(cols = c(data)) %>% drop_na() %>% pull(collection_id) %>% unique() %>% length()
