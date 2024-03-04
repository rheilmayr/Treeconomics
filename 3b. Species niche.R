#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: 1) Characterize climate niche for different species. 
#          2) Standardize annual site data, baseline raster, and CMIP predictions for each species
#
# Input files:
#   merged_ranges.shp:
#   HistoricCWD_AETGrids_Annual.Rdat: Rasters describing historic CWD and AET
#     Generated using historic_cwdraster.R
#   monthlycrubaseline_tas:
#   cmip5_cwdaet_start.Rdat:
#   cmip5_cwdaet_end.Rdat:
#   essentialcwd_data.csv:
#   site_summary.csv:
# 
# Output files:
#   clim_niche.csv: Tabulation of each species' historic climate niche. Parameters used for standardization.
#   sp_clim_predictions.gz: Dataset describing species' standardized historic climate and their predicted climate under CMIP5 across species' full range
#   site_ave_clim.gz: Tables describing each site's species-standardized average historic climate
#   site_an_clim.gz: Tables describing each site's annual species-standardized weather
# 
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
library(terra)
library(readr)
library(tmap)
library(tictoc)
select <- dplyr::select


library(furrr)
n_cores <- 8
future::plan(multisession, workers = n_cores)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Historic climate raster
clim_file <- paste0(wdir, '1_input_processed/climate/HistoricCWD_AETGrids_Annual.Rdat')
load(clim_file)
# cwd_historic <- rast(cwd_historic)
# aet_historic <- rast(aet_historic)

pet_historic <- mean(pet_historic %>% subset(2:80))
cwd_historic <- mean(cwd_historic %>% subset(2:80))
ppt_historic <- mean(ppt_historic %>% subset(2:80))
tmp_historic <- mean(tmp_historic %>% subset(2:80))
petsp_historic <- mean(petsp_historic %>% subset(2:80))
cwdsp_historic <- mean(cwdsp_historic %>% subset(2:80))


# pet_historic <- mean(pet_historic %>% subset(61:80))
# cwd_historic <- mean(cwd_historic %>% subset(61:80))
# ppt_historic <- mean(ppt_historic %>% subset(61:80))
# tmp_historic <- mean(tmp_historic %>% subset(61:80))

names(cwd_historic) = "cwd"
names(pet_historic) = "pet"
names(ppt_historic) = "ppt"
names(tmp_historic) = "tmp"
names(cwdsp_historic) = "cwd_spei"
names(petsp_historic) = "pet_spei"


# Composite into multilayer spatraster
cwd_historic <- rast(cwd_historic)
pet_historic <- rast(pet_historic)
tmp_historic <- rast(tmp_historic)
ppt_historic <- rast(ppt_historic)
cwdsp_historic <- rast(cwdsp_historic)
petsp_historic <- rast(petsp_historic)
clim_historic <- rast(list(cwd_historic, pet_historic, ppt_historic, tmp_historic,cwdsp_historic,petsp_historic))


# # 2. Add terraclimate data for comparison
# cwd_tc <- rast(paste0(wdir,"0_raw/TerraClimate/TerraClimate19611990_def.nc")) %>% 
#   sum()
# pet_tc <- rast(paste0(wdir,"0_raw/TerraClimate/TerraClimate19611990_pet.nc")) %>% 
#   sum()
# ppt_tc <- rast(paste0(wdir,"0_raw/TerraClimate/TerraClimate19611990_ppt.nc")) %>% 
#   sum()
# clim_tc <- rast(list("cwd" = cwd_tc, "pet" = pet_tc, "ppt" = ppt_tc))


# 3. Site-specific historic climate data
site_clim_csv <- paste0(wdir, '1_input_processed/climate/essentialcwd_data.csv')
site_clim_df <- read_csv(site_clim_csv)
site_clim_df <- site_clim_df %>% 
  mutate(site = as.character(site)) %>% 
  rename(collection_id = site)

tc_pet <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_pet.csv"))
tc_cwd <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_def.csv"))
tc_ppt <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_ppt.csv"))
site_clim_df_tc <- tc_pet %>%
  left_join(tc_cwd, by = c("collection_id", "Month", "year")) %>%
  left_join(tc_ppt, by = c("collection_id", "Month", "year")) %>%
  rename(month = Month,
         tc_pet = pet,
         tc_cwd = def,
         tc_ppt = ppt)

site_clim_df <- site_clim_df %>% 
  left_join(site_clim_df_tc, by = c("collection_id", "year", "month"))

# Quick check of agreement
summary(feols(pet ~ tc_pet | collection_id, data = site_clim_df))


# 4. Load species information for sites
site_smry <- read_csv(paste0(wdir, '1_input_processed/dendro/site_summary.csv'))
site_smry <- site_smry %>%
  select(collection_id, sp_id, latitude) %>% 
  mutate(sp_code = tolower(sp_id)) %>% 
  select(-sp_id)


# # NOTE: FIA data not included in replication data repository
# site_smry_fia <- read_csv(paste0(wdir, 'out/dendro/site_summary_fia.csv'))
# site_smry_fia <- site_smry_fia %>%
#   select(collection_id, location_id = plot_cn, sp_id = species_id) %>%
#   mutate(sp_code = tolower(sp_id)) %>%
#   select(-sp_id)
# site_smry <- rbind(site_smry, site_smry_fia)


# 5. Species range maps
range_file <- paste0(wdir, '1_input_processed/species_ranges/merged_ranges_dissolve.shp')
range_sf <- st_read(range_file)

# 6. Climate projections from CMIP5
cmip_end <- load(paste0(wdir, '1_input_processed/climate/cmip5_cwdaet_end.Rdat'))
pet_cmip_end <- aet_raster + cwd_raster
cwd_cmip_end <- cwd_raster
names(cwd_cmip_end) <- NULL # Resetting this due to strange names in file from CMIP processing
rm(cwd_raster)
rm(aet_raster)

cmip_start <- load(paste0(wdir, '1_input_processed/climate/cmip5_cwdaet_start.Rdat'))
pet_cmip_start <- aet_raster + cwd_raster
cwd_cmip_start <- cwd_raster
names(cwd_cmip_start) <- NULL # Resetting this due to strange names in file from CMIP processing
rm(cwd_raster)
rm(aet_raster)




# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Visually inspect data -----------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# tmap_mode("view")
# tm_shape(pet_historic) +
#   tm_raster() +
#   tm_facets(as.layers = TRUE)
# 
# 
# tm_shape(cwd_cmip_end) +
#   tm_raster() +
#   tm_facets(as.layers = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species niches -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull and organize climate distribution for species
pull_clim <- function(spp_code, clim_raster){
  print(spp_code)
  
  # Pull relevant range map
  sp_range <- range_sf %>%
    filter(sp_code == spp_code)
  
  # Pull clim values
  clim_vals <- clim_raster %>% 
    mask(mask = sp_range, touches = TRUE) %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  
  return(clim_vals)
}

pull_clim_base <- partial(.f = pull_clim, clim_raster = clim_historic)

species_list <- range_sf %>%
  pull(sp_code) %>% 
  unique() %>% 
  enframe(name = NULL) %>% 
  select(sp_code = value) %>% 
  arrange(sp_code) %>% 
  drop_na()

clim_df <- species_list %>% 
  mutate(clim_vals = map(sp_code,.f = pull_clim_base))

## Summarize mean and sd of each species' climate
niche_df <- clim_df %>% 
  unnest(clim_vals) %>% 
  group_by(sp_code) %>% 
  summarize(pet_mean = mean(pet),
            pet_sd = sd(pet),
            cwd_mean = mean(cwd),
            cwd_sd = sd(cwd),
            temp_mean = mean(tmp),
            temp_sd = sd(tmp),
            ppt_mean = mean(ppt), 
            ppt_sd = sd(ppt),
            spei_pet_mean = mean(pet_spei),
            spei_pet_sd = sd(pet_spei),
            spei_cwd_mean = mean(cwd_spei),
            spei_cwd_sd = sd(cwd_spei))


## Export species niche description
write_csv(niche_df, paste0(wdir, "2_output/climate/clim_niche.csv"))
niche_df <- read_csv(paste0(wdir, "2_output/climate/clim_niche.csv"))


# ## Repeat for terraclimate data
# pull_clim_tc <- partial(.f = pull_clim, clim_raster = clim_tc)
# 
# clim_df_tc <- species_list %>% 
#   mutate(clim_vals = map(sp_code,.f = pull_clim_tc))
# niche_df_tc <- clim_df_tc %>% 
#   unnest(clim_vals) %>% 
#   group_by(sp_code) %>% 
#   summarize(pet_mean = mean(pet),
#             pet_sd = sd(pet),
#             cwd_mean = mean(cwd),
#             cwd_sd = sd(cwd),
#             ppt_mean = mean(ppt), 
#             ppt_sd = sd(ppt))
# write_csv(niche_df_tc, paste0(wdir, "2_output/climate/clim_niche_tc.csv"))
niche_df_tc <- read_csv(paste0(wdir, "2_output/climate/clim_niche_tc.csv"))
niche_df_tc <- niche_df_tc %>% 
  rename(tc_pet_mean = pet_mean,
         tc_pet_sd = pet_sd,
         tc_cwd_mean = cwd_mean,
         tc_cwd_sd = cwd_sd,
         tc_ppt_mean = ppt_mean,
         tc_ppt_sd = ppt_sd)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Standardize historic climate -------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_standardize <- function(val, sp_mean, sp_sd){
  std_val <- (val - sp_mean) / sp_sd
  return(std_val)
}

sp_std_historic_df <- function(hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd, 
                               temp_mean, temp_sd, ppt_mean, ppt_sd, 
                               tc_pet_mean, tc_pet_sd, tc_cwd_mean, tc_cwd_sd,
                               tc_ppt_mean, tc_ppt_sd, spei_pet_mean, spei_pet_sd,
                               spei_cwd_mean, spei_cwd_sd){
  hist_clim_vals <- hist_clim_vals %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd)) %>% 
    mutate_at(vars(starts_with("temp")), 
              ~sp_standardize(.x, temp_mean, temp_sd)) %>% 
    mutate_at(vars(starts_with("ppt")), 
              ~sp_standardize(.x, ppt_mean, ppt_sd)) %>% 
    mutate_at(vars(starts_with("tc_cwd")), 
              ~sp_standardize(.x, tc_cwd_mean, tc_cwd_sd)) %>% 
    mutate_at(vars(starts_with("tc_pet")), 
              ~sp_standardize(.x, tc_pet_mean, tc_pet_sd)) %>% 
    mutate_at(vars(starts_with("tc_ppt")), 
              ~sp_standardize(.x, tc_ppt_mean, tc_ppt_sd)) %>% 
    mutate_at(vars(starts_with("spei_cwd")), 
              ~sp_standardize(.x, spei_cwd_mean, spei_cwd_sd)) %>% 
    mutate_at(vars(starts_with("spei_pet")), 
              ~sp_standardize(.x, spei_pet_mean, spei_pet_sd))
  return(hist_clim_vals)
}


sp_std_future_df <- function(cmip_df, hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd, temp_mean, temp_sd, ppt_mean, ppt_sd){
  valid_locations <- hist_clim_vals %>% select(x,y)
  cmip_df <- valid_locations %>% 
    left_join(cmip_df, by = c("x", "y"))
  cmip_df <- cmip_df %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd)) %>% 
    mutate_at(vars(starts_with("temp")), 
              ~sp_standardize(.x, temp_mean, temp_sd)) %>% 
    mutate_at(vars(starts_with("ppt")), 
              ~sp_standardize(.x, ppt_mean, ppt_sd))
  return(cmip_df)
}


# clim_df <- clim_df %>% 
#   left_join(niche_df, by = "sp_code")
# 
# clim_df <- clim_df %>% 
#   mutate(clim_historic_sp = future_pmap(list(hist_clim_vals = clim_vals,
#                                              pet_mean = pet_mean,
#                                              pet_sd = pet_sd,
#                                              cwd_mean = cwd_mean,
#                                              cwd_sd = cwd_sd,
#                                              temp_mean = temp_mean,
#                                              temp_sd = temp_sd),
#                                         .f = sp_std_historic_df,
#                                         .options = furrr_options(packages = c( "dplyr"))))
# NOTE: May no longer need this dataframe???


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Apply species standardization to site-level data -----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# niche_df <- niche_df %>% 
#   select(sp_code = sp_code, sp_pet_mean = pet_mean, sp_pet_sd = pet_sd, sp_cwd_mean = cwd_mean, sp_cwd_sd = cwd_sd)


# Shift to water year
site_clim_df <- site_clim_df %>% 
  left_join(site_smry, by = "collection_id") %>% 
  as.data.table()

site_clim_df[,water_year:=year]
site_clim_df[(latitude>=0) & (month>=10),water_year:=year+1] # Northern hemisphere water year is october through september
site_clim_df[(latitude<0) & (month>=7),water_year:=year+1] # Southern hemisphere water year is July through June
site_clim_df <- site_clim_df %>% 
  as_tibble() %>% 
  select(-year) %>% 
  rename(year = water_year)


# Calculate site-level annual climate
site_clim_df = site_clim_df %>%
  group_by(collection_id, year) %>%
  summarise(cwd.an = sum(cwd),
            pet.an = sum(pet),
            ppt.an = sum(ppt),
            temp.an = mean(tmean),
            tc_cwd.an = sum(tc_cwd),
            tc_pet.an = sum(tc_pet),
            tc_ppt.an = sum(tc_ppt),
            spei_pet.an = sum(pet_spei),
            spei_cwd.an = sum(cwd_spei),
            .groups = "drop")


### Calculate site-level, average, historic, relative climate (for second stage)
## TODO: Note - dropping CANA323 because it has null climate data for a few months each year. might want to dig into this
site_clim_df <- site_clim_df %>% 
  filter(collection_id != "CANA323")

ave_site_clim_df <- site_clim_df %>% 
  filter(year>1901, year < 1980) %>% 
  group_by(collection_id) %>% 
  summarise(cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an),
            temp.ave = mean(temp.an),
            ppt.ave = mean(ppt.an),
            tc_cwd.ave = mean(tc_cwd.an, na.rm = TRUE),
            tc_pet.ave = mean(tc_pet.an, na.rm = TRUE),
            tc_ppt.ave = mean(tc_ppt.an, na.rm = TRUE),
            spei_cwd.ave = mean(spei_cwd.an),
            spei_pet.ave = mean(spei_pet.an)) %>% 
  ungroup()

spstd_site_clim_df <- site_smry %>% 
  left_join(ave_site_clim_df, by = "collection_id") %>% 
  group_by(sp_code) %>% 
  nest(data = c(collection_id, cwd.ave, pet.ave, temp.ave, ppt.ave, tc_cwd.ave, tc_pet.ave, tc_ppt.ave, spei_cwd.ave, spei_pet.ave)) %>% 
  left_join(niche_df, by = ("sp_code")) %>%
  left_join(niche_df_tc, by = "sp_code") %>% 
  drop_na() # Dropping some species due to NA niche data

spstd_site_clim_df <- spstd_site_clim_df %>% 
  mutate(site_clim = future_pmap(list(hist_clim_vals = data,
                                      pet_mean = pet_mean,
                                      pet_sd = pet_sd,
                                      cwd_mean = cwd_mean,
                                      cwd_sd = cwd_sd,
                                      temp_mean = temp_mean,
                                      temp_sd = temp_sd,
                                      ppt_mean = ppt_mean,
                                      ppt_sd = ppt_sd,
                                      tc_pet_mean = tc_pet_mean,
                                      tc_pet_sd = tc_pet_sd,
                                      tc_cwd_mean = tc_cwd_mean,
                                      tc_cwd_sd = tc_cwd_sd,
                                      tc_ppt_mean = tc_ppt_mean,
                                      tc_ppt_sd = tc_ppt_sd,
                                      spei_pet_mean = spei_pet_mean,
                                      spei_pet_sd = spei_pet_sd,
                                      spei_cwd_mean = spei_cwd_mean,
                                      spei_cwd_sd = spei_cwd_sd),
                                 .f = sp_std_historic_df,
                                 .options = furrr_options(packages = c( "dplyr"))))

spstd_site_clim_df <- spstd_site_clim_df %>% 
  unnest(site_clim) %>% 
  rename(cwd.spstd = cwd.ave, pet.spstd = pet.ave, temp.spstd = temp.ave, ppt.spstd = ppt.ave, 
         cwd.spstd.tc = tc_cwd.ave, pet.spstd.tc = tc_pet.ave, ppt.spstd.tc = tc_ppt.ave,
         pet.spstd.spei = spei_pet.ave, cwd.spstd.spei = spei_cwd.ave) %>% 
  ungroup() %>% 
  select(collection_id, cwd.spstd, pet.spstd, temp.spstd, ppt.spstd, 
         cwd.spstd.tc, pet.spstd.tc, ppt.spstd.tc, cwd.spstd.spei, pet.spstd.spei)

spstd_site_clim_df <- spstd_site_clim_df %>% 
  left_join(ave_site_clim_df %>% select(collection_id, cwd.ave, pet.ave, temp.ave, ppt.ave, tc_cwd.ave, tc_pet.ave, tc_ppt.ave, spei_pet.ave, spei_cwd.ave), by = "collection_id") %>% 
  rename(cwd.ave.tc = tc_cwd.ave, pet.ave.tc = tc_pet.ave, ppt.ave.tc = tc_ppt.ave, cwd.ave.spei = spei_cwd.ave, pet.ave.spei = spei_pet.ave)

write_rds(spstd_site_clim_df, 
          paste0(wdir, "2_output/climate/site_ave_clim.", compress = "gz"))





### Calculate site-level, annual, historic, relative climate (for first stage) 
an_site_clim_df <- site_smry %>% 
  left_join(site_clim_df, by = "collection_id") %>% 
  group_by(sp_code) %>% 
  nest() %>% 
  left_join(niche_df, by = "sp_code") %>% 
  left_join(niche_df_tc, by = "sp_code") %>% 
  drop_na()

an_site_clim_df <- an_site_clim_df %>% 
  mutate(site_clim = future_pmap(list(hist_clim_vals = data,
                                      pet_mean = pet_mean,
                                      pet_sd = pet_sd,
                                      cwd_mean = cwd_mean,
                                      cwd_sd = cwd_sd,
                                      temp_mean = temp_mean,
                                      temp_sd = temp_sd,
                                      ppt_mean = ppt_mean,
                                      ppt_sd = ppt_sd,
                                      tc_pet_mean = tc_pet_mean,
                                      tc_pet_sd = tc_pet_sd,
                                      tc_cwd_mean = tc_cwd_mean,
                                      tc_cwd_sd = tc_cwd_sd,
                                      tc_ppt_mean = tc_ppt_mean,
                                      tc_ppt_sd = tc_ppt_sd,
                                      spei_cwd_mean = spei_cwd_mean,
                                      spei_cwd_sd = spei_cwd_sd,
                                      spei_pet_mean = spei_pet_mean,
                                      spei_pet_sd = spei_pet_sd),
                                 .f = sp_std_historic_df,
                                 .options = furrr_options(packages = c( "dplyr"))))

an_site_clim_df <- an_site_clim_df %>% 
  ungroup() %>% 
  select(site_clim) %>% 
  unnest(site_clim) %>% 
  rename(cwd.an.spstd = cwd.an, pet.an.spstd = pet.an, temp.an.spstd = temp.an, 
         ppt.an.spstd = ppt.an, cwd.an.spstd.tc = tc_cwd.an, pet.an.spstd.tc = tc_pet.an,
         ppt.an.spstd.tc = tc_ppt.an, cwd.an.spstd.spei = spei_cwd.an, pet.an.spstd.spei = spei_pet.an) %>% 
  ungroup()


write_rds(an_site_clim_df, 
          paste0(wdir, "2_output/climate/site_an_clim.", compress = "gz"))

# Quick check of correlation across datasets
library(fixest)
summary(lm(cwd.spstd ~ cwd.spstd.tc, data = spstd_site_clim_df))
summary(lm(ppt.spstd ~ ppt.spstd.tc, data = spstd_site_clim_df))

summary(lm(cwd.an.spstd ~ cwd.an.spstd.tc, data = an_site_clim_df))
summary(feols(cwd.an.spstd ~ cwd.an.spstd.tc | collection_id, data = an_site_clim_df))
summary(lm(pet.an.spstd ~ pet.an.spstd.tc, data = an_site_clim_df))
summary(feols(pet.an.spstd ~ pet.an.spstd.tc | collection_id, data = an_site_clim_df))
summary(lm(ppt.an.spstd ~ ppt.an.spstd.tc, data = an_site_clim_df))


# ## Exploring source of dropped sites - seems to be entirely driven by sites for species with no range maps
# an_site_clim_df %>% pull(collection_id) %>% unique() %>% length()
# site_clim_df %>% pull(collection_id) %>% unique() %>% length()
# # clim_sites <- clim_df %>% pull(collection_id) %>% unique()
# test_sites <- test %>% pull(collection_id) %>% unique()
# an_site_clim_df %>% unnest(cols = c(data)) %>% pull(collection_id) %>% unique() %>% length()
# an_site_clim_df %>% unnest(cols = c(data)) %>% drop_na() %>% pull(collection_id) %>% unique() %>% length()



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
write_rds(sp_cmip_clim, paste0(wdir, "2_output/climate/sp_clim_predictions.", compress = "gz"))


