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
library(broom)
library(fixest)
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
ss_df <- read_rds(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))

# mod_df <- trim_df
# cwd_mod <- readRDS(paste0(wdir, "out/second_stage/cwd_mod.rds"))
# cwd_vcov <- readRDS(paste0(wdir, "out/second_stage/cwd_mod_vcov.rds"))
# pet_mod <- readRDS(paste0(wdir, "out/second_stage/pet_mod.rds"))
# int_mod <- readRDS(paste0(wdir, "out/second_stage/int_mod.rds"))

constant_sensitivities <- 
  read_rds(paste0(wdir, "out/first_stage/constant_sensitivities.rds"))


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

## Number of species
n_species <- flm_df %>% 
  pull(species_id) %>% 
  unique() %>% 
  length() %>% 
  print()

## Number of sites
n_sites <- flm_df %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length() %>% 
  print()

## Number of trees
n_trees <- flm_df %>% 
  pull(ntrees) %>% 
  sum() %>% 
  print()

## Number of tree-year observations
n_obs <- flm_df %>% 
  pull(nobs) %>% 
  sum() %>% 
  print()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Impacts of water and energy on tree growth  ----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## In the median plot, a one standard deviation increase in CWD was associated with an XX percent decline in annual growth"
# cwd_te_median <- flm_df %>%
#   pull(estimate_cwd.an) %>% 
#   median() %>% 
#   print()
constant_sensitivities$cwd$estimate[1]
constant_sensitivities$cwd$ci_dif[1]



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
# pet_te_median <- flm_df %>%
#   pull(estimate_pet.an) %>%
#   median() %>%
#   print()
constant_sensitivities$pet$estimate[1]
constant_sensitivities$pet$ci_dif[1]



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
# For example, plots sites with a historic, mean CWD exceeding their species’ 
# average historic CWD value experienced across a species’ range experienced an
# 3.3+/-3.7XX percent percent decline in growth in years with in response to a
# one standard deviation increase in CWD (Figure 3, Panel A). 
mod_df <- trim_df %>% 
  filter(cwd.spstd > 0)
cwd_high_mod <- feols(estimate_cwd.an ~ 1,
                       weights = mod_df$cwd_errorweights, data = mod_df,
                       vcov = conley(cutoff = 363, distance = "spherical")) %>% 
  tidy(conf.int = TRUE, conf.level = 0.95) %>% 
  mutate(ci_dif = estimate - conf.low) %>% 
  print()

# In contrast, sites with a lower than average historic CWD (i.e., relatively cooler 
# and wetter regions of a species’ niche) experienced an 8.9+/-3.9 percent decline 
# in growth from a one standard deviation increase in CWD.
mod_df <- trim_df %>% 
  filter(cwd.spstd < 0)
cwd_low_mod <- feols(estimate_cwd.an ~ 1,
                       weights = mod_df$cwd_errorweights, data = mod_df,
                       vcov = conley(cutoff = 363, distance = "spherical")) %>% 
  tidy(conf.int = TRUE, conf.level = 0.95) %>% 
  mutate(ci_dif = estimate - conf.low) %>% 
  print()

# After controlling for the correlation between historic CWD and PET, we estimate that 
# the slope of the line describing the relationship between historic CWD and sensitivity 
# to CWD is XX (Figure 2B). 
ss_df %>%
  summarise(mean_cwd = mean(cwd_cwd),
            lower_cwd = quantile(cwd_cwd, 0.025),
            upper_cwd = quantile(cwd_cwd, 0.975))


# In other words, while a one standard deviation shock to CWD 
# is anticipated to lead to an XX% decline in growth in sites located at a species’ mean 
# CWD, a comparable shock to CWD would lead to only an XX% decline in growth for sites 
# that are historically one standard deviation dryer than the species mean.
pull_marg_fx <- function(at_pet, at_cwd, mod_df){
  cwd_me_predictions <- mod_df$cwd_int + (at_cwd * mod_df$cwd_cwd) + (at_pet * mod_df$cwd_pet)
  cwd_ci_min <- cwd_me_predictions %>% quantile(0.025)
  cwd_ci_max <- cwd_me_predictions %>% quantile(0.975)
  cwd_mean <- cwd_me_predictions %>% mean()
  
  pet_me_predictions <- mod_df$pet_int + (at_cwd * mod_df$pet_cwd) + (at_pet * mod_df$pet_pet)
  pet_ci_min <- pet_me_predictions %>% quantile(0.025)
  pet_ci_max <- pet_me_predictions %>% quantile(0.975)
  pet_mean <- pet_me_predictions %>% mean()
  return(tibble(cwd_mean = cwd_mean, cwd_ci_min = cwd_ci_min, cwd_ci_max = cwd_ci_max,
                pet_mean = pet_mean, pet_ci_min = pet_ci_min, pet_ci_max = pet_ci_max))
}


at_pet <- 0
init_pull_marg_fx = partial(.f = pull_marg_fx, mod_df = ss_df)
cwd_me_df <- tibble(at_cwd = c(-2,-1,0, 1, 2))
cwd_me_df <- cwd_me_df %>%
  mutate(cwd_me = pmap(list(at_pet = at_pet,
                            at_cwd = cwd_me_df$at_cwd),
                       .f = init_pull_marg_fx)) %>% 
  unnest(cwd_me) %>% 
  mutate(cwd_ci_dif = cwd_mean - cwd_ci_min,
         cwd_ci_dif2 = cwd_ci_max - cwd_mean) %>% 
  print()

# For example, 
mod_df <- trim_df %>% 
  filter(pet.spstd < 0)
pet_low_mod <- feols(estimate_pet.an ~ 1,
                      weights = mod_df$pet_errorweights, data = mod_df,
                      vcov = conley(cutoff = 363, distance = "spherical")) %>% 
  tidy(conf.int = TRUE, conf.level = 0.95) %>% 
  mutate(ci_dif = estimate - conf.low) %>% 
  print()

# For example, 
mod_df <- trim_df %>% 
  filter(pet.spstd > 0)
pet_low_mod <- feols(estimate_pet.an ~ 1,
                     weights = mod_df$pet_errorweights, data = mod_df,
                     vcov = conley(cutoff = 363, distance = "spherical")) %>% 
  tidy(conf.int = TRUE, conf.level = 0.95) %>% 
  mutate(ci_dif = estimate - conf.low) %>% 
  print()


# at_cwd <- 0
# init_pull_marg_fx = partial(.f = pull_marg_fx, mod_df = ss_df)
# cwd_me_df <- tibble(at_pet = c(-2,-1, 0, 1, 2))
# cwd_me_df <- cwd_me_df %>%
#   mutate(cwd_me = pmap(list(at_pet = cwd_me_df$at_pet,
#                             at_cwd = at_cwd),
#                        .f = init_pull_marg_fx)) %>% 
#   unnest(cwd_me) %>% 
#   mutate(cwd_ci_dif = cwd_mean - cwd_ci_min,
#          cwd_ci_dif2 = cwd_ci_max - cwd_mean) %>% 
#   print()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Forecasted changes in CWD and PET across species’ climatic niches  --------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gen_dat <- read_csv(paste0(wdir, "out/species_gen_gr.csv")) %>% 
  rename(sp_code = "species_id")

sp_predictions_gen <- sp_predictions %>% 
  left_join(gen_dat)

sp_plot_dat <- sp_predictions_gen %>% 
  group_by(sp_code) %>% 
  mutate(beyond_max_cwd = cwd_cmip_end_mean > max(cwd_hist)) %>% 
  ungroup() %>% 
  mutate(cwd_end_bin = case_when(
    (beyond_max_cwd == TRUE) ~ "Beyond prior max",
    (cwd_cmip_end_mean > 1.5) ~ "1-2 s.d. above prior mean",
    (cwd_cmip_end_mean >= 0) ~ "0-1 s.d. above prior mean",
    (cwd_cmip_end_mean < 0) ~ "Below prior mean"),
    cwd_end_bin = fct_relevel(cwd_end_bin, 
                              "Beyond prior max", "1-2 s.d. above prior mean", 
                              "0-1 s.d. above prior mean", 
                              "Below prior mean"))
sp_plot_dat %>%
  ungroup() %>% 
  group_by(cwd_end_bin) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))


sp_plot_dat %>% 
  mutate(greater_cwd = cwd_cmip_end_mean >= cwd_cmip_start_mean,
         greater_pet = pet_cmip_end_mean >= pet_cmip_start_mean,
         greater_both = greater_cwd * greater_pet) %>% 
  pull(greater_both) %>% 
  mean(na.rm = TRUE)
