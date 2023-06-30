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
rwi_list <- list.files(paste0(wdir, "out/predictions/pred_10000/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))
# sp_predictions <- readRDS(paste0(wdir, "out/predictions/sp_predictions.rds"))

hot_cell_list <- list.files(paste0(wdir, "out/predictions/pred_10000/sp_hot_cells/"), pattern = ".gz", full.names = TRUE)
hot_cells <- do.call('rbind', lapply(hot_cell_list, readRDS))



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


# 8. Change in CWD
cmip_end <- load(paste0(wdir, 'in\\CMIP5 CWD\\cmip5_cwdaet_end.Rdat'))
cwd_cmip_end <- cwd_raster %>% mean()
pet_cmip_end <- (aet_raster + cwd_raster) %>% mean()
rm(cwd_raster)
rm(aet_raster)

cmip_start <- load(paste0(wdir, 'in\\CMIP5 CWD\\cmip5_cwdaet_start.Rdat'))
cwd_cmip_start <- cwd_raster %>% mean()
pet_cmip_end <- (aet_raster + cwd_raster) %>% mean()
rm(cwd_raster)
rm(aet_raster)

cwd_cmip_change <- (cwd_cmip_end - cwd_cmip_start) %>% as.data.frame() %>% drop_na()
pet_cmip_change <- (pet_cmip_end - cwd_cmip_start) %>% as.data.frame() %>% drop_na()



agg_stats <- read_rds(file = paste0(wdir, "out/predictions/sp_rwi_pred_10000/mc_agg_stats.gz"))




block_draw_df <- read_rds(paste0(wdir, "out/second_stage/mc_sample.gz"))

gen_dat <- read_csv(paste0(wdir, "out/species_gen_gr.csv")) %>% 
  rename(sp_code = "species_id")


sp_clim <- read_rds(paste0(wdir, "out/climate/sp_clim_predictions.gz"))


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
# Effects of CWD and PET on tree growth  ----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Specifically, a site that experienced a 1 SD decline in CWD (relative to the distribution of 
# CWD across the species’ full range) experienced a 3.1% decline (-9.9 to -1.0)  in RWI in the 
# same year. The site-level, marginal effect of CWD was significantly negative (p<0.05) for 54% of 
# sites, while only 11% of sites exhibited a significantly positive relationship (Appendix Figure 
# 2, panel A). 
ss_df$cwd_const_sens %>% mean()
ss_df$cwd_const_sens %>%
  quantile(c(0.025, 0.5, 0.975))


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




# In contrast, the relationship between PET and growth varied greatly across sites, 
# with significant, positive relationships in 40% of sites, and significant, negative relationships 
# in 25% of sites (Appendix Figure 2, panel B). 
ss_df$pet_const_sens %>% mean()
ss_df$pet_const_sens %>%
  quantile(c(0.025, 0.5, 0.975))

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
# Drier sites exhibit less sensitivity to drought  --------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# For example, our estimates of site-level sensitivity to CWD indicated that 
# drier-than-average sites experienced a 13% decline (-22% to -7%) in annual growth 
# in response to a 1 SD increase in CWD (Figure 4, A). 
cwd_high_fs_bs <- block_draw_df %>% 
  filter(cwd.spstd > 0) %>% 
  group_by(boot_id) %>% 
  summarise(cwd_coef = mean(cwd_coef))

cwd_high_fs_bs %>%
  pull(cwd_coef) %>% 
  median()

cwd_high_fs_bs %>%
  pull(cwd_coef) %>% 
  quantile(c(0.025, 0.975))


# In contrast, wetter-than-average sites experienced a larger 21% decline (-32% to -11%) 
# in growth from a 1 SD increase in CWD. 
cwd_low_fs_bs <- block_draw_df %>% 
  filter(cwd.spstd < 0) %>% 
  group_by(boot_id) %>% 
  summarise(cwd_coef = mean(cwd_coef))

cwd_low_fs_bs %>%
  pull(cwd_coef) %>% 
  mean() 

cwd_low_fs_bs %>%
  pull(cwd_coef) %>% 
  quantile(c(0.025, 0.975))


# To characterize this variability in CWD sensitivity, we modeled sensitivity as a 
# quadratic function of historic CWD and PET (Figure 4, Panels B and C). The regression 
# results highlight that observed sensitivity declined with increasing CWD (slope at 
# the species’ mean, historic CWD = 0.015, CI = 0.006 to 0.064). 
ss_df %>%
  summarise(mean_cwd = mean(cwd_cwd),
            lower_cwd = quantile(cwd_cwd, 0.025),
            upper_cwd = quantile(cwd_cwd, 0.975))







# # pet_high_fs_bs <- block_draw_df %>% 
# #   filter(pet.spstd > 0) %>% 
# #   group_by(boot_id) %>% 
# #   summarise(pet_coef = mean(pet_coef))
# # 
# # pet_high_fs_bs %>%
# #   pull(pet_coef) %>% 
# #   mean()
# # 
# # pet_high_fs_bs %>%
# #   pull(pet_coef) %>% 
# #   quantile(c(0.025, 0.975))
# # 
# # pet_low_fs_bs <- block_draw_df %>% 
# #   filter(pet.spstd <= 0) %>% 
# #   group_by(boot_id) %>% 
# #   summarise(pet_coef = mean(pet_coef))
# # 
# # pet_low_fs_bs %>%
# #   pull(pet_coef) %>% 
# #   mean() 
# # 
# # pet_low_fs_bs %>%
# #   pull(pet_coef) %>% 
# #   quantile(c(0.025, 0.975))
# # 
# # 
# # 
# # 
# # 
# # mod_df <- trim_df %>% 
# #   filter(cwd.spstd > 0)
# # cwd_high_mod <- feols(estimate_cwd.an ~ 1,
# #                        weights = mod_df$cwd_errorweights, data = mod_df,
# #                        vcov = conley(cutoff = 363, distance = "spherical")) %>% 
# #   tidy(conf.int = TRUE, conf.level = 0.95) %>% 
# #   mutate(ci_dif = estimate - conf.low) %>% 
# #   print()
# 
# # # In contrast, sites with a lower than average historic CWD (i.e., relatively cooler 
# # # and wetter regions of a species’ niche) experienced an 8.5+/-4.2 percent decline 
# # # in growth from a one standard deviation increase in CWD.
# # mod_df <- trim_df %>% 
# #   filter(cwd.spstd < 0)
# # cwd_low_mod <- feols(estimate_cwd.an ~ 1,
# #                        weights = mod_df$cwd_errorweights, data = mod_df,
# #                        vcov = conley(cutoff = 363, distance = "spherical")) %>% 
# #   tidy(conf.int = TRUE, conf.level = 0.95) %>% 
# #   mutate(ci_dif = estimate - conf.low) %>% 
# #   print()
# 
# # After controlling for the correlation between historic CWD and PET, we estimate that 
# # the slope of the line describing the relationship between historic CWD and sensitivity 
# # to CWD is XX (Figure 2B). 
# 
# # In other words, while a one standard deviation shock to CWD 
# # is anticipated to lead to an XX% decline in growth in sites located at a species’ mean 
# # CWD, a comparable shock to CWD would lead to only an XX% decline in growth for sites 
# # that are historically one standard deviation dryer than the species mean.
# pull_marg_fx <- function(at_pet, at_cwd, mod_df){
#   cwd_me_predictions <- mod_df$cwd_int + (at_cwd * mod_df$cwd_cwd) + (at_cwd^2 * mod_df$cwd_cwd2) + (at_pet * mod_df$cwd_pet) + (at_pet^2 * mod_df$cwd_pet2)
#   cwd_ci_min <- cwd_me_predictions %>% quantile(0.025)
#   cwd_ci_max <- cwd_me_predictions %>% quantile(0.975)
#   cwd_mean <- cwd_me_predictions %>% mean()
#   
#   pet_me_predictions <- mod_df$pet_int + (at_cwd * mod_df$pet_cwd) + (at_cwd^2 * mod_df$pet_cwd2) + (at_pet * mod_df$pet_pet) + (at_pet^2 * mod_df$pet_pet2)
#   pet_ci_min <- pet_me_predictions %>% quantile(0.025)
#   pet_ci_max <- pet_me_predictions %>% quantile(0.975)
#   pet_mean <- pet_me_predictions %>% mean()
#   return(tibble(cwd_mean = cwd_mean, cwd_ci_min = cwd_ci_min, cwd_ci_max = cwd_ci_max,
#                 pet_mean = pet_mean, pet_ci_min = pet_ci_min, pet_ci_max = pet_ci_max))
# }
# 
# 
# at_pet <- 0
# init_pull_marg_fx = partial(.f = pull_marg_fx, mod_df = ss_df)
# cwd_me_df <- tibble(at_cwd = c(-2,-1,0, 1, 2))
# cwd_me_df <- cwd_me_df %>%
#   mutate(cwd_me = pmap(list(at_pet = at_pet,
#                             at_cwd = cwd_me_df$at_cwd),
#                        .f = init_pull_marg_fx)) %>% 
#   unnest(cwd_me) %>% 
#   mutate(cwd_ci_dif = cwd_mean - cwd_ci_min,
#          cwd_ci_dif2 = cwd_ci_max - cwd_mean) %>% 
#   print()
# 
# 
# 
# # For example, ITRDB sites located below their species’ mean, historic PET 
# # experienced an average 7.2+/-4.9 percent increase in annual growth in response 
# # to a single-year, one standard deviation increase in PET.
# mod_df <- trim_df %>% 
#   filter(pet.spstd < 0)
# pet_low_mod <- feols(estimate_pet.an ~ 1,
#                       weights = mod_df$pet_errorweights, data = mod_df,
#                       vcov = conley(cutoff = 363, distance = "spherical")) %>% 
#   tidy(conf.int = TRUE, conf.level = 0.95) %>% 
#   mutate(ci_dif = estimate - conf.low) %>% 
#   print()
# 
# # In contrast, increases in PET in sites falling above their species’ historic mean PET 
# # were, on average, not associated with a significant change in growth.
# mod_df <- trim_df %>% 
#   filter(pet.spstd > 0)
# pet_low_mod <- feols(estimate_pet.an ~ 1,
#                      weights = mod_df$pet_errorweights, data = mod_df,
#                      vcov = conley(cutoff = 363, distance = "spherical")) %>% 
#   tidy(conf.int = TRUE, conf.level = 0.95) %>% 
#   mutate(ci_dif = estimate - conf.low) %>% 
#   print()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Climatic changes across species ranges  --------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CMIP5 model ensembles predict that the entirety of every species’ range 
# analyzed will experience increases in both mean PET and CWD 
sp_predictions_gen <- sp_predictions %>% 
  left_join(gen_dat)

sp_plot_dat <- sp_predictions_gen %>% 
  group_by(sp_code) %>% 
  mutate(beyond_max_cwd = cwd_cmip_end_mean > max(cwd_hist),
         beyond_max_pet = pet_cmip_end_mean > max(pet_hist)) %>% 
  ungroup() %>% 
  mutate(cwd_end_bin = case_when(
    (beyond_max_cwd == TRUE) ~ "Beyond prior max",
    (cwd_cmip_end_mean > 1.5) ~ "1-2 s.d. above prior mean",
    (cwd_cmip_end_mean >= 0) ~ "0-1 s.d. above prior mean",
    (cwd_cmip_end_mean < 0) ~ "Below prior mean"),
    cwd_end_bin = fct_relevel(cwd_end_bin, 
                              "Beyond prior max", "1-2 s.d. above prior mean", 
                              "0-1 s.d. above prior mean", 
                              "Below prior mean"),
    pet_end_bin = case_when(
      (beyond_max_pet == TRUE) ~ "Beyond prior max",
      (pet_cmip_end_mean > 1.5) ~ "1-2 s.d. above prior mean",
      (pet_cmip_end_mean >= 0) ~ "0-1 s.d. above prior mean",
      (pet_cmip_end_mean < 0) ~ "Below prior mean"),
    pet_end_bin = fct_relevel(pet_end_bin, 
                              "Beyond prior max", "1-2 s.d. above prior mean", 
                              "0-1 s.d. above prior mean", 
                              "Below prior mean"))
sp_plot_dat %>% 
  mutate(greater_cwd = cwd_cmip_end_mean >= cwd_cmip_start_mean,
         greater_pet = pet_cmip_end_mean >= pet_cmip_start_mean,
         greater_both = greater_cwd * greater_pet) %>% 
  pull(greater_both) %>% 
  mean(na.rm = TRUE)


# On average, the mean PET across a species’ range is anticipated to increase 
# by 1.39 SD (0.45 – 3.04), while CWD is anticipated to increase by 1.41 SD 
# (0.34 – 3.77). 
single_sp_clim <- (sp_clim[1, 8][1] %>% pull())[[1]]

summarise_clim <- function(single_sp_clim){
  single_sp_clim <- single_sp_clim %>% 
    pivot_longer(cols = starts_with("pet_cmip") | starts_with("cwd_cmip"),
                 names_to = c(".value", "mod_number"),
                 names_pattern = "(.*?)(\\d{1,2})$")
  
  single_sp_change <- single_sp_clim %>% 
    mutate(pet_change = pet_cmip_end - pet_cmip_start,
           cwd_change = cwd_cmip_end - cwd_cmip_start) %>% 
    group_by(mod_number) %>% 
    summarise(cwd_change = mean(cwd_change),
              pet_change = mean(pet_change))
  return(single_sp_change)
}

sp_clim <- sp_clim %>% 
  mutate(clim_smry = map(clim_cmip_sp, .f = summarise_clim)) %>% 
  select(sp_code, clim_smry) %>% 
  unnest(clim_smry)


sp_clim %>% 
  pull(cwd_change) %>% 
  mean()

sp_clim %>% 
  pull(cwd_change) %>% 
  quantile(c(0.025, 0.975))

sp_clim %>% 
  pull(pet_change) %>% 
  mean()

sp_clim %>% 
  pull(pet_change) %>% 
  quantile(c(0.025, 0.975))

# Supporting concerns about the emergence of “novel climates” (39, 40), 3.8% of 
# species’ ranges are predicted to face a mean CWD in 2100 that exceeds the
# mean CWD experienced anywhere in that species’ historic, climatic range 
# (Figure 5, D). For some species (e.g., Pinus pinea and Quercus faginea),
# more than half of their current range is projected to be drier than the 
# driest parts of their historic range.
sp_plot_dat %>%
  ungroup() %>% 
  group_by(cwd_end_bin) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  print()

extreme_sp_change <- sp_plot_dat %>%
  ungroup() %>% 
  group_by(sp_code, cwd_end_bin) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  filter(cwd_end_bin == "Beyond prior max") %>% 
  arrange(desc(freq)) %>% 
  print()


# # When compared to a species’ historic climatic distribution, 84 percent of the 
# # area of species’ historic range is forecasted to have a higher average annual 
# # PET than the species’ historic mean PET, while 77 percent of the area is expected 
# # to by dryer than the historic mean CWD. 
# 
# sp_plot_dat %>%
#   ungroup() %>% 
#   group_by(pet_end_bin) %>% 
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n)) %>% 
#   filter(pet_end_bin != "Below prior mean") %>% 
#   pull(freq) %>% 
#   sum() %>% 
#   print()
# 
# sp_plot_dat %>%
#   ungroup() %>% 
#   group_by(cwd_end_bin) %>% 
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n)) %>% 
#   filter(cwd_end_bin != "Below prior mean") %>% 
#   pull(freq) %>% 
#   sum() %>% 
#   print()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Largest growth declines predicted in wet and hot regions  --------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combining our estimates of sensitivity to CWD and PET with predicted changes 
# in climate, we estimate that tree growth will decline by 10.4% (-38.4% to +1.1%) by 2100. 
agg_stats$rwi_pred_change %>% mean()
agg_stats$rwi_pred_change %>% quantile(c(0.025, 0.975))


# 92.5% of grid cells (representing combinations of species and location) have a mean 
# predicted decline in RWI, and 37.2% of grid cells have a distribution of predicted 
# changes in which >97.5% of all predictions are negative. In contrast, 0.4% of grid
# cells have a distribution of predicted RWI changes in which >97.5% of predictions
# are positive. 
changes <- sp_predictions %>%
  mutate(rwi_change_qual = rwi_pred_change_mean < 0) %>% 
  group_by(rwi_change_qual) %>% 
  tally() %>% 
  mutate(proportion = prop.table(n)) %>% 
  print()


changes <- sp_predictions %>%
  mutate(rwi_change_qual = case_when(
    rwi_pred_change_025 > 0 ~ "increase",
    rwi_pred_change_975 < 0 ~ "decrease",
    TRUE ~ "neither"
  )) %>% 
  group_by(rwi_change_qual) %>% 
  tally() %>% 
  mutate(proportion = prop.table(n)) %>% 
  print()


# the wetter but hotter-than-average portions of species’ ranges are projected to experience 
# a 17.3% decline in growth (-53.7% to -0.8%). In contrast, drier but cooler-than-average 
# locations are predicted to experience a 10.6% decline in growth (-33.5% to +0.9%). 

sp_predictions %>%
  filter(cwd_hist<0, pet_hist >0) %>% 
  summarise(mean(rwi_pred_change_mean),
            mean(rwi_pred_change_025),
            mean(rwi_pred_change_975))

sp_predictions %>%
  filter(cwd_hist>0, pet_hist<0) %>% 
  summarise(mean(rwi_pred_change_mean),
            mean(rwi_pred_change_025),
            mean(rwi_pred_change_975))



# In contrast to a neutral model in which predicted sensitivity to PET and CWD 
# are constant across species’ historic, climatic ranges, our model predicts 
# that the wettest 10% of grid cells will experience XX percentage point larger 
# declines in growth (XX-XX), while the dryest 10% of grid cells will experience 
# XX percentage point smaller declines in growth (XX-XX).


pet_band_bounds <- c(.9, 1.1)
subset <- sp_predictions %>%
  filter(pet_hist > pet_band_bounds[1], pet_hist < pet_band_bounds[2])

  
cwd_quantiles <- subset %>% 
  pull(cwd_hist) %>% 
  quantile(c(0.05, 0.95)) %>% 
  print()

subset <- subset %>% 
  mutate(extremes = ifelse(cwd_hist < cwd_quantiles[1], "wet", ifelse(cwd_hist > cwd_quantiles[2], "dry", "neither")))

subset %>%
  filter(extremes != "neither") %>%
  group_by(extremes) %>%
  summarise(rwi_pred_pclim_change_dif_mean = mean(rwi_pred_pclim_change_dif_mean),
            rwi_pred_pclim_change_dif_025 = mean(rwi_pred_pclim_change_dif_025),
            rwi_pred_pclim_change_dif_975 = mean(rwi_pred_pclim_change_dif_975),
            baseline = mean(rwi_pred_change_mean))



subset <- sp_predictions %>%
  filter(pet_hist > 0, cwd_hist < 0)
subset$rwi_pred_pclim_change_dif_mean %>% hist()

subset %>%
  # filter(extremes != "neither") %>% 
  # group_by(extremes) %>% 
  summarise(rwi_pred_pclim_change_dif_mean = mean(rwi_pred_pclim_change_dif_mean),
            rwi_pred_pclim_change_dif_025 = mean(rwi_pred_pclim_change_dif_025),
            rwi_pred_pclim_change_dif_975 = mean(rwi_pred_pclim_change_dif_975),
            baseline = mean(rwi_pred_change_mean))



pet_bins <- seq(-2.5, 2.5, 0.2)
test <- sp_predictions %>% 
  mutate(pet_group = cut(pet_hist, pet_bins))
wet_edge_thresh <- test %>% 
  group_by(sp_code, pet_group) %>% 
  summarize(wet_edge = quantile(cwd_hist, 0.1),
            dry_edge = quantile(cwd_hist, 0.9))

test <- test %>% 
  left_join(wet_edge_thresh, by = c("sp_code", "pet_group"))

subset <- test %>% 
  filter(cwd_hist > dry_edge)

dry_edge_result <- test %>%
  filter(cwd_hist > dry_edge) %>% 
  group_by(pet_group) %>% 
  summarise(rwi_pred_pclim_change_dif_mean = mean(rwi_pred_pclim_change_dif_mean),
            rwi_pred_pclim_change_dif_025 = mean(rwi_pred_pclim_change_dif_025),
            rwi_pred_pclim_change_dif_975 = mean(rwi_pred_pclim_change_dif_975),
            baseline = mean(rwi_pred_change_mean),
            neutral = mean(rwi_pclim_change_mean))

wet_edge_result <- test %>%
  filter(cwd_hist < wet_edge) %>% 
  group_by(pet_group) %>% 
  summarise(rwi_pred_pclim_change_dif_mean = mean(rwi_pred_pclim_change_dif_mean),
            rwi_pred_pclim_change_dif_025 = mean(rwi_pred_pclim_change_dif_025),
            rwi_pred_pclim_change_dif_975 = mean(rwi_pred_pclim_change_dif_975),
            baseline = mean(rwi_pred_change_mean),
            neutral = mean(rwi_pclim_change_mean))

# # pet_band_bounds <- c(0.9, 1.1)
# # cwd_quantiles <- sp_predictions %>% 
# #   filter(pet_hist > pet_band_bounds[1], pet_hist < pet_band_bounds[2]) %>% 
# #   pull(cwd_hist) %>% 
# #   quantile(c(0.01, 0.99)) %>% 
# #   print()
# # 
# # sp_predictions %>% 
# #   filter(pet_hist > pet_band_bounds[1], pet_hist < pet_band_bounds[2]) %>% 
# #   filter(cwd_hist < cwd_quantiles[1]) %>% 
# #   select(rwi_pred_change_mean, rwi_pred_change_025, rwi_pred_change_975) %>% 
# #   summary()
# # 
# # sp_predictions %>% 
# #   filter(pet_hist > pet_band_bounds[1], pet_hist < pet_band_bounds[2]) %>% 
# #   filter(cwd_hist > cwd_quantiles[2]) %>% 
# #   select(rwi_pred_change_mean, rwi_pred_change_025, rwi_pred_change_975) %>% 
# #   summary()
# # 
# # sp_predictions %>% filter(cwd_hist > 1) %>% select(rwi_pred_change_mean) %>% summary()
# # sp_predictions %>% filter(cwd_hist < -0) %>% select(rwi_pred_change_mean) %>% summary()
# # sp_predictions$rwi_pred_change_mean %>% quantile(.5)
# # sp_predictions %>% filter(rwi_pred_change_mean < -0.1) %>% select(cwd_hist) %>% summary()
# # # Proportion of grid cells forecasted to experience significant increase in growth
# 
# 
# # For example, among the grid cells with a historic PET ranging from 0.9 to 1.1, 
# # the wettest 1% of grid cells (historic CWD < -0.41) are predicted to experience 
# # a mean decline in growth of 25.0%. In contrast, the driest 1% of grid cells 
# # (historic CWD > 2.02) are predicted to experience a mean decline in growth of 18.5%.  
# 
# 
hot_cells$pet_hist %>% summary()

cwd_quantiles <- hot_cells %>%
  pull(cwd_hist) %>%
  quantile(c(0.01, 0.99)) %>%
  print()

hot_cells <- hot_cells %>%
  mutate(wet = ifelse(cwd_hist < cwd_quantiles[1], "wet", ifelse(cwd_hist > cwd_quantiles[2], "dry", "neither")))

hot_cell_comparison <- hot_cells %>%
  filter(wet != "neither") %>%
  group_by(iter_idx, wet) %>%
  summarise(rwi_pred_change = mean(rwi_pred_change)) %>%
  pivot_wider(names_from = wet, values_from = rwi_pred_change) %>%
  mutate(difference = dry - wet)

hot_cell_comparison$wet %>% mean()
hot_cell_comparison$wet %>% quantile(c(0.025, 0.975))

hot_cell_comparison$dry %>% mean()
hot_cell_comparison$dry %>% quantile(c(0.025, 0.975))

hot_cell_comparison$difference %>% mean()
hot_cell_comparison$difference %>% quantile(c(0.025, 0.975))

t.test(hot_cell_comparison$wet, hot_cell_comparison)
