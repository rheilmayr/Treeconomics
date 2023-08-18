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

select <- dplyr::select
summarize <- dplyr::summarize

options(scipen=999)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. CMIP climate predictions
sp_clim <- read_rds(paste0(wdir, "2_output/climate/sp_clim_predictions.gz"))

# 2. Site-level regressions
flm_df <- read_csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv"))
flm_df <- flm_df %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

# 3. First stage block bootstraps
block_draw_df <- read_rds(paste0(wdir, "2_output/second_stage/mc_sample.gz"))

# 4. Second stage model
ss_df <- read_rds(paste0(wdir, "2_output/second_stage/ss_bootstrap.rds"))

# 5. Prediction results
rwi_list <- list.files(paste0(wdir, "2_output/predictions/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))
agg_stats <- read_rds(file = paste0(wdir, "out/predictions/mc_agg_stats.gz"))




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


# Trees also exhibited heterogeneity in their responses to annual increases in 
# energy availability (Figure 3, D-F). In the most energy limited portions of 
# species’ range, trees often experienced large increases in growth in response 
# to elevated PET. For example, in the coldest 10% of sites (<1.5 SD below the 
# species’ historic mean PET), growth increased 18.5% (8.4% to 25.0%) in response 
# to a 1 SD increase in PET. However, this positive growth response to elevated 
# annual PET declined across a gradient of historic PET (Figure 3, F). 
# In sites with above-average PET, tree growth was not significantly associated with changes in PET.
pet_q10 <- flm_df$pet.spstd %>% 
  quantile(c(0.1)) %>% 
  print()

pet_low_fs_bs <- block_draw_df %>% 
  filter(pet.spstd < pet_q10) %>%
  group_by(boot_id) %>% 
  summarise(pet_coef = mean(pet_coef))

pet_low_fs_bs %>%
  pull(pet_coef) %>% 
  mean() 

pet_low_fs_bs %>%
  pull(pet_coef) %>% 
  quantile(c(0.025, 0.975))


pet_high_fs_bs <- block_draw_df %>% 
  filter(pet.spstd > 0 ) %>%
  group_by(boot_id) %>% 
  summarise(pet_coef = mean(pet_coef))

pet_high_fs_bs %>%
  pull(pet_coef) %>% 
  mean() 

pet_high_fs_bs %>%
  pull(pet_coef) %>% 
  quantile(c(0.025, 0.975))


# To characterize this variability in CWD sensitivity, we modeled sensitivity as a 
# quadratic function of historic CWD and PET (Figure 4, Panels B and C). The regression 
# results highlight that observed sensitivity declined with increasing CWD (slope at 
# the species’ mean, historic CWD = 0.015, CI = 0.006 to 0.064). 
ss_df %>%
  summarise(mean_cwd = mean(cwd_cwd),
            lower_cwd = quantile(cwd_cwd, 0.025),
            upper_cwd = quantile(cwd_cwd, 0.975))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Climatic changes across species ranges  --------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CMIP5 model ensembles predict that the entirety of every species’ range 
# analyzed will experience increases in both mean PET and CWD 
sp_plot_dat <- sp_predictions %>% 
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
# single_sp_clim <- (sp_clim[1, 10][1] %>% pull())[[1]]

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


extreme_sp_change %>% pull(freq) %>% mean()


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
    rwi_pred_change_050 > 0 ~ "increase",
    rwi_pred_change_950 < 0 ~ "decrease",
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
            mean(rwi_pred_change_975)) %>% 
  print()

sp_predictions %>%
  filter(cwd_hist>0, pet_hist<0) %>% 
  summarise(mean(rwi_pred_change_mean),
            mean(rwi_pred_change_025),
            mean(rwi_pred_change_975)) %>% 
  print()



## Grid cells predicted to experience increases in growth are located almost 
## entirely in the most energy-limited parts of species’ historic ranges.
growth_cells <- sp_predictions %>%
  filter(rwi_pred_change_050>0) 

growth_cells %>% 
  pull(pet_hist) %>% 
  quantile(c(0.25, 0.5, 0.75, 0.8))


# 
# 
# library(janitor)
# species_df <- species_df %>%
#   clean_names() %>%
#   select(-range_map_source)
# 
# species_summary <- species_df %>%
#   mutate(gen_code = substr(species_code, 1, 2)) %>%
#   group_by(gen_code) %>%
#   summarise(n = sum(number_of_sites),
#             mean_pet = weighted.mean(mean_pet, number_of_sites),
#             mean_cwd = weighted.mean(mean_cwd, number_of_sites)) %>%
#   mutate(aet_estimate = mean_pet - mean_cwd,
#          aet_pct = aet_estimate / mean_pet)
# 
# species_summary <- species_summary %>%
#   filter(n > 50) %>%
#   arrange(desc(mean_pet))
# species_summary %>% arrange(desc(mean_cwd))
# 
# 
