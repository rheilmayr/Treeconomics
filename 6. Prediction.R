#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Create predictions of growth impacts from climate change
#
# Input files:
# - ss_mod: R model object saved from Second stage
# - 
#
# ToDo:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(fixest)
# library(raster)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(tidyverse)
library(dtplyr)
library(prediction)
library(tictoc)
library(furrr)
library(snow)
library(profvis)
library(tmap)
library(tidylog)

n_cores <- availableCores() - 6
future::plan(multisession, workers = n_cores)

my_seed <- 5597

n_mc <- 10000


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# Create output directory
out_dir <- paste0(wdir,"out/predictions/sp_rwi_pred_", as.character(n_mc), "/")
dir.create(file.path(out_dir), showWarnings = FALSE)

# 1. Second stage model
mod_df <- read_rds(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))
mod_df <- mod_df %>% 
  rename(iter_idx = boot_id)

# 2. Species-standardized historic and future climate
sp_clim <- read_rds(paste0(wdir, "out/climate/sp_clim_predictions.gz"))
species_list <- sp_clim %>% select(sp_code)

# 3. Constant sensitivites (median of first stage, rather than running 2nd stage)
constant_sensitivities <- read_rds(paste0(wdir, "out/first_stage/constant_sensitivities.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Assign MC coefs and CMIP models  ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Join second stage coefficients to species list
sp_mc <- species_list %>% 
  select(sp_code) %>% 
  crossing(iter_idx = seq(n_mc)) %>% 
  left_join(mod_df, by = "iter_idx")


## Assign specific cmip realization to each MC iteration 
n_cmip_mods <- 47
cmip_assignments <- tibble(iter_idx = seq(1, n_mc)) %>% 
  mutate(cmip_idx = sample(seq(n_cmip_mods), n_mc, replace = TRUE))

## Join cmip model assignments
sp_mc <- sp_mc %>% 
  left_join(cmip_assignments, by = "iter_idx")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define functions ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predict_sens <- function(sppp_code, int_int, int_cwd, int_pet, cwd_int, 
                         cwd_cwd, cwd_pet, pet_int, pet_cwd, pet_pet){
  ## Function used to predict species' sensitivity rasters based on historic 
  ## climate and randomly drawn parameters from second stage model.
  
  select <- dplyr::select
  
  sp_df <- (sp_clim %>% 
    filter(sp_code == sppp_code) %>% 
    pull(clim_historic_sp))[[1]] %>% 
    lazy_dt()
  
  sp_df <- sp_df %>% 
    rename(cwd_hist = cwd,
           pet_hist = pet) %>%
    mutate(cwd_sens = cwd_int + cwd_cwd * cwd_hist + cwd_pet * pet_hist,
           pet_sens = pet_int + pet_cwd * cwd_hist + pet_pet * pet_hist,
           intercept = int_int + int_cwd * cwd_hist + int_pet * pet_hist) %>% 
    select(-cwd_hist, -pet_hist) %>% 
    as_tibble()
  
  return(sp_df)
}


calc_rwi_partials <- function(sppp_code, cmip_id, sensitivity){
  select <- dplyr::select
  
  ## Function used to predict species' RWI rasters based on predicted 
  ## sensitivity raster and assigned CMIP model of future climate. Also
  ## integrates calculations of partialling our climate / sensitivity variations
  
  ## Predict RWI under CMIP scenario
  sp_fut_clim <- sp_clim %>% 
    filter(sp_code == sppp_code) %>% 
    pull(clim_cmip_sp)
  sp_fut_clim <- sp_fut_clim[[1]] %>% 
    lazy_dt() %>% 
    select(x,y,
           cwd_cmip_end = paste0("cwd_cmip_end", as.character(cmip_id)),
           pet_cmip_end = paste0("pet_cmip_end", as.character(cmip_id)),
           cwd_cmip_start = paste0("cwd_cmip_start", as.character(cmip_id)),
           pet_cmip_start = paste0("pet_cmip_start", as.character(cmip_id))) %>% 
    left_join(sensitivity, by = c("x", "y")) %>% 
    mutate(rwi_pred_end = intercept + (pet_cmip_end * pet_sens) + (cwd_cmip_end * cwd_sens),
           rwi_pred_start = intercept + (pet_cmip_start * pet_sens) + (cwd_cmip_start * cwd_sens),
           rwi_pclim_end = constant_sensitivities$int + (pet_cmip_end * constant_sensitivities$pet) + (cwd_cmip_end * constant_sensitivities$cwd),
           rwi_pclim_start = constant_sensitivities$int + (pet_cmip_start * constant_sensitivities$pet) + (cwd_cmip_start * constant_sensitivities$cwd)) %>% 
    select(x,y,
           cwd_sens,
           pet_sens,
           intercept,
           rwi_pred_end,
           rwi_pred_start,
           rwi_pclim_end,
           rwi_pclim_start) %>% 
    as_tibble()
  
  return(sp_fut_clim)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compute sensitivity, RWI for each species by MC combination ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pull_layer <- function(brick, layer_name){
  pulled_layer <- brick %>% subset(layer_name)
}


calc_rwi_quantiles <- function(spp_code, mc_data, parallel = TRUE){
  set.seed(my_seed) # Re-setting seed at start of each iteration to ensure interrupted jobs still produce replicable results
  tic()

  ## Iterating through each species
  print(spp_code)
  
  ## Calculate n_mc versions of species' sensitivity raster
  if (parallel == TRUE) {
    sp_predictions <- mc_data %>% 
      mutate(sensitivity = future_pmap(list(sppp_code = spp_code,
                                            int_int = int_int,
                                            int_cwd = int_cwd,
                                            int_pet = int_pet,
                                            cwd_int = cwd_int,
                                            cwd_cwd = cwd_cwd,
                                            cwd_pet = cwd_pet,
                                            pet_int = pet_int,
                                            pet_cwd = pet_cwd,
                                            pet_pet = pet_pet),
                                       .f = predict_sens,
                                       .options = furrr_options(seed = my_seed, 
                                                                packages = c( "dplyr", "raster", "dtplyr"))))
  } else {
    sp_predictions <- mc_data %>% 
      mutate(sensitivity = pmap(list(sppp_code = spp_code,
                                            int_int = int_int,
                                            int_cwd = int_cwd,
                                            int_pet = int_pet,
                                            cwd_int = cwd_int,
                                            cwd_cwd = cwd_cwd,
                                            cwd_pet = cwd_pet,
                                            pet_int = pet_int,
                                            pet_cwd = pet_cwd,
                                            pet_pet = pet_pet),
                                       .f = predict_sens))
  }
  sp_predictions <- sp_predictions %>% 
    select(iter_idx, cmip_idx, sensitivity)
  gc(verbose = TRUE)

  ## Predict future RWI for each of n_mc run
  if (parallel == TRUE){
    sp_predictions <- sp_predictions %>% 
      mutate(rwi_predictions = future_pmap(list(sppp_code = spp_code,
                                                cmip_id = cmip_idx,
                                                sensitivity = sensitivity),
                                           .f = calc_rwi_partials,
                                           .options = furrr_options(seed = my_seed,
                                                                    packages = c("raster", "dplyr", "dtplyr"))))     
  } else {
    sp_predictions <- sp_predictions %>% 
      mutate(rwi_predictions = pmap(list(sppp_code = spp_code,
                                                cmip_id = cmip_idx,
                                                sensitivity = sensitivity),
                                           .f = calc_rwi_partials)) 
  }
  sp_predictions <- sp_predictions %>% 
    select(iter_idx, rwi_predictions) %>% 
    unnest(rwi_predictions)
  gc(verbose = TRUE)
  
  sp_predictions <- sp_predictions %>% 
    lazy_dt()
  
  # ## Store historic climate to re-join later
  # hist_df <- sp_predictions %>% 
  #   filter(iter_idx == 1) %>% 
  #   select(x, y, cwd_hist, pet_hist)
    
  ## For each species, calculate cell-wise quantiles of variables from n_mc runs
  sp_predictions <- sp_predictions %>% 
    mutate(rwi_pred_change = rwi_pred_end - rwi_pred_start,
           rwi_pclim_change = rwi_pclim_end - rwi_pclim_start) %>% 
    group_by(x, y) %>% 
    summarise(rwi_pred_mean = mean(rwi_pred_end),
              rwi_pred_025 = quantile(rwi_pred_end, 0.025),
              rwi_pred_975 = quantile(rwi_pred_end, 0.975),
              rwi_pclim_mean = mean(rwi_pclim_end),
              rwi_pclim_025 = quantile(rwi_pclim_end, 0.025),
              rwi_pclim_975 = quantile(rwi_pclim_end, 0.975),
              rwi_pred_change_mean = mean(rwi_pred_change),
              rwi_pred_change_025 = quantile(rwi_pred_change, 0.025),
              rwi_pred_change_975 = quantile(rwi_pred_change, 0.975),
              rwi_pclim_change_mean = mean(rwi_pclim_change),
              rwi_pclim_change_025 = quantile(rwi_pclim_change, 0.025),
              rwi_pclim_change_975 = quantile(rwi_pclim_change, 0.975),
              cwd_sens = mean(cwd_sens),
              pet_sens = mean(pet_sens),
              int_sens = mean(intercept),
              # cwd_cmip_start = mean(cwd_cmip_start),
              # pet_cmip_start = mean(pet_cmip_start),
              # cwd_cmip_end = mean(cwd_cmip_end),
              # pet_cmip_end = mean(pet_cmip_end),
              sp_code = spp_code,
              .groups = "drop")

  ## Add back observed climate data
  sp_hist <- (sp_clim %>% 
              filter(sp_code == spp_code) %>% 
              pull(clim_historic_sp))[[1]] %>% 
    lazy_dt() %>% 
    rename(cwd_hist = cwd,
           pet_hist = pet)
  
  sp_predictions <- sp_predictions %>% 
    left_join(sp_hist, by = c("x", "y"))
  
  ## Add observed and predicted climate data
  sp_cmip <- (sp_clim %>% 
    filter(sp_code == spp_code) %>% 
    pull(clim_cmip_sp))[[1]]
  
  sp_cmip <- sp_cmip %>%
    rowwise() %>% 
    mutate(pet_cmip_end_mean = mean(c_across(starts_with("pet_cmip_end"))),
           cwd_cmip_end_mean = mean(c_across(starts_with("cwd_cmip_end"))),
           pet_cmip_start_mean = mean(c_across(starts_with("pet_cmip_start"))),
           cwd_cmip_start_mean = mean(c_across(starts_with("cwd_cmip_start")))) %>% 
    select(x, y, cwd_cmip_start_mean, cwd_cmip_end_mean, pet_cmip_start_mean, pet_cmip_end_mean)

  sp_predictions <- sp_predictions %>% 
    left_join(sp_cmip, by = c("x", "y")) %>% 
    as_tibble()
  
  ## Write out
  sp_predictions %>% 
    write_rds(file = paste0(out_dir, spp_code, ".gz"), compress = "gz")
  
  toc()
  return("done")
}


mc_nests <- sp_mc %>%
  group_by(sp_code) %>%
  nest() %>% 
  drop_na()

# Generally have memory issues with 38 (LAGM), and 93 (PISY) - need to run these with two cores
# large_range_sp <- c("lagm", "pisy")
# spp_code = "abal"
# mc_data = mc_nests %>% filter(sp_code == spp_code) %>% pull(data)
# mc_data = mc_data[[1]]
# parallel = FALSE

# mc_nests_large <- mc_nests %>% 
#   filter((sp_code %in% large_range_sp)) %>% 
#   mutate(predictions = pmap(list(spp_code = sp_code,
#                                  mc_data = data,
#                                  parallel = TRUE),
#                             .f = calc_rwi_quantiles))

mc_nests_small <- mc_nests %>% 
  # filter(!(sp_code %in% large_range_sp)) %>% 
  mutate(predictions = pmap(list(spp_code = sp_code,
                                 mc_data = data,
                                 parallel = TRUE),
                              .f = calc_rwi_quantiles)) 

 



# # Profiling of main function
# spp_code = "juex"
# mc_data = (mc_nests %>% filter(sp_code == spp_code) %>% pull(data))[[1]]
# l = profvis(calc_rwi_quantiles(spp_code, mc_data))






# %>% 
#   select(-data) %>% 
#   unnest(predictions)

# mc_nests %>% 
#   saveRDS(file = paste0(wdir,"out/predictions/rwi_predictions.rds"))


# # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # # Thought experiments - partialling out mechanisms    ---------------
# # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # calc_rwi_partial_sens <- function(cmip_rast, sensitivity){ # NOTE: Should the means be calculated across full range rather than by each species?
# #   mean_fut_cwd <- cmip_rast %>% subset("cwd.spstd") %>% cellStats(stat = "mean")
# #   mean_fut_pet <- cmip_rast %>% subset("pet.spstd") %>% cellStats(stat = "mean")
# #   cwd_sens = sensitivity %>% subset("cwd_sens")
# #   pet_sens = sensitivity %>% subset("pet_sens")
# #   intercept = sensitivity %>% subset("intercept")
# #   rwi_rast <- intercept + (mean_fut_cwd * cwd_sens) + (mean_fut_pet * pet_sens)
# #   names(rwi_rast) = "rwi_psens"
# #   return(rwi_rast)
# # }
# # 
# # calc_rwi_partial_clim <- function(cmip_rast, sensitivity){ # NOTE: Should the means be calculated across full range rather than by each species?
# #   mean_cwd_sens <- sensitivity %>% subset("cwd_sens") %>% cellStats(stat = "mean")
# #   mean_pet_sens <- sensitivity %>% subset("pet_sens") %>% cellStats(stat = "mean")
# #   mean_intercept <-sensitivity %>% subset("intercept") %>% cellStats(stat = "mean") 
# #   rwi_rast <- mean_intercept + (cmip_rast$cwd.spstd * mean_cwd_sens) + (cmip_rast$pet.spstd * mean_pet_sens)
# #   names(rwi_rast) = "rwi_pclim"
# #   return(rwi_rast)
# # }
# # 
# # 
# # sp_predictions <- sp_predictions %>% 
# #   mutate(rwi_predictions_partial_sens = map2(.x = clim_future_sp, .y = sensitivity, calc_rwi_partial_sens),
# #          rwi_predictions_partial_clim = map2(.x = clim_future_sp, .y = sensitivity, calc_rwi_partial_clim))
# 
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Stack rasters into dataframe ------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create_prediction_df <- function(spp_predictions){
#   sp_fut <- (spp_predictions %>% 
#                pull(clim_future_sp))[[1]]
#   names(sp_fut) <- c("cwd.fut", "pet.fut")
#   
#   sp_hist <- (spp_predictions %>% 
#                 pull(clim_historic_sp))[[1]]
#   sp_sens  <- (spp_predictions %>% 
#                  pull(sensitivity))[[1]]
#   sp_rwi  <- (spp_predictions %>% 
#                 pull(rwi_predictions))[[1]]
#   
#   sp_rwi_psens  <- (spp_predictions %>% 
#                       pull(rwi_predictions_partial_sens))[[1]]
#   
#   sp_rwi_pclim  <- (spp_predictions %>% 
#                       pull(rwi_predictions_partial_clim))[[1]]
#   
#   clim_compare <- brick(c(sp_fut, sp_hist, sp_sens, sp_rwi, sp_rwi_psens, sp_rwi_pclim))
#   clim_compare <- clim_compare %>% 
#     as.data.frame(xy = TRUE) %>% 
#     drop_na()
#   return(clim_compare)
# }
# 
# sp_prediction_test <- sp_predictions %>% 
#   group_by(sp_code, iter_idx) %>% 
#   nest() %>% 
#   mutate(pred_df = map(data, create_prediction_df)) %>% 
#   select(-data) %>% 
#   unnest(cols = pred_df) %>% 
#   mutate(cwd_change = cwd.fut - cwd.spstd,
#          pet_change = pet.fut - pet.spstd,
#          rwi_null = cwd.spstd * cwd_sens + pet.spstd * pet_sens + intercept,
#          rwi_change = rwi_pred - rwi_null,
#          rwi_change_psens = rwi_psens - rwi_null,
#          rwi_change_pclim = rwi_pclim - rwi_null)
# 
# 
# sp_predictions %>% 
#   saveRDS(file = paste0(wdir,"out/predictions/sp_predictions.rds") )





# crs_template <- crs(cwd_future)
# cwd_df <- cwd_rast %>% as.data.frame(xy = TRUE) 
# raster_template <- cwd_df %>% select(x,y)
# cwd_df <- cwd_df %>% 
#   drop_na()
# 
# cwd_df2 <- raster_template %>% 
#   left_join(cwd_df, by = c("x", "y"))
# cwd_rast2 <- rasterFromXYZ(cwd_df2, crs = crs)
# 
# tm_raster(cwd_rast2)
# data("World")
# 
# tmap_mode("view")
# tm_shape(cwd_rast) +
#   tm_raster()
# 
# tmap_mode("view")
# tm_shape(cwd_vals) +
#   tm_raster()
# 
# 
# 
# 
# 
# 
# ############ FUNCTION GRAVEYARD
# 
# calc_rwi <- function(sppp_code, cmip_id, sensitivity){
#   ## Function used to predict species' RWI rasters based on predicted 
#   ## sensitivity raster and assigned CMIP model of future climate
#   
#   sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, sppp_code, ".gz"))
#   
#   cmip_rast <- sp_fut_clim %>% 
#     filter(cmip_idx == cmip_id) %>% 
#     pull(clim_future_sp)
#   
#   cwd_sens = sensitivity %>% subset("cwd_sens")
#   pet_sens = sensitivity %>% subset("pet_sens")
#   intercept = sensitivity %>% subset("intercept")
#   rwi_rast <- intercept + (cmip_rast[[1]]$cwd.spstd * cwd_sens) + (cmip_rast[[1]]$pet.spstd * pet_sens)
#   names(rwi_rast) = "rwi_pred"
#   return(rwi_rast)
# }
# 
# quantiles <- function(x){
#   ## Defines quantiles used to summarize MC runs
#   quantile(x, c(0.025, 0.975), na.rm=TRUE)
# }
# 
# # extract_quantiles <- function(rwi_preds){
# #   ## Extracts desired quantiles for MC runs
# #   
# #   rwi_quantiles <- rwi_preds %>% 
# #     pull(rwi_predictions) %>%
# #     brick() %>% 
# #     calc(quantiles)
# #   return(rwi_quantiles)
# # }
# 
# 
# calc_mean_fut_clim <- function(sppp_code){
#   sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, sppp_code, ".gz"))
#   clim_pulls <- sp_fut_clim %>%
#     mutate(cwd.fut = pmap(list(brick = clim_future_sp, layer_name = "cwd.spstd"),
#                           .f = pull_layer),
#            pet.fut = pmap(list(brick = clim_future_sp, layer_name = "pet.spstd"),
#                           .f = pull_layer))
#   
#   cwd.fut.q <- clim_pulls %>%
#     pull(cwd.fut) %>%
#     brick()
#   cwd.fut.q <- raster::mean(cwd.fut.q)
#   names(cwd.fut.q) = "cwd.fut"
#   
#   pet.fut.q <- clim_pulls %>%
#     pull(pet.fut) %>%
#     brick()
#   pet.fut.q <- raster::mean(pet.fut.q)
#   names(pet.fut.q) = "pet.fut"
#   
#   out_brick <- brick(c(cwd.fut.q, pet.fut.q))
#   return(out_brick)
# }
