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
library(raster)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(tidyverse)
library(prediction)
library(tictoc)
library(furrr)
library(snow)

n_cores <- 12
future::plan(multisession, workers = n_cores)

my_seed <- 5597
set.seed(my_seed)

n_mc <- 1000

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'
out_dir <- paste0(wdir,"out/predictions/sp_rwi_pred_", as.character(n_mc), "/")
dir.create(file.path(out_dir), showWarnings = FALSE)


# 1. Second stage model
old_mod_df <- readRDS(paste0(wdir, "out/second_stage/ss_mc_mods.rds"))
boot_ss <- readRDS(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))
 
mod_df <- boot_ss$t
colnames(mod_df) <- c("int_int", "int_cwd", "int_pet",
                       "cwd_int", "cwd_cwd", "cwd_pet",
                       "pet_int", "pet_cwd", "pet_pet")
mod_df <- mod_df %>% 
  as_tibble() %>% 
  mutate(iter_idx = seq(1:10000)) %>% 
  filter(iter_idx <= n_mc)
  
# 2. Species-standardized historic climate
sp_hist_clim <- readRDS(paste0(wdir, "out/climate/sp_clim_historic.rds"))
species_list <- sp_hist_clim %>% select(sp_code)

# 3. Directory of species-standardized future possible climates
sp_fut_clim_dir <- paste0(wdir, "out/climate/sp_clim_predictions/")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Assign MC coefs and CMIP models  ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Join second stage coefficients to species list
sp_mc <- species_list %>% 
  select(sp_code) %>% 
  crossing(iter_idx = seq(n_mc))


# sp_fut_clim_smpl2 <- read_rds(paste0(sp_fut_clim_dir, "abal"))
# write_rds(sp_fut_clim_smpl2, paste0(wdir, "out//climate//sp_clim_predictions//", spp_code, ".rds"))



## Assign specific cmip realization to each MC iteration 
sp_fut_clim_smpl <- readRDS(paste0(sp_fut_clim_dir, "abal.gz"))
n_cmip_mods <- sp_fut_clim_smpl %>% pull(cmip_idx) %>% unique() %>% length()
cmip_assignments <- tibble(iter_idx = seq(1, n_mc)) %>% 
  mutate(cmip_idx = sample(seq(n_cmip_mods), n_mc, replace = TRUE))

## Join cmip model assignments
sp_mc <- sp_mc %>% 
  left_join(cmip_assignments, by = "iter_idx")

## Join second stage coefficients
# mod_df <- mod_df %>% 
#   group_by(iter_idx) %>% 
#   nest() %>% 
#   rename(ss_coefs = data)

sp_mc <- sp_mc %>% 
  left_join(mod_df, by = "iter_idx") 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define functions ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predict_sens <- function(sppp_code, int_int, int_cwd, int_pet, cwd_int, 
                         cwd_cwd, cwd_pet, pet_int, pet_cwd, pet_pet){
  ## Function used to predict species' sensitivity rasters based on historic 
  ## climate and randomly drawn parameters from second stage model.
  
  select <- dplyr::select
  
  clim_historic_sp <- sp_hist_clim %>% 
    filter(sp_code == sppp_code) %>% 
    pull(clim_historic_sp)
  cwd_rast <- clim_historic_sp[[1]] %>% subset("cwd.spstd")
  pet_rast <- clim_historic_sp[[1]] %>% subset("pet.spstd")
  
  cwd_sens <- cwd_int +
              cwd_cwd * cwd_rast +
              cwd_pet * pet_rast
  names(cwd_sens) = "cwd_sens"
  
  pet_sens <- pet_int +
              pet_cwd * cwd_rast +
              pet_pet * pet_rast
  names(pet_sens) = "pet_sens"
  
  intercept <- int_int +
               int_cwd * cwd_rast +
               int_pet * pet_rast
  names(intercept) = "intercept"
  
  sensitivity <- raster::brick(cwd_sens, pet_sens, intercept) 
  return(sensitivity)
}

calc_rwi <- function(sppp_code, cmip_id, sensitivity){
  ## Function used to predict species' RWI rasters based on predicted 
  ## sensitivity raster and assigned CMIP model of future climate
  
  sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, sppp_code, ".gz"))
  
  cmip_rast <- sp_fut_clim %>% 
    filter(cmip_idx == cmip_id) %>% 
    pull(clim_future_sp)
  
  cwd_sens = sensitivity %>% subset("cwd_sens")
  pet_sens = sensitivity %>% subset("pet_sens")
  intercept = sensitivity %>% subset("intercept")
  rwi_rast <- intercept + (cmip_rast[[1]]$cwd.spstd * cwd_sens) + (cmip_rast[[1]]$pet.spstd * pet_sens)
  names(rwi_rast) = "rwi_pred"
  return(rwi_rast)
}



calc_rwi_partials <- function(sppp_code, cmip_id, sensitivity){
  ## Function used to predict species' RWI rasters based on predicted 
  ## sensitivity raster and assigned CMIP model of future climate. Also
  ## integrates calculations of partialling our climate / sensitivity variations
  
  sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, sppp_code, ".gz"))
  
  cmip_rast <- sp_fut_clim %>% 
    filter(cmip_idx == cmip_id) %>% 
    pull(clim_future_sp)
  
  cmip_rast = cmip_rast[[1]]
  
  cwd_sens = sensitivity %>% subset("cwd_sens")
  pet_sens = sensitivity %>% subset("pet_sens")
  intercept = sensitivity %>% subset("intercept")
  rwi_rast <- intercept + (cmip_rast$cwd.spstd * cwd_sens) + (cmip_rast$pet.spstd * pet_sens)
  names(rwi_rast) = "rwi_pred"
  
  # mean_fut_cwd <- cmip_rast %>% subset("cwd.spstd") %>% cellStats(stat = "mean")
  # mean_fut_pet <- cmip_rast %>% subset("pet.spstd") %>% cellStats(stat = "mean")
  # rwi_psens <- intercept + (mean_fut_cwd * cwd_sens) + (mean_fut_pet * pet_sens)
  # names(rwi_psens) = "rwi_psens"
  
  
  mean_cwd_sens <- sensitivity %>% subset("cwd_sens") %>% cellStats(stat = "mean")
  mean_pet_sens <- sensitivity %>% subset("pet_sens") %>% cellStats(stat = "mean")
  mean_intercept <-sensitivity %>% subset("intercept") %>% cellStats(stat = "mean") 
  rwi_pclim <- mean_intercept + (cmip_rast$cwd.spstd * mean_cwd_sens) + (cmip_rast$pet.spstd * mean_pet_sens)
  names(rwi_pclim) = "rwi_pclim"
  
  # return(list("rwi_pred" = rwi_rast, "rwi_psens" = rwi_psens, "rwi_pclim" = rwi_pclim))
  
  return(list("rwi_pred" = c(rwi_rast), "rwi_pclim" = c(rwi_pclim)))
}

quantiles <- function(x){
  ## Defines quantiles used to summarize MC runs
  quantile(x, c(0.025, 0.975), na.rm=TRUE)
}

# extract_quantiles <- function(rwi_preds){
#   ## Extracts desired quantiles for MC runs
#   
#   rwi_quantiles <- rwi_preds %>% 
#     pull(rwi_predictions) %>%
#     brick() %>% 
#     calc(quantiles)
#   return(rwi_quantiles)
# }


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compute sensitivity, RWI for each species by MC combination ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# spp_code = "abal"
# mc_data = mc_nests %>%
#   filter(sp_code == spp_code) %>%
#   pull(data)
# mc_data = mc_data[[1]]

pull_layer <- function(brick, layer_name){
  pulled_layer <- brick %>% subset(layer_name)
}


calc_mean_fut_clim <- function(sppp_code){
  sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, sppp_code, ".gz"))
  clim_pulls <- sp_fut_clim %>%
    mutate(cwd.fut = pmap(list(brick = clim_future_sp, layer_name = "cwd.spstd"),
                          .f = pull_layer),
           pet.fut = pmap(list(brick = clim_future_sp, layer_name = "pet.spstd"),
                          .f = pull_layer))
  
  cwd.fut.q <- clim_pulls %>%
    pull(cwd.fut) %>%
    brick()
  cwd.fut.q <- raster::mean(cwd.fut.q)
  names(cwd.fut.q) = "cwd.fut"
  
  pet.fut.q <- clim_pulls %>%
    pull(pet.fut) %>%
    brick()
  pet.fut.q <- raster::mean(pet.fut.q)
  names(pet.fut.q) = "pet.fut"
  
  out_brick <- brick(c(cwd.fut.q, pet.fut.q))
  return(out_brick)
}


calc_rwi_quantiles <- function(spp_code, mc_data){
  tic()
  ## Iterating through each species
  print(spp_code)
  
  ## Calculate n_mc versions of species' sensitivity raster
  sp_sensitivity <- mc_data %>% 
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
                                                              packages = c( "dplyr", "raster"))))
  
  
  # sp_sensitivity <- mc_data %>% 
  #   mutate(sensitivity = future_pmap(list(sppp_code = spp_code,
  #                                         coefs = ss_coefs),
  #                             .f = predict_sens,
  #                             .options = furrr_options(seed = my_seed, 
  #                                                      packages = c( "dplyr", "raster"))))


  ## Predict future RWI for each of n_mc run
  sp_predictions <- sp_sensitivity %>% 
    mutate(rwi_predictions = future_pmap(list(sppp_code = spp_code,
                                              cmip_id = cmip_idx,
                                              sensitivity = sensitivity),
                                         .f = calc_rwi_partials,
                                         .options = furrr_options(seed = my_seed,
                                                                  packages = "raster"))) 
  
  sp_predictions <- sp_predictions %>% 
    unnest_wider(rwi_predictions) %>% 
    mutate(rwi_pred = rwi_pred[[1]],
           rwi_pclim = rwi_pclim[[1]])
   
  # ## For each species, calculate cell-wise quantiles of rwi from n_mc runs
  # rwi_pred_q <- sp_predictions %>%
  #   pull(rwi_pred) %>%
  #   brick() %>%
  #   calc(quantiles)
  # rwi_psens_q <- sp_predictions %>%
  #   pull(rwi_psens) %>%
  #   brick() %>%
  #   calc(quantiles)
  # rwi_pclim_q <- sp_predictions %>%
  #   pull(rwi_pclim) %>%
  #   brick() %>%
  #   calc(quantiles)
  rwi_pred_q <- sp_predictions %>%
    pull(rwi_pred) %>%
    brick()
  # rwi_psens_q <- sp_predictions %>%
  #   pull(rwi_psens) %>%
  #   brick()
  rwi_pclim_q <- sp_predictions %>%
    pull(rwi_pclim) %>%
    brick()
  
  beginCluster(n = n_cores)
  rwi_pred_mean <- raster::mean(rwi_pred_q)
  rwi_pred_q <- clusterR(rwi_pred_q, calc, args = list(fun = quantiles))
  
  # rwi_psens_q <- clusterR(rwi_psens_q, calc, args = list(fun = quantiles))
  rwi_pclim_mean <- raster::mean(rwi_pclim_q)
  rwi_pclim_q <- clusterR(rwi_pclim_q, calc, args = list(fun = quantiles))
  endCluster()

  names(rwi_pred_q) = c("rwi_pred_025", "rwi_pred_975")
  names(rwi_pred_mean) = c("rwi_pred_mean")
  rwi_pred_q = brick(c(rwi_pred_q, rwi_pred_mean))
  # names(rwi_psens_q) = c("rwi_psens_025", "rwi_psens_975")
  names(rwi_pclim_q) = c("rwi_pclim_025", "rwi_pclim_975")
  names(rwi_pclim_mean) = c("rwi_pclim_mean")
  rwi_pclim_q = brick(c(rwi_pclim_q, rwi_pclim_mean))
  
    
  ## Generate rasters summarizing mean estimate of sensitivity
  sens_pulls <- sp_sensitivity %>%
    mutate(cwd = pmap(list(brick = sensitivity, layer_name = "cwd_sens"),
                      .f = pull_layer),
           pet = pmap(list(brick = sensitivity, layer_name = "pet_sens"),
                      .f = pull_layer),
           intercept = pmap(list(brick = sensitivity, layer_name = "intercept"),
                            .f = pull_layer))
  
  sens_cwd_q <- sens_pulls %>%
    pull(cwd) %>%
    brick()
  sens_cwd_q <- raster::mean(sens_cwd_q)
  names(sens_cwd_q) = "cwd_sens"
  
  
  sens_pet_q <- sens_pulls %>%
    pull(pet) %>%
    brick()
  sens_pet_q <- raster::mean(sens_pet_q)
  names(sens_pet_q) = "pet_sens"
  
  sens_int_q <- sens_pulls %>%
    pull(intercept) %>%
    brick()
  sens_int_q <- raster::mean(sens_int_q)
  names(sens_int_q) = "intercept"
  
  
  ## Pull historic climate
  clim_historic_sp <- sp_hist_clim %>% 
    filter(sp_code == spp_code) %>% 
    pull(clim_historic_sp)
  
  ## Pull mean future climate
  clim_fut_sp <- calc_mean_fut_clim(spp_code)
  
  ## Stack rasters and convert to dataframe
  out_df <- brick(c(clim_historic_sp, 
                    clim_fut_sp,
                    sens_cwd_q, sens_pet_q, sens_int_q, 
                    rwi_pred_q, rwi_pclim_q)) %>% 
    as.data.frame(xy = TRUE) %>% 
    mutate(sp_code = spp_code) %>% 
    drop_na()
  
  ## Write out
  out_df %>% 
    write_rds(file = paste0(out_dir, spp_code, ".gz"), compress = "gz")
  
  ## Clear raster temp files from system
  remove(sp_predictions)
  remove(sens_pulls, sp_sensitivity)
  toc()
  removeTmpFiles(0)
  
  return("done")
}



mc_nests <- sp_mc %>%
  group_by(sp_code) %>%
  nest() %>% 
  # filter(sp_code == "abla") %>%
  drop_na()

mc_nests <- mc_nests %>% 
  mutate(predictions = pmap(list(spp_code = sp_code,
                                 mc_data = data),
                              .f = calc_rwi_quantiles)) 
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




