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
future::plan(multisession, workers = 4)

my_seed <- 5597
set.seed(my_seed)

n_mc <- 5
n_spp <- 4

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Second stage model
mod_df <- readRDS(paste0(wdir, "out\\second_stage\\ss_mc_mods.rds"))

# 2. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family) %>% 
  rename(sp_code = species_id)


# 3. Species-standardized historic climate
sp_hist_clim <- readRDS(paste0(wdir, "out//climate//sp_clim_historic.rds"))
species_list <- sp_hist_clim %>% select(sp_code)

# 4. Directory of species-standardized future possible climates
sp_fut_clim_dir <- paste0(wdir, "out\\climate\\sp_clim_predictions\\")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Assign MC coefs and CMIP models  ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Join second stage coefficients to species list
sp_mc <- species_list %>% 
  select(sp_code) %>% 
  crossing(iter_idx = seq(n_mc))

## Assign specific cmip realization to each MC iteration 
sp_fut_clim_smpl <- readRDS(paste0(sp_fut_clim_dir, "abal"))
n_cmip_mods <- sp_fut_clim_smpl %>% pull(cmip_idx) %>% unique() %>% length()
cmip_assignments <- tibble(iter_idx = seq(1, n_mc)) %>% 
  mutate(cmip_idx = sample(seq(n_cmip_mods), n_mc))

## Join cmip model assignments
sp_mc <- sp_mc %>% 
  left_join(cmip_assignments, by = "iter_idx")

## Join second stage coefficients
mod_df <- mod_df %>% 
  filter(iter_idx %in% seq(n_mc)) %>% 
  group_by(iter_idx) %>% 
  nest() %>% 
  rename(ss_coefs = data)

sp_mc <- sp_mc %>% 
  left_join(mod_df, by = "iter_idx") 

# sp_mc <- sp_mc %>% 
#   group_by(sp_code) %>% 
#   nest()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create climate sensitivity rasters ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predict_sens <- function(spp_code, coefs){
  select <- dplyr::select
  
  clim_historic_sp <- sp_hist_clim %>% 
    filter(sp_code == spp_code) %>% 
    pull(clim_historic_sp)
  cwd_rast <- clim_historic_sp[[1]] %>% subset("cwd.spstd")
  pet_rast <- clim_historic_sp[[1]] %>% subset("pet.spstd")
  
  # coefs <- (coefs %>% deframe())[[1]]
  cwd_coefs <- coefs %>% select(parameter, cwd_coef) %>% deframe()
  pet_coefs <- coefs %>% select(parameter, pet_coef) %>% deframe()
  int_coefs <- coefs %>% select(parameter, int_coef) %>% deframe()
  
  cwd_sens <- cwd_coefs[["(Intercept)"]] +
              cwd_coefs[["cwd.spstd"]] * cwd_rast +
              cwd_coefs[["pet.spstd"]] * pet_rast
  names(cwd_sens) = "cwd_sens"
  
  pet_sens <- pet_coefs[["(Intercept)"]] +
              pet_coefs[["cwd.spstd"]] * cwd_rast +
              pet_coefs[["pet.spstd"]] * pet_rast
  names(pet_sens) = "pet_sens"
  
  intercept <- int_coefs[["(Intercept)"]] +
               int_coefs[["cwd.spstd"]] * cwd_rast +
               int_coefs[["pet.spstd"]] * pet_rast
  names(intercept) = "intercept"
  
  sensitivity <- raster::brick(cwd_sens, pet_sens, intercept) 
  return(sensitivity)
}



## Calculate species by n_mc versions of sensitivity rasters
sp_sensitivity <- sp_mc[1:20,] %>% 
  mutate(sensitivity = future_pmap(list(spp_code = sp_code,
                                        coefs = ss_coefs),
                            .f = predict_sens,
                            .options = furrr_options(seed = 5597, 
                                                     packages = c( "dplyr", "raster"))))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Predict growth under future climate ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calc_rwi <- function(spp_code, cmip_id, sensitivity){
  sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, spp_code))
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


## Predict future RWI
sp_predictions <- sp_sensitivity %>% 
  mutate(rwi_predictions = future_pmap(list(spp_code = sp_code,
                                            cmip_id = cmip_idx,
                                            sensitivity = sensitivity),
                                       .f = calc_rwi,
                                       .options = furrr_options(seed = my_seed,
                                                                packages = "raster")))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize predictions using cell-wise quantiles   ----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quantiles <- function(x){
  quantile(x, c(0.025, 0.5, .975), na.rm=TRUE)
}

extract_quantiles <- function(rwi_preds){
  rwi_quantiles <- rwi_preds %>% 
    pull(rwi_predictions) %>%
    brick() %>% 
    calc(quantiles)
  return(rwi_quantiles)
}

## For each species, calculate cell-wise quantiles from n_mc runs
rwi_quantiles <- sp_predictions %>% 
  select(sp_code, iter_idx, rwi_predictions) %>% 
  group_by(sp_code) %>% 
  nest() %>% 
  mutate(rwi_quantiles = future_map(data, extract_quantiles,
                                    .options = furrr_options(seed = my_seed,
                                                             packages = "raster"))) %>% 
  select(-data)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Thought experiments - partialling out mechanisms    ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calc_rwi_partial_sens <- function(cmip_rast, sensitivity){ # NOTE: Should the means be calculated across full range rather than by each species?
  mean_fut_cwd <- cmip_rast %>% subset("cwd.spstd") %>% cellStats(stat = "mean")
  mean_fut_pet <- cmip_rast %>% subset("pet.spstd") %>% cellStats(stat = "mean")
  cwd_sens = sensitivity %>% subset("cwd_sens")
  pet_sens = sensitivity %>% subset("pet_sens")
  intercept = sensitivity %>% subset("intercept")
  rwi_rast <- intercept + (mean_fut_cwd * cwd_sens) + (mean_fut_pet * pet_sens)
  names(rwi_rast) = "rwi_psens"
  return(rwi_rast)
}

calc_rwi_partial_clim <- function(cmip_rast, sensitivity){ # NOTE: Should the means be calculated across full range rather than by each species?
  mean_cwd_sens <- sensitivity %>% subset("cwd_sens") %>% cellStats(stat = "mean")
  mean_pet_sens <- sensitivity %>% subset("pet_sens") %>% cellStats(stat = "mean")
  mean_intercept <-sensitivity %>% subset("intercept") %>% cellStats(stat = "mean") 
  rwi_rast <- mean_intercept + (cmip_rast$cwd.spstd * mean_cwd_sens) + (cmip_rast$pet.spstd * mean_pet_sens)
  names(rwi_rast) = "rwi_pclim"
  return(rwi_rast)
}


sp_predictions <- sp_predictions %>% 
  mutate(rwi_predictions_partial_sens = map2(.x = clim_future_sp, .y = sensitivity, calc_rwi_partial_sens),
         rwi_predictions_partial_clim = map2(.x = clim_future_sp, .y = sensitivity, calc_rwi_partial_clim))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stack rasters into dataframe ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
create_prediction_df <- function(spp_predictions){
  sp_fut <- (spp_predictions %>% 
               pull(clim_future_sp))[[1]]
  names(sp_fut) <- c("cwd.fut", "pet.fut")
  
  sp_hist <- (spp_predictions %>% 
                pull(clim_historic_sp))[[1]]
  sp_sens  <- (spp_predictions %>% 
                 pull(sensitivity))[[1]]
  sp_rwi  <- (spp_predictions %>% 
                pull(rwi_predictions))[[1]]
  
  sp_rwi_psens  <- (spp_predictions %>% 
                      pull(rwi_predictions_partial_sens))[[1]]
  
  sp_rwi_pclim  <- (spp_predictions %>% 
                      pull(rwi_predictions_partial_clim))[[1]]
  
  clim_compare <- brick(c(sp_fut, sp_hist, sp_sens, sp_rwi, sp_rwi_psens, sp_rwi_pclim))
  clim_compare <- clim_compare %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  return(clim_compare)
}

sp_prediction_test <- sp_predictions %>% 
  group_by(sp_code, iter_idx) %>% 
  nest() %>% 
  mutate(pred_df = map(data, create_prediction_df)) %>% 
  select(-data) %>% 
  unnest(cols = pred_df) %>% 
  mutate(cwd_change = cwd.fut - cwd.spstd,
         pet_change = pet.fut - pet.spstd,
         rwi_null = cwd.spstd * cwd_sens + pet.spstd * pet_sens + intercept,
         rwi_change = rwi_pred - rwi_null,
         rwi_change_psens = rwi_psens - rwi_null,
         rwi_change_pclim = rwi_pclim - rwi_null)


sp_predictions %>% 
  saveRDS(file = paste0(wdir,"out/predictions/sp_predictions.rds") )




