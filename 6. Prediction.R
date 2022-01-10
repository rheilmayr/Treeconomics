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

n_mc <- 10


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Second stage model
mod_df <- readRDS(paste0(wdir, "out\\second_stage\\ss_mc_mods.rds"))

# mod <- readRDS(paste0(wdir, "out\\second_stage\\cwd_mod.rds"))
# pet_mod <- readRDS(paste0(wdir, "out\\second_stage\\pet_mod.rds"))
# int_mod <- readRDS(paste0(wdir, "out\\second_stage\\int_mod.rds"))

# mod <- readRDS(paste0(wdir, "out\\second_stage\\ss_sq_mod.rds"))
# pet_mod <- readRDS(paste0(wdir, "out\\second_stage\\ss_sq_pet_mod.rds"))


# 2. Historic climate raster
clim_file <- paste0(wdir, 'in\\CRUData\\historic_raster\\HistoricCWD_AETGrids.Rdat')
load(clim_file)
cwd_historic <- sum(cwd_historic)
aet_historic <- sum(aet_historic)
pet_historic <- aet_historic + cwd_historic
names(cwd_historic) = "cwd"
names(pet_historic) = "pet"
clim_historic <- raster::brick(list(cwd_historic, pet_historic))

# 3. Climate projections from CMIP5
cmip <- load(paste0(wdir, 'in\\CMIP5 CWD\\cmip5_cwdaet_end.Rdat'))
pet_raster <- aet_raster + cwd_raster
pet_future <- pet_raster
cwd_future <- cwd_raster
rm(pet_raster)
rm(cwd_raster)
rm(aet_raster)

# 4. Species range maps
range_file <- paste0(wdir, 'in\\species_ranges\\merged_ranges.shp')
range_sf <- st_read(range_file)

# 5. Species climate niche
niche <- read.csv(paste0(wdir, "out\\climate\\clim_niche.csv")) %>% 
  select(-X)

# 6. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family) %>% 
  rename(sp_code = species_id)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Organize CMIP models into tibble  ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pull_cmip_model <- function(mod_n){
  cwd_lyr <- cwd_future[[mod_n]]
  names(cwd_lyr) = "cwd"
  pet_lyr <- pet_future[[mod_n]]
  names(pet_lyr) = "pet"
  future_clim <- brick(list(cwd_lyr, pet_lyr))
  return(future_clim)
}

calc_cwd_change <- function(cmip_rast){
  cwd_change <- ((cmip_rast %>% subset("cwd")) - cwd_historic)
  # cwd_change <- (cmip_rast %>% subset("cwd") - cwd_historic) / cwd_historic
  return(cwd_change)
}

n_cmip_mods <- dim(pet_future)[3]
cmip_projections <- tibble(mod_n = 1:n_cmip_mods)
cmip_projections <- cmip_projections %>% 
  mutate(cmip_rast = map(mod_n, pull_cmip_model))

select_cwd <- function(raster_brick){
  cwd <- raster_brick %>% raster::subset("cwd")
  return(cwd)
}

cmip_cwd <- cmip_projections %>% 
  pull(cmip_rast) %>% 
  lapply(select_cwd)

names(cwd_historic) = "cwd.hist"
names(pet_historic) = "pet.hist"
clim_smry <- raster::brick(c(cwd_historic, pet_historic, cmip_cwd))
clim_smry <- clim_smry %>% 
  as.data.frame() %>% 
  drop_na() %>% 
  pivot_longer(c(-cwd.hist, -pet.hist), names_to = "cmip_mod", values_to = "future_cwd") %>% 
  mutate(cwd_change = (future_cwd - cwd.hist),
         cwd_pct_change = 100 * cwd_change / cwd.hist)

# Plot change in cwd
nbins = 26
label_gaps <- 5
label_pattern <- seq(1,nbins,label_gaps)
plot_dat <- clim_smry %>% 
  mutate(cwd.q = as.numeric(ntile(cwd.hist, nbins)),
         pet.q = as.numeric(ntile(pet.hist, nbins)))

cwd.quantiles = quantile(plot_dat$cwd.hist, probs = seq(0, 1, 1/nbins), names = TRUE) %>% 
  round(2) %>% 
  lapply(round)
pet.quantiles = quantile(plot_dat$pet.hist, probs = seq(0, 1, 1/nbins), names = TRUE) %>% 
  round(2) %>% 
  lapply(round)
cwd.breaks = seq(0.5, nbins+0.5, 1)
pet.breaks = seq(0.5, nbins+0.5, 1)

group_dat <- plot_dat %>% 
  filter(cwd.hist>100) %>% 
  group_by(cwd.q, pet.q) %>% 
  dplyr::summarize(cwd_pct_change = mean(cwd_pct_change, na.rm = TRUE),
                   cwd_change = mean(cwd_change, na.rm = TRUE),
                   n = n()) %>% 
  filter(n>1000)

binned_change <- group_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_change)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw(base_size = 22)+
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme(legend.position = "right") +
  labs(fill = "CWD change\nby 2100") +
  scale_x_continuous(labels = cwd.quantiles[label_pattern], breaks = cwd.breaks[label_pattern]) +
  scale_y_continuous(labels = pet.quantiles[label_pattern], breaks = pet.breaks[label_pattern]) +
  ylab(bquote("Historic PET (mmH"[2]*"O)")) +
  xlab(bquote("Historic CWD (mmH"[2]*"O)")) + 
  coord_fixed()

binned_change
ggsave(paste0(wdir, 'figures\\cwd_change.svg'), plot = binned_change)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Rescale CMIP predictions into species standardized climate ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
species_list <- niche %>% 
  select(sp_code) %>% 
  left_join(sp_info, by = "sp_code")

rasterize_spstd <- function(spp_code, clim_raster){
  sp_niche <- niche %>%
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

select_lyr <- function(brick, layer_name){
  lyr <- brick %>% subset(layer_name)
  return(lyr)
}

cmip_projections <- cmip_projections %>% 
  mutate(cwd_rast = map(.x = cmip_rast, subset = "cwd", .f = subset),
         pet_rast = map(.x = cmip_rast, subset = "pet", .f = subset))

cmip_n <- dim(cmip_projections)[1]
cmip_assignments <- tibble(iter_n = seq(1, n_mc)) %>% 
  mutate(cmip_idx = sample(seq(cmip_n), n_mc))

cwd_proj_mean <- brick(cmip_projections$cwd_rast) %>%  
  calc(mean)

pet_proj_mean <- brick(cmip_projections$pet_rast) %>% 
  calc(mean)
  
clim_future <- brick(c(cwd_proj_mean, pet_proj_mean))  # TODO: In future, should incorporate variation across CMIP models rather than just mean. MC based off std in models? Bootstrap individual CMIP models?
names(clim_future) = c("cwd", "pet")
# clim_future <- (cmip_projections[1,2] %>% pull())[[1]]  

sp_predictions <- species_list %>% 
  mutate(clim_future_sp = map(.x = sp_code, clim_raster = clim_future, .f = rasterize_spstd))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate deviation from species' historic range ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_predictions <- sp_predictions %>% 
  mutate(clim_historic_sp = map(.x = sp_code, clim_raster = clim_historic, .f = rasterize_spstd))

# sp_niche <- niche %>%
#   drop_na() %>% 
#   filter(sp_code == spp_code) %>% 
#   select(-sp_code) 
# 
# sp_range <- range_sf %>% 
#   filter(sp_code == spp_code)
# pet_historic_sp <- pet_historic %>%  
#   mask(sp_range)
# pet_historic_spstd <- (pet_historic_sp - sp_niche$pet_mean) / sp_niche$pet_sd
# names(pet_historic_spstd) = "pet.spstd"
# 
# cwd_historic_sp <- cwd_historic %>% 
#   mask(sp_range)
# cwd_historic_spstd <- (cwd_historic_sp - sp_niche$cwd_mean) / sp_niche$cwd_sd
# names(cwd_historic_spstd) = "cwd.spstd"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create climate sensitivity rasters ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# clim_historic_sp <- brick(list(cwd_historic_spstd, pet_historic_spstd))
# pet_sens <- raster::predict(clim_historic_sp, pet_mod)
predict_sens <- function(clim_historic_sp){
  cwd_sens <- raster::predict(clim_historic_sp, mod)
  names(cwd_sens) = "cwd_sens"
  pet_sens <- raster::predict(clim_historic_sp, pet_mod)
  names(pet_sens) = "pet_sens"
  intercept <- raster::predict(clim_historic_sp, int_mod)
  names(intercept) = "intercept"
  sensitivity <- raster::brick(cwd_sens, pet_sens, intercept) 
  return(sensitivity)
}

sp_predictions <- sp_predictions %>% 
  mutate(sensitivity = map(clim_historic_sp, predict_sens))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Predict growth under future climate ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# sp_clip <- function(rast){
#   rast_clip <- mask(rast, sp_range)
#   return(rast_clip)
# }

calc_rwi <- function(cmip_rast, sensitivity){
  cwd_sens = sensitivity %>% subset("cwd_sens")
  pet_sens = sensitivity %>% subset("pet_sens")
  intercept = sensitivity %>% subset("intercept")
  rwi_rast <- intercept + (cmip_rast$cwd.spstd * cwd_sens) + (cmip_rast$pet.spstd * pet_sens)
  names(rwi_rast) = "rwi_pred"
  return(rwi_rast)
}

sp_predictions <- sp_predictions %>% 
  mutate(rwi_predictions = map2(.x = clim_future_sp, .y = sensitivity, calc_rwi))

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

sp_predictions <- sp_predictions %>% 
  group_by(sp_code) %>% 
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























# rwi_cmip_predict <- function(sensitivity){
#   sp_predictions <- cmip_projections %>% 
#     mutate(rwi_rast = map(cmip_rast, calc_rwi, sensitivity = sensitivity)) %>% 
#     pull(rwi_rast)
#   return(sp_predictions)
# }
# 
# sp_predictions <- sp_predictions %>% 
#   mutate(rwi_predictions = map(sensitivity, rwi_cmip_predict))
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Extract to dataframe    ---------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extract_rasters <- function(raster_column){
#   table_column <- raster_column %>% 
#     as.data.frame(xy = T) %>% 
#     drop_na()
#   return(table_column)
# }
# 
# # extract_predictions <- function(clim_historic_sp, rwi_predictions){
# #   predict_rasters <- raster::brick(c(clim_historic_sp, rwi_predictions))
# #   predict_df <- predict_rasters %>% 
# #     as.data.frame() %>% 
# #     drop_na()
# #   return(predict_df)
# # }
# 
# 
# 
# sp_predictions <- sp_predictions %>% 
#   mutate(clim_historic_sp = map(clim_historic_sp, extract_rasters),
#          rwi_predictions = map(rwi_predictions, extract_rasters),
#          rwi_predictions_partial_sens = map(rwi_predictions_partial_sens, extract_rasters),
#          rwi_predictions_partial_clim = map(rwi_predictions_partial_clim, extract_rasters))
# 
# historic_df <- sp_predictions %>% 
#   select(sp_code, clim_historic_sp) %>% 
#   unnest(clim_historic_sp)
# 
# rwi_df <- sp_predictions %>% 
#   select(sp_code, rwi_predictions) %>% 
#   unnest(rwi_predictions)
# 
# rwi_psens_df <- sp_predictions %>% 
#   select(sp_code, rwi_predictions_partial_sens) %>% 
#   unnest(rwi_predictions_partial_sens)
# 
# rwi_pclim_df <- sp_predictions %>% 
#   select(sp_code, rwi_predictions_partial_clim) %>% 
#   unnest(rwi_predictions_partial_clim)
# 
# out_predictions <- historic_df %>% 
#   left_join(rwi_df, by = c("sp_code", "x", "y")) %>% 
#   left_join(rwi_psens_df, by = c("sp_code", "x", "y")) %>% 
#   left_join(rwi_pclim_df, by = c("sp_code", "x", "y"))
# 
# 
# 
# out_predictions %>% 
#   saveRDS(file = paste0(wdir,"out/predictions/sp_predictions.rds") )
# 
# 
# 
# 
# 
# extract_predictions <- function(clim_historic_sp, rwi_predictions, rwi_predictions_partial_sens, rwi_predictions_partial_clim){
#   predict_rasters <- raster::brick(c(clim_historic_sp, rwi_predictions, rwi_predictions_partial_sens, rwi_predictions_partial_clim))
#   predict_df <- predict_rasters %>% 
#     as.data.frame() %>% 
#     drop_na()
#   return(predict_df)
# }
# 
# sp_predictions <- sp_predictions %>% 
#   mutate(predict_df = pmap(list(clim_historic_sp, rwi_predictions, rwi_predictions_partial_sens, rwi_predictions_partial_clim), extract_predictions))
# 
# 
# sp_predictions_df <- sp_predictions %>% 
#   select(sp_code, predict_df) %>% 
#   unnest(predict_df) %>% 
#   pivot_longer(c(-sp_code, -cwd.spstd, -pet.spstd), names_to = "cmip_run", values_to = "rwi")
# 
# 
# sp_predictions_df %>% 
#   filter(sp_code == "pisy") %>% 
#   ggplot(aes(x = cwd.spstd, y = rwi)) +
#   geom_point()
# 
# sp_predictions_df %>% 
#   filter(sp_code == "pila") %>% 
#   ggplot(aes(x = cwd.spstd, y = rwi)) +
#   geom_point()
# 
# 
# ### Binned plot of cwd sensitivity
# nbins = 21
# label_gaps <- 5
# label_pattern <- seq(1,nbins,label_gaps)
# plot_dat <- sp_predictions %>% 
#   # filter(sp_code == "pila") %>% 
#   mutate(cwd.q = as.numeric(ntile(cwd.spstd, nbins)),
#          pet.q = as.numeric(ntile(pet.spstd, nbins)))
# 
# cwd.quantiles = quantile(plot_dat$cwd.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% 
#   round(2) %>% 
#   lapply(round, digits = 1)
# pet.quantiles = quantile(plot_dat$pet.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% 
#   round(2) %>% 
#   lapply(round, digits = 1)
# cwd.breaks = seq(0.5, nbins+0.5, 1)
# pet.breaks = seq(0.5, nbins+0.5, 1)
# 
# group_dat <- plot_dat %>% 
#   group_by(cwd.q, pet.q) %>% 
#   dplyr::summarize(ln_rwi = mean(ln_rwi, na.rm = TRUE),
#                    n = n()) %>% 
#   filter(n>100)
# 
# 
# binned_margins <- group_dat %>% 
#   ggplot(aes(x = cwd.q, y = pet.q, fill = ln_rwi)) +
#   geom_tile() +
#   scale_fill_viridis_c(direction = -1) +
#   theme_bw(base_size = 22)+
#   theme(legend.position = "right") +
#   labs(fill = "Predicted impact\nof climate change\non RWI") +
#   scale_x_continuous(labels = cwd.quantiles[label_pattern], breaks = cwd.breaks[label_pattern]) +
#   scale_y_continuous(labels = pet.quantiles[label_pattern], breaks = pet.breaks[label_pattern]) +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   coord_fixed()
# 
# binned_margins 
# combined_plot <- binned_change / binned_margins
# 
# ggsave(paste0(wdir, 'figures\\predicted_rwi.svg'), plot = binned_margins)
# ggsave(paste0(wdir, 'figures\\climate_change.svg'), plot = combined_plot, width = 7, height = 11, units = "in")
# 
# 
# 
# 
# # calc_mean_rwi <- function(rwi_predictions){
# #   mean_rwi <- raster::brick(rwi_predictions %>% pull(rwi_rast)) %>% mean()
# #   return(mean_rwi)
# # }
# # 
# # sp_predictions <- sp_predictions %>%
# #   mutate(rwi_2100 = map(rwi_predictions, calc_mean_rwi))
# 
# 
# ### Plot drought change by historic cwd and pet
# cwd_historic_sp <- clim_historic_sp %>% subset("cwd.spstd")
# cwd_change <- raster::brick(c(clim_historic_sp, projections$cwd_change))
# cwd_change <- as.data.frame(cwd_change)
# cwd_change <- cwd_change %>% 
#   pivot_longer(-c(cwd.spstd, pet.spstd), names_to = "cmip_mod", values_to = "cwd_change")
# 
# plot_dat <- cwd_change %>% drop_na()
# nbins = 10
# plot_dat <- plot_dat %>% 
#   mutate(cwd.q = as.numeric(ntile(cwd.spstd, nbins)),
#          pet.q = as.numeric(ntile(pet.spstd, nbins)))
# cwd.quantiles = quantile(plot_dat$cwd.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)
# pet.quantiles = quantile(plot_dat$pet.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)
# 
# group_dat <- plot_dat %>% 
#   group_by(cwd.q, pet.q) %>% 
#   summarize(mean = mean(cwd_change),
#             n = n())
# 
# binned_margins <- group_dat %>% 
#   ggplot(aes(x = cwd.q, y = pet.q, fill = mean)) +
#   geom_tile() +
#   # xlim(c(-3, 4)) +
#   #ylim(c(-1.5, 1.5))+
#   # scale_fill_gradientn (colours = c("darkblue","lightblue")) +
#   scale_fill_viridis_c() +
#   # scale_fill_viridis_c() +
#   theme_bw(base_size = 22)+
#   ylab("Deviation from mean PET")+
#   xlab("Deviation from mean CWD")+
#   theme(legend.position = "right") +
#   labs(fill = "Mean predicted\nchange in CWD") +
#   scale_x_continuous(labels = cwd.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   scale_y_continuous(labels = pet.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   coord_fixed()
# 
# binned_margins 
# 
# ### Plot predicted RWI impact by historic cwd and pet
# predict_brick <- raster::brick(c(rwi_change = mean_rwi, cwd_hist = cwd_historic, pet_hist = pet_historic))
# rwi_change_df <- as.data.frame(predict_brick)
# rwi_change_df <- rwi_change_df %>% 
#   drop_na() %>% 
#   mutate(cwd.spstd = (cwd_hist - sp_niche$cwd_mean) / sp_niche$cwd_sd,
#          pet.spstd = (pet_hist - sp_niche$pet_mean) / sp_niche$pet_sd)
# mod <- lm(rwi_change ~ cwd.spstd + pet.spstd, data = rwi_change_df)
# summary(mod)
# 
# xmin = -2.5
# xmax = 2.5
# predictions <- prediction(mod, at = list(cwd.spstd = seq(xmin, xmax, .1)), calculate_se = T) %>% 
#   summary() %>% 
#   rename(cwd.spstd = "at(cwd.spstd)")
# 
# margins_plot <- ggplot(predictions, aes(x = cwd.spstd)) + 
#   geom_line(aes(y = Prediction)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Mean predicted impact on ln(RWI)") + 
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# 
# margins_plot
# 
# 
# plot_dat <- rwi_change_df %>%
#   filter(((abs(cwd.spstd)<3) & (abs(pet.spstd<3)))) %>%
#   drop_na()
# nbins = 10
# plot_dat <- plot_dat %>% 
#   mutate(cwd.q = as.numeric(ntile(cwd.spstd, nbins)),
#          pet.q = as.numeric(ntile(pet.spstd, nbins)))
# 
# cwd.quantiles = quantile(plot_dat$cwd.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)
# pet.quantiles = quantile(plot_dat$pet.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)
# 
# group_dat <- plot_dat %>% 
#   group_by(cwd.q, pet.q) %>% 
#   summarize(mean = mean(rwi_change),
#             n = n()) %>% 
#   filter(n>1)
# 
# binned_margins <- group_dat %>% 
#   ggplot(aes(x = cwd.q, y = pet.q, fill = mean)) +
#   geom_tile() +
#   # xlim(c(-3, 4)) +
#   #ylim(c(-1.5, 1.5))+
#   scale_fill_gradientn (colours = c("darkblue","lightblue")) +
#   # scale_fill_viridis_c() +
#   theme_bw(base_size = 22)+
#   ylab("Deviation from mean PET")+
#   xlab("Deviation from mean CWD")+
#   theme(legend.position = "right") +
#   labs(fill = "Marginal effect\nof CWD") +
#   scale_x_continuous(labels = cwd.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   scale_y_continuous(labels = pet.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   coord_fixed()
# 
# binned_margins 
# 
# 
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Explore species-level climate / sensitivity ---------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Plot cwd deviation against historic cwd
# # Example species: East Coast - tsca; West coast - pipo; Mediterranean -  
# example_species <- "pila" # TODO: Dig into weird CWD censoring for ABAL, TSHE - guessing this is due to 0 CWD values
# sp_fut <- (sp_predictions %>% 
#              filter(sp_code == example_species) %>% 
#              pull(clim_future_sp))[[1]]
# sp_hist <- (sp_predictions %>% 
#               filter(sp_code == example_species) %>% 
#               pull(clim_historic_sp))[[1]]
# sp_sens  <- (sp_predictions %>% 
#                filter(sp_code == example_species) %>% 
#                pull(sensitivity))[[1]]
# 
# names(sp_fut) <- c("cwd.fut", "pet.fut")
# clim_compare <- brick(c(sp_fut, sp_hist, sp_sens))
# clim_compare <- clim_compare %>% 
#   as.data.frame(xy = TRUE)
# 
# clim_compare %>% 
#   ggplot(aes(x = cwd.spstd, y = cwd_sens)) +
#   geom_point() +
#   xlim(c(-2,2)) +
#   ylim(c(-.3,.3)) +
#   theme_bw()
# 
# 
# clim_compare %>% 
#   ggplot(aes(x = cwd.spstd, y = pet_sens)) +
#   geom_point() +
#   xlim(c(-2,2)) +
#   ylim(c(-.3,.3)) +
#   theme_bw()
# 
# clim_compare %>% 
#   ggplot(aes(x = cwd.spstd, y = pet.spstd)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   coord_fixed() +
#   xlim(c(-2,2)) +
#   ylim(c(-2,2)) +
#   theme_bw()
# 
# 
# clim_compare %>% 
#   ggplot(aes(x = cwd.spstd, y = cwd.fut)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   coord_fixed() +
#   xlim(c(-4,4)) +
#   ylim(c(-4,4)) +
#   theme_bw()
# 
# clim_compare %>% 
#   ggplot(aes(x = pet.spstd, y = pet.fut)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   coord_fixed() +
#   xlim(c(-2,2)) +
#   ylim(c(-2,2)) +
#   theme_bw()
# 
# sens  <- (sp_predictions %>% filter(sp_code == example_species) %>% pull(sensitivity))[[1]]
# sens <- sens %>% subset("pet_sens")
# data <- as.data.frame(sens,xy = TRUE)
# data <- data %>% drop_na()
# xlims <- c(round(min(data$x) - 1), round(max(data$x) + 1))
# ylims <- c(round(min(data$y) - 1), round(max(data$y) + 1))
# 
# # Plot species ranges
# world <- ne_coastline(scale = "medium", returnclass = "sf")
# map <- ggplot() +
#   geom_sf(data = world) +
#   geom_raster(data = data, aes(x = x, y = y, fill = pet_sens)) +
#   theme_bw(base_size = 22)+
#   ylab("Latitude")+
#   xlab("Longitude")+
#   scale_fill_viridis_c(direction = -1) +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)
# map
# 
# 
# sens  <- (sp_predictions %>% filter(sp_code == "tsca") %>% pull(sensitivity))[[1]]
# sens <- sens %>% subset("cwd_sens")
# data <- as.data.frame(sens,xy = TRUE)
# data <- data %>% drop_na()
# 
# # Plot species ranges
# world <- ne_coastline(scale = "medium", returnclass = "sf")
# map <- ggplot() +
#   geom_sf(data = world) +
#   geom_raster(data = data, aes(x = x, y = y, fill = cwd_sens)) +
#   theme_bw(base_size = 22)+
#   ylab("Latitude")+
#   xlab("Longitude")+
#   scale_fill_viridis_c(direction = -1) +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)
# map
# 
# 
# 
# sens  <- (sp_predictions %>% filter(sp_code == "tsca") %>% pull(clim_historic_sp))[[1]]
# sens <- sens %>% subset("cwd.spstd")
# data <- as.data.frame(sens,xy = TRUE)
# data <- data %>% drop_na()
# 
# # Plot species ranges
# map <- ggplot() +
#   geom_sf(data = world) +
#   geom_raster(data = data, aes(x = x, y = y, fill = cwd.spstd)) +
#   theme_bw(base_size = 22)+
#   ylab("Latitude")+
#   xlab("Longitude")+
#   scale_fill_viridis_c(direction = -1) +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)
# map
# 
# 
# sens  <- (sp_predictions %>% filter(sp_code == "tsca") %>% pull(clim_historic_sp))[[1]]
# sens <- sens %>% subset("pet.spstd")
# data <- as.data.frame(sens,xy = TRUE)
# data <- data %>% drop_na()
# 
# # Plot species ranges
# world <- ne_coastline(scale = "medium", returnclass = "sf")
# map <- ggplot() +
#   geom_sf(data = world) +
#   geom_raster(data = data, aes(x = x, y = y, fill = pet.spstd)) +
#   theme_bw(base_size = 22)+
#   ylab("Latitude")+
#   xlab("Longitude")+
#   scale_fill_viridis_c(direction = -1) +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)
# map
# 
# 
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Map historic climate and sensitivity ---------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# xlims <- c(-132, -100)
# ylims <- c(30, 56)
# cwd_sensitivity_df <- cwd_sens %>% 
#   as.data.frame(xy = TRUE) %>% 
#   drop_na() %>% 
#   filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
#   rename(sensitivity = layer)
# 
# cwd_historic_df <- cwd_historic_spstd %>% 
#   as.data.frame(xy = TRUE) %>% 
#   filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>%
#   drop_na()
# 
# pet_historic_df <- pet_historic_spstd %>% 
#   as.data.frame(xy = TRUE) %>% 
#   filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
#   drop_na()
# 
# 
# 
# world <- ne_coastline(scale = "medium", returnclass = "sf")
# sens_map <- ggplot() +
#   geom_raster(data = cwd_sensitivity_df, aes(x = x, y = y, fill = sensitivity)) +
#   # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
#   scale_fill_viridis_c(direction = -1,
#                        guide = guide_legend(
#                          label.theme = element_text(angle = 90),
#                          label.position = "bottom"
#                        )) +
#   theme_bw(base_size = 22)+
#   theme(legend.position = "bottom") +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   geom_sf(data = world, color = 'black') +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US
# 
# pet_map <- ggplot() +
#   geom_raster(data = pet_historic_df, aes(x = x, y = y, fill = pet.spstd)) +
#   # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
#   scale_fill_viridis_c(guide = guide_legend(
#     label.theme = element_text(angle = 90),
#     label.position = "bottom"
#   )) +
#   theme_bw(base_size = 22)+
#   theme(legend.position = "bottom") +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   geom_sf(data = world, color = 'black') +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US
# 
# cwd_map <- ggplot() +
#   geom_raster(data = cwd_historic_df, aes(x = x, y = y, fill = cwd.spstd)) +
#   # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
#   scale_fill_viridis_c(guide = guide_legend(
#     label.theme = element_text(angle = 90),
#     label.position = "bottom"
#   )) +
#   theme_bw(base_size = 22) +
#   theme(legend.position = "bottom") +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   geom_sf(data = world, color = 'black') +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US
# 
# sens_map | cwd_map
# 
# 
# rwi_rast <- mean_rwi
# 
# # rwi_rast <- projections %>% 
# #   filter(mod_n == 10) %>% 
# #   pull(rwi_rast)
# 
# 
# 
# rwi_df <- rwi_rast[[1]] %>% 
#   as.data.frame(xy = TRUE) %>% 
#   drop_na() %>% 
#   filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
#   rename(log_rwi = layer)
# 
# rwi_map <- ggplot() +
#   geom_raster(data = rwi_df, aes(x = x, y = y, fill = log_rwi)) +
#   # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
#   scale_fill_viridis_c(direction = -1,
#                        guide = guide_legend(
#                          label.theme = element_text(angle = 90),
#                          label.position = "bottom"
#                        )) +
#   theme_bw(base_size = 22)+
#   theme(legend.position = "bottom") +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   geom_sf(data = world, color = 'black') +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US
# 
# rwi_map
# 
# 
# 
# 
# 
# cwd_prj <- projections %>% 
#   filter(mod_n==1) %>% 
#   pull(cmip_rast)
# 
# cwd_change <- cwd_prj[[1]]$cwd - cwd_historic
# 
# change_df <- cwd_change %>% 
#   as.data.frame(xy = TRUE) %>% 
#   drop_na() %>% 
#   filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
#   rename(change = layer)
# 
# change_map <- ggplot() +
#   geom_raster(data = change_df, aes(x = x, y = y, fill = change)) +
#   # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
#   scale_fill_viridis_c(guide = guide_legend(
#                          label.theme = element_text(angle = 90),
#                          label.position = "bottom"
#                        )) +
#   theme_bw(base_size = 22)+
#   theme(legend.position = "bottom") +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   geom_sf(data = world, color = 'black') +
#   coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US
# 
# change_map
