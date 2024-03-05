#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(MASS)
library(tidyverse)
library(broom)
library(purrr)
library(margins)
library(tidylog)
library(fixest)
library(gstat)
library(sf)
library(units)
library(dtplyr)
library(marginaleffects)

set.seed(5597)

select <- dplyr::select

n_mc <- 10000


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_std.csv'))
flm_df_tc <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_tc_std.csv'))
flm_df_spei <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_spei_std.csv'))
flm_df_ppt <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_ppt_std.csv'))
flm_df_tc <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_tc_std.csv'))


flm_df <- flm_df_tc %>%
  select(collection_id, estimate_pet.an.spstd.tc, std.error_pet.an.spstd.tc, estimate_cwd.an.spstd.tc, std.error_cwd.an.spstd.tc) %>% 
  inner_join(flm_df, by = "collection_id")

flm_df <- flm_df_spei %>% 
  select(collection_id, estimate_pet.an.spstd.spei, std.error_pet.an.spstd.spei, estimate_cwd.an.spstd.spei, std.error_cwd.an.spstd.spei) %>% 
  inner_join(flm_df, by = "collection_id")

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz"))
flm_df <- flm_df %>% 
  left_join(ave_site_clim, by = c("collection_id"))

# 3. Site information
site_df <- read_csv(paste0(wdir, '1_input_processed/dendro/site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# # 4. Species information
# sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
# sp_info <- sp_info %>% 
#   select(species_id, genus, gymno_angio, family)
# site_df <- site_df %>% 
#   left_join(sp_info, by = "species_id")

# Merge back into main flm_df
flm_df <- flm_df %>% 
  left_join(site_df, by = "collection_id")






# Add weighting based on inverse of first stage variance
flm_df <- flm_df %>% 
  mutate(cwd_errorweights = 1 / (std.error_cwd.an),
         cwd_errorweights_tc = 1 / (std.error_cwd.an.spstd.tc),
         cwd_errorweights_spei = 1 / (std.error_cwd.an.spstd.spei),
         errorweights2 = sqrt(ntrees),
         pet_errorweights = 1 / (std.error_pet.an.spstd.tc),
         int_errorweights = 1 / (std.error_intercept))

# Identify and trim extreme outliers
cwd_est_bounds.spei = quantile(flm_df$estimate_cwd.an.spstd.spei, c(0.01, 0.99),na.rm=T)
pet_est_bounds.spei = quantile(flm_df$estimate_pet.an.spstd.spei, c(0.01, 0.99),na.rm=T)
cwd_est_bounds.tc = quantile(flm_df$estimate_cwd.an.spstd.tc, c(0.01, 0.99),na.rm=T)
pet_est_bounds.tc = quantile(flm_df$estimate_pet.an.spstd.tc, c(0.01, 0.99),na.rm=T)
cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds = quantile(flm_df$estimate_pet.an, c(0.01, 0.99),na.rm=T)
cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)
cwd_spstd_bounds_tc = quantile(flm_df$cwd.spstd.tc, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds_tc = quantile(flm_df$pet.spstd.tc, c(0.01, 0.99), na.rm = T)
cwd_spstd_bounds_spei = quantile(flm_df$cwd.spstd.spei, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds_spei = quantile(flm_df$pet.spstd.spei, c(0.01, 0.99), na.rm = T)

# flm_df <- flm_df %>%
#   mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
#            (estimate_cwd.an>cwd_est_bounds[2]) |
#            (estimate_pet.an<pet_est_bounds[1]) |
#            (estimate_pet.an>pet_est_bounds[2]) |
#            (cwd.spstd<cwd_spstd_bounds[1]) |
#            (cwd.spstd>cwd_spstd_bounds[2]) |
#            (pet.spstd<pet_spstd_bounds[1]) |
#            (pet.spstd>pet_spstd_bounds[2]))


flm_df <- flm_df %>%
  mutate(outlier = 
           (estimate_cwd.an.spstd.spei<cwd_est_bounds.spei[1]) |
           (estimate_cwd.an.spstd.spei>cwd_est_bounds.spei[2]) |
           (estimate_pet.an.spstd.spei<pet_est_bounds.spei[1]) |
           (estimate_pet.an.spstd.spei>pet_est_bounds.spei[2]) |
           (cwd.spstd.spei<cwd_spstd_bounds_spei[1]) |
           (cwd.spstd.spei>cwd_spstd_bounds_spei[2]) |
           (pet.spstd.spei<pet_spstd_bounds_spei[1]) |
           (pet.spstd.spei>pet_spstd_bounds_spei[2]))

flm_df <- flm_df %>%
  mutate(outlier = 
           (estimate_cwd.an.spstd.tc<cwd_est_bounds.tc[1]) |
           (estimate_cwd.an.spstd.tc>cwd_est_bounds.tc[2]) |
           (estimate_pet.an.spstd.tc<pet_est_bounds.tc[1]) |
           (estimate_pet.an.spstd.tc>pet_est_bounds.tc[2]))


# flm_df <- flm_df %>%
#   mutate(outlier = 
#            (estimate_cwd.an<cwd_est_bounds[1]) |
#            (estimate_cwd.an>cwd_est_bounds[2]) |
#            (estimate_pet.an<pet_est_bounds[1]) |
#            (estimate_pet.an>pet_est_bounds[2]) |
#            (estimate_cwd.an.spstd.tc<cwd_est_bounds.tc[1]) |
#            (estimate_cwd.an.spstd.tc>cwd_est_bounds.tc[2]) |
#            (estimate_pet.an.spstd.tc<pet_est_bounds.tc[1]) |
#            (estimate_pet.an.spstd.tc>pet_est_bounds.tc[2]) |
#            (estimate_cwd.an.spstd.spei<cwd_est_bounds.spei[1]) |
#            (estimate_cwd.an.spstd.spei>cwd_est_bounds.spei[2]) |
#            (estimate_pet.an.spstd.spei<pet_est_bounds.spei[1]) |
#            (estimate_pet.an.spstd.spei>pet_est_bounds.spei[2]) |
#            (cwd.spstd<cwd_spstd_bounds[1]) |
#            (cwd.spstd>cwd_spstd_bounds[2]) |
#            (pet.spstd<pet_spstd_bounds[1]) |
#            (pet.spstd>pet_spstd_bounds[2]) |
#            (cwd.spstd.tc<cwd_spstd_bounds_tc[1]) |
#            (cwd.spstd.tc>cwd_spstd_bounds_tc[2]) |
#            (pet.spstd.tc<pet_spstd_bounds_tc[1]) |
#            (pet.spstd.tc>pet_spstd_bounds_tc[2]))

# # Save out full flm_df to simplify downstream scripts and ensure consistency
# flm_df %>% write.csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv"))

# Trim outliers
trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

site_points=st_as_sf(trim_df,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

vg <-variogram(estimate_cwd.an.spstd.tc~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
plot(vg, vg.fit)
# print(paste0("Range before hitting sill (km): "), as.character(vg.fit[2,3]))

vg.range = vg.fit[2,3] * 1000
vg.range


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Quick test of primary regression ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formula = as.formula("estimate_cwd.an ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)

formula = as.formula("estimate_cwd.an.spstd.tc ~ cwd.spstd.tc + (cwd.spstd.tc^2) + pet.spstd.tc + (pet.spstd.tc^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights_tc,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)

marg_fx_df <- function(mod){
  inc <- 0.1
  min <- -2.5
  max <- 2.5
  cwd_pred <- predictions(mod, newdata = datagrid(pet.spstd.tc = 0, cwd.spstd.tc = seq(min,max,inc))) %>% 
    mutate(variation = "cwd")
  pet_pred <- predictions(mod, newdata = datagrid(pet.spstd.tc = seq(min,max,inc), cwd.spstd.tc = 0)) %>% 
    mutate(variation = "pet")
  return(rbind(cwd_pred, pet_pred))
}


preds <- marg_fx_df(cwd_mod)

cwd_mfx_plot <- preds %>% 
  filter(variation == "cwd") %>% 
  ggplot(aes(x = cwd.spstd.tc)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  xlab("Historic standardized CWD") +
  ylab("Estimated sensitivity to CWD") +
  theme_bw()
cwd_mfx_plot






formula = as.formula("estimate_cwd.an ~ cwd.spstd.tc + (cwd.spstd.tc^2) + pet.spstd.tc + (pet.spstd.tc^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)

formula = as.formula("estimate_cwd.an.spstd.tc ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)")
mod_data <- trim_df
cwd_mod_tc <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights_tc,
                    vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod_tc)

formula = as.formula("estimate_cwd.an.spstd.spei ~ cwd.spstd.spei + (cwd.spstd.spei^2) + pet.spstd.spei + (pet.spstd.spei^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights_spei,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)


formula = as.formula("estimate_cwd.an.spstd.tc ~ cwd.spstd.tc + (cwd.spstd.tc^2) + pet.spstd.tc + (pet.spstd.tc^2)")
mod_data <- trim_df
cwd_mod_tc <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights_tc,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod_tc)

summary(lm(estimate_cwd.an.spstd.tc ~ estimate_cwd.an, data = trim_df))
trim_df %>% 
  ggplot(aes(x = estimate_cwd.an, y = estimate_cwd.an.spstd.tc)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth()


summary(lm(cwd.spstd.tc ~ cwd.spstd, data = trim_df))
trim_df %>% 
  ggplot(aes(x = cwd.spstd, y = cwd.spstd.tc)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0)


summary(lm(std.error_cwd.an.spstd.tc ~ std.error_cwd.an, data = trim_df))
trim_df %>% 
  ggplot(aes(x = std.error_cwd.an, y = std.error_cwd.an.spstd.tc)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0,.3) +
  ylim(0, .3)


trim_df %>% 
  ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) +
  geom_point() +
  geom_smooth()

trim_df %>% 
  ggplot(aes(x = cwd.spstd.tc, y = estimate_cwd.an.spstd.tc)) +
  geom_point() +
  geom_smooth()


marg_fx_df <- function(mod){
  inc <- 0.1
  min <- -2.5
  max <- 2.5
  ppt_pred <- predictions(mod, newdata = datagrid(pet.spstd = 0, ppt.spstd = seq(min,max,inc))) %>% 
    mutate(variation = "ppt")
  pet_pred <- predictions(mod, newdata = datagrid(pet.spstd = seq(min,max,inc), ppt.spstd = 0)) %>% 
    mutate(variation = "pet")
  return(rbind(ppt_pred, pet_pred))
}


preds <- marg_fx_df(cwd_mod)

cwd_mfx_plot <- preds %>% 
  filter(variation == "ppt") %>% 
  ggplot(aes(x = ppt.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
cwd_mfx_plot

pet_mfx_plot <- preds %>% 
  filter(variation == "pet") %>% 
  ggplot(aes(x = pet.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
pet_mfx_plot

formula = as.formula("estimate_pet.an ~ ppt.spstd + pet.spstd + (ppt.spstd^2) + (pet.spstd^2)")
pet_mod <- feols(formula, weights = mod_data$pet_errorweights, data = mod_data,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(pet_mod)
preds <- marg_fx_df(pet_mod)

cwd_mfx_plot <- preds %>% 
  filter(variation == "ppt") %>% 
  ggplot(aes(x = ppt.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
cwd_mfx_plot

pet_mfx_plot <- preds %>% 
  filter(variation == "pet") %>% 
  ggplot(aes(x = pet.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
pet_mfx_plot
