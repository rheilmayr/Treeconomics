#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run regressions to explore impact of historical climate on weather sensitivity
#
# Input files:
#   site_ave_clim.gz:
#   site_pet_cwd_std.csv:
#   site_summary.csv:
# 
# Output files:
#   ss_bootstrap.rds
#
# ToDo:
# - Update / finalize genus analyses
# - Rebuild robustness tests based on final baseline model
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(MASS)
library(car)
library(tidyverse)
library(broom)
library(purrr)
library(patchwork)
library(dbplyr)
library(RSQLite)
library(lmtest)
library(sandwich)
library(prediction)
library(Hmisc)
library(hexbin)
library(ggpubr)
library(ggiraphExtra)
library(modi)
library(margins)
library(tidylog)
library(fixest)
library(biglm)
library(gstat)
library(sf)
library(units)
library(dtplyr)
library(marginaleffects)
library(tidylog)

set.seed(5597)

select <- dplyr::select



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'out/first_stage/fia_pet_cwd_std.csv'))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz")) %>% 
  distinct()
flm_df <- flm_df %>% 
  left_join(ave_site_clim, by = c("collection_id"))

# 3. Site information
site_df <- read_csv(paste0(wdir, 'out/dendro/site_summary_fia.csv'))
site_df <- site_df %>%
  select(collection_id, species_id, latitude, longitude) %>% 
  distinct()
site_df <- site_df %>%
  mutate(species_id = str_to_lower(species_id))


# Merge back into main flm_df
flm_df <- flm_df %>%
  left_join(site_df, by = "collection_id")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep and trim data -----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add weighting based on inverse of first stage variance
flm_df <- flm_df %>% 
  mutate(cwd_errorweights = 1 / (std.error_cwd.an),
         errorweights2 = sqrt(ntrees),
         pet_errorweights = 1 / (std.error_pet.an),
         int_errorweights = 1 / (std.error_intercept))

# Identify and trim extreme outliers
cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds = quantile(flm_df$estimate_pet.an, c(0.01, 0.99),na.rm=T)
cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)

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
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]))


# Trim outliers
trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spatial autocorrelation of trim_df ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_points=st_as_sf(trim_df,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
plot(vg, vg.fit)
# print(paste0("Range before hitting sill (km): "), as.character(vg.fit[2,3]))

vg.range = vg.fit[2,3] * 1000


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Quick test of primary regression ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)


cwd_median <- trim_df$cwd.spstd %>% median()
mod_slopes <- slopes(cwd_mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = cwd_median))
mod_slopes <- mod_slopes %>% 
  as_tibble() %>% 
  filter(term == "cwd.spstd")
mod_slopes


marg_fx_df <- function(mod){
  inc <- 0.1
  min <- -2.5
  max <- 2.5
  cwd_pred <- predictions(mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = seq(min,max,inc))) %>% 
    mutate(variation = "cwd")
  pet_pred <- predictions(mod, newdata = datagrid(pet.spstd = seq(min,max,inc), cwd.spstd = 0)) %>% 
    mutate(variation = "pet")
  return(rbind(cwd_pred, pet_pred))
}


preds <- marg_fx_df(cwd_mod)

cwd_mfx_plot <- preds %>% 
  filter(variation == "cwd") %>% 
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
cwd_mfx_plot

pet_mfx_plot <- preds %>% 
  filter(variation == "pet") %>% 
  ggplot(aes(x = pet.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
pet_mfx_plot

formula = as.formula("estimate_pet.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
pet_mod <- feols(formula, weights = mod_data$pet_errorweights, data = mod_data,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(pet_mod)
preds <- marg_fx_df(pet_mod)

cwd_mfx_plot <- preds %>% 
  filter(variation == "cwd") %>% 
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
cwd_mfx_plot

pet_mfx_plot <- preds %>% 
  filter(variation == "pet") %>% 
  ggplot(aes(x = pet.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
pet_mfx_plot



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compare to ITRDB ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Site-level regressions
flm_itrdb <- read_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std.csv'))

# 2. Historic site-level climate
flm_itrdb <- flm_itrdb %>% 
  left_join(ave_site_clim, by = c("collection_id"))

# 3. Site information
site_itrdb <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_itrdb <- site_itrdb %>% 
  select(collection_id, sp_id, latitude, longitude)
site_itrdb <- site_itrdb %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# Merge back into main flm_df
flm_itrdb <- flm_itrdb %>% 
  left_join(site_itrdb, by = "collection_id")

# Add weighting based on inverse of first stage variance
flm_itrdb <- flm_itrdb %>% 
  mutate(cwd_errorweights = 1 / (std.error_cwd.an),
         errorweights2 = sqrt(ntrees),
         pet_errorweights = 1 / (std.error_pet.an),
         int_errorweights = 1 / (std.error_intercept))

# Identify and trim extreme outliers
cwd_est_bounds = quantile(flm_itrdb$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds = quantile(flm_itrdb$estimate_pet.an, c(0.01, 0.99),na.rm=T)
cwd_spstd_bounds = quantile(flm_itrdb$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(flm_itrdb$pet.spstd, c(0.01, 0.99), na.rm = T)

# flm_df <- flm_df %>%
#   mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
#            (estimate_cwd.an>cwd_est_bounds[2]) |
#            (estimate_pet.an<pet_est_bounds[1]) |
#            (estimate_pet.an>pet_est_bounds[2]) |
#            (cwd.spstd<cwd_spstd_bounds[1]) |
#            (cwd.spstd>cwd_spstd_bounds[2]) |
#            (pet.spstd<pet_spstd_bounds[1]) |
#            (pet.spstd>pet_spstd_bounds[2]))


flm_itrdb <- flm_itrdb %>%
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]))




test_itrdb <- flm_itrdb %>%
  filter(species_id %in% c("pipo", "psme", "pied", "pifl"),
         outlier == 0)

test_itrdb %>% select(estimate_cwd.an, std.error_cwd.an, ntrees) %>% summary()

test_fia <- trim_df %>%
  filter(species_id %in% c("pipo", "psme", "pied", "pifl"))
test_fia %>% select(estimate_cwd.an, std.error_cwd.an, ntrees) %>% summary()

cwd_range_fia <- flm_df %>% 
  filter(species_id %in% c("pipo", "psme", "pied", "pifl")) %>% 
  pull(cwd.spstd) %>% 
  quantile(c(0.025, 0.975))
cwd_range_itrdb <- flm_itrdb %>% 
  filter(species_id %in% c("pipo", "psme", "pied", "pifl")) %>% 
  pull(cwd.spstd) %>% 
  quantile(c(0.025, 0.975))
pet_range_fia <- flm_df %>% 
  filter(species_id %in% c("pipo", "psme", "pied", "pifl")) %>% 
  pull(pet.spstd) %>% 
  quantile(c(0.025, 0.975))
pet_range_itrdb <- flm_itrdb %>% 
  filter(species_id %in% c("pipo", "psme", "pied", "pifl")) %>% 
  pull(pet.spstd) %>% 
  quantile(c(0.025, 0.975))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spatial autocorrelation of datasets ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_points=st_as_sf(test_itrdb,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
vg.range_itrdb = vg.fit[2,3]


site_points=st_as_sf(test_fia,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
vg.range_fia = vg.fit[2,3]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A: ME of CWD by cwd.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
mod_itrdb <- feols(formula, data = test_itrdb, weights = test_itrdb$cwd_errorweights,
                 vcov = conley(cutoff = vg.range_itrdb, distance = "spherical"))
mod_fia <- feols(formula, data = test_fia, weights = test_fia$cwd_errorweights,
                   vcov = conley(cutoff = vg.range_fia, distance = "spherical"))

summary(mod_itrdb)
summary(mod_fia)

at_pet <- 0
seq_inc <- 0.1

### Marginal effects of CWD - ITRDB
pred_fia <- predictions(mod_fia, newdata = datagrid(pet.spstd = at_pet, cwd.spstd = seq(cwd_range_fia[1],cwd_range_fia[2],seq_inc))) %>% 
  mutate(dataset = "fia")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(mod_itrdb, newdata = datagrid(pet.spstd = at_pet, cwd.spstd = seq(cwd_range_itrdb[1],cwd_range_itrdb[2],seq_inc))) %>% 
  mutate(dataset = "itrdb")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
cwd_cwdstd_plot <- pred_all %>% 
  ggplot(aes(x = cwd.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = dataset, fill = dataset)) + 
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme_bw() +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD")  +
  scale_color_manual(values = c("itrdb" = "darkred",
                                "fia"="darkblue")) +
  scale_fill_manual(values = c("itrdb" = "darkred",
                               "fia"="darkblue"))
cwd_cwdstd_plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# B: ME of CWD by pet.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Marginal effects of CWD - ITRDB
pred_fia <- predictions(mod_fia, newdata = datagrid(pet.spstd = seq(pet_range_fia[1],pet_range_fia[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(dataset = "fia")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(mod_itrdb, newdata = datagrid(pet.spstd = seq(pet_range_itrdb[1],pet_range_itrdb[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(dataset = "itrdb")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
cwd_petstd_plot <- pred_all %>% 
  ggplot(aes(x = pet.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = dataset, fill = dataset)) + 
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme_bw() +  
  xlab("Historic PET\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD")  +
  scale_color_manual(values = c("itrdb" = "darkred",
                                "fia"="darkblue")) +
  scale_fill_manual(values = c("itrdb" = "darkred",
                                "fia"="darkblue"))
cwd_petstd_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# C: ME of PET by cwd.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formula = as.formula("estimate_pet.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
mod_itrdb <- feols(formula, data = test_itrdb, weights = test_itrdb$pet_errorweights,
                   vcov = conley(cutoff = vg.range_itrdb, distance = "spherical"))
mod_fia <- feols(formula, data = test_fia, weights = test_fia$pet_errorweights,
                 vcov = conley(cutoff = vg.range_itrdb, distance = "spherical"))

summary(mod_itrdb)
summary(mod_fia)

### Marginal effects of CWD - ITRDB
pred_fia <- predictions(mod_fia, newdata = datagrid(pet.spstd = 0, cwd.spstd = seq(cwd_range_fia[1],cwd_range_fia[2],seq_inc))) %>% 
  mutate(dataset = "fia")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(mod_itrdb, newdata = datagrid(pet.spstd = at_pet, cwd.spstd = seq(cwd_range_itrdb[1],cwd_range_itrdb[2],seq_inc))) %>% 
  mutate(dataset = "itrdb")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
pet_cwdstd_plot <- pred_all %>% 
  ggplot(aes(x = cwd.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = dataset, fill = dataset)) + 
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme_bw()  +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to PET")  +
  scale_color_manual(values = c("itrdb" = "darkred",
                                "fia"="darkblue")) +
  scale_fill_manual(values = c("itrdb" = "darkred",
                               "fia"="darkblue"))
pet_cwdstd_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# D: ME of PET by pet.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Marginal effects of CWD - ITRDB
pred_fia <- predictions(mod_fia, newdata = datagrid(pet.spstd = seq(pet_range_fia[1],pet_range_fia[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(dataset = "fia")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(mod_itrdb, newdata = datagrid(pet.spstd = seq(pet_range_itrdb[1],pet_range_itrdb[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(dataset = "itrdb")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
pet_petstd_plot <- pred_all %>% 
  ggplot(aes(x = pet.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = dataset, fill = dataset)) + 
  geom_line(colors = c("darkred", "darkblue")) +
  geom_ribbon(alpha=0.2) +
  theme_bw() +
  xlab("Historic PET\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to PET")  +
  scale_color_manual(values = c("itrdb" = "darkred",
                                "fia"="darkblue")) +
  scale_fill_manual(values = c("itrdb" = "darkred",
                               "fia"="darkblue"))
pet_petstd_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine plots ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combined_plot <- (cwd_cwdstd_plot / cwd_petstd_plot) | (pet_cwdstd_plot / pet_petstd_plot)  

combined_plot <- combined_plot +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

combined_plot



flm_itrdb %>% pull(ntrees) %>% summary()
flm_df %>% pull(ntrees) %>% summary()
test_itrdb %>% pull(species_id) %>% unique()
test_fia %>% pull(species_id) %>% unique()
test_fia %>% pull(plo)