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

out_dir <- paste0(wdir, '3_results/figures/comments/')

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_std.csv'))
flm_df_tc <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_tc_std.csv'))
flm_df_spei <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_spei_std.csv'))
flm_df_cru <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_cru_std.csv'))
flm_df_ppt <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_ppt_std.csv'))
flm_df_tc <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_tc_std.csv'))
flm_df_58 <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_std_58.csv'))


flm_df <- flm_df_tc %>%
  select(collection_id, estimate_pet.an.spstd.tc, std.error_pet.an.spstd.tc, estimate_cwd.an.spstd.tc, std.error_cwd.an.spstd.tc) %>% 
  inner_join(flm_df, by = "collection_id")

flm_df <- flm_df_spei %>% 
  select(collection_id, estimate_pet.an.spstd.spei, std.error_pet.an.spstd.spei, estimate_cwd.an.spstd.spei, std.error_cwd.an.spstd.spei) %>% 
  inner_join(flm_df, by = "collection_id")

flm_df <- flm_df_cru %>% 
  select(collection_id, estimate_pet.an.spstd.cru, std.error_pet.an.spstd.cru, estimate_cwd.an.spstd.cru, std.error_cwd.an.spstd.cru) %>% 
  inner_join(flm_df, by = "collection_id")

flm_df <- flm_df_ppt %>% 
  select(collection_id, estimate_pet.an.ppt = estimate_pet.an, std.error_pet.an.ppt = std.error_pet.an, estimate_ppt.an, std.error_ppt.an) %>% 
  inner_join(flm_df, by = "collection_id")

flm_df <- flm_df_58 %>% 
  select(collection_id, estimate_pet.an.spstd.58 = estimate_pet.an, std.error_pet.an.spstd.58  = std.error_pet.an, 
         estimate_cwd.an.spstd.58 = estimate_cwd.an, std.error_cwd.an.spstd.58 = std.error_cwd.an) %>% 
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
         cwd_errorweights_cru = 1 / (std.error_cwd.an.spstd.cru),
         cwd_errorweights_58 = 1 / (std.error_cwd.an.spstd.58),
         ppt_errorweights = 1 / (std.error_ppt.an),
         errorweights2 = sqrt(ntrees),
         pet_errorweights = 1 / (std.error_pet.an.spstd.tc),
         int_errorweights = 1 / (std.error_intercept))

# Identify and trim extreme outliers
cwd_est_bounds.cru = quantile(flm_df$estimate_cwd.an.spstd.cru, c(0.01, 0.99),na.rm=T)
pet_est_bounds.cru = quantile(flm_df$estimate_pet.an.spstd.cru, c(0.01, 0.99),na.rm=T)
cwd_est_bounds.spei = quantile(flm_df$estimate_cwd.an.spstd.spei, c(0.01, 0.99),na.rm=T)
pet_est_bounds.spei = quantile(flm_df$estimate_pet.an.spstd.spei, c(0.01, 0.99),na.rm=T)
cwd_est_bounds.tc = quantile(flm_df$estimate_cwd.an.spstd.tc, c(0.01, 0.99),na.rm=T)
pet_est_bounds.tc = quantile(flm_df$estimate_pet.an.spstd.tc, c(0.01, 0.99),na.rm=T)
cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds = quantile(flm_df$estimate_pet.an, c(0.01, 0.99),na.rm=T)
cwd_est_bounds.58 = quantile(flm_df$estimate_cwd.an.spstd.58, c(0.01, 0.99),na.rm=T)
pet_est_bounds.58 = quantile(flm_df$estimate_pet.an.spstd.58, c(0.01, 0.99),na.rm=T)
ppt_est_bounds.ppt = quantile(flm_df$estimate_ppt.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds.ppt = quantile(flm_df$estimate_pet.an.ppt, c(0.01, 0.99),na.rm=T)
cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)
cwd_spstd_bounds_tc = quantile(flm_df$cwd.spstd.tc, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds_tc = quantile(flm_df$pet.spstd.tc, c(0.01, 0.99), na.rm = T)
cwd_spstd_bounds_spei = quantile(flm_df$cwd.spstd.spei, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds_spei = quantile(flm_df$pet.spstd.spei, c(0.01, 0.99), na.rm = T)
cwd_spstd_bounds_cru = quantile(flm_df$cwd.spstd.cru, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds_cru = quantile(flm_df$pet.spstd.cru, c(0.01, 0.99), na.rm = T)

# flm_df <- flm_df %>%
#   mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
#            (estimate_cwd.an>cwd_est_bounds[2]) |
#            (estimate_pet.an<pet_est_bounds[1]) |
#            (estimate_pet.an>pet_est_bounds[2]) |
#            (cwd.spstd<cwd_spstd_bounds[1]) |
#            (cwd.spstd>cwd_spstd_bounds[2]) |
#            (pet.spstd<pet_spstd_bounds[1]) |
#            (pet.spstd>pet_spstd_bounds[2]))
# 
# flm_df <- flm_df %>%
#   mutate(outlier = 
#            (estimate_cwd.an.spstd.cru<cwd_est_bounds.cru[1]) |
#            (estimate_cwd.an.spstd.cru>cwd_est_bounds.cru[2]) |
#            (estimate_pet.an.spstd.cru<pet_est_bounds.cru[1]) |
#            (estimate_pet.an.spstd.cru>pet_est_bounds.cru[2]))
# # |
# #            (cwd.spstd.spei<cwd_spstd_bounds_spei[1]) |
# #            (cwd.spstd.spei>cwd_spstd_bounds_spei[2]) |
# #            (pet.spstd.spei<pet_spstd_bounds_spei[1]) |
# #            (pet.spstd.spei>pet_spstd_bounds_spei[2]))
# 
# 
# flm_df <- flm_df %>%
#   mutate(outlier = 
#            (estimate_cwd.an.spstd.tc<cwd_est_bounds.tc[1]) |
#            (estimate_cwd.an.spstd.tc>cwd_est_bounds.tc[2]) |
#            (estimate_pet.an.spstd.tc<pet_est_bounds.tc[1]) |
#            (estimate_pet.an.spstd.tc>pet_est_bounds.tc[2]))
# 
# 
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
# 
# # Save out full flm_df to simplify downstream scripts and ensure consistency
# flm_df %>% write.csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv"))
# 
# # Trim outliers
# trim_df <- flm_df %>% 
#   filter(outlier==0) %>% 
#   drop_na()

site_points=st_as_sf(trim_df,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

vg <-variogram(estimate_cwd.an.spstd.tc~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
plot(vg, vg.fit)
# print(paste0("Range before hitting sill (km): "), as.character(vg.fit[2,3]))

vg.range = vg.fit[2,3] * 1000
vg.range


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main specification ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_data <- flm_df %>%
  mutate(outlier =
           (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2])) %>% 
  filter(outlier==0)
formula = as.formula("estimate_cwd.an ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)")
cwd_mod <- feols(formula, data = cwd_data, weights = cwd_data$cwd_errorweights,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)


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
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  xlab("Historic standardized CWD") +
  ylab("Estimated sensitivity to CWD") +
  theme_bw()
cwd_mfx_plot
ggsave(paste0(out_dir, "cwd_mfx.png"), cwd_mfx_plot, height = 4, width = 6)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compare to terraclimate ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tc_data <- flm_df %>%
  mutate(outlier = 
           (estimate_cwd.an.spstd.tc<cwd_est_bounds.tc[1]) |
           (estimate_cwd.an.spstd.tc>cwd_est_bounds.tc[2]) |
           (estimate_pet.an.spstd.tc<pet_est_bounds.tc[1]) |
           (estimate_pet.an.spstd.tc>pet_est_bounds.tc[2])) %>% 
  filter(outlier == 0)

formula = as.formula("estimate_cwd.an.spstd.tc ~ cwd.spstd.tc + (cwd.spstd.tc^2) + pet.spstd.tc + (pet.spstd.tc^2)")
tc_mod <- feols(formula, data = tc_data, weights = tc_data$cwd_errorweights_tc,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(tc_mod)

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
preds <- marg_fx_df(tc_mod)
tc_mfx_plot <- preds %>% 
  filter(variation == "cwd") %>% 
  ggplot(aes(x = cwd.spstd.tc)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  xlab("Historic standardized CWD") +
  ylab("Estimated sensitivity to CWD") +
  theme_bw()
tc_mfx_plot
ggsave(paste0(out_dir, "tc_mfx.png"), tc_mfx_plot, height = 4, width = 6)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Restrict primary analysis to post-58 ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd58_data <- flm_df %>%
  mutate(outlier =
           (estimate_cwd.an.spstd.58<cwd_est_bounds.58[1]) |
           (estimate_cwd.an.spstd.58>cwd_est_bounds.58[2]) |
           (estimate_pet.an.spstd.58<pet_est_bounds.58[1]) |
           (estimate_pet.an.spstd.58>pet_est_bounds.58[2])) %>% 
  filter(outlier==0)
formula = as.formula("estimate_cwd.an.spstd.58 ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)")
cwd58_mod <- feols(formula, data = cwd58_data, weights = cwd58_data$cwd_errorweights_58,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd58_mod)


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
preds <- marg_fx_df(cwd58_mod)
cwd58_mfx_plot <- preds %>% 
  filter(variation == "cwd") %>% 
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  xlab("Historic standardized CWD") +
  ylab("Estimated sensitivity to CWD") +
  theme_bw()
cwd58_mfx_plot

ggsave(paste0(out_dir, "cwd58_mfx.png"), cwd58_mfx_plot, height = 4, width = 6)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compare to CRU  ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cru_data <- flm_df %>%
  mutate(outlier = 
           (estimate_cwd.an.spstd.cru<cwd_est_bounds.cru[1]) |
           (estimate_cwd.an.spstd.cru>cwd_est_bounds.cru[2]) |
           (estimate_pet.an.spstd.cru<pet_est_bounds.cru[1]) |
           (estimate_pet.an.spstd.cru>pet_est_bounds.cru[2])) %>% 
  filter(outlier == 0)

formula = as.formula("estimate_cwd.an.spstd.cru ~ cwd.spstd.cru + (cwd.spstd.cru^2) + pet.spstd.cru + (pet.spstd.cru^2)")
cru_mod <- feols(formula, data = cru_data, weights = cru_data$cwd_errorweights_cru,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cru_mod)

marg_fx_df <- function(mod){
  inc <- 0.1
  min <- -2.5
  max <- 2.5
  cwd_pred <- predictions(mod, newdata = datagrid(pet.spstd.cru = 0, cwd.spstd.cru = seq(min,max,inc))) %>% 
    mutate(variation = "cwd")
  pet_pred <- predictions(mod, newdata = datagrid(pet.spstd.cru = seq(min,max,inc), cwd.spstd.cru = 0)) %>% 
    mutate(variation = "pet")
  return(rbind(cwd_pred, pet_pred))
}
preds <- marg_fx_df(cru_mod)
cru_mfx_plot <- preds %>% 
  filter(variation == "cwd") %>% 
  ggplot(aes(x = cwd.spstd.cru)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  xlab("Historic standardized CWD") +
  ylab("Estimated sensitivity to CWD") +
  theme_bw()
cru_mfx_plot
ggsave(paste0(out_dir, "cru_mfx.png"), cru_mfx_plot, height = 4, width = 6)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PPT-focused analysis  ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppt_data <- flm_df %>%
  mutate(outlier = 
           (estimate_ppt.an<ppt_est_bounds.ppt[1]) |
           (estimate_ppt.an>ppt_est_bounds.ppt[2]) |
           (estimate_pet.an.ppt<pet_est_bounds.ppt[1]) |
           (estimate_pet.an.ppt>pet_est_bounds.ppt[2])) %>% 
  filter(outlier == 0)

formula = as.formula("estimate_ppt.an ~ ppt.spstd + (ppt.spstd^2) + pet.spstd + (pet.spstd^2)")
ppt_mod <- feols(formula, data = ppt_data, weights = ppt_data$ppt_errorweights,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(ppt_mod)

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
preds <- marg_fx_df(ppt_mod)
ppt_mfx_plot <- preds %>% 
  filter(variation == "ppt") %>% 
  ggplot(aes(x = ppt.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  xlab("Historic standardized PPT") +
  ylab("Estimated sensitivity to PPT") +
  theme_bw()
ppt_mfx_plot
ggsave(paste0(out_dir, "ppt_mfx.png"), ppt_mfx_plot, height = 4, width = 6)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SPEI-based PET calculations  ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spei_data <- flm_df %>%
  mutate(outlier =
           (estimate_cwd.an.spstd.spei<cwd_est_bounds.spei[1]) |
           (estimate_cwd.an.spstd.spei>cwd_est_bounds.spei[2]) |
           (estimate_pet.an.spstd.spei<pet_est_bounds.spei[1]) |
           (estimate_pet.an.spstd.spei>pet_est_bounds.spei[2]))

formula = as.formula("estimate_cwd.an.spstd.spei ~ cwd.spstd.spei + (cwd.spstd.spei^2) + pet.spstd.spei + (pet.spstd.spei^2)")
spei_mod <- feols(formula, data = spei_data, weights = spei_data$cwd_errorweights_spei,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(spei_mod)

marg_fx_df <- function(mod){
  inc <- 0.1
  min <- -2.5
  max <- 2.5
  cwd_pred <- predictions(mod, newdata = datagrid(pet.spstd.spei = 0, cwd.spstd.spei = seq(min,max,inc))) %>% 
    mutate(variation = "cwd")
  pet_pred <- predictions(mod, newdata = datagrid(pet.spstd.spei = seq(min,max,inc), cwd.spstd.spei = 0)) %>% 
    mutate(variation = "pet")
  return(rbind(cwd_pred, pet_pred))
}
preds <- marg_fx_df(spei_mod)
spei_mfx_plot <- preds %>% 
  filter(variation == "cwd") %>% 
  ggplot(aes(x = cwd.spstd.spei)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  xlab("Historic standardized CWD") +
  ylab("Estimated sensitivity to CWD") +
  theme_bw()
spei_mfx_plot
# ggsave(paste0(out_dir, "cru_mfx.png"), cru_mfx_plot, height = 4, width = 6)





# feols(estimate_cwd.an ~ estimate_cwd.an.spstd.tc, data = flm_df)
# flm_df %>% 
#   ggplot(aes(x = estimate_cwd.an, y = estimate_cwd.an.spstd.tc)) + 
#   geom_point() +
#   geom_smooth()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compare raw climate data  ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Site-specific historic climate data
site_clim_csv <- paste0(wdir, '1_input_processed/climate/essentialcwd_data.csv')
site_clim_df <- read_csv(site_clim_csv)
site_clim_df %>% filter(is.na(cwd)) %>% pull(site) %>% unique()  ## Note: these are largely NAN due to missing SWC or topo data from previous scripts
site_clim_df <- site_clim_df %>% 
  mutate(site = as.character(site)) %>% 
  rename(location_id = site)

tc_pet <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_pet.csv"))
tc_cwd <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_def.csv"))
tc_ppt <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_ppt.csv"))
site_clim_df_tc <- tc_pet %>%
  left_join(tc_cwd, by = c("collection_id", "Month", "year")) %>%
  left_join(tc_ppt, by = c("collection_id", "Month", "year")) %>%
  rename(location_id = collection_id, 
         month = Month,
         pet_tc = pet,
         cwd_tc = def,
         ppt_tc = ppt)

site_clim_df <- site_clim_df %>% 
  left_join(site_clim_df_tc, by = c("location_id", "year", "month"))

clim_df <- site_clim_df %>% 
  drop_na() %>% # Mainly dropping FIA sites since we haven't pulled TC data for those
  rename(collection_id = location_id) %>% 
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
            cru_pet.an = sum(pet_cru),
            cru_cwd.an = sum(cwd_cru),
            .groups = "drop")


pet_compare <- clim_df %>%
  sample_n(20000) %>%
  ggplot(aes(y = tc_pet.an, x = pet.an)) +
  geom_point(alpha = .1) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red") +
  theme_bw() +
  xlab("Corrected annual PET") +
  ylab("Terraclimate annual PET")
pet_compare
ggsave(paste0(out_dir, "pet_compare.png"), pet_compare, height = 4, width = 6)


cwd_compare <- clim_df %>%
  sample_n(20000) %>%
  ggplot(aes(y = tc_cwd.an, x = cwd.an)) +
  geom_point(alpha = .1) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red") +
  theme_bw() +
  xlab("Corrected annual CWD") +
  ylab("Terraclimate annual CWD")
cwd_compare
ggsave(paste0(out_dir, "cwd_compare.png"), cwd_compare, height = 4, width = 6)


ppt_compare <- clim_df %>%
  sample_n(10000) %>%
  ggplot(aes(y = tc_ppt.an, x = ppt.an)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red") +
  theme_bw() +
  xlab("Corrected annual CWD") +
  ylab("Terraclimate annual CWD")
ppt_compare


mod <- lm(cwd_tc ~ cwd, clim_df)
summary(mod)

mod <- feols(cwd_tc ~ cwd | collection_id, clim_df)
summary(mod)

mod <- lm(pet_tc ~ pet, clim_df)
summary(mod)

mod <- feols(pet_tc ~ pet | collection_id, clim_df)
summary(mod)

mod <- lm(pet ~ pet_spei, clim_df)
summary(mod)

mod <- feols(pet ~ pet_spei | collection_id, clim_df)
summary(mod)

mod <- lm(pet_tc ~ pet_cru, clim_df)
summary(mod)

mod <- lm(pet ~ pet_cru, clim_df)
summary(mod)







clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = pet_spei, x = pet)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red") +
  ylim(0, 250)

clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = pet_tc, x = pet_spei)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")


clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = cwd_tc, x = cwd)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")


clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = cwd_cru, x = cwd)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")


















formula = as.formula("estimate_cwd.an ~ cwd.spstd.tc + (cwd.spstd.tc^2) + pet.spstd.tc + (pet.spstd.tc^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)



formula = as.formula("estimate_cwd.an.spstd.spei ~ cwd.spstd.spei + (cwd.spstd.spei^2) + pet.spstd.spei + (pet.spstd.spei^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights_spei,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)

formula = as.formula("estimate_cwd.an.spstd.tc ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)")
mod_data <- trim_df
cwd_mod_tc <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights_tc,
                    vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod_tc)




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
