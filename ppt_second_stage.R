# Add weighting based on inverse of first stage variance
flm_df <- flm_df %>% 
  mutate(ppt_errorweights = 1 / (std.error_ppt.an),
         errorweights2 = sqrt(ntrees),
         pet_errorweights = 1 / (std.error_pet.an),
         int_errorweights = 1 / (std.error_intercept))

# Identify and trim extreme outliers
ppt_est_bounds = quantile(flm_df$estimate_ppt.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds = quantile(flm_df$estimate_pet.an, c(0.01, 0.99),na.rm=T)
ppt_spstd_bounds = quantile(flm_df$ppt.spstd, c(0.01, 0.99), na.rm = T)
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
  mutate(outlier = (estimate_ppt.an<ppt_est_bounds[1]) |
           (estimate_ppt.an>ppt_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]))

# Save out full flm_df to simplify downstream scripts and ensure consistency
flm_df %>% write.csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv"))

# Trim outliers
trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

site_points=st_as_sf(trim_df,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

vg <-variogram(estimate_ppt.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
plot(vg, vg.fit)
# print(paste0("Range before hitting sill (km): "), as.character(vg.fit[2,3]))

vg.range = vg.fit[2,3] * 1000



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Quick test of primary regression ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formula = as.formula("estimate_ppt.an ~ ppt.spstd + (ppt.spstd^2) + pet.spstd + (pet.spstd^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$ppt_errorweights,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
summary(cwd_mod)

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
