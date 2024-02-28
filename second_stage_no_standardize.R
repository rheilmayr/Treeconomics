#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load packages --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("f_cwd_function_new.R")
library(dplyr)
library(naniar)
library(data.table)
library(tidyverse)
library(geosphere)
library(zoo)
library(tictoc)

library(tidyverse)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic()
# Define path
wdir <- 'remote/'

site_clim_csv <- paste0(wdir, '1_input_processed/climate/essentialcwd_data.csv')
site_clim_df <- read_csv(site_clim_csv)
site_clim_df <- site_clim_df %>% 
  mutate(site = as.character(site)) %>% 
  rename(collection_id = site)
site_clim_df <- site_clim_df %>% 
  group_by(collection_id, year) %>% 
  summarise(pet = sum(pet),
            cwd = sum(cwd))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz"))

# 3. Site information
dendro_dir <- paste0(wdir, "1_input_processed/dendro/")
site_df <- read_csv(paste0(dendro_dir, 'site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# 4. Species information
sp_info <- read_csv(paste0(wdir, '1_input_processed/species_ranges/species_metadata.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_df <- site_df %>% 
  left_join(sp_info, by = "species_id")


# 5. Dendrochronologies - used for single-stage model comparison
dendro_df <- read_csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            .groups = "drop") %>% 
  as_tibble()


## Merge data
site_clim_df <- site_clim_df %>% 
  inner_join(ave_site_clim, by = "collection_id")

dendro_df <- site_clim_df %>% 
  inner_join(dendro_df, by = c("collection_id", "year"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define first stage regression model  -------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_mod <- function(site_data, outcome = "rwi", water_var = "cwd.an", energy_var = "pet.an", mod_type = "lm"){
  failed <- F
  reg_error <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
  no_cwd_var <- (site_data %>% select(all_of(water_var)) %>% n_distinct() == 1)
  no_pet_var <- (site_data %>% select(all_of(energy_var)) %>% n_distinct() == 1)
  
  if (no_cwd_var | no_pet_var) {
    message(paste0("Site has no variation in ", water_var, " or ", energy_var))
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / pet data 
    tryCatch(
      expr = {
        formula <- as.formula(paste0(outcome, " ~ ", energy_var, " + ", water_var))
        if (mod_type == "lm"){
          mod <- lm(formula, data = site_data)
        }
        if (mod_type == "lme"){
          mod <- nlme::lme(formula,
                           data=site_data, method="REML",
                           random = ~ 1 | tree,
                           correlation = nlme::corAR1(form=~year|tree))
        }
        
        mod_sum <- summary(mod)
        mod_vcov <- vcov(mod)
        # cov <- list(int_cwd = mod_vcov[1, 2], 
        #             int_pet = mod_vcov[1, 3], 
        #             pet_cwd = mod_vcov[2, 3])
        nobs <- nobs(mod)
        mod <- tidy(mod) %>%
          mutate(term = term %>% str_replace("\\(Intercept\\)", "intercept")) %>% 
          filter(term %in% c('intercept', water_var, energy_var)) %>% 
          pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
        # mod <- mod %>% 
        #   rename_all(funs(stringr::str_replace_all(., energy_var, 'energy.an')))
        cov_var_name <- paste0("cov_int_", water_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("(Intercept)"), c(water_var)]
        cov_var_name <- paste0("cov_int_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("(Intercept)"), c(energy_var)]
        cov_var_name <- paste0("cov_", water_var %>% str_replace(".an", ""), "_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c(water_var), c(energy_var)]
        mod$r2 = mod_sum$r.squared
      },
      error = function(e){ 
        message("Returned regression error")
        reg_error <<- e[1]
        failed <<- T
      }
    )    
  }
  if (failed){
    return(NA)
  }
  return(tibble(mod = list(mod), nobs = nobs, ntrees = ntrees, error = reg_error))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run site-level regressions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_df <- dendro_df %>%
  # filter(year>1957) %>% 
  # drop_na() %>% 
  rename(cwd.an = cwd,
         pet.an = pet) %>%
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  # filter(nobs>10) %>% 
  nest()


fs_mod_bl <- partial(fs_mod, outcome = "rwi", water_var = "cwd.an", energy_var = "pet.an", mod_type = "lm")

site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = fs_mod_bl))

data_df <- site_df %>% 
  select(collection_id,data)

## Primary first stage results using CWD and PET
fs_df <- site_df %>% 
  select(collection_id, fs_result) %>% 
  unnest(fs_result)
fs_df <- fs_df[which(!(fs_df %>% pull(mod) %>% is.na())),]
fs_df <- fs_df %>% 
  unnest(mod) %>% 
  select(-error)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep and trim data for second stage -----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flm_df <- fs_df %>% 
  left_join(ave_site_clim, by = "collection_id")

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

flm_df <- flm_df %>%
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]) |
           (cwd.spstd<cwd_spstd_bounds[1]) |
           (cwd.spstd>cwd_spstd_bounds[2]) |
           (pet.spstd<pet_spstd_bounds[1]) |
           (pet.spstd>pet_spstd_bounds[2]))


flm_df <- flm_df %>%
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]))

# Save out full flm_df to simplify downstream scripts and ensure consistency
flm_df %>% write.csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv"))

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
formula = as.formula("estimate_cwd.an ~ cwd.ave + (cwd.ave^2) + pet.ave + (pet.ave^2)")
mod_data <- trim_df
cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights)
summary(cwd_mod)

marg_fx_df <- function(mod){
  inc <- 100
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
