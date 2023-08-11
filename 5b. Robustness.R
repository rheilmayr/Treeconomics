#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run regressions to explore impact of historical climate on weather sensitivity
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using species_niche.R
# - tree_ring_data_V2.db: Compiled database of ITRDB observations 
# - essentialcwd_data.csv: File detailing plot-level weather history
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
source("f_spec_chart_function.R")
library(lme4)
library(multcomp)

select <- dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
fs_spl <- read_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std.csv'))
fs_nb <- read_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std_nb.csv'))
fs_ar <- read_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std_ar.csv'))

fs_re <- read_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std_re.csv')) # Uses lme with AR-correlated standard errors
fs_temp <- read_csv(paste0(wdir, 'out/first_stage/site_temp_cwd_std.csv'))
fs_cum <- read_rds(paste0(wdir, 'out/first_stage/dnlm_cum_effects.rds'))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))

# 3. Site information
site_df <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_df <- site_df %>% 
  left_join(sp_info, by = "species_id")


niche_df <- read.csv(paste0(wdir, "out/climate/clim_niche.csv")) %>%
  select(-X)


# 5. Dendrochronologies - used for single-stage model comparison
dendro_dir <- paste0(wdir, "out/dendro/")
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


# 6. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "out/climate/site_an_clim.gz"))
dendro_df <- dendro_df %>% 
  left_join(an_site_clim, by = c("collection_id", "year"))

dendro_df <- dendro_df %>% 
  left_join(ave_site_clim, by = "collection_id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep and trim data -----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_fs <- function(flm_df){
  flm_df <- flm_df %>% 
    left_join(ave_site_clim, by = c("collection_id"))
  
  # Add weighting based on inverse of first stage variance
  flm_df <- flm_df %>% 
    mutate(cwd_errorweights = 1 / (std.error_cwd.an),
           tree_errorweights = sqrt(ntrees),
           pet_errorweights = 1 / (std.error_pet.an),
           int_errorweights = 1 / (std.error_intercept),
           equalweights = 1)
  
  # Identify and trim extreme outliers
  cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
  pet_est_bounds = quantile(flm_df$estimate_pet.an, c(0.01, 0.99),na.rm=T)
  cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
  pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)
  
  
  flm_df <- flm_df %>% 
    mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
             (estimate_cwd.an>cwd_est_bounds[2]) |
             (estimate_pet.an<pet_est_bounds[1]) |
             (estimate_pet.an>pet_est_bounds[2]))
  
  flm_df <- flm_df %>% 
    left_join(site_df, by = "collection_id")
  
  flm_df <- flm_df %>% 
    mutate(trim_y = !outlier,
           trim_x = cwd.spstd>cwd_spstd_bounds[1] & cwd.spstd<cwd_spstd_bounds[2] & 
             pet.spstd>pet_spstd_bounds[1] & pet.spstd<pet_spstd_bounds[2])

  return(flm_df)
}

process_fs_temp <- function(flm_df){
  flm_df <- flm_df %>% 
    left_join(ave_site_clim, by = c("collection_id"))
  
  # Add weighting based on inverse of first stage variance
  flm_df <- flm_df %>% 
    mutate(cwd_errorweights = 1 / (std.error_cwd.an),
           tree_errorweights = sqrt(ntrees),
           temp_errorweights = 1 / (std.error_temp.an),
           int_errorweights = 1 / (std.error_intercept),
           equalweights = 1)
  
  # Identify and trim extreme outliers
  cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
  temp_est_bounds = quantile(flm_df$estimate_temp.an, c(0.01, 0.99),na.rm=T)
  cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
  temp_spstd_bounds = quantile(flm_df$temp.spstd, c(0.01, 0.99), na.rm = T)
  
  
  flm_df <- flm_df %>% 
    mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
             (estimate_cwd.an>cwd_est_bounds[2]) |
             (estimate_temp.an<temp_est_bounds[1]) |
             (estimate_temp.an>temp_est_bounds[2]))
  
  flm_df <- flm_df %>% 
    left_join(site_df, by = "collection_id")
  
  flm_df <- flm_df %>% 
    mutate(trim_y = !outlier,
           trim_x = cwd.spstd>cwd_spstd_bounds[1] & cwd.spstd<cwd_spstd_bounds[2] & 
             temp.spstd>temp_spstd_bounds[1] & temp.spstd<temp_spstd_bounds[2])
  
  
  return(flm_df)
}

fs_spl <- fs_spl %>% 
  process_fs()
fs_nb <- fs_nb %>% 
  process_fs()
fs_ar <- fs_ar %>% 
  process_fs()
fs_re <- fs_re %>% 
  process_fs()
fs_temp <- fs_temp %>% 
  process_fs_temp()

fs_cum <- fs_cum %>%  # Has already been trimmed and had weights calculated in dnlm script
  # left_join(ave_site_clim, by = c("collection_id")) %>% 
  left_join(site_df, by = "collection_id") %>% 
  rename(estimate_cwd.an = cwd_cum, cwd_errorweights = errorweights)
fs_cum <- fs_cum %>% 
  mutate(trim_y = 1,
         trim_x = 1)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spatial autocorrelation of trim_df ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_points=st_as_sf(fs_spl %>% filter(outlier == FALSE),coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
plot(vg, vg.fit)
vg.range = vg.fit[2,3] * 1000



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Robustness tests --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameters to vary
#  - Trim X
#  - Trim y
#  - Weighting
#  - Linear
#  - Species control
#  - Detrend method
#  - Energy control (temp in place of pet)
#  - Single-stage RE model

robustness_test <- function(params){
  data <- params$despline_data
  formula <- params$formula
  if (params$t_x) {data <- data %>% filter(trim_x==TRUE)}
  if (params$t_y) {data <- data %>% filter(trim_y==TRUE)}
  cwd_median <- data %>% select(cwd.spstd) %>% drop_na() %>% pull(cwd.spstd) %>% median()
  error_weights <- data[[params$weights]]
  if (params$mod_type == "fe") {
    mod <- feols(formula, weights = error_weights, data = data,
                 vcov = conley(cutoff = vg.range/1000, distance = "spherical"))

    if (formula == as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd")) {
      test_str <- "cwd.spstd = 0"
    } else{
      test_str <- paste0("cwd.spstd + (2 * \`I(cwd.spstd^2)\` * ", as.character(cwd_median), ") = 0")
    }
    }
  if (params$mod_type == "re") {
    mod <- lmer(formula, data = data)
    test_str <- paste0("cwd.an.spstd:cwd.spstd + (2 * cwd.an.spstd:I(cwd.spstd^2) * ", as.character(cwd_median), ") = 0")
  }
  lincom <- glht(mod, linfct = c(test_str))
  lincom <- summary(lincom)
  estimate <- lincom$test$coefficients 
  std.error <- lincom$test$sigma
  mod_slopes <- tibble("estimate" = estimate, "std.error" = std.error)
  return(mod_slopes)
}


specs <- data.frame(coef=NaN, 
                    se=NaN,
                    detrend_spline = NaN,
                    detrend_nb = NaN,
                    detrend_ar = NaN,
                    pet = NaN,
                    temp = NaN,
                    contemp_rwi = NaN,
                    cum_dnlm = NaN,
                    two_stage = NaN,
                    single_stage = NaN,
                    squared_term = NaN,
                    linear_mod = NaN,
                    tmprl = NaN,
                    species_control = NaN, 
                    weight_se = NaN,
                    trim_y = NaN,
                    trim_x = NaN) 


## Baseline model
params <- list(despline_data = fs_spl,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        pet = TRUE,
                        temp = FALSE,
                      tmprl = FALSE,
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                        species_control = FALSE,
                        squared_term = TRUE,
                        linear_mod = FALSE,
                        detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        cum_dnlm = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)








## Test alternate desplining
# Negative binomial
params <- list(despline_data = fs_nb,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = FALSE, 
                        detrend_nb = TRUE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)


# AR
params <- list(despline_data = fs_ar,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = FALSE, 
                        detrend_nb = FALSE,
                        detrend_ar = TRUE)
specs <- rbind(specs, new_row)

## Temp in place of PET
params <- list(despline_data = fs_temp,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + temp.spstd + (cwd.spstd^2) + (temp.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                      se = mod_slopes %>% pull("std.error"),
                      pet = FALSE,
                      temp = TRUE,
                      tmprl = FALSE,
                      trim_x = params$t_x,
                      trim_y = params$t_y,
                      weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                      species_control = FALSE,
                      squared_term = TRUE,
                      linear_mod = FALSE,
                      detrend_spline = TRUE, 
                      detrend_nb = FALSE,
                      cum_dnlm = FALSE,
                      detrend_ar = FALSE)
specs <- rbind(specs, new_row)


## Cumulative impact from dynamic lag
params <- list(despline_data = fs_cum,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                      se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                      trim_x = params$t_x,
                      trim_y = params$t_y,
                      weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = FALSE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                      species_control = FALSE,
                      squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = TRUE,
                      detrend_spline = TRUE, 
                      detrend_nb = FALSE,
                      detrend_ar = FALSE)
specs <- rbind(specs, new_row)


## Test alternate model structures
## Single stage RE model
params <- list(despline_data = dendro_df,
               mod_type = "re",
               formula = as.formula("rwi ~ cwd.an.spstd:pet.spstd + cwd.an.spstd:I(pet.spstd**2) +
               cwd.an.spstd:cwd.spstd + cwd.an.spstd:I(cwd.spstd**2) +
               pet.an.spstd:pet.spstd + pet.an.spstd:I(pet.spstd**2) +
               pet.an.spstd:cwd.spstd + pet.an.spstd:I(cwd.spstd**2) +
                                    (1 + cwd.an.spstd + pet.an.spstd | collection_id)"),
               t_x = FALSE,
               t_y = FALSE,
               weights = "equalweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                      se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                      trim_x = params$t_x,
                      trim_y = params$t_y,
                      weight_se = params$weights == "cwd_errorweights",
                      species_control = FALSE,
                      contemp_rwi = TRUE,
                      two_stage = FALSE,
                      single_stage = TRUE,
                      squared_term = TRUE,
                      linear_mod = FALSE,
                      detrend_spline = TRUE, 
                      detrend_nb = FALSE,
                      cum_dnlm = FALSE,
                      detrend_ar = FALSE)
specs <- rbind(specs, new_row)


## First stage controlling for temporal autocorrelation
params <- list(despline_data = fs_re,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                      se = mod_slopes %>% pull("std.error"),
                      pet = FALSE,
                      temp = FALSE,
                      tmprl = TRUE,
                      trim_x = params$t_x,
                      trim_y = params$t_y,
                      weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                      species_control = FALSE,
                      squared_term = TRUE,
                      linear_mod = FALSE,
                      detrend_spline = TRUE, 
                      detrend_nb = FALSE,
                      cum_dnlm = FALSE,
                      detrend_ar = FALSE)
specs <- rbind(specs, new_row)


# Linear model
params <- list(despline_data = fs_spl,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                        species_control = FALSE,
                        squared_term = FALSE,
                        linear_mod = TRUE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)


# # Species fixed effects
# params <- list(despline_data = fs_spl,
#                mod_type = "fe",
#                formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2) | species_id"),
#                t_x = FALSE,
#                t_y = TRUE,
#                weights = "cwd_errorweights")
# mod_slopes <- robustness_test(params)
# new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
#                         se = mod_slopes %>% pull("std.error"),
#                       pet = TRUE,
#                       temp = FALSE,
#                       tmprl = FALSE,
#                         trim_x = params$t_x,
#                         trim_y = params$t_y,
#                         weight_se = params$weights == "cwd_errorweights",
#                       contemp_rwi = TRUE,
#                       two_stage = TRUE,
#                       single_stage = FALSE,
#                         species_control = TRUE,
#                         squared_term = TRUE,
#                       linear_mod = FALSE,
#                       cum_dnlm = FALSE,
#                       detrend_spline = TRUE, 
#                         detrend_nb = FALSE,
#                         detrend_ar = FALSE)
# specs <- rbind(specs, new_row)



## Test alternate trimming or weighting
# Don't drop any outliers
params <- list(despline_data = fs_spl,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = FALSE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)

# Don't weight observations
params <- list(despline_data = fs_spl,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "equalweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)

# Drop outliers in X
params <- list(despline_data = fs_spl,
               mod_type = "fe",
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = TRUE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                      se = mod_slopes %>% pull("std.error"),
                      pet = TRUE,
                      temp = FALSE,
                      tmprl = FALSE,
                      trim_x = params$t_x,
                      trim_y = params$t_y,
                      weight_se = params$weights == "cwd_errorweights",
                      contemp_rwi = TRUE,
                      two_stage = TRUE,
                      single_stage = FALSE,
                      species_control = FALSE,
                      squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                      detrend_nb = FALSE,
                      detrend_ar = FALSE)
specs <- rbind(specs, new_row)

# ## Test model controlling for site-level SD of cwd and pet
# params <- list(despline_data = fs_spl,
#                formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2) + pet.sd + cwd.sd"),
#                t_x = FALSE,
#                t_y = TRUE,
#                weights = "cwd_errorweights")
# mod_slopes <- robustness_test(params)
# new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
#                       se = mod_slopes %>% pull("std.error"),
#                       trim_x = params$t_x,
#                       trim_y = params$t_y,
#                       weight_se = params$weights == "cwd_errorweights",
#                       species_control = FALSE,
#                       squared_term = TRUE,
#                       linear_mod = FALSE,
#                       detrend_spline = TRUE, 
#                       detrend_nb = FALSE,
#                       cum_dnlm = FALSE,
#                       detrend_ar = FALSE)
# specs <- rbind(specs, new_row)
# 
# 
# formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2) + cwd.sd + pet.sd")
# data <- fs_spl %>% filter(outlier == 0)
# mod <- feols(formula, weights = data$cwd_errorweights, data = data,
#              vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
# summary(mod)
# 
# fs_spl %>% ggplot(aes(x = cwd.spstd, y = cwd.sd)) +
#   geom_point() +
#   geom_smooth() +
#   xlim(-10,10)
# 
# 
# formula = as.formula("estimate_cwd.an ~ cwd.spstd + temp.spstd + (cwd.spstd^2) + (temp.spstd^2)")
# data <- fs_temp %>% filter(outlier == 0)
# mod <- feols(formula, weights = data$cwd_errorweights, data = data,
#              vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
# summary(mod)

## Create figure
highlight_n <- 1

specs <- specs %>% 
  select(-species_control) %>% 
  drop_na()

labels <- list("Detrending" = c("Spline", "Negative binomial", "Autoregressive"),
               "Energy control" = c("PET", "Mean temperature"),
               "First stage" = c("Annual RWI", "Cumulative\ndynamic lag"),
               "Model structure" = c("Two stage", "Single stage", "Squared\nCWD and PET", "Linear\nCWD and PET", "Temporally\ncorrelated error"),
               "Trimming and weighting" = c("Weight by\ninverse of s.e.", "Trim\noutliers in y", "Trim\noutliers in X"))


svg(paste0(wdir, 'figures\\a4_robustness.svg'), width = 9, height = 14)
specs <- specs %>% drop_na()
robustness_fig <- schart(specs, labels, highlight=highlight_n, order = "asis", 
                         heights = c(.4,.6),
                         n=c(1, 2, 1, 1, 3, 3), ci = c(.95), 
                         ylab = "Slope of line relating sites' historic\nCWD to CWD's impact on growth\n(evaluated at median historic CWD)",
                         col.est=c("grey80", "dodgerblue4"),
                         col.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         bg.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         lwd.symbol=1)
dev.off()

## Contrast Conley CI to bootstrapped confidence interval
bs_df <- read_rds(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))
beta_1 <- bs_df %>% pull(cwd_cwd)
beta_2 <- bs_df %>% pull(cwd_cwd2)
cwd_median <- fs_spl %>% 
  filter(outlier == 0) %>% 
  pull(cwd.spstd) %>% 
  median()
bs_df$margfx <- beta_1 + (2 * beta_2 * cwd_median)
bs_df$dummy = 1
bs_df %>% 
  ggplot(aes(x = dummy, y = margfx)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0, .1)) 

bs_df %>% 
  ggplot(aes(x = dummy, y = margfx)) +
  geom_violin()
  
  
### CIs for response to reviewers
cis <- specs %>% 
  mutate(ci = se * 1.96) 

# baseline
cis[1,]

# temp (Comment 1.5)
cis %>% 
  filter(temp==1) %>% 
  select(coef, ci)

# temporal autocorrelation (Comment 3.2)
cis %>% 
  filter(tmprl==1) %>% 
  select(coef, ci)

# Single integrated model (Comment 3.3)
cis %>% 
  filter(single_stage==1) %>% 
  select(coef, ci)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate variation by genus  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trim_df <- fs_spl %>% 
  filter(outlier==0) %>% 
  drop_na()

run_ss_conley <- function(data, outcome = "cwd_coef"){
  if (outcome == "cwd_coef") {
    error_weights = data$cwd_errorweights
  } else if (outcome == "pet_coef") {
    error_weights = data$pet_errorweights
  } else if (outcome == "int_coef") {
    error_weights = data$int_errorweights
  }
  formula <- as.formula(paste(outcome, " ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)"))
  mod <- feols(formula, data=data, weights = error_weights,
               vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
  
  # formula <- as.formula(paste(outcome, " ~ cwd.spstd + I(cwd.spstd**2) + pet.spstd + I(pet.spstd**2)"))
  # mod <- lm(formula, data=data, weights = error_weights)
  return(mod)
}

# run_ss_conley <- function(data, outcome = "estimate_cwd.an"){
#   formula <- as.formula(paste(outcome, " ~ cwd.spstd + pet.spstd")) ## ADD BACK WEIGHTING HERE?
#   mod <- feols(formula, data = data, weights = data$cwd_errorweights, 
#         vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
#   return(mod)
# }

genus_lookup <- trim_df %>% select(collection_id, genus)
genus_freq <- trim_df %>% 
  group_by(genus) %>% 
  summarise(n_collections = n_distinct(collection_id),
            min_cwd = min(cwd.spstd),
            max_cwd = max(cwd.spstd),
            range_cwd = max_cwd - min_cwd) %>% 
  arrange(desc(n_collections))

gymno_key <- sp_info %>% 
  select(genus, gymno_angio) %>% 
  unique()



# N = 10: 16 genera; 26: 12 genera; 50: 9 genera
genus_keep <- genus_freq %>% 
  filter(n_collections>50) %>%
  pull(genus)

genus_df <- trim_df %>% 
  mutate(int_coef = estimate_intercept,
         pet_coef = estimate_pet.an, 
         cwd_coef = estimate_cwd.an) %>% 
  filter(genus %in% genus_keep) %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(model_estimates = map(data, run_ss_conley))

genus_df$model_estimates

write_rds(genus_df, paste0(wdir, "out/second_stage/ss_conley_genus.rds"))


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Investigate variation by species  ----------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# sp_lookup <- trim_df %>% select(collection_id, species_id)
# sp_freq <- trim_df %>% 
#   group_by(species_id) %>% 
#   summarise(n_collections = n_distinct(collection_id),
#             min_cwd = min(cwd.spstd),
#             max_cwd = max(cwd.spstd),
#             range_cwd = max_cwd - min_cwd) %>% 
#   arrange(desc(n_collections))
# 
# gymno_key <- sp_info %>% 
#   select(species_id, gymno_angio) %>% 
#   unique()
# 
# sp_keep <- sp_freq %>% 
#   filter(n_collections>50) %>%
#   pull(species_id)
# 
# sp_df <- trim_df %>% 
#   mutate(int_coef = estimate_intercept,
#          pet_coef = estimate_pet.an, 
#          cwd_coef = estimate_cwd.an) %>% 
#   filter(species_id %in% sp_keep) %>% 
#   group_by(species_id) %>% 
#   nest() %>% 
#   mutate(model_estimates = map(data, run_ss_conley))
# 
# sp_df$model_estimates
# 
# write_rds(sp_df, paste0(wdir, "out/second_stage/ss_conley_species.rds"))
# 
# 
# 
# 
# ## REVIEW COMMENT 1.7 - assess association between gamma and species climate niche
# ## Seems like species in high PET 
# pull_me_median <- function(gen_mod, gen_dat){
#   median_cwd <- gen_dat %>% pull(cwd.spstd) %>% median()
#   lht <- linearHypothesis(gen_mod, c(paste0('cwd.spstd + ', as.character(median_cwd), ' * I(cwd.spstd^2) = 0')))
#   pvalue <- lht$`Pr(>Chisq)`[2] %>% round(digits = 3)
#   coefs = gen_mod$coefficients
#   me <- (coefs['cwd.spstd'] + median_cwd * 2 * coefs['I(cwd.spstd^2)']) %>% round(digits = 3)
#   return(me)
# }
# 
# pull_sig <- function(gen_mod, gen_dat){
#   median_cwd <- gen_dat %>% pull(cwd.spstd) %>% median()
#   lht <- linearHypothesis(gen_mod, c(paste0('cwd.spstd + ', as.character(median_cwd), ' * I(cwd.spstd^2) = 0')))
#   pvalue <- lht$`Pr(>Chisq)`[2] %>% round(digits = 3)
#   return(pvalue)
# }
# 
# me_median <- sp_df %>%
#   mutate(me_median = map2(model_estimates, data, pull_coefs)) %>%
#   mutate(me_sig = map2(model_estimates, data, pull_sig)) %>%
#   select(species_id, me_median, me_sig) %>%
#   unnest(me_median) %>%
#   unnest(me_sig) %>% 
#   arrange(species_id)
# 
# me_median <- me_median %>%
#   filter(me_sig < 0.05) %>%
#   left_join(niche_df, by = c("species_id" = "sp_code"))
# 
# 
# me_median %>% ggplot(aes(x = pet_mean, y = me_median)) + geom_point()
# me_median %>% ggplot(aes(x = cwd_mean, y = me_median)) + geom_point()
# 
# mod <- lm(me_median ~ cwd_mean + pet_mean, data = me_median)
# summary(mod)

# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Old version --------------------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
# trim_y <- list(TRUE, FALSE)
# trim_x <- list(TRUE, FALSE)
# 
# # cluster_level <- list(list(TRUE, "species_id"), list(FALSE, "collection_id"))
# cluster_level <- list(list(FALSE, "collection_id"))
# weighting <- list(list(TRUE, FALSE, "cwd_errorweights"), list(FALSE, FALSE, "equalweights"))
# #  list(FALSE, TRUE, "tree_errorweights"), 
# 
# # gymno_angio <- list(list(TRUE, FALSE, "gymno"), list(FALSE, TRUE, "angio"), list(TRUE, TRUE, c("gymno", "angio")))
# gymno_angio <- list(list(TRUE, TRUE, c("gymno", "angio")))
# equations <- list(list(FALSE, FALSE, "estimate_cwd.an ~ cwd.spstd + pet.spstd"), 
#                   list(FALSE, TRUE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"))
# # list(TRUE, FALSE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + factor(species_id)"),
# # list(TRUE, TRUE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd**2) + I(pet.spstd**2) + factor(species_id)"))
# detrend_method <- list(list(TRUE, FALSE, FALSE, fs_spl),
#                        list(FALSE, TRUE, FALSE, fs_nb),
#                        list(FALSE, FALSE, TRUE, fs_ar))
# 
# specs <- data.frame(coef=NaN, 
#                     se=NaN, 
#                     trim_x = NaN, trim_y = NaN, 
#                     weight_se = NaN, weight_n = NaN,
#                     include_gymno = NaN, include_angio = NaN,
#                     species_control = NaN, squared_term = NaN,
#                     species_cluster = NaN,
#                     detrend_spline = NaN,
#                     detrend_nb = NaN,
#                     detrend_ar = NaN)
# i = 1
# for (t_y in trim_y){
#   for (t_x in trim_x){
#     for(wght in weighting){
#       for (g_a in gymno_angio){
#         for (eqn in equations){
#           for(c_lvl in cluster_level){
#             for(dtrnd in detrend_method){
#               data <- dtrnd[[4]] %>% filter(gymno_angio %in% g_a[[3]])
#               cwd_median <- data$cwd.spstd %>% median()
#               if (t_x) {data <- data %>% filter(trim_x==TRUE)}
#               if (t_y) {data <- data %>% filter(trim_y==TRUE)}
#               formula <- eqn[[3]] %>% as.formula()
#               error_weights <- data[[wght[[3]]]]
#               # mod <- feols(formula, weights = data[[wght[[3]]]], data = data)
#               mod <- feols(formula, weights = error_weights, data = data,
#                            vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
#               # mod <- lm(eqn[[3]], weights = data[[wght[[3]]]], data=data)
#               # cluster_vcov <- vcovCL(mod, cluster = data[[c_lvl[[2]]]])
#               # mod <- coeftest(mod, vcov = vcovCL, cluster = data[[c_lvl[[2]]]])
#               print(mod)
#               mod_slopes <- slopes(mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = cwd_median))
#               mod_slopes <- mod_slopes %>% 
#                 as_tibble() %>% 
#                 filter(term == "cwd.spstd")
#               
#               specs[i,] <- data_frame(coef = mod_slopes %>% pull("estimate"),
#                                       se = mod_slopes %>% pull("std.error"),
#                                       trim_x = t_x,
#                                       trim_y = t_y,
#                                       weight_se = wght[[1]],
#                                       weight_n = wght[[2]],
#                                       include_gymno = g_a[[1]],
#                                       include_angio = g_a[[2]],
#                                       species_control = eqn[[1]],
#                                       squared_term = eqn[[2]],
#                                       species_cluster = c_lvl[[1]],
#                                       detrend_spline = dtrnd[[1]], 
#                                       detrend_nb = dtrnd[[2]],
#                                       detrend_ar = dtrnd[[3]])
#               i <- i+1
#               
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# 
# specs <- specs %>% drop_na()
# specs <- specs %>% 
#   select(-include_gymno, -include_angio)
# 
# specs <- specs %>%
#   select(-species_cluster)
# 
# highlight_n <- which(specs$trim_y == TRUE &
#                        specs$trim_x == FALSE &
#                        specs$weight_se == TRUE &
#                        specs$weight_n == FALSE &
#                        # specs$species_control == FALSE &
#                        specs$squared_term == TRUE &
#                        specs$detrend_spline == TRUE)
# 
# 
# schart(specs, highlight=highlight_n, order = "increasing", ci = c(.95, .99))
# 
# 
# labels <- list("Trimming:" = c("Trim outliers"),
#                "Weighting:" = c("Inverse s.e.", "Square root of n"),
#                "Controls" = c("Squared CWD and PET", "Species fixed effects"))
# 
# schart(specs, labels, highlight=highlight_n, order = "increasing", ci = c(.95, .99), adj=c(0,0), offset=c(5.5,5), ylab = "Slope of line relating\nsites' historic CWD to\nCWD's impact on growth")
# 
# 
# 
# despline_data <- fs_spl
# formula <- as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
# t_x <- TRUE
# t_y <- FALSE
# weights <- "cwd_errorweights"
# params <- list(despline_data = fs_spl,
#                formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
#                t_x = TRUE,
#                t_y = FALSE,
#                weights = "cwd_errorweights")


