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
source("spec_chart_function.R")

select <- dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote\\'

# 1. Site-level regressions
fs_spl <- read_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv'))
fs_nb <- read_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std_nb.csv'))
fs_ar <- read_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std_ar.csv'))

fs_cum <- read_rds(paste0(wdir, 'out/first_stage/dnlm_cum_effects.rds'))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_ave_clim.gz"))

# 3. Site information
site_df <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
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

fs_spl <- fs_spl %>% 
  process_fs()
fs_nb <- fs_nb %>% 
  process_fs()
fs_ar <- fs_ar %>% 
  process_fs()

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
# print(paste0("Range before hitting sill (km): "), vg.fit[2,3])

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

robustness_test <- function(params){
  data <- params$despline_data
  formula <- params$formula
  if (params$t_x) {data <- data %>% filter(trim_x==TRUE)}
  if (params$t_y) {data <- data %>% filter(trim_y==TRUE)}
  cwd_median <- data$cwd.spstd %>% median()
  error_weights <- data[[params$weights]]
  # mod <- feols(formula, weights = data[[wght[[3]]]], data = data)
  mod <- feols(formula, weights = error_weights, data = data,
               vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
  mod_slopes <- slopes(mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = cwd_median))
  mod_slopes <- mod_slopes %>% 
    as_tibble() %>% 
    filter(term == "cwd.spstd")
  return(mod_slopes)
}


# Parameters to vary
#  - Trim X
#  - Trim y
#  - Weighting
#  - Linear
#  - Species control
#  - Detrend method

specs <- data.frame(coef=NaN, 
                    se=NaN,
                    detrend_spline = NaN,
                    detrend_nb = NaN,
                    detrend_ar = NaN,
                    cum_dnlm = NaN,
                    squared_term = NaN,
                    linear_mod = NaN,
                    species_control = NaN, 
                    weight_se = NaN,
                    trim_y = NaN,
                    trim_x = NaN) 


## Baseline model
params <- list(despline_data = fs_spl,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = FALSE,
                        squared_term = TRUE,
                        linear_mod = FALSE,
                        detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        cum_dnlm = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)


## Test alternate desplining
params <- list(despline_data = fs_nb,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = FALSE, 
                        detrend_nb = TRUE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)


params <- list(despline_data = fs_ar,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = FALSE, 
                        detrend_nb = FALSE,
                        detrend_ar = TRUE)
specs <- rbind(specs, new_row)


## Cumulative impact from dynamic lag
params <- list(despline_data = fs_cum,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                      se = mod_slopes %>% pull("std.error"),
                      trim_x = params$t_x,
                      trim_y = params$t_y,
                      weight_se = params$weights == "cwd_errorweights",
                      species_control = FALSE,
                      squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = TRUE,
                      detrend_spline = TRUE, 
                      detrend_nb = FALSE,
                      detrend_ar = FALSE)
specs <- rbind(specs, new_row)


## Test alternate model structures
params <- list(despline_data = fs_ar,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = FALSE,
                        squared_term = FALSE,
                        linear_mod = TRUE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)



params <- list(despline_data = fs_nb,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2) | species_id"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = TRUE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)


## Test alternate trimming or weighting
params <- list(despline_data = fs_spl,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = TRUE,
               t_y = TRUE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)


params <- list(despline_data = fs_spl,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = FALSE,
               weights = "cwd_errorweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)


params <- list(despline_data = fs_spl,
               formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)"),
               t_x = FALSE,
               t_y = TRUE,
               weights = "equalweights")
mod_slopes <- robustness_test(params)
new_row <- data_frame(coef = mod_slopes %>% pull("estimate"),
                        se = mod_slopes %>% pull("std.error"),
                        trim_x = params$t_x,
                        trim_y = params$t_y,
                        weight_se = params$weights == "cwd_errorweights",
                        species_control = FALSE,
                        squared_term = TRUE,
                      linear_mod = FALSE,
                      cum_dnlm = FALSE,
                      detrend_spline = TRUE, 
                        detrend_nb = FALSE,
                        detrend_ar = FALSE)
specs <- rbind(specs, new_row)

specs <- specs %>% 
  drop_na()

## Create figure
highlight_n <- which(specs$trim_y == TRUE &
                       specs$trim_x == FALSE &
                       specs$weight_se == TRUE &
                       specs$species_control == FALSE &
                       specs$squared_term == TRUE &
                       specs$detrend_spline == TRUE &
                       specs$cum_dnlm == FALSE)


# chart <- schart(specs, highlight=highlight_n, ci = c(.95), n=c(1, 3,2,2))


labels <- list("Detrending" = c("Spline", "Negative binomial", "Autoregressive"),
               "First stage" = "Cumulative\ndynamic lag",
               "Model structure" = c("Squared\nCWD and PET", "Linear\nCWD and PET", "Species\nfixed effects"),
               "Trimming and\nweighting" = c("Weight by\ninverse of s.e.", "Trim\noutliers in y", "Trim\noutliers in X"))


svg(paste0(wdir, 'figures\\a4_robustness.svg'), width = 7, height = 11)

robustness_fig <- schart(specs, labels, highlight=highlight_n, order = "asis", 
                         n=c(1, 2, 1, 2,3), ci = c(.95), 
                         ylab = "Slope of line relating sites' historic\nCWD to CWD's impact on growth\n(evaluated at median historic CWD)",
                         col.est=c("grey80", "dodgerblue4"),
                         col.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         bg.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         lwd.symbol=1)
dev.off()

# adj=c(0,0), offset=c(5.5,5), 




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





