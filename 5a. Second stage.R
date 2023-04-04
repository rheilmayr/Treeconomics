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

set.seed(5597)

select <- dplyr::select

n_mc <- 10000


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote\\'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv'))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_ave_clim.gz"))
flm_df <- flm_df %>% 
  left_join(ave_site_clim, by = c("collection_id"))

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
#            (cwd.spstd<cwd_spstd_bounds[1]) |
#            (cwd.spstd>cwd_spstd_bounds[2]) |
#            (pet.spstd<pet_spstd_bounds[1]) |
#            (pet.spstd>pet_spstd_bounds[2]))
           

flm_df <- flm_df %>% 
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
                   (estimate_cwd.an>cwd_est_bounds[2]) |
                   (estimate_pet.an<pet_est_bounds[1]) |
                   (estimate_pet.an>pet_est_bounds[2]))

# Save out full flm_df to simplify downstream scripts and ensure consistency
flm_df %>% write.csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

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
print(paste0("Range before hitting sill (km): "), vg.fit[2,3])

vg.range = vg.fit[2,3] * 1000



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Identify spatially proximate blocks of sites ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_list <- trim_df %>%
  pull(collection_id) %>%
  unique()
n_sites <- length(site_list)

site_dist=st_distance(site_points)
rownames(site_dist)=site_points$collection_id
colnames(site_dist)=site_points$collection_id
# save(site_dist,file=paste0(wdir,"out/site_distances.Rdat"))
# load(paste0(wdir,"out/site_distances.Rdat"))

dist_df <- as_tibble(site_dist) %>% 
  drop_units() 

dist_df <- dist_df %>%
  lazy_dt() %>% 
  mutate(collection_id = names(dist_df)) %>% 
  # select(collection_id, site_list) %>% 
  # filter(collection_id %in% site_list) %>% 
  mutate(across(.cols = !collection_id, ~(.x < vg.range))) %>% 
  # mutate(across(.cols = !collection_id, ~ifelse((.x < range), collection_id, "DROP"))) %>% 
  as_tibble()

block_list <- c()
for (site in site_list){
  block_sites <- dist_df %>% 
    filter(get(site) == TRUE) %>% 
    pull(collection_id)
  block_list[site] <- list(block_sites)
}
save(block_list,file=paste0(wdir,"out/spatial_blocks.Rdat"))
load(file=paste0(wdir,"out/spatial_blocks.Rdat"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create block bootstrap draws  ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
draw_blocks <- function(site_sample){
  samp <- site_sample %>% pull(samp)
  blocked_draw <- (block_list[samp] %>% unlist() %>% unname())[1:n_sites]
  return(blocked_draw)  
}

n_obs = n_mc * n_sites
block_draw_df <- tibble(boot_id = rep(1:n_mc, each = n_sites)) %>% 
  lazy_dt() %>% 
  mutate(samp = sample(site_list,size=n_obs,replace=TRUE)) %>%
  group_by(boot_id) %>% 
  nest()

block_draw_df <- block_draw_df %>% 
  mutate(sites = map(.x = data, .f = draw_blocks)) %>% # COULD PARALLELIZE HERE?
  select(boot_id, sites) %>% 
  as_tibble() %>% 
  unnest(sites) 

block_draw_df <- block_draw_df %>% 
  rename(collection_id = sites)

## Identify number of draws needed for each site
n_draws <- block_draw_df %>% 
  group_by(collection_id) %>% 
  tally() %>% 
  rename(n_draw = n)

trim_df <- trim_df %>% 
  left_join(n_draws, by = "collection_id")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Random draws of coefs from first stage ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Function to create random draws of first stage coefficients
draw_coefs <- function(n, cwd_est, pet_est, int_est, cwd_ste, pet_ste, int_ste, 
                       cwd_pet_cov, pet_int_cov, cwd_int_cov){
  mu <- c("cwd_coef" = cwd_est, 
          "pet_coef" = pet_est, 
          "int_coef" = int_est)
  vcov <- matrix(data = NA, nrow = 3, ncol = 3)
  vcov[1,1] <- cwd_ste^2
  vcov[2,2] <- pet_ste^2
  vcov[3,3] <- int_ste^2
  vcov[1,2] <- cwd_pet_cov
  vcov[2,1] <- cwd_pet_cov
  vcov[1,3] <- cwd_int_cov
  vcov[3,1] <- cwd_int_cov
  vcov[2,3] <- pet_int_cov
  vcov[3,2] <- pet_int_cov
  
  draw <- mvrnorm(n, mu, vcov)
  draw <- as_tibble(draw)
  draw$iter_idx <- seq(1,n)
  draw <- draw %>% select(iter_idx, cwd_coef, pet_coef, int_coef)
  return(draw)
}

## Create needed number (n_draw) of random draws of first stage coefficients for each site
trim_df <- trim_df %>% 
  drop_na()

mc_df <- trim_df %>%
  mutate(coef_draws = pmap(list(n = trim_df$n_draw + 1, 
                                cwd_est = trim_df$estimate_cwd.an, 
                                pet_est = trim_df$estimate_pet.an,
                                int_est = trim_df$estimate_intercept, 
                                cwd_ste = trim_df$std.error_cwd.an,
                                pet_ste = trim_df$std.error_pet.an, 
                                int_ste = trim_df$std.error_intercept,
                                cwd_pet_cov = trim_df$cov_cwd_pet, 
                                cwd_int_cov = trim_df$cov_int_cwd,
                                pet_int_cov = trim_df$cov_int_pet), 
                           draw_coefs))


## Unnest to create dataframe of n_site X n_draw coefficient estimates
mc_df <- mc_df %>% 
  unnest(coef_draws) %>% 
  select(collection_id, iter_idx, cwd_coef, pet_coef, int_coef, cwd.spstd, 
         pet.spstd, latitude, longitude, cwd_errorweights, pet_errorweights, int_errorweights)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge first stage draws back to bootstrap dataframe -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
block_draw_df <- block_draw_df %>% 
  group_by(collection_id) %>% 
  mutate(iter_idx = 1:n())

block_draw_df <- block_draw_df %>% 
  left_join(mc_df, by = c("collection_id", "iter_idx"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export first stage draws to pull summary stats -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
block_draw_df %>% 
  # select(boot_id, collection_id, cwd_coef, pet_coef, int_coef, cwd.spstd, pet.spstd) %>% 
  write_rds(paste0(wdir, "out/second_stage/mc_sample.gz"), compress = "gz")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run bootstrap estimation of second stage model -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Function defining main second stage model
run_ss <- function(data, outcome = "cwd_coef"){
  if (outcome == "cwd_coef") {
    error_weights = data$cwd_errorweights
  } else if (outcome == "pet_coef") {
    error_weights = data$pet_errorweights
  } else if (outcome == "int_coef") {
    error_weights = data$int_errorweights
  }
  formula <- as.formula(paste(outcome, " ~ cwd.spstd + I(cwd.spstd^2) + pet.spstd + I(pet.spstd^2)"))
  mod <- lm(formula, data=data, weights = error_weights)
  # mod <- lm(formula, data=data)
  coefs <- mod %>% 
    tidy() %>% 
    pull(estimate) 
  return(coefs)
}


## Function to run second stage model for first stage CWD, 
## PET and Intercept terms and organize coefficients
bs_ss <- function(data){
  # data <- data %>% # Needed to add this since block bootstrap is returning nested tibble
  #   unnest(cols = c(data))
  cwd_mod <- data %>% 
    run_ss(outcome = "cwd_coef")
  pet_mod <- data %>% 
    run_ss(outcome = "pet_coef")
  int_mod <- data %>% 
    run_ss(outcome = "int_coef")
  # return(list("int_mod" = c(int_mod),
  #             "cwd_mod" = c(cwd_mod),
  #             "pet_mot" = c(pet_mod)))
  return(list(int_int = int_mod[1],
              int_cwd = int_mod[2],
              int_cwd2 = int_mod[3],
              int_pet = int_mod[4],
              int_pet2 = int_mod[5],
              cwd_int = cwd_mod[1],
              cwd_cwd = cwd_mod[2],
              cwd_cwd2 = cwd_mod[3],
              cwd_pet = cwd_mod[4],
              cwd_pet2 = cwd_mod[5],
              pet_int = pet_mod[1],
              pet_cwd = pet_mod[2],
              pet_cwd2 = pet_mod[3],
              pet_pet = pet_mod[4],
              pet_pet2 = pet_mod[5]))
}

bs_constant <- function(data){
  # data <- data %>% # Needed to add this since block bootstrap is returning nested tibble
  #   unnest(cols = c(data))
  const_sens <- data %>% 
    summarise(cwd_const_sens = weighted.mean(cwd_coef, cwd_errorweights),
          pet_const_sens = weighted.mean(pet_coef, pet_errorweights),
          int_const_sens = weighted.mean(int_coef, pet_errorweights))
    # summarise(cwd_const_sens = mean(cwd_coef),
    #           pet_const_sens = mean(pet_coef),
    #           int_const_sens = mean(pet_coef))
    
  return(const_sens)
}


## Create dataframe holding bootstrap samples
boot_df <- block_draw_df %>%
  select(-iter_idx) %>% 
  group_by(boot_id) %>% 
  nest()

# # n_sites <- boot_df %>% pull(collection_id) %>% unique() %>% length()
# # n_obs = n_mc * n_sites
# 
# 
# 
# # boot_df <- boot_df %>% 
# #   ## TODO: Here's where we probably want to introduce block bootstrapping design...
# #   sample_n(size = n_obs, replace = TRUE, weight = cwd_errorweights) %>% 
# #   ## TODO: Should weighting be done separately for CWD, PET and Intercept terms? 
# #   ## Each has a different standard error in first stage...
# #   mutate(iter_idx = rep(1:n_mc, each=n_sites)) %>% 
# #   relocate(iter_idx) %>% 
# #   group_by(iter_idx) %>% 
# #   nest() 
#   ## Note: goal at end of bootstrap is to have a dataframe with 10,000 nested datasets 
#   ## that can be used for second stage model...
# 
# source("block_bootstrapping.R")
# ncores=detectCores()-2
# my.cluster=makeCluster(ncores);registerDoParallel(cl=my.cluster)
# boot_df=blockbootstrap_func(boot_df)
# saveRDS(boot_df,file=paste0(wdir,"out/blockbootstrap_samples.rds"))
# 
# boot_df <- readRDS(paste0(wdir,"out/blockbootstrap_samples.rds"))


## Estimate second stage models
boot_df <- boot_df %>% 
  mutate(const_sens = map(.x = data, .f = bs_constant)) %>%
  unnest_wider(const_sens) %>%
  mutate(estimates = map(.x = data, .f = bs_ss)) %>% 
  unnest_wider(estimates) %>% 
  select(-data) %>% 
  ungroup()


## Summarize bootstrap results
bs_coefs <- apply(boot_df %>% select(-boot_id), 2, mean) %>% print()
bs_ste <- apply(boot_df %>% select(-boot_id), 2, sd) %>% print()


## Save out bootstrapped coefficients that reflect uncertainty from both first and second stage models
write_rds(boot_df, paste0(wdir, "out/second_stage/ss_bootstrap.rds"))









# ## TODO: Should we be worried that bootstrap se very similar (but smaller?) than simple regression result?
# ## Bootstrap is currently returning s.e. of 0.001849, compare to 0.002054 for simple regression:
# mod <- feols(estimate_cwd.an ~ cwd.spstd + pet.spstd, 
#              weights = trim_df$cwd_errorweights, data = trim_df)
# summary(mod) 
# 
# library(marginaleffects)
# mod_df <- trim_df
# mod_df <- mod_df %>% 
#   mutate(int_coef = estimate_intercept,
#          pet_coef = estimate_pet.an, 
#          cwd_coef = estimate_cwd.an)
# 
# run_ss(mod_df, outcome = "cwd_coef")
# mod <- feols(estimate_cwd.an ~ cwd.spstd + cwd.spstd**2 + pet.spstd + pet.spstd**2, 
#              weights = mod_df$cwd_errorweights, data = mod_df)
# summary(mod) 
# plot_cap(mod, condition = "pet.spstd") + xlim(-3,3) + ylim(-.2, 0)
# plot_cap(mod, condition = "cwd.spstd") + xlim(-3,3) + ylim(-.2, 0)
# 
# newdata <- tibble("cwd.spstd" = 1, 
#                   "pet.spstd" = 0)
# predict(mod, newdata)
# 
# 
# mod <- feols(estimate_cwd.an ~ 1, 
#              weights = mod_df$cwd_errorweights, data = mod_df)
# mod_df %>% summarise(test = weighted.mean(estimate_cwd.an, cwd_errorweights))
# 
# 
# mod <- feols(estimate_pet.an ~ cwd.spstd + cwd.spstd**2 + pet.spstd + pet.spstd**2,
#              weights = mod_df$pet_errorweights, data = mod_df)
# summary(mod) 
# plot_cap(mod, condition = "pet.spstd") + xlim(-3,3) + ylim(-.1, .3)
# plot_cap(mod, condition = "cwd.spstd") + xlim(-3,3) + ylim(-.2, .2)
# 
# newdata <- tibble("cwd.spstd" = 0, 
#                   "pet.spstd" = 0)
# predict(mod, newdata)
# 
# 
# 
# mod <- feols(estimate_cwd.an ~ cwd.spstd + cwd.spstd**2 + pet.spstd + pet.spstd**2, 
#              data = mod_df)
# summary(mod) 
# plot_cap(mod, condition = "pet.spstd") + xlim(-3,3) + ylim(-.4, .1)
# plot_cap(mod, condition = "cwd.spstd") + xlim(-3,3) + ylim(-.4, .1)
# 
# newdata <- tibble("cwd.spstd" = 0, 
#                   "pet.spstd" = 1)
# predict(mod, newdata)
# 
# 
# mod <- feols(estimate_pet.an ~ cwd.spstd + cwd.spstd**2 + pet.spstd + pet.spstd**2,
#              data = mod_df)
# summary(mod) 
# plot_cap(mod, condition = "pet.spstd") + xlim(-3,3) + ylim(-.4, .5)
# plot_cap(mod, condition = "cwd.spstd") + xlim(-3,3) + ylim(-.4, .5)
# 
# newdata <- tibble("cwd.spstd" = -2, 
#                   "pet.spstd" = -2)
# predict(mod, newdata)
# 
# 
# 
# mod <- feols(estimate_intercept ~ cwd.spstd + cwd.spstd**2 + pet.spstd + pet.spstd**2, 
#              weights = mod_df$int_errorweights, data = mod_df)
# summary(mod) 
# plot_cap(mod, condition = "pet.spstd") + xlim(-3,3) + ylim(.6,1.3)
# plot_cap(mod, condition = "cwd.spstd") + xlim(-3,3) + ylim(.6,1.3)
# mod_df %>% summarise(test = weighted.mean(estimate_intercept, int_errorweights))
# 
# ## OLD VERSION OF BOOTSTRAPPING USING BOOT PACKAGE
# # mod_df <- trim_df %>% filter(genus == "Araucaria")
# # mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data = mod_df)
# # summary(mod) 
# mod <- feols(estimate_cwd.an ~ cwd.spstd + pet.spstd,
#              weights = trim_df$cwd_errorweights, data = trim_df,
#              vcov = conley(cutoff = 363, distance = "spherical"))
# summary(mod)
# 
# mod <- feols(estimate_pet.an ~ cwd.spstd + pet.spstd,
#              weights = trim_df$cwd_errorweights, data = trim_df,
#              vcov = conley(cutoff = 363, distance = "spherical"))
# summary(mod)
# 
# 
# # 
# # 
# # 
# # mod <- feols(estimate_cwd.an ~ cwd.spstd + pet.spstd, data = trim_df,
# #           vcov = conley(distance = "spherical"))
# # summary(mod) ## TODO: Is bootstrap se too similar to simple regression result? Revisit...
# # 
# # 
# # boot_df <- tibble(idx = seq(1, n_obs, 1), n_samples = n_sites) %>% 
# #   mutate(data = pmap(.l = list(tbl = mc_df, 
# #                                replace = TRUE, 
# #                                size = n_samples), 
# #                      .f = sample_n ))
# # 
# # apply_boot <- function(mc_df){
# #   n_obs <- dim(mc_df)[1]
# #   n_sites <- trim_df %>% pull(collection_id) %>% unique() %>% length()
# #   
# #   
# #   boot_ss_coefs <- boot(mc_df, bs_ss, R = n_mc,
# #                         sim = "parametric",  
# #                         ran.gen=function(d,p) d[sample(1:n_obs, n_sites), ])
# #   # boot_ss_coefs <- boot_ss_coefs$t
# #   # colnames(boot_ss_coefs) <- c("int_int", "int_cwd", "int_pet",
# #   #                              "cwd_int", "cwd_cwd", "cwd_pet",
# #   #                              "pet_int", "pet_cwd", "pet_pet")
# #   # boot_ss_coefs <- boot_ss_coefs %>% 
# #   #   as_tibble() %>% 
# #   #   mutate(iter_idx = seq(1:n_mc))%>% 
# #   #   relocate(iter_idx)
# #   
# #   return(boot_ss_coefs)
# # }
# # 
# # boot_ss_coefs <- apply_boot(mc_df)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate variation by genus  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



# summarize_models <- function(mod, dat) {
#   median_cwd <- dat %>% pull(cwd.spstd) %>% median()
#   lht <- linearHypothesis(mod, c(paste0('cwd.spstd + ', as.character(median_cwd), ' * I(cwd.spstd^2) = 0')))
#   pvalue <- lht$`Pr(>Chisq)`[2]
#   coefs = mod$coefficients
#   me <- coefs['cwd.spstd'] + median_cwd * 2 * coefs['I(cwd.spstd^2)']
#   return(tibble("median_cwd" = median_cwd, "pvalue" = pvalue, "me" = me))
# }
# 
# genus_df <- genus_df %>% 
#   mutate(test_stats = map2(model_estimates, data, .f = summarize_models))
# genus_df <- genus_df %>% 
#   unnest(test_stats)





# 
# genus_id = "Juniperus"
# plot_df <- genus_df %>% 
#   filter(genus == genus_id) %>% 
#   pull(data)
# 
# plot_df = plot_df[[1]] %>% 
#   mutate(predicted = estimate_cwd.an)
# # plot_df %>% 
# #   ggplot() +
# #   geom_point(data = plot_df %>% filter(genus == genus_id), aes(x = cwd.spstd, y = estimate_cwd.an) +
# #   xlim(-2.5, 2.5)
# # 
# 
# 
# gen_plot <- ggplot() +
#   geom_line(data = genus_predictions %>% filter(genus == genus_id), aes(x = cwd.spstd, y = predicted)) +
#   geom_ribbon(data = genus_predictions  %>% filter(genus == genus_id), aes(x = cwd.spstd, ymin=conf.low, ymax=conf.high), alpha=0.2) +
#   geom_point(data = plot_df, aes(x = cwd.spstd, y = predicted)) +
#   # geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkblue") +
#   theme_bw(base_size = 22) + 
#   # theme(
#   #   # panel.grid.major = element_blank(),
#   #   # panel.grid.minor = element_blank(),
#   #   strip.background = element_blank(),
#   #   panel.border = element_rect(colour = "black", fill = NA)) +
#   # geom_line(aes(y = conf.low), linetype = 3) +
#   # geom_line(aes(y = conf.high), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") +
#   xlim(c(-3, 3))
# 
# 
# gen_plot


# mc_df <- mc_df %>% 
#   left_join(genus_lookup, by = "collection_id")
# 
# genus_df <- mc_df %>% 
#   filter(genus %in% genus_keep) %>% 
#   group_by(genus) %>% 
#   nest() %>% 
#   mutate(model_estimates = map(data, apply_boot))
# 
# genus_df <- genus_df %>% 
#   left_join(genus_freq, by = "genus") %>% 
#   left_join(genus_key, by = "genus")
# 
# genus_df <- genus_df %>% 
#   select(-data) %>% 
#   unnest(model_estimates)
# 




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Primary second stage model --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run regression and cluster s.e.

# trim_df <- trim_df %>%
#   filter(gymno_angio=="gymno")

mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = cwd_errorweights, data=flm_df)
cluster_vcov <- vcovCL(mod, cluster = flm_df$collection_id)
coeftest(mod, vcov = vcovCL, cluster = flm_df$collection_id)

mod_df = trim_df
mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = cwd_errorweights, data=mod_df)
cluster_vcov <- vcovCL(mod, cluster = mod_df$species_id)
coeftest(mod, vcov = vcovCL, cluster = mod_df$species_id)
saveRDS(mod, paste0(wdir, "out\\second_stage\\cwd_mod.rds"))
saveRDS(cluster_vcov, paste0(wdir, "out\\second_stage\\cwd_mod_vcov.rds"))
# saveRDS(sq_pet_mod, paste0(wdir, "out\\second_stage\\ss_sq_pet_mod.rds"))

pet_mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd, weights = pet_errorweights, data=mod_df)
pet_cluster_vcov <- vcovCL(pet_mod, cluster = mod_df$collection_id)
coeftest(pet_mod, cluster = mod_df$collection_id)
saveRDS(pet_mod, paste0(wdir, "out\\second_stage\\pet_mod.rds"))
saveRDS(pet_cluster_vcov, paste0(wdir, "out\\second_stage\\pet_mod_vcov.rds"))


int_mod <- lm(estimate_intercept ~ cwd.spstd + pet.spstd, weights = int_errorweights, data = mod_df)
int_cluster_vcov <- vcovCL(int_mod, cluster = mod_df$collection_id)
coeftest(int_mod, cluster = mod_df$collection_id)
saveRDS(int_mod, paste0(wdir, "out\\second_stage\\int_mod.rds"))
saveRDS(int_cluster_vcov, paste0(wdir, "out\\second_stage\\int_mod_vcov.rds"))


# Squared terms
sq_cwd_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd^2) + I(pet.spstd^2), weights = errorweights, data=mod_df)
coeftest(sq_cwd_mod, vcov = vcovCL, cluster = trim_df$collection_id)
sq_cwd_cluster_vcov <- vcovCL(sq_cwd_mod, cluster = trim_df$collection_id)
saveRDS(sq_cwd_mod, paste0(wdir, "out\\second_stage\\sq_cwd_mod.rds"))
saveRDS(sq_cwd_cluster_vcov, paste0(wdir, "out\\second_stage\\sq_cwd_mod_vcov.rds"))

sq_pet_mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd + I(cwd.spstd^2) + I(pet.spstd^2), weights = pet_errorweights, data=mod_df)
sq_pet_cluster_vcov <- vcovCL(sq_pet_mod, cluster = trim_df$collection_id)
saveRDS(sq_pet_mod, paste0(wdir, "out\\second_stage\\sq_pet_mod.rds"))
saveRDS(sq_pet_cluster_vcov, paste0(wdir, "out\\second_stage\\sq_pet_mod_vcov.rds"))

sq_int_mod <- lm(estimate_intercept ~ cwd.spstd + pet.spstd + I(cwd.spstd^2) + I(pet.spstd^2), weights = int_errorweights, data=mod_df)
sq_int_cluster_vcov <- vcovCL(sq_int_mod, cluster = trim_df$collection_id)
saveRDS(sq_int_mod, paste0(wdir, "out\\second_stage\\sq_int_mod.rds"))
saveRDS(sq_int_cluster_vcov, paste0(wdir, "out\\second_stage\\sq_int_mod_vcov.rds"))


mod_df %>%
  mutate(pred_rwi = estimate_intercept + (estimate_pet.an * pet.spstd) + (estimate_cwd.an * cwd.spstd)) %>%
  select(pred_rwi) %>%
  summary() #CAUTION: Shouldn't mean here be much closer to 1? Investigate in first stage script....







# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Margins plots --------------------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# trim_df <- trim_df %>% 
#   mutate(significant =   p.value_cwd.an<0.05)
# trim_df %>% 
#   ggplot(aes(x = cwd.spstd, y = estimate_cwd.an, color = std.error_cwd.an)) +
#   stat_smooth(method = "loess") +
#   theme_bw()
# 
# xmin = -2
# xmax = 2
# sq_predictions <- prediction(sq_mod, at = list(cwd.spstd = seq(xmin, xmax, .1)), vcov = sq_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(cwd.spstd = "at(cwd.spstd)")
# margins_plot <- ggplot(sq_predictions, aes(x = cwd.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto CWD") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot
# 
# 
# sq_predictions <- prediction(sq_mod, at = list(pet.spstd = seq(xmin, xmax, .1)), vcov = sq_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(pet.spstd = "at(pet.spstd)")
# margins_plot <- ggplot(sq_predictions, aes(x = pet.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic PET\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto CWD") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot
# 
# sq_pet_predictions <- prediction(sq_pet_mod, at = list(pet.spstd = seq(xmin, xmax, .1)), vcov = sq_pet_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(pet.spstd = "at(pet.spstd)")
# 
# margins_plot <- ggplot(sq_pet_predictions, aes(x = pet.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic PET\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto PET") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot
# 
# sq_pet_predictions <- prediction(sq_pet_mod, at = list(cwd.spstd = seq(xmin, xmax, .1)), vcov = sq_pet_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(cwd.spstd = "at(cwd.spstd)")
# 
# margins_plot <- ggplot(sq_pet_predictions, aes(x = cwd.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto PET") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Robustness tests --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# No weights
flm_df <- flm_df %>% 
  mutate(equal_weights = 1)

# Trimming
flm_df <- flm_df %>% 
  mutate(trim_y = estimate_cwd.an>cwd_est_bounds[1] & estimate_cwd.an<cwd_est_bounds[2],
         trim_x = cwd.spstd>cwd_spstd_bounds[1] & cwd.spstd<cwd_spstd_bounds[2] & 
                  pet.spstd>pet_spstd_bounds[1] & pet.spstd<pet_spstd_bounds[2])

trim_y <- list(TRUE)
trim_x <- list(TRUE, FALSE)
# cluster_level <- list(list(TRUE, "species_id"), list(FALSE, "collection_id"))
cluster_level <- list(list(FALSE, "collection_id"))
weighting <- list(list(TRUE, FALSE, "errorweights"), list(FALSE, TRUE, "errorweights2"), list(FALSE, FALSE, "equalweights"))
# gymno_angio <- list(list(TRUE, FALSE, "gymno"), list(FALSE, TRUE, "angio"), list(TRUE, TRUE, c("gymno", "angio")))
gymno_angio <- list(list(TRUE, TRUE, c("gymno", "angio")))
equations <- list(list(FALSE, FALSE, "estimate_cwd.an ~ cwd.spstd + pet.spstd"), 
                  list(FALSE, TRUE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd**2) + I(pet.spstd**2)"), 
                  list(TRUE, FALSE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + factor(species_id)"),
                  list(TRUE, TRUE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd**2) + I(pet.spstd**2) + factor(species_id)"))


specs <- data.frame(coef=NaN, 
                    se=NaN, 
                    trim_x = NaN, trim_y = NaN, 
                    weight_se = NaN, weight_n = NaN,
                    include_gymno = NaN, include_angio = NaN,
                    species_control = NaN, squared_term = NaN,
                    species_cluster = NaN)
i = 1
for (t_y in trim_y){
  for (t_x in trim_x){
    for(wght in weighting){
      for (g_a in gymno_angio){
        for (eqn in equations){
          for(c_lvl in cluster_level){
            data <- flm_df %>% filter(gymno_angio %in% g_a[[3]])
            if (t_x) {data <- data %>% filter(trim_x==TRUE)}
            if (t_y) {data <- data %>% filter(trim_y==TRUE)}
            mod <- lm(eqn[[3]], weights = data[[wght[[3]]]], data=data)
            cluster_vcov <- vcovCL(mod, cluster = data[[c_lvl[[2]]]])
            mod <- coeftest(mod, vcov = vcovCL, cluster = data[[c_lvl[[2]]]])
            specs[i,] <- data_frame(coef = mod[2,1],
                                    se = mod[2,2],
                                    trim_x = t_x,
                                    trim_y = t_y,
                                    weight_se = wght[[1]],
                                    weight_n = wght[[2]],
                                    include_gymno = g_a[[1]],
                                    include_angio = g_a[[2]],
                                    species_control = eqn[[1]],
                                    squared_term = eqn[[2]],
                                    species_cluster = c_lvl[[1]])
            i <- i+1
          }
        }
      }
    }
  }
}

specs <- specs %>% drop_na()
specs <- specs %>% 
  select(-include_gymno, -include_angio, -trim_y)

specs <- specs %>% 
  select(-species_cluster)

highlight_n <- which(specs$trim_x == TRUE  &
                     # specs$species_cluster == FALSE &
                     # specs$trim_y == TRUE &
                     specs$weight_se == TRUE &
                     specs$weight_n == FALSE &
                     specs$species_control == FALSE &
                     specs$squared_term == FALSE)


schart(specs, highlight=highlight_n, order = "increasing", ci = c(.95, .99))


labels <- list("Trimming:" = c("Trim outliers"),
               "Weighting:" = c("Inverse s.e.", "Square root of n"),
               "Controls" = c("Squared CWD and PET", "Species fixed effects"))

schart(specs, labels, highlight=highlight_n, order = "increasing", ci = c(.95, .99), adj=c(0,0), offset=c(5.5,5), ylab = "Slope of line relating\nsites' historic CWD to\nCWD's impact on growth")







# 
# sq_all_predictions = sq_pet_predictions %>% rename(cwd.spstd = pet.spstd) + sq_predictions
# 
# flm_df %>% filter(species_id == "tsme") %>%  ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) + geom_point() + ylim(c(-.5,.5))





# Drop PET
nopet_mod <- lm(estimate_cwd.an ~ cwd.spstd, weights = errorweights, data=trim_df)
nopet_cluster_vcov <- vcovCL(nopet_mod, cluster = trim_df$collection_id)
coeftest(nopet_mod, vcov = nopet_cluster_vcov)

# Don't trim dataset
allobs_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=flm_df)
allobs_cluster_vcov <- vcovCL(allobs_mod, cluster = flm_df$collection_id)
coeftest(allobs_mod, vcov = allobs_cluster_vcov)

# Limit to gymnosperms
subset_df <- trim_df %>% filter(gymno_angio=="gymno")
# subset_df <- trim_df %>% filter(gymno_angio=="gymno") %>% filter(collection_id %in% old_sites)
subset_df <- trim_df %>% filter(genus=="Pinus")
# subset_df <- trim_df %>% filter(genus=="Pinus") %>% filter(collection_id %in% old_sites)
subset_df %>% dim()
subset_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=subset_df)
subset_cluster_vcov <- vcovCL(subset_mod, cluster = subset_df$collection_id)
coeftest(subset_mod, vcov = subset_cluster_vcov)


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Investigate mechanism - tree or site --------------------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Focus on trees from genera with effect, and for which we have full weather history
# young_df <- trim_df %>%
#   filter(young==TRUE) %>%
#   drop_na()
# 
# mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trstd + pet.spstd + pet.trstd, weights = errorweights, data=young_df)
# cluster_vcov <- vcovCL(mod, cluster = young_df$collection_id)
# coeftest(mod, vcov = vcovCL, cluster = young_df$collection_id)
# 
# # Could include all trees from sites with any variation in tree-level drought history (including trees born <1901)
# young_sites <- trim_df %>%
#   filter(young == T) %>%
#   pull(collection_id) %>%
#   unique()
# youngsite_df <- trim_df %>%
#   filter(collection_id %in% young_sites) %>%
#   drop_na()
# mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trstd + pet.spstd + pet.trstd, weights = errorweights, data=youngsite_df)
# coeftest(mod, vcov = vcovCL, cluster = youngsite_df$collection_id)
# 
# # Include collection fixed effects
# mod <- feols(estimate_cwd.an ~ cwd.trstd + pet.trstd | collection_id, weights = trim_df$errorweights, data=trim_df)
# summary(mod)
# 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate variation by genus - old version  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# predict_cwd_effect <- function(trim_df){
#   trim_df <- trim_df %>%
#     filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
#            abs(cwd.spstd)<5,
#            abs(pet.spstd)<5) %>%
#     drop_na()
#   gen_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=trim_df)
#   cluster_vcov <- vcovCL(gen_mod, cluster = trim_df$collection_id)
#   print(coeftest(gen_mod, vcov = vcovCL, cluster = trim_df$collection_id))
#   predictions <- prediction(gen_mod, at = list(cwd.spstd = seq(-3, 3, .5)), vcov = cluster_vcov, calculate_se = T) %>% 
#     summary() %>% 
#     rename(cwd.spstd = "at(cwd.spstd)") 
#   return(predictions)
# }
# # ## Plot marginal effects by gymno / angio
# grp_freq <- flm_df %>% 
#   group_by(gymno_angio) %>% 
#   summarise(n_collections = n_distinct(collection_id)) %>% 
#   arrange(desc(n_collections))
# 
# grp_keep <- grp_freq %>% 
#   filter(n_collections>50) %>% 
#   pull(gymno_angio)
# 
# grp_df <- flm_df %>%
#   filter(gymno_angio %in% grp_keep) %>% 
#   group_by(gymno_angio) %>% 
#   nest() %>%
#   mutate(prediction = map(data, predict_cwd_effect))
# 
# grp_df <- grp_df %>% 
#   unnest(prediction) %>% 
#   select(-data)
# 
# grp_df$gymno_angio <- grp_df$gymno_angio %>% recode("angio" = "Angiosperm", "gymno" = "Gymnosperm")
# 
# margins_plot <- ggplot(grp_df, aes(x = cwd.spstd)) + 
#   geom_line(aes(y = Prediction)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   facet_wrap(~gymno_angio, scales = "free", ncol = 1) +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") + 
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::comma)
# margins_plot
# ggsave(paste0(wdir, 'figures\\gymno_angio.svg'), plot = margins_plot)
# 
# ## Plot marginal effects by gymno / angio
# grp_freq <- flm_df %>% 
#   group_by(gymno_angio) %>% 
#   summarise(n_collections = n_distinct(collection_id)) %>% 
#   arrange(desc(n_collections))
# 
# grp_keep <- grp_freq %>% 
#   filter(n_collections>50) %>% 
#   pull(gymno_angio)
# 
# grp_df <- flm_df %>%
#   filter(gymno_angio %in% grp_keep) %>% 
#   group_by(gymno_angio) %>% 
#   nest() %>%
#   mutate(prediction = map(data, predict_cwd_effect))
# 
# grp_df <- grp_df %>% 
#   unnest(prediction) %>% 
#   select(-data)
# 
# grp_df$gymno_angio <- grp_df$gymno_angio %>% recode("angio" = "Angiosperm", "gymno" = "Gymnosperm")
# 
# margins_plot <- ggplot(grp_df, aes(x = cwd.spstd)) + 
#   geom_line(aes(y = Prediction)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   facet_wrap(~gymno_angio, scales = "free", ncol = 1) +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") + 
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::comma)
# margins_plot
# ggsave(paste0(wdir, 'figures\\gymno_angio.svg'), plot = margins_plot)


old_mod_df <- read_csv(paste0(wdir, "out/first_stage/archive/site_pet_cwd_std_augmented.csv"))
old_mod_df <- old_mod_df %>% 
  filter(outlier==0)


predict_cwd_effect <- function(mod_df){
  mod_df <- mod_df %>%
    drop_na()
  gen_mod <- feols(estimate_cwd.an ~ cwd.spstd + pet.spstd, data=mod_df,
                   weights = mod_df$errorweights, vcov = "conley")
  mod_cl <- tidy(coeftest(gen_mod))
  pred_min <- trim_df$cwd.spstd %>% quantile(0.025)
  pred_max <- trim_df$cwd.spstd %>% quantile(0.975)
  predictions <- prediction(gen_mod, at = list(cwd.spstd = seq(pred_min, pred_max, 0.1)), vcov = cluster_vcov, calculate_se = T) %>% 
    summary() %>% 
    rename(cwd.spstd = "at(cwd.spstd)")
  out <- tibble(predictions = list(predictions), model = list(mod_cl))
  return(out)
}




## Plot marginal effects by genus
genus_freq <- trim_df %>% 
  group_by(genus) %>% 
  summarise(n_collections = n_distinct(collection_id),
            min_cwd = min(cwd.spstd),
            max_cwd = max(cwd.spstd),
            range_cwd = max_cwd - min_cwd) %>% 
  arrange(desc(n_collections))

genus_keep <- genus_freq %>% 
  filter(n_collections>10) %>%
  pull(genus)


genus_key <- sp_info %>% 
  select(genus, gymno_angio) %>% 
  unique()

genus_df <- trim_df %>% 
  filter(genus %in% genus_keep) %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(prediction = map(data, predict_cwd_effect)) %>% 
  unnest(prediction) %>% 
  left_join(genus_key, by = "genus")

genus_coefs <- genus_df %>%
  unnest(model) %>% 
  filter(term == "cwd.spstd") %>% 
  select(genus, estimate, p.value)

genus_coefs <- genus_coefs %>% 
  mutate(lab = paste0("Slope: ", as.character(format(estimate, digits = 3, scientific = F)), 
                      "\nP value: ", as.character(format(p.value, digits = 3, scientific = T)))) %>% 
  left_join(genus_freq, by = "genus") %>% 
  arrange(desc(n_collections)) %>%
  print()

genus_predictions <- genus_df %>% 
  unnest(predictions) %>% 
  select(-data, -model) %>% 
  left_join(genus_freq, by = "genus") %>% 
  left_join(genus_coefs, by = "genus")

saveRDS(genus_predictions, paste0(wdir, "out/second_stage/genus_mods.rds"))

# histogram <- ggplot(flm_df %>% filter(genus %in% genus_keep), aes(x = cwd.spstd)) + 
#   facet_wrap(~genus) +
#   geom_histogram(aes(fill = gymno_angio), bins = 40, alpha=0.4) +
#   xlim(c(-3, 3)) +
#   theme_bw(base_size = 22) + 
#   ylab("# sites") +
#   xlab("Historic CWD\n(Deviation from species mean)")
# histogram


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Deprecated code below this  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
margins(mod, at  = list(pet.spstd = -2:2))
cplot(mod, 'pet.spstd', what = "prediction")
persp(mod, 'cwd.spstd', 'pet.spstd')
img <- image(mod, 'cwd.spstd', 'pet.spstd')




group_dat <- plot_dat %>% 
  group_by(cwd.q) %>% 
  summarize(wvar = weighted.var(estimate_cwd.an, errorweights),
            wsd = sqrt(wvar),
            wmean = weighted.mean(estimate_cwd.an, errorweights),
            n = n(),
            error = qt(0.975, df = n-1)*wsd/sqrt(n),
            lower = wmean - error,
            upper = wmean + error)


ggplot(group_dat, aes(x=cwd.q, y=wmean)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) 


pet_dat <- plot_dat %>% 
  group_by(pet.q) %>% 
  summarize(wvar.cwd = weighted.var(estimate_cwd.an, errorweights),
            wsd.cwd = sqrt(wvar),
            wmean = weighted.mean(estimate_cwd.an, errorweights),
            n = n(),
            error = qt(0.975, df = n-1)*wsd/sqrt(n),
            lower = wmean - error,
            upper = wmean + error)

ggplot(group_dat, aes(x=pet.q, y=wmean)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) 


p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p1



## Plot estimates by historic cwd / pet - would be awesome to combine these and make prettier
plot_dat <- flm_df
coef_plot1 <- 
  plot_dat %>% 
  ggplot(aes(x = cwd.spstd, y = pet.spstd, 
             z = estimate_pet.an, weight = errorweights)) +
  stat_summary_hex() +
  xlim(c(-3, 4)) +
  #ylim(c(-1.5, 1.5))+
  scale_fill_gradientn (colours = c("darkblue","lightblue", "lightgrey")) +
  theme_bw()+
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme(legend.title = element_blank(),
        legend.position = "bottom")

coef_plot1 

coef_plot2 <- plot_dat %>% ggplot(aes(x = cwd.spstd, y = estimate_cwd.an, weight = errorweights)) +
  stat_summary_bin(bins = 20)+
  theme_bw()+
  xlim(-3, 4)+
  ylab("Estimate CWD")+
  xlab("Deviation from mean CWD")


coef_plot2

coef_plot3 <- plot_dat %>% ggplot(aes(x = pet.spstd, y = estimate_cwd.an, weight = nobs)) +
  stat_summary_bin(bins = 20)+
  xlim(-3,4)+
  xlab("Deviation from mean PET")+
  ylab("Estimate CWD")+
  theme_bw()+
  rotate()

coef_plot3
coef_plot1/coef_plot2


ggarrange(coef_plot2, NULL,coef_plot1, coef_plot3, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = T)




## Plot prediction along one axis
cwd.x <- seq(-2, 2, length.out=1000)
pred <- predict(cwd.mod, data.frame(cwd.spstd=cwd.x, pet.spstd = 0),
                interval="prediction",level=0.90)
pred <- cbind(pred, cwd.x) %>%
  as_tibble()
plot(cwd.x, pred$fit, type="l", xlab="Historic CWD",ylab="Estimated effect of CWD",las=1,lwd=2,col="#0dcfca",
     ylim = c(pred$lwr %>% min(), pred$upr %>% max()))
polygon(c(cwd.x,rev(cwd.x)),c(pred$lwr,rev(pred$upr)),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)


## Plot prediction along two axes
cwd.x <- seq(-2, 2, length.out=50)
pet.x <- seq(-2, 2, length.out=50)
grid.x <- expand.grid(cwd.x, pet.x)

pred <- predict(cwd.mod,data.frame(cwd.spstd=grid.x$Var1, pet.spstd = grid.x$Var2),
                interval="prediction",level=0.90)
pred <- cbind(pred,grid.x) %>% 
  as_tibble()

p <- pred %>% ggplot(aes(x = Var1, y = Var2, fill = fit)) +
  geom_tile()
p



## Bootstrap distributions from first stage
N <- 100

draw_est_smpl <- function(dat){
  var1 <- as.numeric(dat[which(names(dat)=="std.error_cwd.an")])^2
  var2 <- as.numeric(dat[which(names(dat)=="std.error_pet.an")])^2
  cov <- as.numeric(dat[which(names(dat)=="pet_cwd_cov")])^2
  vcov <- matrix(c(var1, cov, cov, var2), 2)
  mu <- c(data$estimate_cwd.an, data$estimate_pet.an)
  smpl <- mvrnorm(N, mu = mu, Sigma = vcov )
  colnames(smpl)=c("cwd.est.smpl","pet.est.smpl")
  smpl <- smpl %>% as_tibble()
  return(smpl)
}

draw_df <- flm_df %>% 
  rowwise %>% 
  do(draw = draw_est_smpl(.))
draw_df <- cbind(flm_df, draw_df) %>% 
  unnest(draw)
  

## Boxplots of effects by cwd / pet bins


sig_df <- flm_df %>% 
  filter(p.value_cwd.an<0.05)

plot_dat <- flm_df %>%
  # filter(genus %in% c("pi", "ps", "ju")) %>% 
  group_by("species_id") %>%
  mutate(cwd.q = as.factor(ntile(cwd.ave, 5)),
         pet.q = as.factor(ntile(pet.ave, 5)),
         cwd.high = as.factor(ntile(cwd.ave, 2)-1),
         pet.high = as.factor(ntile(pet.ave, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high)) %>%
  ungroup()



plot_dat %>% 
  mutate(inv_std = 1/ std.error_cwd.an) %>% 
  group_by(cwd.q) %>% 
  summarize(wmean = weighted.mean(estimate_cwd.an, inv_std, na.rm = TRUE))

plot_dat %>% 
  group_by(cwd.q) %>% 
  summarize(mean = mean(estimate_cwd.an, na.rm = TRUE))


p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p2 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p3 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p4 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)

(p1 | p2) / (p3 | p4)





sp <- "pipo"
sp_data <- flm_df %>% filter(species_id == sp)
sp_data %>% filter(cwd.spstd<1) %>% 
  ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) + 
  geom_point() + 
  geom_point(data = sp_data %>% filter(collection_id=="CO559"), color = 'red', alpha = 1, size = 5) + 
  geom_smooth(method = lm) +
  theme_bw(base_size = 25) +
  ylab("Sensitivity to CWD (\U03B2)")+
  xlab("Historic CWD\n(Deviation from species mean)") +
  geom_hline(yintercept = 0)

sp_data %>% arrange(estimate_cwd.an) %>% select(collection_id, estimate_cwd.an, p.value_cwd.an, cwd.ave, cwd.spstd)
(sp_data %>% arrange(desc(estimate_cwd.an)) %>% select(collection_id, estimate_cwd.an, p.value_cwd.an, cwd.ave, cwd.spstd))[10:20,]
