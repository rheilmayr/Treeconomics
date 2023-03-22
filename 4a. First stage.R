#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Run plot-level regressions of RWI sensitivity to annual weather variability
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using species_niche.R
# - tree_ring_data_V2.db: Compiled database of ITRDB observations 
# - essentialcwd_data.csv: File detailing plot-level weather history
# - dendro_dir: Directory containing processed RWI data from "1. Dendro preprocess.R"
#
# ToDo:
# - fix joins to prevent duplicate species_id
# - think through how to deal with CWD outliers
# - track down lost observations - currently dropping a lot due to NAN or failed RWI generation
# - Incorporate code from tree-level analysis script to generate DLNM plots? 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyr)
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(lfe)
library(broom)
library(purrr)
library(fixest)
library(dtplyr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# # 2. Site-specific weather history
# cwd_csv <- paste0(wdir, 'out\\climate\\essentialcwd_data.csv')
# cwd_df <- read_csv(cwd_csv)
# cwd_df <- cwd_df %>% 
#   mutate("site_id" = as.character(site))

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out\\dendro\\")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            rwl = mean(rwl),
            rwi_ar = mean(rwi_ar),
            rwi_nb = mean(rwi_nb),
            .groups = "drop") %>% 
  as_tibble()


# 2. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_an_clim.gz"))
dendro_df <- dendro_df %>% 
  left_join(an_site_clim, by = c("collection_id", "year"))


# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)

dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')

# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))
psme_list <- site_smry %>% 
  filter(species_id == "psme") %>% 
  pull(collection_id)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export example sites for presentations  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ex_sites <- c("CO559", "CA585")
dendro_df %>% 
  filter(collection_id %in% ex_sites) %>% 
  write.csv(paste0(wdir, "out\\dendro\\example_sites.csv"))

# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Explore variance by site ------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Calculate site-level historic climate
# hist_clim_df <- clim_df %>%
#   group_by(collection_id) %>%
#   filter(year<1980) %>%
#   summarise(aet.ave = mean(aet.an),
#             cwd.ave = mean(cwd.an),
#             pet.ave = mean(pet.an),
#             cwd.min = min(cwd.an)) %>% 
#   left_join(site_smry, by = "collection_id") %>% 
#   select(-genus, -gymno_angio, -species_id)
# 
# 
# hist_clim_df <- hist_clim_df %>% 
#   left_join(niche_df, by = c("species_id" = "sp_code")) %>% 
#   mutate(pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
#          cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)
# 
# 
# site_stats <- dendro_df %>% 
#   group_by(collection_id) %>% 
#   summarise(rwi_var = var(rwi),
#             rwi_mean = mean(rwi)) %>% 
#   left_join(hist_clim_df, by = c("collection_id"))
# 
# site_stats %>% ggplot(aes(x = cwd.spstd, y = rwi_var)) +
#   geom_point(shape = 23) +
#   geom_smooth(method=lm) +
#   theme_bw(base_size = 25) +
#   ylab("Variance in RWI")+
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   ylim(0, 0.5) +
#   xlim(-3, 3)
# 
# site_stats %>% ggplot(aes(x = cwd.spstd, y = rwi_mean)) +
#   geom_point(shape = 23) +
#   geom_smooth(method=lm) +
#   theme_bw(base_size = 25) +
#   ylab("Mean RWI")+
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   ylim(0.9, 1.1) +
#   xlim(-3, 3)
# 
# 
# var_mod <- lm(rwi_var ~ cwd.spstd + pet.spstd, data = site_stats)
# mean_mod <- lm(rwi_mean ~ cwd.spstd + pet.spstd, data = site_stats)








# flm_df <- flm_df %>% 
#   left_join(hist_clim_df, by = c("collection_id"))
# 
# trim_df %>% 
#   filter(errorweights>5300) %>% 
#   ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) +
#   geom_point(shape = 23) +
#   geom_smooth(method=lm) +
#   theme_bw(base_size = 25) +
#   ylab("Sensitivity to CWD")+
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   ylim(-0.001, 0.0005) +
#   xlim(-3, 3)
# 
# flm_df <- flm_df %>% 
#   mutate(errorweights = 1 / (std.error_cwd.an),
#          pet_errorweights = 1 / (std.error_pet.an))
# 
# flm_df <- flm_df %>%
#   group_by(species_id) %>%
#   mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
#          cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T)) %>%
#   ungroup()
# 
# trim_df <- flm_df %>%
#   filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
#          abs(cwd.spstd)<5,
#          abs(pet.spstd)<5,
#          abs(std.error_cwd.an) < 50 * abs(estimate_cwd.an),
#          abs(std.error_cwd.an) < 50 * abs(estimate_cwd.an)) %>% 
#   drop_na()
# coef_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, data = trim_df, weights = errorweights)
# coef_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, data = trim_df %>% filter(errorweights>5300))
# summary(coef_mod)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define regression model and helper funcs -------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# dendro_lm <- function(s_id, sp_id){
#   dendro_file <- paste0('sid-', s_id, '_spid-', sp_id, '.csv')
#   dendro_dat <- read.csv(paste0(dendro_dir, dendro_file)) %>% 
#     mutate(site_id = as.character(s_id))
#   dendro_dat <- dendro_dat %>% 
#     inner_join(clim_df, by = c("site_id", "year"))
#   
#   failed <- F
#   # Try to run felm. Typically fails if missing cwd / aet data 
#   tryCatch(
#     expr = {
#       dendro_dat <- dendro_dat %>%
#         mutate(ln_rwi = log(rwi))
#       mod <- felm(ln_rwi ~ pet.an + cwd.an |tree_id|0|tree_id+year, data = dendro_dat)
#     },
#     error = function(e){ 
#       message("Returned error on site ",s_id, " and species ", sp_id)
#       print(e)
#       failed <<- T
#     }
#   )
#   if (failed){
#     return(NULL)
#   }
#   return(mod)
# }
# 
# 
# getcov <- function(felm_mod, cov_vars = c("pet.an", "cwd.an")){
#   failed <- F
#   tryCatch(
#     expr = {
#       vcov <- felm_mod$vcv
#       pet_cwd_cov <- vcov %>% 
#         subset(rownames(vcov) == cov_vars[1]) %>% 
#         as_tibble() %>% 
#         pull(cov_vars[2])
#     },
#     error = function(e){
#       failed <<- T
#     }
#   )
#   if (failed){
#     return(NULL)
#   }
#   return(pet_cwd_cov)
# }


getcov <- function(felm_mod){
  failed <- F
  tryCatch(
    expr = {
      vcov <- felm_mod$vcv
      pet_cwd_cov <- vcov %>% 
        subset(rownames(vcov) == "cwd.an") %>% 
        as_tibble() %>% 
        pull("pet.an")
    },
    error = function(e){
      failed <<- T
    }
  )
  if (failed){
    return(NULL)
  }
  return(pet_cwd_cov)
}


getnobs <- function(felm_mod){
  failed <- F
  tryCatch(
    expr = {
      nobs <- felm_mod$N
    },
    error = function(e){
      failed <<- T
    }
  )
  if (failed){
    return(NULL)
  }
  return(nobs)
}


ihsTransform <- function(y) {log(y + (y ^ 2 + 1) ^ 0.5)}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run site-level regressions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_mod <- function(site_data, outcome = "rwi"){
  rwl_mean <- site_data$rwl %>% mean(na.rm = TRUE)
  rwl_sd <- site_data$rwl %>% sd(na.rm = TRUE)
  
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
  # ncores <- site_data %>% select(tree, core) %>%  n_distinct()
  no_cwd_var <- (site_data %>% select(cwd.an) %>% n_distinct() == 1)
  no_pet_var <- (site_data %>% select(pet.an) %>% n_distinct() == 1)
  
  # TODO: Update to integrate lags over three years of cwd data?
  # site_data <- site_data %>% 
  #   arrange(tree, core, year) %>% 
  #   group_by(tree, core) %>% 
  #   mutate(lag1_cwd = lag(cwd.an, order_by = year),
  #          lag2_cwd = lag(cwd.an, n = 2, order_by = year))
  # femod <- feols(width ~ cwd.an + lag1_cwd + lag2_cwd + pet.an | tree, site_data)
  # 
  
  if (no_cwd_var | no_pet_var) {
    message("Site has no variation in CWD or PET")
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / aet data 
    tryCatch(
      expr = {
        # TODO: 7-9-21; Should we switch back to tree-level FE model? Complicates intercepts...
        formula <- as.formula(paste0(outcome, " ~ pet.an + cwd.an"))
        mod <- lm(formula, data = site_data)

        # mod <- lm(rwi ~ pet.an + cwd.an, data = site_data)
        # if (ntrees==1){
        #   mod <- lm(rwi ~ cwd.an + pet.an, site_data)
        # } else{
        #   mod <- lm(rwi ~ cwd.an + pet.an, site_data)
        # }
        mod_sum <- summary(mod)
        mod_vcov <- vcov(mod)
        cov <- list(int_cwd = mod_vcov[1, 2], 
                    int_pet = mod_vcov[1, 3], 
                    pet_cwd = mod_vcov[2, 3])
        nobs <- nobs(mod)
        mod <- tidy(mod) %>%
          mutate(term = term %>% str_replace("\\(Intercept\\)", "intercept")) %>% 
          filter(term %in% c('intercept', 'cwd.an', 'pet.an')) %>% 
          pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
        mod$cov_int_cwd = mod_vcov[c("(Intercept)"), c("cwd.an")]
        mod$cov_int_pet = mod_vcov[c("(Intercept)"), c("pet.an")]
        mod$cov_cwd_pet = mod_vcov[c("cwd.an"), c("pet.an")]
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
  return(tibble(mod = list(mod), nobs = nobs, ntrees = ntrees, rwl_mean = rwl_mean, rwl_sd = rwl_sd, error = reg_error))
}
site_df <- dendro_df %>% 
  drop_na() %>% 
  mutate(cwd.an = cwd.an.spstd,
         pet.an = pet.an.spstd) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  filter(nobs>10) %>% 
  nest()

fs_mod_nb <- partial(fs_mod, outcome = "rwi_nb")
fs_mod_ar <- partial(fs_mod, outcome = "rwi_ar")

site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = fs_mod),
         fs_result_nb = map(data, .f = fs_mod_nb),
         fs_result_ar = map(data, .f = fs_mod_ar))


data_df <- site_df %>% 
  select(collection_id,data)

fs_df <- site_df %>% 
  select(collection_id, fs_result) %>% 
  unnest(fs_result)

fs_df <- fs_df[which(!(fs_df %>% pull(mod) %>% is.na())),]
fs_df <- fs_df %>% 
  unnest(mod)

fs_df <- fs_df %>% 
  select(-error)

fs_df %>% write_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv'))


## Repeat using results from nb detrended data
fs_nb <- site_df %>% 
  select(collection_id, fs_result_nb) %>% 
  unnest(fs_result_nb)
fs_nb <- fs_nb[which(!(fs_nb %>% pull(mod) %>% is.na())),]
fs_nb <- fs_nb %>% 
  unnest(mod) %>% 
  select(-error)
fs_nb %>% write_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std_nb.csv'))


## Repeat using results from ar detrended data
fs_ar <- site_df %>% 
  select(collection_id, fs_result_ar) %>% 
  unnest(fs_result_ar)
fs_ar <- fs_ar[which(!(fs_ar %>% pull(mod) %>% is.na())),]
fs_ar <- fs_ar %>% 
  unnest(mod) %>% 
  select(-error)
fs_ar %>% write_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std_ar.csv'))




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Test for site-level concavity --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nonlinear_fs_mod <- function(site_data){
  rwl_mean <- site_data$rwl %>% mean(na.rm = TRUE)
  rwl_sd <- site_data$rwl %>% sd(na.rm = TRUE)
  
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
  # ncores <- site_data %>% select(tree, core) %>%  n_distinct()
  no_cwd_var <- (site_data %>% select(cwd.an) %>% n_distinct() == 1)
  no_pet_var <- (site_data %>% select(pet.an) %>% n_distinct() == 1)
  
  # TODO: Update to integrate lags over three years of cwd data?
  # site_data <- site_data %>% 
  #   arrange(tree, core, year) %>% 
  #   group_by(tree, core) %>% 
  #   mutate(lag1_cwd = lag(cwd.an, order_by = year),
  #          lag2_cwd = lag(cwd.an, n = 2, order_by = year))
  # femod <- feols(width ~ cwd.an + lag1_cwd + lag2_cwd + pet.an | tree, site_data)
  # 
  
  if (no_cwd_var | no_pet_var) {
    message("Site has no variation in CWD or PET")
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / aet data 
    tryCatch(
      expr = {
        # TODO: 7-9-21; Should we switch back to tree-level FE model? Complicates intercepts...
        mod <- lm(rwi ~ pet.an + I(pet.an**2) + cwd.an + I(cwd.an**2), data = site_data)
        # if (ntrees==1){
        #   mod <- lm(rwi ~ cwd.an + pet.an, site_data)
        # } else{
        #   mod <- lm(rwi ~ cwd.an + pet.an, site_data)
        # }
        mod <- tidy(mod) %>%
          pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
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
  return(tibble(mod = list(mod), error = reg_error))
}

site_df <- dendro_df %>% 
  drop_na() %>% 
  mutate(cwd.an = cwd.an.spstd,
         pet.an = pet.an.spstd) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  filter(nobs>10) %>% 
  nest()


library(marginaleffects)
i = 2408
test_data <- site_df[i,2][[1]][[1]]
test_mod <- lm(rwi ~ pet.an + I(pet.an**2) + cwd.an + I(cwd.an**2), data = test_data)
summary(test_mod)
plot_cap(test_mod, condition = "cwd.an") + 
  xlim(-2,2)
test_mod <- nonlinear_fs_mod(test_data)
test_mod$mod


site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = nonlinear_fs_mod))

pull_dist <- function(data){
  cwd_stats = data$cwd.an.spstd %>% quantile(c(0.01, 0.1, 0.5, 0.9, .99)) %>% 
    as_tibble() %>%
    mutate(quantile = c('q0.01', 'q0.1', 'q0.5', 'q0.9', 'q0.99'))
  cwd_stats = pivot_wider(cwd_stats, names_from = quantile, values_from = value)
  return(tibble(cwd_stats))
}

site_df <- site_df %>%
  mutate(cwd_dist = map(data, .f = pull_dist))

site_df <- site_df %>% 
  select(collection_id, fs_result, cwd_dist) %>% 
  unnest(cwd_dist) %>%
  unnest(fs_result)

site_df <- site_df[which(!(site_df %>% pull(mod) %>% is.na())),]
site_df <- site_df %>% 
  unnest(mod)

site_df <- site_df %>% 
  select(-error)

site_df <- site_df %>%
  mutate(slope_median = estimate_cwd.an + 2 * `estimate_I(cwd.an^2)` * `q0.5`,
         slope_extreme = estimate_cwd.an + 2 * `estimate_I(cwd.an^2)` * `q0.9`,
         slope_difs = slope_extreme - slope_median)

site_df %>%
  select(slope_median, slope_extreme, slope_difs) %>%
  summary()

site_df %>%
  filter(slope_median < 0) %>%
  select(slope_median, slope_extreme, slope_difs) %>%
  summary()



site_df %>%
  ggplot(aes(x = slope_median, y = slope_extreme)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlim(-.3, .11) +
  ylim(-.3, .1)

site_df %>% 
  mutate(neg = estimate_cwd.an < 0,
         concave = `estimate_I(cwd.an^2)` < 0,
         sig = `p.value_I(cwd.an^2)` < 0.05) %>%
  group_by(neg, sig, concave) %>%
  tally()


# # test some results
# cid <- "CA585"
# 
# coefs <- site_df %>% 
#   filter(collection_id==cid)
# data <- (data_df %>% 
#   filter(collection_id==cid) %>% 
#   pull(data))[[1]]
# 
# data <- data %>% 
#   mutate(pred = coefs$estimate_intercept + (coefs$estimate_cwd.an * cwd.an) + (coefs$estimate_pet.an * pet.an))
# 
# data %>% 
#   ggplot(aes(x = rwi, y = pred)) +
#   geom_point()


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Identify tree-level weather history ------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# tree_df <- dendro_df %>%
#   select(collection_id, year, tree, core, rwi, cwd.an.spstd, pet.an.spstd) %>%
#   drop_na() %>% 
#   mutate(cwd.an = cwd.an.spstd,
#          pet.an = pet.an.spstd,
#          tree_id = paste(collection_id, tree, sep = "_")) %>% 
#   distinct()
# 
# tree_clim <- tree_df %>%
#   filter(year<1980) %>%
#   group_by(tree_id) %>%
#   summarize(first_year = min(year),
#             last_year = max(year),
#             young = first_year > 1901,
#             cwd.trmean = mean(cwd.an.spstd, na.rm = T),
#             cwd.trsd = sd(cwd.an.spstd, na.rm = T),
#             pet.trmean = mean(pet.an.spstd, na.rm = T),
#             pet.trsd = sd(pet.an.spstd, na.rm = T))
# 
# 
# # tree_clim %>% 
# #   filter(young==T) %>% 
# #   select(collection_id, tree, cwd.trmean) %>% 
# #   unique() %>% 
# #   group_by(collection_id) %>% 
# #   tally() %>% 
# #   filter(n>1)
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Create tree-level regressions ------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# tree_df <- tree_df %>% 
#   group_by(tree_id) %>%
#   add_tally(name = 'nobs') %>% 
#   filter(nobs>10) %>% 
#   nest()
# 
# tree_fs <- tree_df %>% 
#   mutate(fs_result = map(data, fs_mod))
# 
# tree_data <- tree_fs %>% 
#   select(tree_id, data)
# 
# tree_fs <- tree_fs %>% 
#   select(tree_id, fs_result) %>% 
#   unnest(fs_result)
# 
# tree_fs <- tree_fs[which(!(tree_fs %>% pull(mod) %>% is.na())),]
# 
# tree_fs <- tree_fs %>% 
#   unnest(cols = c(mod))
# 
# tree_fs <- tree_fs %>% 
#   select(-error) %>% 
#   ungroup() %>% 
#   left_join(tree_clim, by = "tree_id")
# 
# tree_fs <- tree_fs %>% 
#   select(-fs_result, -rwl_mean, -rwl_sd)
#   
# 
# tree_fs %>% write.csv(paste0(wdir, 'out\\first_stage\\tree_pet_cwd_std.csv'))
# 
# 
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Create illustrative figure --------------------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# high_sens = "CO559"
# low_sens = "CA585"
# site_data <- dendro_df %>% 
#   filter(collection_id == high_sens)
# 
# site_data %>% ggplot(aes(x = cwd.an, y = rwi)) +
#   geom_point(shape = 23) +
#   geom_smooth(method=lm) +
#   theme_bw(base_size = 25) +
#   ylab("Ring width index")+
#   xlab("Climatic water deficit")

  

## TODO: add error reports to site_summary document?

# # sites <- cwd_dendro_sites[1:10,]
# # i = 2
# # s_id <- cwd_dendro_sites[i,] %>%
# #   pull(site_id)
# # sp_id <- cwd_dendro_sites[i,] %>%
# #   pull(species_id)
# 
# 
# 
# site_lm <- cwd_dendro_sites %>% 
#   group_by(site_id, species_id) %>% 
#   nest() %>% 
#   mutate(mod = map2(site_id, species_id, dendro_lm),
#          cov = map(mod, getcov),
#          nobs = map(mod, getnobs),
#          mod = map(mod, tidy))
# site_lm <- site_lm %>% 
#   select(-data) %>% 
#   unnest(c(mod, cov, nobs)) %>%
#   filter(term %in% c('cwd.an', 'pet.an'))
# 
# siteCoef <- site_lm %>%
#   pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
# 
# siteCoef %>% write.csv(paste0(wdir, 'first_stage\\', 'log_cwd_pet.csv'))
# 
# 
# 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run tree-level regressions ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# load_obs <- function(s_id, sp_id){
#   dendro_file <- paste0('sid-', s_id, '_spid-', sp_id, '.csv')
#   dendro_dat <- read.csv(paste0(dendro_dir, dendro_file)) %>% 
#     mutate(site_id = as.character(s_id)) %>% 
#     select(-X, -site_id) %>% 
#     drop_na()
#   return(dendro_dat)
# }
# 
# fs_mod <- function(tree_data){
#   failed <- F
#   reg_error <- NA
#   femod <- NA
#   pet_cwd_cov <- NA
#   nobs <- NA
#   ncores <- tree_data %>% select(core) %>%  n_distinct()
#   no_cwd_var <- (tree_data %>% select(cwd.an) %>% n_distinct() == 1)
#   no_pet_var <- (tree_data %>% select(pet.an) %>% n_distinct() == 1)
#   if (no_cwd_var | no_pet_var) {
#     message("Site has no variation in CWD or PET")
#     failed <- T
#   } else{
#     # Try to run felm. Typically fails if missing cwd / aet data 
#     tryCatch(
#       expr = {
#         # mod <- lm(ln_rwi ~ pet.an + cwd.an, data = tree_data)
#         if (ncores==1){
#           femod <- feols(ln_rwi ~ cwd.an + pet.an, tree_data)
#         } else{
#           femod <- feols(ln_rwi ~ cwd.an + pet.an | core, tree_data)
#         }
#         vcov <- femod$cov.unscaled
#         pet_cwd_cov <- vcov %>% 
#           subset(rownames(vcov) == "cwd.an") %>% 
#           as_tibble() %>% 
#           pull("pet.an")
#         nobs <- femod$nobs
#         femod <- tidy(femod) %>%
#           filter(term %in% c('cwd.an', 'pet.an')) %>% 
#           pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
#       },
#       error = function(e){ 
#         message("Returned regression error")
#         reg_error <<- e[1]
#         failed <<- T
#       }
#     )    
#   }
#   if (failed){
#     return(tibble(mod = list(femod), cov = pet_cwd_cov, nobs = nobs, ncores = ncores, error = reg_error))
#   }
#   return(tibble(mod = list(femod), cov = pet_cwd_cov, nobs = nobs, ncores = ncores, error = reg_error))
# }
# 
# 
# tree_df <- dendro_df %>% 
#   drop_na() %>% 
#   mutate(ln_rwi = log(rwi)) %>% 
#   group_by(collection_id, tree) %>%
#   add_tally(name = 'nobs') %>% 
#   filter(nobs>10) %>% 
#   nest()
# 
# tree_df <- tree_df %>% 
#   mutate(fs_result = map(data, fs_mod)) 
# 
# tree_df <- tree_df %>% 
#   select(collection_id, tree, fs_result, data) %>% 
#   unnest(fs_result) 
# 
# tree_df <- tree_df[which(!(tree_df %>% pull(mod) %>% is.na())),]
# tree_df <- tree_df %>% 
#   filter(map_lgl(error, is.null)) %>% 
#   unnest(mod)
# 
# tree_df <- tree_df %>% 
#   select(-error, -data)
# 
# tree_df <- tree_df %>% 
#   left_join(tree_clim, by = c("collection_id", "tree"))
# 
# # tree_df %>% filter(p.value_pet.an<0.05) %>% select(estimate_cwd.an, estimate_pet.an) %>% summary()
# 
# tree_df %>% write.csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv'))
# 
# 
# 
# 




# 
# # select_sites = cwd_dendro_sites[1:20,]
# 
# tree_df <- cwd_dendro_sites %>% 
#   group_by(site_id, species_id) %>% 
#   nest() %>% 
#   mutate(data = map2(site_id, species_id, load_obs)) %>% 
#   unnest(data)
# 
# 
# tree_df <- tree_df %>% 
#   left_join(clim_df, by = c("site_id", "year")) %>% 
#   drop_na()
# 
# tree_df <- tree_df %>% 
#   mutate(ln_rwi = log(rwi)) %>% 
#   drop_na()
# 
# 
# tree_summary <- tree_df %>% 
#   select(tree_id, species_id, site_id) %>% 
#   distinct()
# 
# tree_lm <- tree_df %>% 
#   group_by(tree_id) %>%
#   nest() %>% 
#   mutate(mod = map(data, fs_mod),
#          mod = map(mod, tidy))
# 
# tree_lm <- tree_lm %>% 
#   unnest(mod) %>%
#   select(-data) %>% 
#   filter(term %in% c('cwd.an', 'pet.an'))
# 
# tree_coef <- tree_lm %>%
#   pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
# 
# tree_coef <- tree_coef %>% 
#   left_join(tree_summary, by = c("tree_id"))


