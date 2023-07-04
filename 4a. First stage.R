#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Run plot-level regressions of RWI sensitivity to annual weather variability
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using "3b. Species_niche.R"
# - site_an_clim.gz: File detailing site-level weather history. Generated using "3b. Species_niche.R"
# - rwi_long.csv: Directory containing processed RWI data from "1b. Parse ITRDB.R"
# - species_gen_gr.csv: Annotated data about species.
# - site_summary.csv: Generated using "1b. Parse ITRDB.R"
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
wdir <- 'remote/'


# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out/dendro/")
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
an_site_clim <- read_rds(paste0(wdir, "out/climate/site_an_clim.gz"))
dendro_df <- dendro_df %>% 
  left_join(an_site_clim, by = c("collection_id", "year"))


# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
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
fs_mod <- function(site_data, outcome = "rwi", energy_var = "pet.an"){
  rwl_mean <- site_data$rwl %>% mean(na.rm = TRUE)
  rwl_sd <- site_data$rwl %>% sd(na.rm = TRUE)
  
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
  no_cwd_var <- (site_data %>% select(cwd.an) %>% n_distinct() == 1)
  no_pet_var <- (site_data %>% select(energy_var) %>% n_distinct() == 1)
  
  if (no_cwd_var | no_pet_var) {
    message(paste0("Site has no variation in cwd.an or ", energy_var))
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / pet data 
    tryCatch(
      expr = {
        formula <- as.formula(paste0(outcome, " ~ ", energy_var, " + cwd.an"))
        mod <- lm(formula, data = site_data)

        mod_sum <- summary(mod)
        mod_vcov <- vcov(mod)
        # cov <- list(int_cwd = mod_vcov[1, 2], 
        #             int_pet = mod_vcov[1, 3], 
        #             pet_cwd = mod_vcov[2, 3])
        nobs <- nobs(mod)
        mod <- tidy(mod) %>%
          mutate(term = term %>% str_replace("\\(Intercept\\)", "intercept")) %>% 
          filter(term %in% c('intercept', 'cwd.an', energy_var)) %>% 
          pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
        # mod <- mod %>% 
        #   rename_all(funs(stringr::str_replace_all(., energy_var, 'energy.an')))
        mod$cov_int_cwd = mod_vcov[c("(Intercept)"), c("cwd.an")]
        cov_var_name <- paste0("cov_int_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("(Intercept)"), c(energy_var)]
        cov_var_name <- paste0("cov_cwd_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("cwd.an"), c(energy_var)]
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
  rename(cwd.an = cwd.an.spstd,
         pet.an = pet.an.spstd,
         temp.an = temp.an.spstd) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  filter(nobs>10) %>% 
  nest()


# i = 200
# test_data <- site_df[i,2][[1]][[1]]
# test_data
# site_data = test_data
# outcome = "rwi"
# energy_var = "pet.an"
# mod <- fs_mod(test_data, outcome = outcome, energy_var = energy_var)


# comb_mod <- lm(rwi ~ cwd.an.spstd + pet.an.spstd + temp.an.spstd, data = test_data)
# summary(comb_mod)
# temp_mod <- lm(rwi ~ cwd.an.spstd + temp.an.spstd, data = test_data)
# summary(temp_mod)
# pet_mod <- lm(rwi ~ cwd.an.spstd + pet.an.spstd, data = test_data)
# summary(pet_mod)


# #### REVIEWER COMMENTS - TESTING ALTERNATE MODELS
# i = 3500
# data_df <- site_df[i, 2]
# data_df <- data_df[[1]][[1]]
# mod <- lm(rwi ~ cwd.an + pet.an, data_df)
# summary(mod)
# 
# data_df <- data_df %>%
#   mutate(tree_id = as.factor(tree))
# data_df$tree %>% unique() %>% length()
# mod_re <- lmer(rwi ~ 1 + cwd.an + pet.an + (1 | tree_id), data_df)
# summary(mod_re)
# ## Identical coefficients to main specification, but singular / boundary warnings
# 
# data_df <- data_df %>%
#   mutate(tree_id = as.factor(tree))
# data_df$tree %>% unique() %>% length()
# mod_fe <- feols(rwi ~ cwd.an + pet.an | tree_id, data_df)
# summary(mod_fe)
# ## Very similar coefficients to main specification. Could include this in robustness plot, but doing predictions with this would be a huge headache
# 
# mod_bayes <- stan_lmer(rwi ~ cwd.an + pet.an + (1 | tree_id), data_df)
# summary(mod_bayes)
# 
# ave_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_ave_clim.gz"))
# dendro_df <- dendro_df %>%
#   left_join(ave_site_clim, by = "collection_id")
# 
# mod_full_fe <- lm(rwi ~ cwd.an.spstd * cwd.spstd * I(cwd.spstd**2)  + cwd.an.spstd * pet.spstd * I(pet.spstd**2) + pet.an.spstd * cwd.spstd * I(cwd.spstd**2) + pet.an.spstd * pet.spstd * I(pet.spstd**2), data = dendro_df)
# summary(mod_full_fe)
# 
# library(lme4)
# mod_full_re <- lmer(rwi ~ cwd.an.spstd * cwd.spstd * I(cwd.spstd**2)  + cwd.an.spstd * pet.spstd * I(pet.spstd**2) + pet.an.spstd * cwd.spstd * I(cwd.spstd**2) + pet.an.spstd * pet.spstd * I(pet.spstd**2) + (1 | collection_id), data = dendro_df)
# summary(mod_full_fe)
# ## Yields fairly similar coefficient estimate to two-stage model. Generate ME plot?
# 
# #### REVIEWER COMMENTS - TESTING ALTERNATE MODELS
# 
# 
# #### REVIEWER COMMENTS - DRIVEN BY VARIABILITY?
# sd_df <- dendro_df %>% 
#   group_by(collection_id) %>% 
#   summarise(cwd_sd = sd(cwd.an.spstd, na.rm = TRUE),
#             cwd_mean = mean(cwd.an.spstd, na.rm = TRUE))
# 
# mod <- lm(cwd_mean ~ cwd_sd, sd_df)
# summary(mod)
# sd_df %>% ggplot(aes(x = cwd_mean, y = cwd_sd)) +
#   geom_point()
# 
# 
# #### REVIEWER COMMENTS - DRIVEN BY VARIABILITY?



fs_mod_bl <- partial(fs_mod, outcome = "rwi", energy_var = "pet.an")
fs_mod_nb <- partial(fs_mod, outcome = "rwi_nb", energy_var = "pet.an")
fs_mod_ar <- partial(fs_mod, outcome = "rwi_ar", energy_var = "pet.an")
fs_mod_temp <- partial(fs_mod, outcome = "rwi", energy_var = "temp.an")

site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = fs_mod_bl),
         fs_result_nb = map(data, .f = fs_mod_nb),
         fs_result_ar = map(data, .f = fs_mod_ar),
         fs_result_temp = map(data, .f = fs_mod_temp))


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

# fs_df %>%
#   # filter(p.value_cwd.an<0.05) %>% 
#   select(estimate_cwd.an) %>% 
#   summary()
# 
# fs_df %>%
#   filter(p.value_temp.an<0.05) %>% 
#   select(estimate_temp.an) %>% 
#   summary()
# 
# pet_fs_df %>%
#   filter(p.value_pet.an<0.05) %>% 
#   select(collection_id, estimate_pet.an) %>% 
#   summary()
# 
# pet_fs_df %>%
#   # filter(p.value_cwd.an<0.05) %>% 
#   select(collection_id, estimate_cwd.an) %>% 
#   summary()

# fs_df %>% write_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv'))
fs_df %>% write_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std.csv'))


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


## Repeat using results from temp model
fs_temp <- site_df %>% 
  select(collection_id, fs_result_temp) %>% 
  unnest(fs_result_temp)
fs_temp <- fs_temp[which(!(fs_temp %>% pull(mod) %>% is.na())),]
fs_temp <- fs_temp %>% 
  unnest(mod) %>% 
  select(-error)
fs_temp %>% write_csv(paste0(wdir, 'out\\first_stage\\site_temp_cwd_std.csv'))




