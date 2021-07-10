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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'out\\climate\\essentialcwd_data.csv')
cwd_df <- read_csv(cwd_csv)
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site))

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out\\dendro\\")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)
# old_dendro_df <- read_csv(paste0(dendro_dir, "rwi_long_old.csv"))
# old_dendro_df <- old_dendro_df %>% 
#   filter(year>1900) %>% 
#   select(-core_id)
# dendro_sites <- read.csv(paste0(dendro_dir, "2_valid_sites.csv")) %>% 
#   select(-X) %>% 
#   mutate(file_name = paste0('sid-', site_id, '_spid-', species_id, '.csv'))



# cwd_sites <- cwd_df %>% 
#   select(site) %>% 
#   distinct()
# cwd_dendro_sites <- dendro_sites %>% 
#   inner_join(cwd_sites, by = c("site_id" = "site"))

# 4. Site information
site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)

dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')

# 5. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))
psme_list <- site_smry %>% 
  filter(species_id == "psme") %>% 
  pull(collection_id)


# 6. Historic species niche data
niche_csv <- paste0(wdir, 'out/climate/clim_niche.csv')
niche_df <- read_csv(niche_csv) %>% 
  select(-X1)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site historic climate ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  rename(collection_id = site_id) %>% 
  group_by(collection_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum((aet+cwd)))
            # temp.an = mean(tmean),
            # ppt.an = sum(ppt))

dendro_df <- dendro_df %>% 
  left_join(clim_df, by = c("collection_id", "year"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Standardize annual weather by species niche  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niche_df <- niche_df %>% 
  select(species_id = sp_code, sp_pet_mean = pet_mean, sp_pet_sd = pet_sd, sp_cwd_mean = cwd_mean, sp_cwd_sd = cwd_sd)

dendro_df <- dendro_df %>% 
  left_join(niche_df, by = c("species_id")) %>% 
  mutate(pet.an.spstd = (pet.an - sp_pet_mean) / sp_pet_sd,
         cwd.an.spstd = (cwd.an - sp_cwd_mean) / sp_cwd_sd)

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


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Identify tree-level weather history ------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# tree_clim <- dendro_df %>% 
#   select(collection_id, year, tree, aet.an, cwd.an, pet.an) %>% 
#   distinct()
# 
# tree_clim <- tree_clim %>%
#   filter(year<1980) %>% 
#   group_by(collection_id, tree) %>% 
#   summarize(first_year = min(year),
#             young = first_year > 1901,
#             cwd.trmean = mean(cwd.an),
#             cwd.trsd = sd(cwd.an),
#             pet.trmean = mean(pet.an),
#             pet.trsd = sd(pet.an))
# 
# # tree_clim %>% 
# #   filter(young == T) %>% 
# #   group_by(collection_id) %>% 
# #   tally() %>% 
# #   filter(n>1)


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
fs_mod <- function(site_data){
  rwl_mean <- site_data$rwl %>% mean(na.rm = TRUE)
  rwl_sd <- site_data$rwl %>% sd(na.rm = TRUE)
  
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
  ncores <- site_data %>% select(tree, core) %>%  n_distinct()
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
        mod <- lm(rwi ~ pet.an + cwd.an, data = site_data)
        # if (ntrees==1){
        #   mod <- lm(rwi ~ cwd.an + pet.an, site_data)
        # } else{
        #   mod <- lm(rwi ~ cwd.an + pet.an, site_data)
        # }
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
      },
      error = function(e){ 
        message("Returned regression error")
        reg_error <<- e[1]
        failed <<- T
      }
    )    
  }
  if (failed){
    return(tibble(mod = list(mod), nobs = nobs, ncores = ncores, ntrees = ntrees, rwl_mean = rwl_mean, rwl_sd = rwl_sd, error = reg_error))
  }
  return(tibble(mod = list(mod), nobs = nobs, ncores = ncores, ntrees = ntrees, rwl_mean = rwl_mean, rwl_sd = rwl_sd, error = reg_error))
}

site_df <- dendro_df %>% 
  drop_na() %>% 
  mutate(cwd.an = cwd.an.spstd,
         pet.an = pet.an.spstd) %>% 
  mutate(ln_rwi = log(rwi)) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  filter(nobs>10) %>% 
  nest()


site_df <- site_df %>% 
  mutate(fs_result = map(data, fs_mod))

data_df <- site_df %>% 
  select(collection_id,data)

site_df <- site_df %>% 
  select(collection_id, fs_result) %>% 
  unnest(fs_result)

site_df <- site_df[which(!(site_df %>% pull(mod) %>% is.na())),]
site_df <- site_df %>% 
  unnest(mod)

site_df <- site_df %>% 
  select(-error)

site_df %>% write.csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv'))



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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create illustrative figure --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
high_sens = "CO559"
low_sens = "CA585"
site_data <- dendro_df %>% 
  filter(collection_id == high_sens)

site_data %>% ggplot(aes(x = cwd.an, y = rwi)) +
  geom_point(shape = 23) +
  geom_smooth(method=lm) +
  theme_bw(base_size = 25) +
  ylab("Ring width index")+
  xlab("Climatic water deficit")

  

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

fs_mod <- function(tree_data){
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ncores <- tree_data %>% select(core) %>%  n_distinct()
  no_cwd_var <- (tree_data %>% select(cwd.an) %>% n_distinct() == 1)
  no_pet_var <- (tree_data %>% select(pet.an) %>% n_distinct() == 1)
  if (no_cwd_var | no_pet_var) {
    message("Site has no variation in CWD or PET")
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / aet data 
    tryCatch(
      expr = {
        # mod <- lm(ln_rwi ~ pet.an + cwd.an, data = tree_data)
        if (ncores==1){
          femod <- feols(ln_rwi ~ cwd.an + pet.an, tree_data)
        } else{
          femod <- feols(ln_rwi ~ cwd.an + pet.an | core, tree_data)
        }
        vcov <- femod$cov.unscaled
        pet_cwd_cov <- vcov %>% 
          subset(rownames(vcov) == "cwd.an") %>% 
          as_tibble() %>% 
          pull("pet.an")
        nobs <- femod$nobs
        femod <- tidy(femod) %>%
          filter(term %in% c('cwd.an', 'pet.an')) %>% 
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
    return(tibble(mod = list(femod), cov = pet_cwd_cov, nobs = nobs, ncores = ncores, error = reg_error))
  }
  return(tibble(mod = list(femod), cov = pet_cwd_cov, nobs = nobs, ncores = ncores, error = reg_error))
}


tree_df <- dendro_df %>% 
  drop_na() %>% 
  mutate(ln_rwi = log(rwi)) %>% 
  group_by(collection_id, tree) %>%
  add_tally(name = 'nobs') %>% 
  filter(nobs>10) %>% 
  nest()

tree_df <- tree_df %>% 
  mutate(fs_result = map(data, fs_mod)) 

tree_df <- tree_df %>% 
  select(collection_id, tree, fs_result, data) %>% 
  unnest(fs_result) 

tree_df <- tree_df[which(!(tree_df %>% pull(mod) %>% is.na())),]
tree_df <- tree_df %>% 
  filter(map_lgl(error, is.null)) %>% 
  unnest(mod)

tree_df <- tree_df %>% 
  select(-error, -data)

tree_df <- tree_df %>% 
  left_join(tree_clim, by = c("collection_id", "tree"))

# tree_df %>% filter(p.value_pet.an<0.05) %>% select(estimate_cwd.an, estimate_pet.an) %>% summary()

tree_df %>% write.csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv'))




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


