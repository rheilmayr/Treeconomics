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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "rwi_data\\")
dendro_sites <- read.csv(paste0(dendro_dir, "2_valid_sites.csv")) %>% 
  select(-X) %>% 
  mutate(file_name = paste0('sid-', site_id, '_spid-', species_id, '.csv'))

# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site))

cwd_sites <- cwd_df %>% 
  select(site) %>% 
  distinct()
cwd_dendro_sites <- dendro_sites %>% 
  inner_join(cwd_sites, by = c("site_id" = "site"))



ihsTransform <- function(y) {log(y + (y ^ 2 + 1) ^ 0.5)}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site historic climate ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  group_by(site_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum(petm),
            temp.an = mean(tmean),
            ppt.an = sum(ppt))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define regression model and helper funcs -------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dendro_lm <- function(s_id, sp_id){
  dendro_file <- paste0('sid-', s_id, '_spid-', sp_id, '.csv')
  dendro_dat <- read.csv(paste0(dendro_dir, dendro_file)) %>% 
    mutate(site_id = as.character(s_id))
  dendro_dat <- dendro_dat %>% 
    inner_join(clim_df, by = c("site_id", "year"))
  
  failed <- F
  # Try to run felm. Typically fails if missing cwd / aet data 
  tryCatch(
    expr = {
      dendro_dat <- dendro_dat %>%
        mutate(ln_rwi = log(rwi))
      mod <- felm(ln_rwi ~ pet.an + cwd.an |tree_id|0|tree_id+year, data = dendro_dat)
    },
    error = function(e){ 
      message("Returned error on site ",s_id, " and species ", sp_id)
      print(e)
      failed <<- T
    }
  )
  if (failed){
    return(NULL)
  }
  return(mod)
}


getcov <- function(felm_mod, cov_vars = c("pet.an", "cwd.an")){
  failed <- F
  tryCatch(
    expr = {
      vcov <- felm_mod$vcv
      pet_cwd_cov <- vcov %>% 
        subset(rownames(vcov) == cov_vars[1]) %>% 
        as_tibble() %>% 
        pull(cov_vars[2])
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run site-level regressions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# sites <- cwd_dendro_sites[1:10,]
# i = 2
# s_id <- cwd_dendro_sites[i,] %>%
#   pull(site_id)
# sp_id <- cwd_dendro_sites[i,] %>%
#   pull(species_id)



site_lm <- cwd_dendro_sites %>% 
  group_by(site_id, species_id) %>% 
  nest() %>% 
  mutate(mod = map2(site_id, species_id, dendro_lm),
         cov = map(mod, getcov),
         nobs = map(mod, getnobs),
         mod = map(mod, tidy))
site_lm <- site_lm %>% 
  select(-data) %>% 
  unnest(c(mod, cov, nobs)) %>%
  filter(term %in% c('cwd.an', 'pet.an'))

siteCoef <- site_lm %>%
  pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))

siteCoef %>% write.csv(paste0(wdir, 'first_stage\\', 'log_log_pet_cwd.csv'))








#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run tree-level regressions ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_mod <- function(tree_df){
  failed <- F
  # Try to run felm. Typically fails if missing cwd / aet data 
  tryCatch(
    expr = {
      mod <- lm(ln_rwi ~ pet.an + cwd.an, data = tree_df)
    },
    error = function(e){ 
      message("Returned error")
      print(e)
      failed <<- T
    }
  )
  if (failed){
    return(NULL)
  }
  return(mod)
}


tree_summary <- tree_df %>% 
  select(tree_id, species_id, site_id) %>% 
  distinct()

tree_lm <- tree_df %>% 
  group_by(tree_id) %>%
  nest() %>% 
  mutate(mod = map(data, fs_mod),
         mod = map(mod, tidy))

tree_lm <- tree_lm %>% 
  unnest(mod) %>%
  select(-data) %>% 
  filter(term %in% c('cwd.an', 'pet.an'))

tree_coef <- tree_lm %>%
  pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))

tree_coef <- tree_coef %>% 
  left_join(tree_summary, by = c("tree_id"))

# Merge historic climate
tree_coef <- tree_coef %>% 
  left_join(hist_clim_df, by = c("site_id")) %>% 
  left_join(sp_clim_df, by = c("species_id")) %>% 
  mutate(pet.spstd = (pet.ave - pet.sp.ave) / pet.sp.sd,
         cwd.spstd = (cwd.ave - cwd.sp.ave) / cwd.sp.sd)

tree_coef <- tree_coef %>% 
  mutate(errorweights = 1/(std.error_cwd.an^2)) %>% 
  drop_na()
mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = tree_coef$errorweights, data = tree_coef)
summary(mod)