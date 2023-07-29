#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Run plot-level regressions of RWI sensitivity to annual weather variability
#
# Input files:
# - site_an_clim.gz: File detailing site-level weather history. Generated using "3b. Species_niche.R"
# - rwi_long.csv: Data containing processed RWI data. Generated using "1b. Parse ITRDB.R"
# - site_summary.csv: Summary data about each site. Generated using "1b. Parse ITRDB.R"
#
# Output files:
# - example_sites.csv: Dendrochronologies for two example sites. Used for methods summary figure.
# - site_pet_cwd_std.csv: Table of first stage regression parameters for baseline specification.
# - site_pet_cwd_std_nb.csv: Table of first stage regression parameters for robustness model using NB desplining.
# - site_pet_cwd_std_ar.csv: Table of first stage regression parameters for robustness model using AR desplining.
# - site_temp_cwd_std.csv:Table of first stage regression parameters for robustness model using temperature in place of PET.
#
#
# ToDo:
# - fix joins to prevent duplicate species_id
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyr)
library(tidyverse)
# library(tidylog)
library(dbplyr)
library(broom.mixed)
library(broom)
library(purrr)
library(fixest)
library(dtplyr)
library(furrr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import and integrate data --------------------------------------------------------
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export example sites for presentations  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ex_sites <- c("CO559", "CA585")
dendro_df %>% 
  filter(collection_id %in% ex_sites) %>% 
  write.csv(paste0(wdir, "out\\dendro\\example_sites.csv"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define regression model  -------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_mod <- function(site_data, outcome = "rwi", energy_var = "pet.an", mod_type = "lm"){
  failed <- F
  reg_error <- NA
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
  return(tibble(mod = list(mod), nobs = nobs, ntrees = ntrees, error = reg_error))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run site-level regressions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_df <- dendro_df %>% 
  drop_na() %>% 
  rename(cwd.an = cwd.an.spstd,
         pet.an = pet.an.spstd,
         temp.an = temp.an.spstd) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  filter(nobs>10) %>% 
  nest()


fs_mod_bl <- partial(fs_mod, outcome = "rwi", energy_var = "pet.an", mod_type = "lm")
fs_mod_nb <- partial(fs_mod, outcome = "rwi_nb", energy_var = "pet.an", mod_type = "lm")
fs_mod_ar <- partial(fs_mod, outcome = "rwi_ar", energy_var = "pet.an", mod_type = "lm")
fs_mod_temp <- partial(fs_mod, outcome = "rwi", energy_var = "temp.an", mod_type = "lm")
fs_mod_re <- partial(fs_mod, outcome = "rwi", energy_var = "pet.an", mod_type = "lme")

site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = fs_mod_bl),
         fs_result_nb = map(data, .f = fs_mod_nb),
         fs_result_ar = map(data, .f = fs_mod_ar),
         fs_result_temp = map(data, .f = fs_mod_temp),
         fs_result_re = map(data, .f = fs_mod_re))


site_df <- site_df %>% 
  mutate(fs_result_re = map(data, .f = fs_mod_re))


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

fs_df %>% write_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std.csv'))


## Repeat using results from nb detrended data
fs_nb <- site_df %>% 
  select(collection_id, fs_result_nb) %>% 
  unnest(fs_result_nb)
fs_nb <- fs_nb[which(!(fs_nb %>% pull(mod) %>% is.na())),]
fs_nb <- fs_nb %>% 
  unnest(mod) %>% 
  select(-error)
fs_nb %>% write_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std_nb.csv'))


## Repeat using results from ar detrended data
fs_ar <- site_df %>% 
  select(collection_id, fs_result_ar) %>% 
  unnest(fs_result_ar)
fs_ar <- fs_ar[which(!(fs_ar %>% pull(mod) %>% is.na())),]
fs_ar <- fs_ar %>% 
  unnest(mod) %>% 
  select(-error)
fs_ar %>% write_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std_ar.csv'))


## Repeat using results from temp model
fs_temp <- site_df %>% 
  select(collection_id, fs_result_temp) %>% 
  unnest(fs_result_temp)
fs_temp <- fs_temp[which(!(fs_temp %>% pull(mod) %>% is.na())),]
fs_temp <- fs_temp %>% 
  unnest(mod) %>% 
  select(-error)
fs_temp %>% write_csv(paste0(wdir, 'out/first_stage/site_temp_cwd_std.csv'))


## Repeat using results from re model
fs_re <- site_df %>% 
  select(collection_id, fs_result_re) %>% 
  unnest(fs_result_re)
fs_re <- fs_re[which(!(fs_re %>% pull(mod) %>% is.na())),]
fs_re <- fs_re %>% 
  unnest(mod) %>% 
  select(-error)
fs_re %>% write_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std_re.csv'))
