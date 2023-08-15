#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/11/2023
# Review comment 1.4: 
  # Including an individual tree random effect in equation 1 would better 
  # constrain the parameters and should be done to account for repeated measured of 
  # individual trees/nonindependence of observations. 
#
# Approach:
# - Re-run first stage model using RE, Bayesian RE and FE structures. Compare coefficients.
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
library(rstanarm)
library(broom.mixed)
library(furrr)
library(patchwork)
library(nlme) ##TODO: Shift to nlme lme model?



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data (Same as 4a. First stage imports ---------------------------
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
# Compare alternate models (Reviewer comment) ----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_mod_compare <- function(site_data){
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
  no_cwd_var <- (site_data %>% select(cwd.an) %>% n_distinct() == 1)
  no_pet_var <- (site_data %>% select(pet.an) %>% n_distinct() == 1)
  
  if (no_cwd_var | no_pet_var) {
    message(paste0("Site has no variation in cwd.an or pet.an"))
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / pet data 
    tryCatch(
      expr = {
        site_data <- site_data %>%
          mutate(tree_id = as.factor(tree))
        base_mod <- lm(rwi ~ pet.an + cwd.an, data = site_data)
        formula = as.formula("rwi ~pet.an + cwd.an")
        # re_mod <- nlme:lme(formula, random = ~ 1 | tree, data = site_data)
        fe_mod <- feols(rwi ~ pet.an + cwd.an | tree_id, site_data )
        bayes_mod <- stan_lmer(rwi ~ cwd.an + pet.an + (1 | tree_id), site_data)
        
        base_coef <- broom::tidy(base_mod) %>% filter(term == "cwd.an") %>% select(estimate, std.error) %>% mutate(mod = "base")
        # re_coef <- broom.mixed::tidy(re_mod) %>% filter(term == "cwd.an") %>% select(estimate, std.error) %>% mutate(mod = "re")
        fe_coef <- broom::tidy(fe_mod) %>% filter(term == "cwd.an") %>% select(estimate, std.error) %>% mutate(mod = "fe")
        bayes_coef <- broom.mixed::tidy(bayes_mod) %>% filter(term == "cwd.an") %>% select(estimate, std.error) %>% mutate(mod = "bayes_re")
        
        out_df <- rbind(base_coef, fe_coef, bayes_coef)
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
  return(out_df)
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

n_cores <- availableCores() - 6
future::plan(multisession, workers = n_cores)
my_seed <- 5597

sample_df <- site_df %>% ungroup() %>% slice_sample(n = 100)
start.time <- Sys.time()
sample_df <- sample_df %>%
  mutate(coefficients = future_map(data,
                                   .f = fs_mod_compare,
                                   .options = furrr_options(seed = my_seed,
                                                            packages = c( "nlme", "broom", "broom.mixed"))))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

sample_df <- sample_df %>% 
  unnest(coefficients) %>% 
  select(-data) %>% 
  drop_na() %>%
  pivot_wider(id_cols = "collection_id", names_from = "mod", values_from = c("estimate", "std.error"))



mod <- lm(estimate_base ~ estimate_fe, data = sample_df)
r2_1 <- (mod %>% summary())$r.squared
summary(mod)
plot_1 <- sample_df %>% 
  ggplot(aes(x = estimate_fe, y = estimate_base)) +
  geom_point(size = 2) +
  geom_smooth(method='lm', formula= y~x, linetype="dashed") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  xlab("CWD coefficients from first-stage model with tree fixed effects") +
  ylab("CWD coefficients from baseline first-stage model") +
  geom_text(aes(label = paste0("R2 = ", as.character(round(r2_1, 4))), x = 0.1, y = 0.9),
            hjust = 0, vjust = 0)
plot_1


mod <- lm(estimate_base ~ estimate_bayes_re, data = sample_df)
r2_2 <- (mod %>% summary())$r.squared
summary(mod)
plot_2 <- sample_df %>% 
  ggplot(aes(x = estimate_bayes_re, y = estimate_base)) +
  geom_point(size = 2) +
  geom_smooth(method='lm', formula= y~x, linetype="dashed") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  xlab("CWD coefficients from Bayesian first-stage model with tree random effects") +
  ylab("CWD coefficients from baseline first-stage model") +
  geom_text(aes(label = paste0("R2 = ", as.character(round(r2_2, 4))), x = 0.1, y = 0.9),
            hjust = 0, vjust = 0) +
  ggtitle("Correlation between coefficients from baseline specification\nand Bayesian random effects model")
plot_2


mod <- lm(std.error_base ~ std.error_fe, data = sample_df)
r2_3 <- (mod %>% summary())$r.squared
summary(mod)
plot_3 <- sample_df %>% 
  ggplot(aes(x = std.error_fe, y = std.error_base)) +
  geom_point(size = 2) +
  geom_smooth(method='lm', formula= y~x, linetype="dashed") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  xlab("CWD standard error from first-stage model with tree fixed effects") +
  ylab("CWD standard error from baseline first-stage model") +
  geom_text(aes(label = paste0("R2 = ", as.character(round(r2_3, 4))), x = 0.1, y = 0.9),
            hjust = 0, vjust = 0)
plot_3

mod <- lm(std.error_base ~ std.error_bayes_re, data = sample_df)
r2_4 <- (mod %>% summary())$r.squared
summary(mod)
plot_4 <- sample_df %>% 
  ggplot(aes(x = std.error_bayes_re, y = std.error_base)) +
  geom_point(size = 2) +
  geom_smooth(method='lm', formula= y~x, linetype="dashed") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  xlab("CWD standard error from Bayesian first-stage model with tree random effects") +
  ylab("CWD standard error from baseline first-stage model") +
  geom_text(aes(label = paste0("R2 = ", as.character(round(r2_4, 4))), x = 0.43, y = 0.65),
            hjust = 0, vjust = 0) +
  ggtitle("Correlation between coefficients from baseline specification\nand Bayesian random effects model")
plot_4

plot_2 | plot_4
