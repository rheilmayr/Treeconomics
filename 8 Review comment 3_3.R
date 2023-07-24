#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/11/2023
# Review comment 3.3: 
#   the state-of-the-art way to analyze equation 1 and 2 from a statistical 
#   perspective would be to form a single hierarchical Bayesian model that 
#   would then allow model-based uncertainty quantification and inference 
#   without having to resort to the block bootstrap –which is notoriously 
#   sensitive to specification of blocks and in highly nonstationary (spatial 
#   and/or temporal) settings is quite suspect.  
# 
#
# Approach:
# - Run one integrated model?
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
library(marginaleffects)
library(prediction)
library(margins)
library(lme4)
library(multcomp)



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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add site-level averages from second stage ---------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))

dendro_df <- dendro_df %>% 
  left_join(ave_site_clim, by = "collection_id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run model ---------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

mod_re <- lmer(rwi ~ cwd.an.spstd:pet.spstd + cwd.an.spstd:I(pet.spstd**2) + 
              cwd.an.spstd:cwd.spstd + cwd.an.spstd:I(cwd.spstd**2) +
              pet.an.spstd:pet.spstd + pet.an.spstd:I(pet.spstd**2) + 
              pet.an.spstd:cwd.spstd + pet.an.spstd:I(cwd.spstd**2) + 
              (1 + cwd.an.spstd + pet.an.spstd | collection_id), data = dendro_df)
summary(mod_re)


cwd_median <- ss_df %>% dplyr::select(cwd.spstd) %>% drop_na() %>% pull(cwd.spstd) %>% median()
test_str <- paste0("cwd.an.spstd:cwd.spstd + (2 * cwd.an.spstd:I(cwd.spstd^2) * ", as.character(cwd_median), ") = 0")
# test_str <- paste0("cwd.an.spstd:cwd.spstd = 0")
lincom <- glht(mod_re, linfct = c(test_str))
lincom <- summary(lincom)
coef <- lincom$test$coefficients
se <- lincom$test$sigma

confint(lincom)




## Contrast to primary specification confidence interval
ss_bs <- read_rds(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))
ss_bs %>% pull(cwd_cwd) %>% quantile(c(0.025, 0.975))



mod_fe <- feols(rwi ~ cwd.an.spstd:cwd.spstd | 
                  collection_id, data = dendro_df)




# mod_df <- dendro_df %>% 
#   sample_n(100000)
# 
# ## FE model
# mod_fe <- feols(rwi ~ cwd.an.spstd + cwd.an.spstd:pet.spstd + cwd.an.spstd:pet.spstd**2 + 
#                cwd.an.spstd:cwd.spstd + cwd.an.spstd:cwd.spstd**2 +
#                pet.an.spstd + pet.an.spstd:pet.spstd + pet.an.spstd:pet.spstd**2 + 
#                pet.an.spstd:cwd.spstd + pet.an.spstd:cwd.spstd**2 +
#                cwd.an.spstd * collection_id + pet.an.spstd * collection_id | 
#                collection_id, data = mod_df)
# summary(mod_fe)

ss_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))
ss_modifications <- ss_df %>% 
  select(collection_id, outlier, cwd_errorweights)

dendro_df <- dendro_df %>% 
  left_join(ss_modifications, by = "collection_id")

## FE model - recreate primary results
mod_df <- dendro_df %>% 
  filter(outlier == 0)

site_counts <- mod_df %>% 
  group_by(collection_id) %>% 
  tally()

mod_df <- mod_df %>% 
  left_join(site_counts, by = "collection_id") %>% 
  mutate(weights = 1 / n)

mod <- feols(rwi ~ cwd.an.spstd + cwd.an.spstd:pet.spstd + cwd.an.spstd:pet.spstd**2 + 
               cwd.an.spstd:cwd.spstd + cwd.an.spstd:cwd.spstd**2 +
               pet.an.spstd + pet.an.spstd:pet.spstd + pet.an.spstd:pet.spstd**2 + 
               pet.an.spstd:cwd.spstd + pet.an.spstd:cwd.spstd**2 | 
               collection_id, data = mod_df)
summary(mod)

mod <- feols(rwi ~ cwd.an.spstd + cwd.an.spstd:pet.spstd + 
               cwd.an.spstd:cwd.spstd +
               pet.an.spstd + pet.an.spstd:pet.spstd + 
               pet.an.spstd:cwd.spstd | 
               collection_id, data = mod_df)
summary(mod)

mod <- feols(rwi ~ cwd.an.spstd + pet.an.spstd | collection_id, 
             data = mod_df)
summary(mod)
ss_df %>% filter(outlier==0) %>% select(estimate_pet.an, estimate_cwd.an) %>%  summary()

ss_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, data = ss_df %>% filter(outlier == 0))
summary(ss_mod)

## RE model
mod <- lmer(rwi ~ cwd.an.spstd + cwd.an.spstd:pet.spstd + cwd.an.spstd:I(pet.spstd**2) + 
              cwd.an.spstd:cwd.spstd + cwd.an.spstd:I(cwd.spstd**2) +
               pet.an.spstd + pet.an.spstd:pet.spstd + pet.an.spstd:I(pet.spstd**2) + 
              pet.an.spstd:cwd.spstd + pet.an.spstd:I(cwd.spstd**2) + 
              (1 | collection_id), data = dendro_df)
summary(mod)





mod_slopes <- slopes(mod, newdata = datagrid(cwd.an.spstd = 1, pet.an.spstd = 0, 
                                             pet.spstd = 0, cwd.spstd = cwd_median))



slope_plot <- plot_slopes(mod, variables = "cwd.an.spstd", condition = "cwd.spstd") 
slope_plot +
  xlim(-2.5, 2.5) +
  theme_bw() +
  ylim(-0.25, 0)

slope_plot <- plot_slopes(mod, variables = "cwd.an.spstd", condition = "pet.spstd") 
slope_plot +
  xlim(-2.5, 2.5) +
  theme_bw() +
  ylim(-0.25, 0)

slope_plot <- plot_slopes(mod, variables = "pet.an.spstd", condition = "cwd.spstd") 
slope_plot +
  xlim(-2.5, 2.5) +
  theme_bw() +
  ylim(0.4, 0.8)

slope_plot <- plot_slopes(mod, variables = "pet.an.spstd", condition = "pet.spstd") 
slope_plot +
  xlim(-2.5, 2.5) +
  theme_bw() +
  ylim(0.5, 0.75)



marg_fx_df <- function(mod){
  inc <- 0.1
  min <- -2.5
  max <- 2.5
  cwd_pred_0 <- predictions(mod, newdata = datagrid(cwd.an.spstd = 0, pet.an.spstd = 0, pet.spstd = 0, cwd.spstd = 0)) %>% 
    mutate(variation = "cwd")
  cwd_pred_1 <- predictions(mod, newdata = datagrid(cwd.an.spstd = 1, pet.an.spstd = 0, pet.spstd = 0, cwd.spstd = seq(min,max,inc)))
  return(rbind(cwd_pred, pet_pred))
}

cwd_pred_0 <- predictions(mod, newdata = datagrid(cwd.an.spstd = 0, pet.an.spstd = 0, pet.spstd = 0, cwd.spstd = seq(min,max,inc))) %>% 
  mutate(variation = "cwd")
cwd_pred_1 <- predictions(mod, newdata = datagrid(cwd.an.spstd = 1, pet.an.spstd = 0, pet.spstd = 0, cwd.spstd = seq(min,max,inc)))

adjustment <- cwd_pred_0$estimate[1]
cwd_pred_1$adj_estimate <- cwd_pred_1$estimate - adjustment
cwd_pred_1$adj_conf.high <- cwd_pred_1$conf.high - adjustment
cwd_pred_1$adj_conf.low <- cwd_pred_1$conf.low - adjustment

cwd_mfx_plot <- cwd_pred_1 %>% 
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = adj_estimate)) +
  geom_ribbon(aes(ymin=adj_conf.low, ymax=adj_conf.high), alpha=0.2)
cwd_mfx_plot




# #### REVIEWER COMMENTS - Single model?
# library(lme4)
# library(rstanarm)
# library(broom.mixed)
# 
# ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))
# all_dendro <- dendro_df %>% 
#   left_join(ave_site_clim, by = c("collection_id"))
# 
# all_dendro <- dendro_df %>% 
#   drop_na() %>% 
#   rename(cwd.an = cwd.an.spstd,
#          pet.an = pet.an.spstd,
#          temp.an = temp.an.spstd)
# 
# mod_df <- all_dendro
# 
# mod <- lmer(rwi ~ cwd.an * pet.spstd * cwd.spstd + pet.an + (1 | collection_id), data = mod_df)
# summary(mod)
# 
# mod <- feols(rwi ~ cwd.an + cwd.an:pet.spstd + cwd.an:(pet.spstd**2) + cwd.an:cwd.spstd + cwd.an:(cwd.spstd**2) + 
#                pet.an + pet.an:pet.spstd + pet.an:(pet.spstd**2) + pet.an:cwd.spstd + pet.an:(cwd.spstd**2)| collection_id, data = mod_df)
# summary(mod)
# 
# #### REVIEWER COMMENTS - Single model?
