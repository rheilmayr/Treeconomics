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
mod <- feols(rwi ~ cwd.an.spstd + cwd.an.spstd:pet.spstd + cwd.an.spstd:pet.spstd**2 + cwd.an.spstd:cwd.spstd + cwd.an.spstd:cwd.spstd**2 +
               pet.an.spstd + pet.an.spstd:pet.spstd + pet.an.spstd:pet.spstd**2 + pet.an.spstd:cwd.spstd + pet.an.spstd:cwd.spstd**2 | collection_id + species_id, data = dendro_df)
# mod <- lmer(rwi ~ cwd.an.spstd + cwd.an.spstd:pet.spstd + cwd.an.spstd:pet.spstd**2 + cwd.an.spstd:cwd.spstd + cwd.an.spstd:cwd.spstd**2 +
#                pet.an.spstd + pet.an.spstd:pet.spstd + pet.an.spstd:pet.spstd**2 + pet.an.spstd:cwd.spstd + pet.an.spstd:cwd.spstd**2 + (1 | collection_id), data = dendro_df)

summary(mod)


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

preds <- marg_fx_df(mod)
adjustment <- cwd_pred_0$estimate[1]
cwd_pred_1$adj_estimate <- cwd_pred_1$estimate - adjustment
cwd_pred_1$adj_conf.high <- cwd_pred_1$conf.high - adjustment
cwd_pred_1$adj_conf.low <- cwd_pred_1$conf.low - adjustment

cwd_mfx_plot <- cwd_pred_1 %>% 
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = adj_estimate)) +
  geom_ribbon(aes(ymin=adj_conf.low, ymax=adj_conf.high), alpha=0.2)
cwd_mfx_plot
