#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 9/21/20
# Purpose: Assess tree-level evolution of drought sensitivity
#
# Input files:
# - 
#
# ToDo:
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
library(dlnm)
library(tidyverse)
library(data.table)
library(tidylog)
library(dplyr)
library(data.table)
library(stringr)
library(dtplyr)
library(ggplotify)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out/dendro/")
dendro_df <- read_csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            rwl = mean(rwl),
            .groups = "drop") %>% 
  mutate(tree_id = paste0(collection_id, "_", tree)) %>% 
  as_tibble()


# 2. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "out/climate/site_an_clim.gz"))
# dendro_df <- dendro_df %>% 
#   left_join(an_site_clim, by = c("collection_id", "year"))


# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)


# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))
dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')


# 5. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create lagged cwd  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlags=15
shock = 1

keep_trees <- dendro_df %>%  # Need to remove trees with too few observations to fit crossbasis
  group_by(tree_id) %>% 
  tally() %>% 
  filter(n>20) %>% 
  pull(tree_id)

dendro_df <- dendro_df %>% 
  filter(tree_id %in% keep_trees)

cwd_lag_names <- str_c("cwd_L", str_pad(0:nlags, width = 2, pad = 0))
pet_lag_names <- str_c("pet_L", str_pad(0:nlags, width = 2, pad = 0))
clim_lagged <- an_site_clim %>% 
  group_by(collection_id) %>% 
  arrange(collection_id, year) %>% 
  do(data.frame(., setNames(shift(.$cwd.an.spstd, 0:nlags), cwd_lag_names))) %>%
  do(data.frame(., setNames(shift(.$pet.an.spstd, 0:nlags), pet_lag_names)))

# Drop NAs from missing lags
clim_lagged <- clim_lagged %>% 
  ungroup() %>% 
  filter(complete.cases(clim_lagged[,grep("cwd_L", colnames(clim_lagged))]),
         complete.cases(clim_lagged[,grep("pet_L", colnames(clim_lagged))]))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create cross-basis matrices  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_lagged = data.frame(clim_lagged[,grep("cwd_L",colnames(clim_lagged))])
pet_lagged = data.frame(clim_lagged[,grep("pet_L",colnames(clim_lagged))])

an_site_clim <- an_site_clim %>% drop_na()
cb_knots=logknots(nlags,3)
cwd_cb = crossbasis(an_site_clim$cwd.an.spstd,
                    lag=c(0,nlags),
                    group = an_site_clim$collection_id,
                    argvar=list(fun = "lin"),
                    arglag=list(fun = "ns",
                                knots=cb_knots,
                                Boundary.knots=c(0,nlags)))
pet_cb = crossbasis(an_site_clim$pet.an.spstd,
                    lag=c(0,nlags),
                    group = an_site_clim$collection_id,
                    argvar=list(fun = "lin"),
                    arglag=list(fun = "ns",
                                knots=cb_knots,
                                Boundary.knots=c(0,nlags)))


cwd_cb_dat <- cwd_cb %>% 
  as_tibble() %>% 
  rename_with(~ stringr::str_replace(.x, 
                                     pattern = "v1", 
                                     replacement = "cwd_cbv1"), 
              matches("v1")) 
pet_cb_dat <- pet_cb %>% 
  as_tibble()  %>% 
  rename_with(~ stringr::str_replace(.x, 
                                     pattern = "v1", 
                                     replacement = "pet_cbv1"), 
              matches("v1")) 

cb_df <- an_site_clim %>% 
  cbind(cwd_cb_dat, pet_cb_dat)


cb_df <- dendro_df %>% 
  inner_join(cb_df, by = c("collection_id", "year")) %>% 
  left_join(ave_site_clim, by = "collection_id") 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prep dnlm regression ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_vars <- dim(cwd_cb_dat)[2]
bylag = 0.1

petnam <- paste("pet_cbv1.l", 1:n_vars, sep="")
cwdnam <- paste("cwd_cbv1.l", 1:n_vars, sep="")
xnam <- append(cwdnam, petnam)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Site-level DNLM  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_dnlm <- function(mod_df){
  failed <- F
  tryCatch(
    expr = {
      fmla <- as.formula(paste("rwi ~ ", paste(xnam, collapse= "+")))
      lagmod=feols(fmla, data=mod_df)
      
      cwd_cp = crosspred(cwd_cb, lagmod, cen = 0, at = -1:1, bylag = bylag, cumul = TRUE)
      pet_cp = crosspred(pet_cb, lagmod, cen = 0, at = -1:1, bylag = bylag, cumul = TRUE)
      
      lag_effects <- tibble(lag = seq(0,nlags,bylag),
                            cwd_effect = cwd_cp$matfit[3,],
                            pet_effect = pet_cp$matfit[3,],
                            cwd_cum =  cwd_cp$cumfit[3,nlags + 1],
                            pet_cum =  pet_cp$cumfit[3,nlags + 1],
                            cwd_cum_se = cwd_cp$cumse[3,nlags + 1])
    },
    error = function(e) {
      message("Returned dlnm error")
      reg_error <<- e[1]
      failed <<- T
    }
  )
  if (failed){
    return(NA)
  }
  return(lag_effects)
}

mod_df <- cb_df %>% 
  group_by(collection_id) %>% 
  nest()

mod_df <- mod_df %>% 
  mutate(dnlm = map(data, run_dnlm))


lagged_effects <- mod_df %>%
  select(collection_id, dnlm) %>%
  unnest(dnlm) %>%
  left_join(ave_site_clim, by = "collection_id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export final cumulative effects (for robustness table run)  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cum_effects <- lagged_effects %>%
  select(collection_id, lag, cwd_cum, pet_cum, cwd_cum_se, cwd.spstd, pet.spstd) %>%
  distinct()  %>% 
  filter(lag == 15)

y_lims <- cum_effects$cwd_cum %>% quantile(c(0.025, 0.975))

cum_mod_df <- cum_effects %>% 
  filter(cwd_cum > y_lims[1],
         cwd_cum < y_lims[2]) %>%
  mutate(errorweights = 1 / cwd_cum_se)

cum_mod_df %>% 
  write_rds(paste0(wdir, "out/first_stage/dnlm_cum_effects.rds"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export lagged effects (for appendix sub-plot)  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lagged_effects <- lagged_effects %>%
  # filter(cwd.spstd>0) %>% 
  select(-cwd_cum, -pet_cum) %>%
  drop_na()

lagged_effects$cwd_quantile <- ntile(lagged_effects$cwd.spstd, 3) %>% as.factor()
lagged_effects$pet_quantile = ntile(lagged_effects$pet.spstd, 3) %>% as.factor()

lagged_effects %>% 
  write_rds(paste0(wdir, "out/first_stage/dnlm_lagged_effects.rmd"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Illustrate appendix figure ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_median_effect <- lagged_effects %>%
  # filter(cwd_quantile==4) %>%
  group_by(lag) %>%
  summarise(med_effect = median(cwd_effect),
            upper = quantile(cwd_effect, 0.66, na.rm = T),
            lower = quantile(cwd_effect, 0.33, na.rm = T))
cwd_plot <- cwd_median_effect %>%
  ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.2) +
  theme_bw() +
  xlim(0, 10) +
  ylim(-0.3, 0.1)
cwd_plot +
  xlab("Lag (years)") +
  ylab("Median effect of CWD=1 shock on RWI")


pet_median_effect <- lagged_effects %>%
  # filter(pet_quantile==4) %>%
  group_by(lag) %>%
  summarise(med_effect = median(pet_effect),
            upper = quantile(pet_effect, 0.66, na.rm = T),
            lower = quantile(pet_effect, 0.33, na.rm = T))
pet_plot <- pet_median_effect %>%
  ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.2) +
  theme_bw() +
  xlim(0, 10) +
  ylim(-0.15, 0.3)

pet_plot +
  xlab("Lag (years)") +
  ylab("Median effect of PET=1 shock on RWI")




