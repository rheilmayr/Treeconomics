#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 4/2/23
# Purpose: Generate intro hypotheses figure
#
# Input files:
#
# ToDo:
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wdir <- 'remote/'


# 1. Prediction rasters
rwi_list <- list.files(paste0(wdir, "out/predictions/pred_10000/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep dataframe --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins <- seq(-2.5,2.5,0.2)
labels <- seq(-2.4,2.4,0.2)
bin_length <- length(bins)
dn_sens <- seq(-0.1, 0, (0.1) / (bin_length - 2))
cons_sens <- rep(-0.05, times = bin_length - 1)
res_sens <- seq(0, -0.1, (-0.1) / (bin_length - 2))

exposure_df <- sp_predictions %>% 
  filter(abs(cwd_hist) < 2.5) %>% 
  mutate(cwd_bin = cut(cwd_hist, bins)) %>% 
  mutate(cwd_change = cwd_cmip_end_mean - cwd_cmip_start_mean) %>% 
  group_by(cwd_bin) %>% 
  summarise(cwd_change = mean(cwd_change),
            cwd_end = mean(cwd_cmip_end_mean),
            cwd_start = mean(cwd_cmip_start_mean)) %>% 
  mutate(labels = labels)

exposure_df <- exposure_df %>% 
  mutate(dn_sens = dn_sens,
         cons_sens = cons_sens,
         res_sens = res_sens)

# Sensitivity
sens_df <- exposure_df  %>% 
  pivot_longer(cols = c("dn_sens", "cons_sens", "res_sens")) %>% 
  select(cwd_bin, labels, name, sens = value)

# Exposure
exp_df <- exposure_df %>% 
  select(cwd_bin, labels, cwd_change)

# Vulnerability
vuln_df <- exposure_df  %>% 
  mutate(dn_vuln = (dn_sens * cwd_end) - (dn_sens * cwd_start),
         cons_vuln = (cons_sens * cwd_end) - (cons_sens * cwd_start),
         res_vuln = (res_sens * cwd_end) - (res_sens * cwd_start)) %>% 
  pivot_longer(cols = c("dn_vuln", "cons_vuln", "res_vuln")) %>% 
  select(cwd_bin, labels, name, rwi_change = value)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sens_df %>% 
  ggplot(aes(x = labels, y = sens, fill = name, group = name, color = name)) +
  geom_line() +
  theme_bw()

exp_df %>% 
  ggplot(aes(x = labels, y = cwd_change)) +
  geom_col() +
  theme_bw()

vuln_df %>% 
  ggplot(aes(x = labels, y = rwi_change, fill = name, group = name, color = name)) +
  geom_line() +
  theme_bw()

vuln_df %>%
  ggplot(aes(x = labels, y = rwi_change, fill = name, group = name)) +
  geom_col(position = "dodge")


