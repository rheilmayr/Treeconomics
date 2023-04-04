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

# 1. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_an_clim.gz"))


# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_ave_clim.gz"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Characterize recent changes --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recent_clim <- an_site_clim %>% 
  filter(year >= 2000) %>% 
  group_by(collection_id) %>% 
  summarise(cwd.recent = mean(cwd.an.spstd),
            pet.recent = mean(pet.an.spstd))

compare_df <- recent_clim %>% 
  left_join(ave_site_clim, by = "collection_id")

compare_df <- compare_df %>% 
  mutate(cwd_change = cwd.recent - cwd.spstd,
         pet_change = pet.recent - pet.spstd)

compare_df %>% 
  ggplot(aes(x = pet.spstd, y = pet.recent)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  xlim(-10,10) +
  ylim(-10,15)


compare_df %>% 
  ggplot(aes(x = cwd.spstd, y = cwd.recent)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  xlim(-10,10) +
  ylim(-10,15)


compare_df %>% 
  ggplot(aes(x = cwd.spstd, y = cwd_change)) +
  geom_point() +
  geom_smooth(method='lm') +
  ylim(-2.5, 2.5)


compare_df %>% 
  ggplot(aes(x = pet.spstd, y = cwd_change)) +
  geom_point() +
  geom_smooth(method='lm') +
  ylim(-2.5, 2.5)


