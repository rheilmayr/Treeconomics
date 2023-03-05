#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Create predictions of growth impacts from climate change
#
# Input files:
# - ss_mod: R model object saved from Second stage
# - 
#
# ToDo:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(fixest)
library(patchwork)
library(tidyverse)
library(dtplyr)
library(prediction)
library(tictoc)
library(tmap)
library(tidylog)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv')) %>%
  select(-X1)

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_ave_clim.gz"))


# 3. Site information
site_df <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

flm_df <- flm_df %>% 
  left_join(ave_site_clim, by = c("collection_id")) %>%
  left_join(site_df, by = c("collection_id"))

# 4. Spatially proximate blocks of sites
load(file=paste0(wdir,"out/spatial_blocks.Rdat"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot sensitivity --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_df <- flm_df %>%
  filter(species_id == "pipo")


plot_df %>%
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
  geom_point()

## Look at specific site and its nearby sites
site_id <- "AZ036"
site_id <- "NM512"
site_id <- "CA525"

proximate_ids <- block_list[site_id] %>% unlist(use.names = FALSE)

plot_df <- plot_df %>%
  filter(collection_id %in% proximate_ids)

plot_df %>%
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
  geom_point()
