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
wdir <- 'NewRemote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std.csv')) #%>%
#select(-X1)

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))


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
species_df <- flm_df %>%
  filter(species_id == "quma")

sort(unique(flm_df$species_id))

#d'orangville study: filter(species_id %in%c("quve", "lasi", "asch", "litu", "qual", "thoc", "pire", "pcgl"))
#BEPA
#PCGL
#Larix sibirica
##Pinus sibirica
meanpet_all=round(mean(species_df$pet.spstd), digits = 3)
meancwd_all=round(mean(species_df$cwd.spstd), digits=3)
species_df %>%
  #filter(estimate_cwd.an>-.75) %>% 
  filter(estimate_cwd.an<2) %>% 
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
  geom_point()+
  guides(color=F)+
  geom_smooth(method="lm")+
  ggtitle(paste("PET", meanpet_all, species_df$species_id, "CWD", meancwd_all))


unique(species_df$collection_id)


## Look at specific site and its nearby sites
site_id <- "ID006"
species=unique(species_df$species_id)

proximate_ids <- block_list[site_id] %>% unlist(use.names = FALSE)

local_df <- species_df %>%
  filter(collection_id %in% proximate_ids)

meanpet=round(mean(local_df$pet.spstd), digits = 3)
meancwd=round(mean(local_df$pet.spstd), digits = 3)

local_df %>%
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd, color=cwd.spstd)) +
  geom_point(size=2)+
  guides(color=F)+
  geom_smooth(method="lm")+
  ggtitle(paste("PET", meanpet, site_id, species, "CWD", meancwd))

