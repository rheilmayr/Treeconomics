#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Create appendix table providing some species-level summary
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
library(broom)
library(purrr)
library(sf)
library(gstat)
library(units)
library(stringr)
library(readxl)

select <- dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. ITRDB data
itrdb_df <- read_csv(paste0(wdir, "/out/dendro/site_summary.csv"))


# 2. Species climates
niche_df <- read_csv(paste0(wdir, "out/climate/clim_niche.csv")) %>% 
  select(-X1)


# # 3. Species ranges
# range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
# range_sf <- st_read(range_file)

# 4. Range metadata
source_df <- read_excel(paste0(wdir, "in/species_ranges/species_summary.xlsx"), sheet = "sources")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize ITRDB sites and merge to range metadata  --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
species_df <- itrdb_df %>% 
  mutate(sp_code = tolower(sp_id)) %>% 
  group_by(sp_code) %>%
  summarize(n_sites = n(),
            n_trees = sum(n_trees, na.rm = TRUE))

species_df <- species_df %>% 
  left_join(niche_df, by = "sp_code")

species_df <- species_df %>% 
  arrange(desc(n_sites))

source_df <- source_df %>% 
  rename(sp_code = species_id) %>% 
  mutate(source = case_when(
    source == "little" ~ "Little, 1971",
    source == "caudullo et al" ~ "Caudullo et al., 2017",
    source == "malyshev" ~ "Malyshev, 2008",
    source == "leslie" ~ "Sundaram et al., 2019" 
  ))

species_df <- species_df %>% 
  left_join(source_df, by = "sp_code")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species without range data for methods narrative --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
species_df <- species_df %>% 
  mutate(generic_sp = str_sub(sp_code,-2,-1)=="sp") %>% 
  filter(generic_sp == FALSE)

species_df <- species_df %>% 
  mutate(missing_range = is.na(source))

## Number of missing species
sp_included = species_df %>%
  group_by(missing_range) %>% 
  summarise(n_sites = n())
sp_included[1,2] / sum(sp_included$n_sites)


## Share of sites from species with missing ranges
n_included = species_df %>%
  group_by(missing_range) %>% 
  summarise(n_sites = sum(n_sites))
  
n_included[1,2] / sum(n_included$n_sites)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Format appendix table --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
species_df <- species_df %>% 
  filter(n_sites > 2) %>% 
  drop_na()

species_df <- species_df %>%
  arrange(sp_code) %>% 
  mutate(pet_mean = round(pet_mean),
         pet_sd = round(pet_sd),
         cwd_mean = round(cwd_mean),
         cwd_sd = round(cwd_sd)) %>% 
  select("Species code" = sp_code,
         "Species name" = spp,
         "Number of sites" = n_sites,
         "Range map source" = source,
         "Mean CWD" = cwd_mean,
         "S.d. CWD" = cwd_sd,
         "Mean PET" = pet_mean,
         "S.d. PET" = pet_sd)

species_df %>% 
  write_csv(paste0(wdir, "tables/a1_sp_summary.csv"))
               