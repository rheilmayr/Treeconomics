#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/3/2023
# Purpose: Process and detrend FIA ring widths from Klesse et al. 2008
#
# Input files:
#   FIA_TreeRingMeta_Klesse2018.txt: Core metadata from Klesse et al 2018. Shared by John Shaw.
#   FIA_TreeRingRW_Klesse2018.txt: Confidential RW measurements from Klesse et al 2018. Shared by John Shaw
#   FIA Plot tables: Accessed from FIA Datamart
#
# Output:
#   rwi_long_fia.csv: RWI series for each FIA core since 1900
#   site_summary_fia.csv: Metadata (locations) for each FIA collection 
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(purrr)
library(dplR)
library(tidylog)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define paths
# wdir <- 'remote/'
wdir <- 'G:/.shortcut-targets-by-id/10TtqG9P3BY70rcYp-WACmO38J5zBeflA/Treeconomics/Data/'

out_dir <- paste0('remote/1_input_processed/dendro/')

# Read in Klesse metadata
meta <- read_csv(paste0(wdir, "in/klesse_2018/FIA_TreeRingMeta_Klesse2018.txt"), col_types = cols(.default = "d", CN = "character", PLT_CN = "character"))
# meta <- read_csv('G:/.shortcut-targets-by-id/10TtqG9P3BY70rcYp-WACmO38J5zBeflA/Treeconomics/Data/in/klesse_2018/FIA_TreeRingMeta_Klesse2018.txt', col_types = cols(.default = "d", CN = "character", PLT_CN = "character"))

# Read in tree ring widths
rwl_df <- read_csv(paste0(wdir, "in/klesse_2018/FIA_TreeRingRW_Klesse2018.txt"), col_types = cols(.default = "d", CN = "character"))
# rwl_df <- read_csv('G:/.shortcut-targets-by-id/10TtqG9P3BY70rcYp-WACmO38J5zBeflA/Treeconomics/Data/in/klesse_2018/FIA_TreeRingRW_Klesse2018.txt', col_types = cols(.default = "d", CN = "character"))

# Read in FIA plot locations
plot_locs <- tibble("states" = c("CO", "ID", "UT", "WY", "MT"))
get_fia_plot_loc <- function(state){
  # cond_data <- read_csv(paste0(wdir, "in/fia/", state, "_COND.csv"), col_types = cols(.default = "?", CN = "character", PLT_CN = "character"))
  # cond_data <- cond_data %>% 
  #   select(COND_CN = CN, PLT_CN)
  plot_data <- read_csv(paste0(wdir, "in/fia/", state, "_PLOT.csv"), col_types = cols(.default = "?", "CN" = "character", 'LAT' = "d", 'LON' = "d"))
  plot_locs <- plot_data %>% 
    select(PLT_CN = CN, LAT, LON)
  return(plot_locs)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Join datasets --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_locs <- plot_locs %>% 
  mutate(locs = map(states, get_fia_plot_loc)) %>% 
  unnest(locs)

meta <- meta %>% 
  inner_join(plot_locs, by = "PLT_CN")

sp_lookup <- tibble("SPCD" = c(202, 113, 106, 122), "species_id" = c("PSME", "PIFL", "PIED", "PIPO"))

meta <- meta %>%
  left_join(sp_lookup, by = "SPCD")

rwl_df <- rwl_df %>% 
  left_join(meta, by = "CN")

rwl_df <- rwl_df %>% 
  select(RW, Year, CN, PLT_CN, LAT, LON, species_id) %>% 
  drop_na()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Convert to RWI --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format_rwl <- function(rwl_dat){
  rwl_dat <- rwl_dat %>% distinct() # Some CNs are duplicated in the Klesse et al data
  rwl_dat <- rwl_dat %>% 
    pivot_wider(id_cols = Year,
                names_from = CN,
                values_from = RW)
  
  rwl_dat <- rwl_dat %>% 
    arrange(Year) %>% 
    column_to_rownames("Year") %>% 
    as.data.frame()
  return(rwl_dat)
}


detrend_rwl <- function(rwl_dat, method) {
  # Purpose:
  #   Apply dplR to convert ring widths (RWL) to ring width index (RWI)
  # Inputs:
  #   rwl_dat: data.table
  #     Generated from pull_rwl()
  # Outputs:
  #   rwi_dat: data.table
  #     Table of de-trended ring width indices
  rwi_dat <- rwl_dat %>%
    detrend(method = method, make.plot = FALSE, verbose = FALSE) # uses a spline that is 0.67 the series length
  return(rwi_dat)
}

reformat_long <- function(rwi_dat) {
  rwi_dat$year = row.names(rwi_dat)
  rwi_dat <- rwi_dat %>% 
    pivot_longer(names_to = "core_cn", values_to = "rwi", cols = -year)
  return(rwi_dat)
}

rwl_nested <- rwl_df %>%
  select(PLT_CN, species_id, CN, Year, RW) %>% 
  group_by(PLT_CN, species_id) %>% 
  nest()

rwl_nested <- rwl_nested %>% 
  mutate(rwl = map(data, format_rwl),
         rwi = map2(rwl, "Spline", detrend_rwl))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clean up formatting -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rwi_df <- rwl_nested %>% 
  mutate(rwi = map(rwi, reformat_long)) %>% 
  select(-rwl, -data) %>% 
  unnest(rwi) %>% 
  select(core_cn, plot_cn = PLT_CN, species_id, year, rwi)

meta <- meta %>% 
  mutate(tree_id = paste(PLT_CN, as.character(SUBP), as.character(TREE), sep = '_'),
         collection_id = paste("fia", PLT_CN, species_id, sep = "_")) 

tree_ids <- meta %>% 
  select(core_cn = CN, collection_id, tree_id)
  
rwi_df <- rwi_df %>% 
  left_join(tree_ids, by = "core_cn") %>% 
  select(core_cn, plot_cn, collection_id, tree_id, species_id, year, rwi) %>% 
  filter(year > 1900)

fia_site_smry <- meta %>% 
  select(collection_id, plot_cn = PLT_CN, species_id, latitude = LAT, longitude = LON) %>% 
  unique()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Write data -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_csv(rwi_df, paste0(out_dir, "rwi_long_fia.csv"))
write_csv(fia_site_smry, paste0(out_dir, "site_summary_fia.csv"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Optional - visualize data -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(sf)
library(tmap)
tmap_mode("view")

fia_site_smry <- fia_site_smry %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
tm_shape(fia_site_smry) +
  tm_dots("species_id")
