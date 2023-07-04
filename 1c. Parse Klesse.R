#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/6/20
# Purpose: Clean ITRDB rwl files
#
# Input files:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#devtools::install_github("tidyverse/tidyr") #if necessary for pivot_wider function
library(tidyverse)
library(purrr)
library(dplR)
library(tidylog)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# Read in Klesse metadata
meta <- read_csv(paste0(wdir, "in/klesse_2018/FIA_TreeRingMeta_Klesse2018.txt"), col_types = cols(.default = "d", CN = "character", PLT_CN = "character"))

# Read in tree ring widths
rwl_df <- read_csv(paste0(wdir, "in/klesse_2018/FIA_TreeRingRW_Klesse2018.txt"), col_types = cols(.default = "d", CN = "character"))

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


rwl_nested <- rwl_df %>%
  select(PLT_CN, species_id, CN, Year, RW) %>% 
  group_by(PLT_CN, species_id) %>% 
  nest()

rwl_nested <- rwl_nested %>% 
  mutate(rwl = map(data, format_rwl),
         rwi = map2(rwl, "Spline", detrend_rwl))



# PLT_CN = "3036494010690", species_id = "PSME"

test <- rwl_nested %>% filter(PLT_CN == "3036494010690") %>% pull(data)
test <- test[[1]]
format_rwl(test)
