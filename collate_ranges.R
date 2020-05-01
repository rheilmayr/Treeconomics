#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Research support: Melanie Leung
# Project: Treeconomics
# Date: 4/29/19
# Purpose: Collate individual species range maps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set-up --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Packages
library(raster)
library(sf)
library(tidyverse)

# Create file list
# data_dir <- 'remote\\'
ranges_dir <- 'D:\\cloud\\Dropbox\\collaborations\\treeconomics\\species_ranges\\processed_data\\'
species_list <- list.files(ranges_dir)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load and merge data ----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_crs <- 4326
load_shp <- function(sp_name){
  file_name <- paste0(ranges_dir, "\\", sp_name, "\\", sp_name, ".shp")
  sp_sf <- st_read(file_name)
  sp_sf <- sp_sf %>%
    st_transform(new_crs) %>% 
    select(geometry) %>% 
    mutate(sp_code = sp_name)  
  return(sp_sf)
}

sp_sf_list <- lapply(species_list, load_shp)
merged_sf <- do.call(rbind, sp_sf_list)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export new file --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st_write(merged_sf, paste0(ranges_dir, '..\\', 'merged_ranges.shp'), update = TRUE)
         