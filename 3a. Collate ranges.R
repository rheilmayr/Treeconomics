#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Research support: Melanie Leung
# Project: Treeconomics
# Date: 5/27/20
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
wdir <- 'remote\\in\\species_ranges\\'
# wdir <- 'D:\\cloud\\Dropbox\\collaborations\\treeconomics\\species_ranges\\'
ranges_dir <- paste0(wdir, 'processed_data\\')
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
# Add Leslie ranges --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leslie_shp <- paste0(wdir, 'raw_data\\conifers_of_the_world\\ranges\\coniferranges.shp')
leslie_idx <- paste0(wdir, 'leslie_index.csv')
leslie_idx <- read.csv(leslie_idx) %>% 
  filter(source=="leslie")

leslie_list <- leslie_idx %>% 
  pull(conifersoftheworld_taxon)

leslie_ranges <- st_read(leslie_shp)
select_ranges <- leslie_ranges %>% 
  left_join(leslie_idx, by = c("species" = "conifersoftheworld_taxon"))

print(paste0("All records matched: ", dim(leslie_idx)[1]==dim(select_ranges)[1]))

select_ranges <- select_ranges %>% 
  rename(sp_code = ï..species_id) %>% 
  select(sp_code, geometry)

select_ranges <- select_ranges %>% 
  st_transform(crs(merged_sf))

merged_sf <- rbind(merged_sf, select_ranges)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export new file --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st_write(merged_sf, paste0(wdir, 'merged_ranges.shp'), overwrite = TRUE)

# merged_sf %>% st_set_geometry(NULL) %>% select(sp_code) %>% distinct()
