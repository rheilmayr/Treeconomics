#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Research support: Melanie Leung
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Collate individual species range maps
# 
# Input files:
#   /species_ranges/: Directory of individual species range maps created from lit review.
# 
# Output files:
#   merged_ranges.shp: Compiled shapefile of range maps for all species included in this study
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load packages --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Packages
library(raster)
library(sf)
library(readxl)
library(tidyverse)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load and merge data ----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wdir <- 'remote/in/species_ranges/'

## 1. Range metadata
source_df <- read_excel(paste0(wdir, "species_summary.xlsx"), sheet = "sources")

## 2. Genus labels
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))

source_df <- source_df %>% 
  select(-spp) %>% 
  left_join(sp_info, by = "species_id")

## 3. Load range shapefiles
# Create file list
ranges_dir <- paste0(wdir, 'processed_data/')
species_list <- list.files(ranges_dir, full.names = FALSE)
species_list <- species_list[which(species_list!="desktop.ini")]

new_crs <- 4326
load_shp <- function(sp_name){
  file_name <- paste0(ranges_dir, "/", sp_name, "/", sp_name, ".shp")
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
leslie_idx <- read_csv(leslie_idx) %>% 
  filter(source=="leslie")

leslie_list <- leslie_idx %>% 
  pull(conifersoftheworld_taxon)

leslie_ranges <- st_read(leslie_shp)
select_ranges <- leslie_ranges %>% 
  left_join(leslie_idx, by = c("species" = "conifersoftheworld_taxon"))

select_ranges <- select_ranges %>% 
  drop_na()

print(paste0("All records matched: ", dim(leslie_idx)[1]==dim(select_ranges)[1]))

select_ranges <- select_ranges %>% 
  rename(sp_code = species_id) %>% 
  select(sp_code, geometry)

select_ranges <- select_ranges %>% 
  st_transform(crs(merged_sf))

merged_sf <- rbind(merged_sf, select_ranges)

merged_sf <- merged_sf %>% 
  left_join(source_df, by = c("sp_code" = "species_id"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export new file --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st_write(merged_sf, paste0(wdir, 'merged_ranges.shp'), overwrite = TRUE, append = FALSE)
## Note: Using ArcGIS to dissolve due to sf bugs

source_df %>% write_csv(paste0(wdir, 'species_metadata.csv'))
