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
library(tidyr)
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(lfe)
library(broom)
library(purrr)
library(dplR)
library(stringi)
library(varhandle)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

data_dir <- paste0(wdir, 'raw_in\\itrdb\\rwi\\')
header_files <- list.files(data_dir, pattern = 'noaa.rwl')
# filenames <- list.files(data_dir, pattern = '.rwl')
# header_bool <- stri_detect_fixed(filenames, "-noaa")
# header_files <- filenames[which(header_bool)]
# rwl_files <- filenames[which(!header_bool)]

out_dir <- paste0(wdir, 'out\\dendro\\')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parse headers --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
findVar <- function (header, search_string){
  row_id <- which(startsWith(header, search_string))
  row <- header[row_id]
  split <- strsplit(row, ": ")
  value <- split[[1]][2]
  return(value)
}

idMetric <- function(data_id){
  last_chr <- str_sub(data_id,-1)
  if (varhandle::check.numeric(last_chr)) { # For datasets with no suffix, data is total ring width
    metric <- "total_ring_width"
  }
  else {
    if (last_chr=='e') {
      metric <- "early_ring_width"
    } else if (last_chr=='l') {
      metric <- "late_ring_width"
    } else {
      metric <- "other"
    }
  }
  return(metric)
}


site_list <- list()

for (file in header_files) {
  print(file)
  txt <- readLines(paste0(data_dir, file))
  last_row <- which(txt == "# Variables ")
  header <- txt[1:(last_row-1)]
  search_list <- list(nLat = "# Northernmost_Latitude:",
                      sLat = "# Southernmost_Latitude:",
                      wLon = "# Westernmost_Longitude:",
                      eLon = "# Easternmost_Longitude:",
                      location = "# Location:",
                      elevation = "# Elevation:",
                      site_name = "# Site_Name:",
                      collection_id = "# Collection_Name:",
                      start_year = "# Earliest_Year:",
                      end_year = "# Most_Recent_Year:",
                      time_unit = "# Time_Unit:",
                      sp_scient = "# Species_Name:",
                      sp_common = "# Common_Name:",
                      sp_id = "# Tree_Species_Code:")
  parsed <- lapply(search_list, function(x) {findVar(header, x)})
  parsed$filename <- file
  parsed <- as_tibble(parsed)
  site_list[[file]] <- parsed
}

site_data <- bind_rows(site_list)

site_data <- site_data %>% 
  mutate(nLat = as.numeric(nLat),
         sLat = as.numeric(sLat),
         wLon = as.numeric(wLon),
         eLon = as.numeric(eLon),
         elevation = as.numeric(elevation),
         start_year = as.numeric(start_year),
         end_year = as.numeric(end_year),
         site_name = str_squish(site_name),
         collection_id = str_squish(collection_id))

site_data <- site_data %>% 
  rowwise() %>% 
  mutate(latitude = mean(c(nLat, sLat), na.rm = TRUE),
         longitude = mean(c(eLon, wLon), na.rm = TRUE),
         filename = str_remove(filename, "-noaa.rwl"),
         metric = idMetric(filename))

total_sites <- site_data %>% 
  filter(metric == "total_ring_width") %>% 
  pull(collection_id)

early_sites <- site_data %>% 
  filter(metric == "early_ring_width") %>% 
  pull(collection_id)

late_sites <- site_data %>% 
  filter(metric == "late_ring_width") %>% 
  pull(collection_id)

early_sites <- early_sites[which(!early_sites %in% total_sites)]
late_sites <- late_sites[which(!late_sites %in% total_sites)]
sum_sites <- early_sites[which(early_sites %in% late_sites)]
keep_sites <- c(total_sites, sum_sites)

site_summary <- site_data %>% 
  select(collection_id, site_name, location, latitude, longitude, elevation, start_year, end_year, time_unit, sp_id, sp_scient, sp_common) %>% 
  distinct() %>% 
  filter(collection_id %in% keep_sites)
write_csv(site_summary, paste0(out_dir, "site_summary.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Process rwl files --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_to_filename <- function(site, suffix){
  site <- str_to_lower(site)
  site <- paste0(site, suffix, '.rwl')
}


open_rwl <- function(file){
  caught_error <- NA
  caught_warning <- NA
  rwl <- NA
  ids <- NA
  tryCatch(
    expr = {
      rwl <- read.tucson(paste0(data_dir, file))
      ids <- autoread.ids(rwl)
    },
    error = function(e){ 
      message("Returned error on site ", file)
      print(e)
      caught_error <<- e
    },
    warning = function(w){
      message("Returned warning on file ", file)
      print(w)
      caught_warning <<- w
    }
  )
  out <- tibble(rwl = list(rwl), ids = list(ids), error = caught_error[1], warning = caught_warning[1])
  return(out)
}


# read_ids <- function(rwl){
#   caught_error <- NA
#   caught_warning <- NA
#   ids <- NA
#   tryCatch(
#     expr = {
#       ids <- autoread.ids(rwl)
#     },
#     error = function(e){ 
#       message("Returned error on site ", file)
#       print(e)
#       caught_error <<- e
#     },
#     warning = function(w){
#       message("Returned warning on file ", file)
#       print(w)
#       caught_warning <<- w
#     }
#   )
#   out <- tibble(ids = list(ids), error = caught_error[1], warning = caught_warning[1])
#   return(out)
# }


rwl_files <- sapply(total_sites, site_to_filename, suffix = '')
rwl_data <- tibble::enframe(rwl_files, name = "collection_id", value = "file")
rwl_data <- rwl_data[1:50,]
rwl_data$rwl <- map(rwl_data$file, open_rwl) 
rwl_data <- rwl_data %>% unnest(rwl)
warnings <- rwl_data %>% 
  unnest(warning) %>% 
  select(collection_id, warning)
errors <- rwl_data %>% 
  unnest(error) %>% 
  select(collection_id, error)

error_collections <- errors %>% pull(collection_id)
warning_collections <- warnings %>% pull(collection_id)
error_collections <- c(error_collections, warning_collections)

clean_data <- rwl_data %>% 
  select(-error, -warning) %>% 
  filter(!(collection_id %in% error_collections))

ids <- clean_data[1,4] %>% pull()
rwi <- clean_data[1,3] %>% pull()
rwi <- rwi[[1]]
rwi <- rwi %>% 
  rownames_to_column("year")
rwi <- rwi %>% 
  pivot_longer(cols = -year, names_to = "core_id", values_to = "width") %>% 
  arrange(core_id, year)

ids <- ids[[1]]
ids <- ids %>% 
  rownames_to_column("core_id")

rwi <- rwi %>% 
  left_join(ids, by = "core_id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOTES:
# NOAA header files  missing species data block and collection name. Added manually
#   cana575-cana586; paki041
#
# Collections ARGE and ME recoded as ARGE999 and ME999
# 
# 
# Updates over original dataset -
#   Added data from more recent ITRDB uploads
#   Correctly dealing with early / latewood observations
#   Separating out multiple cores from a single tree
#   Fewer parsing errors?
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
