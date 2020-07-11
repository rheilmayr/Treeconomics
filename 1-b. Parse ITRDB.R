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
library(testthat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

data_dir <- paste0(wdir, 'in\\itrdb\\rwi\\')
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
# Process total rwl files --------------------------------------------------------
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

separate_errors <- function(rwl_data){
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
  errors <- rwl_data %>% 
    select(-rwl, -ids) %>% 
    filter(collection_id %in% error_collections) %>% 
    unnest(c(warning, error))
  return(list(clean_data = clean_data, errors = errors))
}

rwl_files <- sapply(total_sites, site_to_filename, suffix = '')
rwl_data <- tibble::enframe(rwl_files, name = "collection_id", value = "file")
rwl_data$rwl <- map(rwl_data$file, open_rwl) 
rwl_data <- rwl_data %>% unnest(rwl)
rwl_results <- separate_errors(rwl_data)
clean_data <- rwl_results$clean_data


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Process earlywood/latewood files ---------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_rwl <- function(e_dat, l_dat){
  failed <- F
  error <- NA
  tryCatch(
    expr = {
      expect_identical(rownames(e_dat), rownames(l_dat))
      expect_identical(colnames(e_dat), colnames(l_dat))
    },
    error = function(e){ 
      failed <<- T
      error <<- e[1]
    }
  )
  if (failed){
    return(tibble(rwl = list(NULL), error = error))
    # return(list(rwl = NULL, error = error))
  } 
  else {
    e_na <- is.na(e_dat)
    l_na <- is.na(l_dat)
    both_na <- (e_na & l_na)
    e_dat[is.na(e_dat)] <- 0
    l_dat[is.na(l_dat)] <- 0
    sum_dat <- e_dat + l_dat
    sum_dat[both_na] <- NA
    rwl <- as.rwl(sum_dat)
    return(tibble(rwl = list(rwl), error = error))
    }
}

erwl_files <- sapply(sum_sites, site_to_filename, suffix = 'e')
erwl_data <- tibble::enframe(erwl_files, name = "collection_id", value = "file")
erwl_data$rwl <- map(erwl_data$file, open_rwl) 
erwl_data <- erwl_data %>% unnest(rwl)
erwl_results <- separate_errors(erwl_data)
e_clean_data <- erwl_results$clean_data
e_clean_data <- e_clean_data %>% 
  rename(erwl = rwl)

lrwl_files <- sapply(sum_sites, site_to_filename, suffix = 'l')
lrwl_data <- tibble::enframe(lrwl_files, name = "collection_id", value = "file")
lrwl_data$rwl <- map(lrwl_data$file, open_rwl) 
lrwl_data <- lrwl_data %>% unnest(rwl)
lrwl_results <- separate_errors(lrwl_data)
l_clean_data <- lrwl_results$clean_data
l_clean_data <- l_clean_data %>% 
  rename(lrwl = rwl)

s_clean_data <- e_clean_data %>% 
  select(-ids, -file) %>% 
  inner_join(l_clean_data, by = "collection_id")

s_clean_data$rwl <- map2(s_clean_data$erwl, s_clean_data$lrwl, sum_rwl)
s_clean_data <- s_clean_data %>% 
  unnest(rwl) %>% 
  filter(map_lgl(error, is.na)) %>% #currently all files are merging correctly
  select(collection_id, file, rwl, ids)

clean_data <- rbind(clean_data, s_clean_data)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Detrend rwl to create rwi --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detrend_rwl <- function(rwl_dat) {
  # Purpose:
  #   Apply dplR to convert ring widths (RWL) to ring width index (RWI)
  # Inputs:
  #   rwl_dat: data.table
  #     Generated from pull_rwl()
  # Outputs:
  #   rwi_dat: data.table
  #     Table of de-trended ring width indices
  rwi_dat <- rwl_dat %>%
    detrend(method = "Spline", make.plot = FALSE, verbose = FALSE) # uses a spline that is 0.67 the series length
  return(rwi_dat)
}

clean_data$rwi <- map(clean_data$rwl, detrend_rwl)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pivot to long --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pivot_rw <- function(rw, ids){
  rw <- rw %>% as.data.frame()
  rw <- rw %>% 
    rownames_to_column("year")
  rw <- rw %>% 
    pivot_longer(cols = -year, names_to = "core_id", values_to = "width") %>% 
    arrange(core_id, year)
  
  ids <- ids %>% as.data.frame()
  ids <- ids %>% 
    rownames_to_column("core_id")
  
  rw <- rw %>% 
    left_join(ids, by = "core_id") %>% 
    drop_na()
  
  return(rw)
}

clean_data$rwi_long <- map2(clean_data$rwi, clean_data$ids, pivot_rw)
clean_data <- clean_data %>% 
  select(collection_id, rwi_long) %>% 
  unnest(rwi_long) %>% 
  mutate(site = replace_na(1)) %>% 
  select(collection_id, core_id, year, tree, core, site, width)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export results and error log --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_csv(clean_data, paste0(out_dir, "rwi_long.csv"))

error_log <- rbind(rwl_results$errors, erwl_results$errors, lrwl_results$errors)
write_csv(clean_data, paste0(out_dir, "itrdb_error_log.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize number of observations ---------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_obs <- clean_data %>% dim()
n_trees <- clean_data %>% 
  select(collection_id, tree) %>% 
  distinct() %>% 
  dim()
n_collections <- clean_data %>% 
  select(collection_id) %>% 
  distinct() %>% 
  dim()
n_usable_obs <- clean_data %>% 
  filter(year>1900) %>% 
  select(year, collection_id, tree) %>% 
  distinct() %>% 
  dim()


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
#   Summing early / latewood observations
#   Separating out multiple cores from a single tree
#   Fewer parsing errors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
