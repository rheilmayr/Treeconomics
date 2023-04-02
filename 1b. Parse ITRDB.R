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

### Pre-remove sites with spans outside of study period?
# site_data <- site_data %>% 
#   filter((end_year<1900) | (end_year>3000))


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
  tryCatch(
    expr = {
      rwl <- read.tucson(paste0(data_dir, file)) # NOTE: Any way to extract RWL units at this point?
    },
    error = function(e){
      tryCatch(
        expr = {
          path_nohead <- remove_head(file)
          rwl <<- read.tucson(paste0(path_nohead))
        },
        error = function(e){
          message("Returned error on site ", file)
          print(e)
          caught_error <<- e
        }
      )
    }
    # ,
    # warning = function(w){
    #   message("Returned warning on file ", file)
    #   print(w)
    #   caught_warning <<- w
    # }
  )
  # out <- tibble(rwl = list(rwl), read_error = list(caught_error), read_warning = list(caught_warning))
  out <- tibble(rwl = list(rwl), read_error = list(caught_error))
  return(out)
}


remove_head <- function(file){
  path <- paste0(data_dir, file)
  temp_path <- paste0(data_dir, "temp_nohead.rwl")
  txt <- readLines(path)
  txt <- txt[-(1:3)]
  fileConn<-file(temp_path)
  writeLines(txt, fileConn)
  close(fileConn)
  return(temp_path)
}

read_ids <- function(rwl){
  # rwl <- rwl[[1]]
  caught_error <- NA
  caught_warning <- NA
  ids <- NA
  tryCatch(
    expr = {
      ids <- autoread.ids(rwl)
    },
    error = function(e){ 
      message("Returned error on site ", file)
      print(e)
      caught_error <<- e
    }
    # ,
    # warning = function(w){
    #   message("Returned warning on file ", file)
    #   print(w)
    #   caught_warning <<- w
    #   ids <<- autoread.ids(rwl)
    # }
  )
  # out <- tibble(ids = list(ids), id_warning = list(caught_warning), id_error = list(caught_error))
  out <- tibble(ids = list(ids), id_error = list(caught_error))
  return(out)
}


separate_errors <- function(rwl_data){
  errors <- rwl_data %>% 
    filter((!is.na(read_error)) | (!is.na(id_error)))

  error_collections <- errors %>% pull(collection_id)
  
  clean_data <- rwl_data %>% 
    select(-c(read_error, id_error)) %>% 
    filter(!(collection_id %in% error_collections))
  errors <- rwl_data %>% 
    select(-rwl, -ids) %>% 
    filter(collection_id %in% error_collections)
  return(list(clean_data = clean_data, errors = errors))
}

rwl_files <- sapply(total_sites, site_to_filename, suffix = '')
rwl_data <- tibble::enframe(rwl_files, name = "collection_id", value = "file")
rwl_data$rwl <- map(rwl_data$file, open_rwl) 
rwl_data <- rwl_data %>% unnest(rwl)
rwl_data <- rwl_data %>% 
  mutate(ids = map(rwl_data$rwl, read_ids))
rwl_data <-rwl_data %>% unnest(ids)
rwl_results <- separate_errors(rwl_data)
rwl_clean_data <- rwl_results$clean_data

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
erwl_data <- erwl_data %>% 
  mutate(ids = map(erwl_data$rwl, read_ids))
erwl_data <- erwl_data %>% unnest(ids)
erwl_results <- separate_errors(erwl_data)
e_clean_data <- erwl_results$clean_data
e_clean_data <- e_clean_data %>% 
  rename(erwl = rwl)

lrwl_files <- sapply(sum_sites, site_to_filename, suffix = 'l')
lrwl_data <- tibble::enframe(lrwl_files, name = "collection_id", value = "file")
lrwl_data$rwl <- map(lrwl_data$file, open_rwl) 
lrwl_data <- lrwl_data %>% unnest(rwl)
lrwl_data <- lrwl_data %>% 
  mutate(ids = map(lrwl_data$rwl, read_ids))
lrwl_data <- lrwl_data %>% unnest(ids)
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

clean_data <- rbind(rwl_clean_data, s_clean_data)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Detrend rwl to create rwi --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

clean_data$rwi <- map2(clean_data$rwl, "Spline", detrend_rwl)
clean_data$rwi_nb <- map2(clean_data$rwl, "ModNegExp", detrend_rwl)
clean_data$rwi_ar <- map2(clean_data$rwl, "Ar", detrend_rwl)


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

clean_data$rwl_long <- map2(clean_data$rwl, clean_data$ids, pivot_rw)
clean_data$rwi_long <- map2(clean_data$rwi, clean_data$ids, pivot_rw)
clean_data$rwi_nb_long <- map2(clean_data$rwi_nb, clean_data$ids, pivot_rw)
clean_data$rwi_ar_long <- map2(clean_data$rwi_ar, clean_data$ids, pivot_rw)

clean_data_ar <- clean_data %>% 
  select(collection_id, rwi_ar_long) %>% 
  unnest(c(rwi_ar_long), names_sep = '_') %>% 
  select(collection_id, year = rwi_ar_long_year, tree = rwi_ar_long_tree, core = rwi_ar_long_core, 
         core_id = rwi_ar_long_core_id, rwi_ar = rwi_ar_long_width)
# %>% 
#   mutate(site = replace_na(1))

clean_data <- clean_data %>% 
  select(collection_id, rwi_long, rwl_long, rwi_nb_long) %>% 
  unnest(c(rwi_long, rwl_long, rwi_nb_long), names_sep = '_') %>% 
  select(collection_id, year = rwi_long_year, tree = rwi_long_tree, core = rwi_long_core, 
         core_id = rwi_long_core_id, rwi = rwi_long_width, rwl = rwl_long_width, 
         rwi_nb = rwi_nb_long_width) 
# %>% 
#   mutate(site = replace_na(1))

clean_data <- clean_data %>% 
  left_join(clean_data_ar, by = c("collection_id", "year", "tree", "core", 'core_id'))

# test_that("No collections contain multiple sites",{
#   expect_equal(clean_data %>% select(site) %>% unique() %>% pull(), 1)
# })

clean_data <- clean_data  %>% 
  mutate(tree = vctrs::vec_group_id(tree),
         core = vctrs::vec_group_id(core),
         year = as.integer(year)) %>% 
  select(collection_id, core_id, tree, core, year, rwi, rwl, rwi_ar, rwi_nb)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tabulate sites that are dropped from analysis  -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_smry <- site_data %>% 
  select(collection_id, location, latitude, longitude, elevation, sp_id) %>% 
  unique()

l_parse_errors <- lrwl_results$errors %>% pull(collection_id)
e_parse_errors <- erwl_results$errors %>% pull(collection_id)
total_parse_errors <- rwl_results$errors %>% pull(collection_id)

data_smry <- clean_data %>% 
  group_by(collection_id) %>% 
  summarise(obs_start_year = min(year),
            obs_end_year = max(year),
            n_trees = n_distinct(tree))

site_smry <- site_smry %>% 
  left_join(data_smry, by = "collection_id")

site_smry <- site_smry %>% 
  mutate(valid_dates = obs_end_year >= 1901 | obs_end_year > 2020,
         valid_rwl = collection_id %in% c(sum_sites, total_sites),
         valid_parse = !(collection_id %in% c(l_parse_errors, e_parse_errors, total_parse_errors)))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Filter to usable data  --------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keep_sites <- site_smry %>% 
  filter(valid_dates == TRUE) %>% 
  pull(collection_id)

export_data <- clean_data %>% 
  filter(collection_id %in% keep_sites,
         year >= 1901)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Final data quality checks  --------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_check <- export_data %>% 
  group_by(collection_id, tree, core, year) %>% 
  summarise(n = n())
errors <- n_check %>%
  filter(n>1) %>% 
  ungroup() %>% 
  select(collection_id, tree) %>% 
  unique()

# test_that("No core has multiple observations in a year",{
#   expect_equal(dim(errors)[1], 0)
# })

# Drop trees (n=22) that have multiple observations in a single year
error_list <- errors %>%   
  mutate(temp_id = paste0(collection_id, "_", tree)) %>% 
  pull(temp_id)
error_collections <- errors %>% select(collection_id) %>% unique()

export_data <- export_data %>% 
  mutate(temp_id = paste0(collection_id, "_", tree)) %>%
  filter(!(temp_id %in% error_list)) %>% 
  select(-temp_id)

test_that("Core ids have been assigned",{
  expect_equal(is.na(export_data$core_code) %>% sum(), 0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export data and site summary  ------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_csv(export_data, paste0(out_dir, "rwi_long.csv"))

write_csv(site_smry, paste0(out_dir, "site_summary.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize number of observations ---------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# number of original files
site_data$collection_id %>% unique() %>% length()

# number of sites dropped due to "other" rwl files
site_smry %>% 
  filter(valid_rwl == FALSE) %>% 
  pull(collection_id) %>% 
  length()

# number of sites dropped due to parsing errors
site_smry %>% 
  filter(valid_parse == FALSE) %>% 
  pull(collection_id) %>% 
  length()

# number of sites dropped due to no data in period of interest
site_smry %>% 
  filter(valid_dates == FALSE) %>% 
  pull(collection_id) %>% 
  length()

# number of sites dropped due to some error
site_smry %>% 
  filter(!valid_dates | !valid_parse | !valid_rwl) %>% 
  pull(collection_id) %>% 
  length()

# number of observations
clean_data %>% 
  pull(collection_id) %>% 
  length()

# number of trees
clean_data %>% 
  select(collection_id, tree) %>% 
  distinct() %>%
  pull(collection_id) %>% 
  length()

# number of sites in final dataset
clean_data %>% 
  select(collection_id) %>% 
  distinct() %>%
  pull(collection_id) %>% 
  length()
  





# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOTES:
# NOAA header files  missing species data block and collection name. Added manually
#   cana575-cana586; paki041
#
# Series with inconsistent names edited manually
#   
#
# Collections ARGE and ME recoded as ARGE999 and ME999
#   AK020; CV05SW: After 1860 shifted to lowercase
#   AK077; ZC-5-6W: After 1760 was missing first digit (Z) of code
#   AK109; FTR25WAS: Typo corrected after 1810
#   AR073; PB07C: First row missing first digit (O)
#   AR077; p160t046: Last row has typo
#   BRIT10; all series: Removed 0s after precision indicators
#   CA527; 041031: Added precision
#   CA599; TPR17B: Added precision
#   CA621; KER17C: Corrected typo, Shouldve been KER17B
#   CA645; B2707A
#   CANA101: Added space to first line
#   CANA151: Corrected lots of typos
#   CANA161; 1s012a
#   CANA260; p0504b
#   CANA273; bs496b2: Added precision
#   CANA275: Corrected typos
#   CANA521: Corrected typo
#   
# Collections with duplicate series IDs - removed one tree (when possible, tree that didn't overlap with period of interest)
# Some have invalidly long ids - removed those
#   AK027
#   AK030
#   AK048
#   AK050

## CAN BE USED TO DEBUG RWI FILES WITH ERRORS
# i = 4
# file = list[i]
# file
# read.tucson(paste0(data_dir, file, ".rwl"))

# Updates over original dataset -
#   Added data from more recent ITRDB uploads
#   Summing early / latewood observations
#   Separating out multiple cores from a single tree
#   Fewer parsing errors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
