#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/6/20
# Purpose: Clean ITRDB rwl files
#
# Input files:
#   - ITRDBmetadata12January2022.xlsx: Metadata catalog of the ITRDB v8. Used 
#       to index sites
#   - /itrdb-v713-cleaned-rwl/: Directory holding the reformatted ITRDB (rITRDB)
#       Downloaded from https://www.ncei.noaa.gov/access/paleo-search/study/25570
#   - /rwi/: Directory holding ITRDB noaa.rwl files from ITRDB v7.22.
#       Downloaded jsing 1a. Pull ITRDB.R
# 
# Output files:
#   - rwi_long.csv: Long dataset (collection by tree by year) of all detrended ITRDB data 
#   - site_summary.csv: Summary of each ITRDB site included in rwi_long
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
library(readxl)
library(tidylog)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

## ITRDB metadata (v8)
itrdb_meta <- read_excel(paste0(wdir, 'in/itrdb/ITRDBmetadata12January2022.xlsx'))

## Revisited ITRDB (v7.13)
ritrdb_data_dir <- paste0(wdir, 'in/itrdb_zhao_corrected/AppendixS1/Cleaned datasets/itrdb-v713-cleaned-rwl/')
ritrdb_meta <- read_csv(paste0(ritrdb_data_dir, "rwl_metadata.csv"))

## Original ITRDB data (v7.22)
itrdb_data_dir <- paste0(wdir, 'in/itrdb/rwi/')
header_files <- list.files(itrdb_data_dir, pattern = 'noaa.rwl')

## Define output directory
out_dir <- paste0(wdir, 'out\\dendro\\')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Narrow based on metadata --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itrdb_df <- itrdb_meta %>% 
  filter(LastYear > 1901,
         !is.na(Species)) %>% 
  select(collection_id = ITRDB_Code, latitude = Latitude, longitude = Longitude, elevation_itrdb = Elevation, species_id = Species) %>% 
  drop_na()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ID sites in revisited ITRDB --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ritrdb_df <- ritrdb_meta %>% 
  filter(type %in% c("Ring Width","Earlywood Width", "Latewood Width"),
         missing_file == 0) 

ritrdb_df <- ritrdb_df %>% 
  mutate(collection_id = gsub('e$|l$|w$', '', id),
         collection_id = str_to_upper(collection_id)) %>% 
  select(collection_id, region, id, type) %>% 
  mutate(filepath = paste0(ritrdb_data_dir, region, "/", id, ".rwl"))

total_sites <- ritrdb_df %>% 
  filter(type == "Ring Width") %>% 
  pull(collection_id)

early_sites <- ritrdb_df %>% 
  filter(type == "Earlywood Width") %>% 
  pull(collection_id)

late_sites <- ritrdb_df %>% 
  filter(type == "Latewood Width") %>% 
  pull(collection_id)

early_sites <- early_sites[which(!early_sites %in% total_sites)]
late_sites <- late_sites[which(!late_sites %in% total_sites)]
sum_sites <- early_sites[which(early_sites %in% late_sites)]
keep_sites <- c(total_sites, sum_sites)

problem_ritrdb_ids <- c()

total_sites_ritrdb <- ritrdb_df %>% 
  filter(collection_id %in% total_sites,
         type == "Ring Width") 

## Drop some duplicated and incorrect files
total_sites_ritrdb <- total_sites_ritrdb %>% 
  filter(!((collection_id %in% c("GEOR002", "GEOR003", "GEOR004", "GEOR005", "GEOR006", "GEOR007", "GEOR008", "GEOR009", "GEOR010")) & (region == "europe"))) %>% 
  filter(!(id %in% c("ausl038w", "ausl039w", "ausl043w" )))

total_sites_ritrdb <- total_sites_ritrdb %>% 
  select(collection_id, t_rwl = filepath) %>% 
  mutate(e_rwl = NaN,
         l_rwl = NaN)

sum_sites_ritrdb <- ritrdb_df %>% 
  filter(collection_id %in% sum_sites,
         type %in% c("Earlywood Width", "Latewood Width"))

sum_sites_ritrdb <- sum_sites_ritrdb %>% 
  mutate(type = ifelse(type == "Earlywood Width", "e_rwl", "l_rwl")) %>% 
  pivot_wider(id_cols = collection_id, names_from = type, values_from = filepath) %>% 
  mutate(t_rwl = NaN) %>% 
  select(collection_id, t_rwl, e_rwl, l_rwl)

ritrdb_df <- total_sites_ritrdb %>% 
  rbind(sum_sites_ritrdb) %>% 
  mutate(datasource = 'ritrdb_7.13')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Look for remaining sites in ITRDB download  ----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ritrdb_sites <- ritrdb_df %>% pull(collection_id)
itrdb_sites <- itrdb_df %>% pull(collection_id)
missing_collections <- itrdb_sites[which(!(itrdb_sites %in% ritrdb_sites))]
  
header_files <- list.files(itrdb_data_dir, pattern = 'noaa.rwl')

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
  txt <- readLines(paste0(itrdb_data_dir, file))
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

## Focus on sites that still need data
site_data <- site_data %>% 
  filter(collection_id %in% missing_collections)

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
test_that("No sum sites are left after total sites are already accounted for",
          expect_equal(length(sum_sites),0))
keep_sites <- total_sites

site_data <- site_data %>% 
  filter(collection_id %in% keep_sites,
         metric == "total_ring_width") %>% 
  select(collection_id, metric, filename)

site_data <- site_data %>% 
  mutate(datasource = "itrdb_7.22")

site_data <- site_data %>% 
  mutate(e_rwl = NaN, l_rwl = NaN) %>% 
  mutate(t_rwl = paste0(itrdb_data_dir, filename, ".rwl")) %>% 
  select(collection_id, t_rwl, e_rwl, l_rwl, datasource)

site_data[site_data$collection_id %>% duplicated(),]

itrdb722_sites <- site_data %>% pull(collection_id)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine datasets --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("No duplicates across two data sources",
          expect_equal(itrdb722_sites[which((itrdb722_sites %in% ritrdb_sites))] %>% length(), 0))
site_files <- site_data %>% 
  rbind(ritrdb_df)

itrdb_df <- itrdb_df %>% 
  left_join(site_files, by = "collection_id")  ## Note - dropping some sites since they don't have series spanning into study period 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Process total rwl files --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# site_to_filename <- function(site, suffix){
#   site <- str_to_lower(site)
#   site <- paste0(site, suffix, '.rwl')
# }

open_rwl <- function(file){
  caught_error <- NA
  caught_warning <- NA
  rwl <- NA
  tryCatch(
    expr = {
      rwl <- read.tucson(file)
    },
    error = function(e){
      tryCatch(
        expr = {
          path_nohead <- remove_head(file)
          rwl <<- read.tucson(paste0(path_nohead))
          file.remove(path_nohead)
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
  path <- paste0(file)
  temp_path <- paste0("temp_nohead.rwl")
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

rwl_data <- itrdb_df %>% select(collection_id, file = t_rwl) %>% drop_na()
# rwl_files <- sapply(total_sites, site_to_filename, suffix = '')
# rwl_data <- tibble::enframe(rwl_files, name = "collection_id", value = "file")
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

erwl_data <- itrdb_df %>% select(collection_id, file = e_rwl) %>% filter(file != "NaN")
# erwl_files <- sapply(sum_sites, site_to_filename, suffix = 'e')
# erwl_data <- tibble::enframe(erwl_files, name = "collection_id", value = "file")
erwl_data$rwl <- map(erwl_data$file, open_rwl) 
erwl_data <- erwl_data %>% unnest(rwl)
erwl_data <- erwl_data %>% 
  mutate(ids = map(erwl_data$rwl, read_ids))
erwl_data <- erwl_data %>% unnest(ids)
erwl_results <- separate_errors(erwl_data)
e_clean_data <- erwl_results$clean_data
e_clean_data <- e_clean_data %>% 
  rename(erwl = rwl)

lrwl_data <- itrdb_df %>% select(collection_id, file = l_rwl) %>% filter(file != "NaN")
# lrwl_files <- sapply(sum_sites, site_to_filename, suffix = 'l')
# lrwl_data <- tibble::enframe(lrwl_files, name = "collection_id", value = "file")
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
site_smry <- itrdb_df %>% 
  select(collection_id, latitude, longitude, elevation_itrdb, sp_id = species_id, datasource) %>% 
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
         valid_dates = replace_na(valid_dates, FALSE),
         valid_rwl = collection_id %in% c(sum_sites, total_sites, sum_sites_ritrdb %>% pull(collection_id), total_sites_ritrdb %>% pull(collection_id)),
         valid_parse = if_else(collection_id == "MD007", TRUE, if_else(valid_rwl==FALSE, FALSE, !(collection_id %in% c(l_parse_errors, e_parse_errors, total_parse_errors))))) ## Note MD007 has id parsing error, but produces ok rwi data

site_smry %>% 
  group_by(datasource, valid_rwl, valid_parse, valid_dates) %>% 
  tally() %>% 
  print()


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
itrdb_df$collection_id %>% unique() %>% length()

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
export_data %>% 
  select(collection_id, tree) %>% 
  distinct() %>%
  pull(collection_id) %>% 
  length()

# number of sites in final dataset
export_data %>% 
  select(collection_id) %>% 
  distinct() %>%
  pull(collection_id) %>% 
  length()
  
site_smry %>% 
  filter(valid_dates == TRUE, valid_rwl == TRUE, valid_parse == TRUE) %>% 
  pull(collection_id) %>% 
  length()

