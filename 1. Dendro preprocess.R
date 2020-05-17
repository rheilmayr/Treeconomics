#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Process dendrochronologies to output plot-level growth deviations
#
# Input files:
# - tree_ring_data_V2.db: Compiled database of ITRDB observations 
# - cwd_csv: File detailing plot-level weather history
#
# ToDo: 
# - Shift to tree-level output based on Pederson et al 2020 paper?
# - Parallelize code using furrr and future packages?
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
# library(future)
# library(furrr)
# library("tidylog", warn.conflicts = FALSE)
# future::plan(multiprocess)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'
# 'D:/cloud/Google Drive/Treeconomics/Data/'
# #for Fran
# wdir="C:/Users/fmoore/Google Drive/Treeconomics/Data/"
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
cwd_csv = paste0(wdir, 'essentialcwd_data.csv')

# Connect to database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)

# Open tables
tree_db = as.data.frame(tbl(conn, 'trees'))
spp_db = as.data.frame(tbl(conn,'species'))
site_db = as.data.frame(tbl(conn, "sites"))
obs_db = tbl(conn, 'observations_new')

# Create index of site-species combinations
sites <- tree_db %>%
  select(species_id, site_id) %>%
  distinct() %>%
  collect()

# Create unique tree ids that incorporate site id
tree_db <- tree_db %>% 
  mutate(tree_id = paste0("s", site_id, "_t", tree_id))
obs_db <- obs_db %>% 
  mutate(tree_id = paste0("s", site_id, "_t", tree_id))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define functions for dendro pipeline  ----------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pull_rwl <- function(s_id, sp_id){
  # Purpose:
  #   Pulls ring width data for a specific site and species
  # Parameters:
  #   s_id: str
  #     Site id to pull from
  #   sp_id: str
  #     Species id to pull from
  # 
  # Outputs:
  #   rwl_dat: data.frame
  #     Data frame of annual ring widths for site-species combination
  
  print(paste0(s_id, " - ", sp_id))    
  # Identify trees in site / species combination
  tree_ids = tree_db %>%
    filter(species_id == sp_id,
           site_id == s_id) %>%
    select('tree_id')
  
  # Pull observations of identified trees
  obs = obs_db %>%
    filter(tree_id %in% local(tree_ids$tree_id)) %>%
    arrange(tree_id, desc(year)) %>%
    collect()
  
  obs <- obs %>%  
    arrange(year) %>%
    mutate(year = as.character(year)) %>%
    select(tree_id, year, ring_width) 
  
  # Check for invalid duplicates in site-species-year
  any_duplicates <- obs %>%
    select(tree_id, year) %>%
    duplicated() %>%
    any()
  if (any_duplicates) {
    print(paste0("Duplicate tree-year combinations; skipping site-species: ", s_id, " - ", sp_id))
    e = "Duplicate tree-year observations"
    return(NULL)
  }

  obs <- obs %>%
    pivot_wider(names_from = tree_id, values_from = ring_width) %>%
    column_to_rownames("year")
  
  failed <- F
  
  # Check for invalid data that can't be converted to rwl
  tryCatch(
    expr = {
      rwl_dat <- obs %>% as.rwl()
    },
    error = function(e){ 
      message("Returned error on site ", s_id)
      print(e)
      failed <<- T
    }
  )
  if (failed){
    return(NULL)
  }
  
  # Check for invalid data that doesn't allow for rwl report
  tryCatch(
    expr = {
      rwl_report <- rwl.report(rwl_dat)
      internal_na <- rwl.report(rwl_dat)
      internal_na <- internal_na[14][1]$internalNAs
    },
    error = function(e){
      message("Warning on site ", s_id)
      print(e)
      failed <<- T
    }
  )
  if (failed){
    return(NULL)
  }
  
  # For series with internal NAs, drop old observations prior to NA
  if (length(internal_na)>0){
    na_trees <- internal_na %>% names()
    for (t in 1:length(na_trees)) {
      na_tree <- na_trees[t]
      last_na <- max(unlist(internal_na[t]))
      na_years <- (rwl_dat %>% row.names()) <= last_na
      rwl_dat <- rwl_dat %>%
        mutate(!!na_tree := na_if(na_years, na_tree))
    }
    rwl_dat <- rwl_dat %>% as.rwl()
  }
  return(rwl_dat)
}

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

create_crn <- function(rwi_dat){
  # Purpose:
  #   Apply dplR to convert tree-level ring widths indices (RWL) to site-level chronology (crn)
  # Inputs:
  #   rwi_dat: data.table
  #     Generated from detrend_rwl()
  # Outputs:
  #   crn_dat: data.table
  #     Table of site chronologies
  crn_dat <- rwi_dat %>% 
    chron(prefix = "CRN", prewhiten = TRUE) %>%
    # select(CRNstd, CRNres) %>%
    rownames_to_column("year") %>%
    as_tibble()
  return(crn_dat)
}


process_dendro <- function(s_id, sp_id) {
  # Purpose:
  #   Run full dendro workflow
  # Inputs:
  #   s_id: str
  #     Site id
  #   sp_id: str
  #     Species id
  # Outputs:
  #   crn_dat: data.table
  #     Table of site chronologies
  pb$tick()$print()
  rwl_dat <- pull_rwl(s_id, sp_id)
  if (rwl_dat %>% is.null()) {
    return(NaN)
  }
  rwi_dat <- rwl_dat %>% detrend_rwl()
  crn_dat <- rwi_dat %>% create_crn()
  # Diagnostic plots
  # rwl_dat %>% rwl.report()
  # rwl_dat %>% plot.rwl()
  # crn_dat %>% plot()
  # print(paste0("Successfully processed chronology for: ", s_id, sp_id))
  return(crn_dat)  
}



export_rwi <- function(s_id, sp_id) {
  # Purpose:
  #   Run full dendro workflow
  # Inputs:
  #   s_id: str
  #     Site id
  #   sp_id: str
  #     Species id
  # Outputs:
  #   crn_dat: data.table
  #     Table of site chronologies
  pb$tick()$print()
  rwl_dat <- pull_rwl(s_id, sp_id)
  if (rwl_dat %>% is.null()) {
    return(NaN)
  }
  rwi_dat <- rwl_dat %>% detrend_rwl()
  rwi_dat <- rwi_dat %>% 
    rownames_to_column("year") %>% 
    pivot_longer(cols = -year,
                 names_to = "tree_id",
                 values_to = "rwi")

  write.csv(rwi_dat, paste0(wdir, 'rwi_data\\sid-', s_id, '_spid-', sp_id, '.csv'))
  # Diagnostic plots
  # rwl_dat %>% rwl.report()
  # rwl_dat %>% plot.rwl()
  # crn_dat %>% plot()
  # print(paste0("Successfully processed chronology for: ", s_id, sp_id))
  # return(crn_dat)  
  return(dim(rwi_dat)[1])
  }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run chronology generation  ----------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# sites <- sites[1:10,]
# site <- sites[2,]
# sp_id <- site %>% pull(species_id)
# s_id <- site %>% pull(site_id)

ptm <- proc.time()
pb <- progress_estimated(dim(sites)[1])
sites$nobs <- map2(sites$site_id, sites$species_id, export_rwi) 
ptm <- (proc.time() - ptm)
ptm <- ptm[3] / 60

invalid_sites <- sites %>%
  filter(nobs %>% is.na()) %>%
  as_tibble() %>% 
  select(species_id, site_id)
n_invalid = dim(invalid_sites)[1]
invalid_sites %>% write.csv(paste0(wdir, 'rwi_data\\1_invalid_sites.csv'))

valid_sites <- sites %>% 
  filter(!nobs %>% is.na()) %>%
  unnest('nobs') %>% 
  as_tibble()
n_valid = dim(valid_sites)[1]
valid_sites %>% write.csv(paste0(wdir, 'rwi_data\\2_valid_sites.csv'))

nobs <- valid_sites %>% 
  pull(nobs) %>% 
  sum()

fileConn <- file(paste0(wdir, 'rwi_data\\3_dendro_summary.txt'))
writeLines(c(paste0("Dendro processing completed with ",
                    as.character(n_valid),
                    " valid site/species combinations and ",
                    as.character(n_invalid),
                    " invalid site/species combinations."),
             paste0("In aggregate, processed dataset includes ",
                    as.character(nobs),
                    " tree-year observations"),
             paste0("Total processing time was ",
                    as.character(ptm),
                    " minutes.")), fileConn)
close(fileConn)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generate and export cleaned chronology  ----------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pb <- progress_estimated(dim(sites)[1])
sites$crn <- map2(sites$site_id, sites$species_id, process_dendro) 
invalid_sites <- sites %>%
  filter(crn %>% is.na()) %>%
  as_tibble()

valid_sites <- sites %>% 
  filter(!crn %>% is.na()) %>%
  as_tibble()

valid_data <- valid_sites %>%
  unnest() %>%
  drop_na()

write.csv(valid_data, paste0(wdir, "clean_crn.csv"))

