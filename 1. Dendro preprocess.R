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

future::plan(multiprocess)

# Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
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

# sites <- sites[1:20,]

# ### DATA QUALITY NOTES:
# ### Several sites incorrectly parsed. Part of early year year variable has been added to tree_id
# sites <- sites %>%
#   filter(site_id %in% c("101 1", "kr", "rm0std", "gm0", "l18", "nl", "yks"))
# 
# s_id <- "101 1"
# tree_db %>% filter(site_id == s_id)
# sp_id <- "abal"
# obs <- pull_obs(s_id, sp_id)
# 
# 
# #### Deal with incorrect parsing of year / tree_id (101 1)
# tree_err <- function(x) {}
#   (" 190" %in% x)
# 
# str_to_shift = " 190"
# obs$errors <- lapply(obs$tree_id, function(x) str_detect(x, str_to_shift))
# obs <- obs %>% 
#   mutate(tree_id = if_else(errors==TRUE, str_replace(tree_id, " 190", ""), tree_id),
#          tree_id = str_trim(tree_id, "both"),
#          year = if_else(errors==TRUE, paste0("190", year), year))
# ## Check if tree_ids don't match tree_db
# db_trees <- tree_db %>%
#   filter(site_id == s_id, species_id == sp_id) %>%
#   distinct() %>%
#   pull(tree_id)
# n_obs <- dim(obs)[1]
# obs_trees <- obs %>%
#   distinct(tree_id) %>%
#   pull(tree_id)
# all(unlist(lapply(obs_trees, function(x) x %in% db_trees)))
# 
# 
# pull_obs <- function(s_id, sp_id){
#   # Identify trees in site / species combination
#   tree_ids = tree_db %>%
#     filter(species_id == sp_id,
#            site_id == s_id) %>%
#     select('tree_id')
#   
#   # Pull observations of identified trees
#   obs = obs_db %>%
#     filter(site_id == s_id,
#            tree_id %in% local(tree_ids$tree_id)) %>%
#     arrange(tree_id, desc(year)) %>%
#     collect()
#   
#   obs <- obs %>%  
#     arrange(year) %>%
#     mutate(year = as.character(year)) %>%
#     select(tree_id, year, ring_width) 
# }
# 
# check_obs <- function(obs, s_id, sp_id){
#   ## Check duplicates
#   any_duplicates <- obs %>%
#     select(tree_id, year) %>%
#     duplicated() %>%
#     any()
#   if (any_duplicates) {
#     return(NaN, "Error: Multiple obs for single tree_id / year combination")
#   }
#   
#   ## Check for unusual ring widths
#   if (any(obs$ring_width >100)) {
#     return(NaN, "Error: Ring widths unusually large")
#   }
# 
#   if (any(obs$ring_width <- 0.001)) {
#     return(NaN, "Error: Ring widths unusually large")
#   }  
# }


pull_rwl <- function(s_id, sp_id){
  print(paste0(s_id, " - ", sp_id))    
  # Identify trees in site / species combination
  tree_ids = tree_db %>%
    filter(species_id == sp_id,
           site_id == s_id) %>%
    select('tree_id')
  
  # Pull observations of identified trees
  obs = obs_db %>%
    filter(site_id == s_id,
           tree_id %in% local(tree_ids$tree_id)) %>%
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
    return(NaN)
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
    return(NaN)
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
    return(NaN)
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
  rwi_dat <- rwl_dat %>%
    detrend(method = "Spline", make.plot = FALSE, verbose = FALSE) # uses a spline that is 0.67 the series length
  return(rwi_dat)
}

create_crn <- function(rwi_dat){
  crn_dat <- rwi_dat %>% 
    chron(prefix = "CRN", prewhiten = TRUE) %>%
    # select(CRNstd, CRNres) %>%
    rownames_to_column("year") %>%
    as_tibble()
  return(crn_dat)
}

process_dendro <- function(s_id, sp_id) {
  pb$tick()$print()
  rwl_dat <- pull_rwl(s_id, sp_id)
  if (rwl_dat %>% is.na()) {
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

pb <- progress_estimated(dim(sites)[1])
sites$crn <- map2(sites$site_id, sites$species_id, process_dendro) 
invalid_sites <- sites %>%
  filter(crn%>% is.na()) %>%
  as_tibble()

valid_sites <- sites %>% 
  filter(!crn %>% is.na()) %>%
  as_tibble()

valid_data <- valid_sites %>%
  unnest() %>%
  drop_na()

write.csv(valid_data, paste0(wdir, "clean_crn.csv"))

# s_id <- sites %>% pull(site_id) %>% nth(i)
# sp_id <- sites %>% pull(species_id) %>% nth(i)
# crn_dat <- process_dendro(s_id, sp_id)
# 
# sites <- sites %>%
#   filter(!site_id %in% c("101 1", "kr", "rm0std", "gm0", "l18", "nl", "cibola", "rah", "yks", 
#                          "obb", "q 3qr0"))
# test_sites$crn <- future_map2(test_sites$site_id, test_sites$species_id, 
#                               process_dendro)

# 
#   # check for sufficient trees
#   n_trees <- tree_ids %>% n_distinct
#   if (n_trees<5){
#     print(paste0("Too few trees for site-species combination: ", s_id, sp_id))
#     next
#   }
#   
#   # check for 
#   
# }
#   
# 
# }  
#   ,
#               crn = map(data, possibly(create_crn, otherwise = NaN))) ### NOTE: Several sites can't generate chronologies because of lack of observations, or interrupted time series
#   
#   # crn_dat %>% plot(add.spline=TRUE, nyrs=20)
#   if (dim(obs_crn)[1]==0) {
#     print(paste0("No valid sites for: ", sp_id))
#     next
#   }
# 
# 
# # Define regression model for first stage
# fit_mod <- function(d) {
#   felm(ring_width ~ cwd.an + pet.an + age+I(age^2)+I(age^3)+I(age^4)|tree_id|0|tree_id+year, data = d )
# }
# 
# getcov <- function(m){
#   return(m$vcv[2,1])
# }
# 
# 
# fit_mod <- function(d) {
#   felm(crn ~ cwd.an + pet.an |0|0|year, data = d )
# }
# 
# 
# # Load site climate data
# cwd_df <- read.csv(cwd_csv, sep=',')
# cwd_df <- cwd_df %>% 
#   mutate("site_id" = as.character(site)) %>%
#   select(-site)
# 
# clim_df = cwd_df %>%
#   group_by(site_id, year) %>%
#   summarise(aet.an = sum(aet),
#             cwd.an = sum(cwd),
#             pet.an = sum(petm),
#             temp.an = mean(tmean),
#             ppt.an = sum(ppt))
# hist_clim_df <- clim_df %>%
#   group_by(site_id) %>%
#   filter(year<1980) %>%
#   summarise(aet.ave = mean(aet.an),
#             cwd.ave = mean(cwd.an),
#             pet.ave = mean(pet.an),
#             temp.ave = mean(temp.an),
#             ppt.ave = mean(ppt.an))
# 
# 
#   
#   crn_list <- list()
#   n_dropped = 0
#   for (i in 1:dim(obs_crn)[1]){
#     sid <- obs_crn$site_id[[i]]
#     crn_df <- obs_crn %>% 
#       pull(crn) %>% 
#       nth(i)
#     if (!crn_df %>% is.list()){
#       n_dropped <- n_dropped + 1
#       print(paste0(sid, ": site dropped due to unprocessed chronology"))
#     } else {
#       crn_df <- crn_df %>%
#         select(crn = CRNstd) %>%
#         rownames_to_column("year") %>%
#         mutate(year = as.integer(year),
#                "site_id" = sid)
#       crn_df <- crn_df %>%
#         inner_join(clim_df, by = c("site_id", "year"))
#       crn_list[[i]] <- crn_df
#     }
#   }
#   print(paste0("Total sites dropped due to unprocessed chronology: ", n_dropped))
#   crn_df <- bind_rows(crn_list)
# 
#   if (dim(crn_df)[1]==0) {
#     print(paste0("No valid sites for: ", sp_id))
#     next
#   }
#   
#   site_lm <- crn_df %>%
#     group_by(site_id) %>%
#     nest() %>%
#     mutate(mod = map(data, possibly(fit_mod, otherwise = NaN))) %>%
#     filter(!is.na(mod)) %>% 
#     mutate(cwd.pet.cov = map(mod, possibly(getcov, otherwise = NaN)),
#            mod = map(mod, tidy)) %>%
#     unnest(mod,cwd.pet.cov) %>%
#     filter(term %in% c('cwd.an', 'pet.an'))
# 
#   if (dim(site_lm)[1]==0) {
#     print(paste0("No valid sites for: ", sp_id))
#     next
#   }
#     
#   siteCoef <- site_lm %>%
#     pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value")) %>%
#     select(-data)
#   
#   # Attach climate, nobs and ntrees for each site
#   site_df <- obs %>%
#     select(c("site_id")) %>%
#     distinct() %>%
#     left_join(hist_clim_df, by = c("site_id")) %>%
#     left_join(siteCoef, by = c("site_id")) 
#   # %>%
#     # left_join(nobs, by = c("site_id", "species_id")) %>%
#     # left_join(ntrees, by = c("site_id", "species_id"))
# 
#   site_df$species_id <- sp_id
#   
#   # Calculate standardized historic climate relative to species mean / std
#   site_df <- site_df %>%
#     mutate(cwd.spstd = scale(cwd.ave)[,1],
#            aet.spstd = scale(aet.ave)[,1],
#            pet.spstd = scale(pet.ave)[,1],
#            ppt.spstd = scale(ppt.ave)[,1],
#            temp.spstd = scale(temp.ave)[,1])
# 
#   df_list[[i]] <- site_df
# }
# 
# full_df = bind_rows(df_list)
# write.csv(full_df, paste0(wdir, "first_stage.csv"))
# 
# full_df
# mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, data=full_df)
# mod %>% summary()
