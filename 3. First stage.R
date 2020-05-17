#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Combine species niche and dendro data to assess climate responsiveness
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using species_niche.R
# - tree_ring_data_V2.db: Compiled database of ITRDB observations 
# - essentialcwd_data.csv: File detailing plot-level weather history
#
# ToDo:
# - add nobs
# - fix joins to prevent duplicate species_id
# - think through how to deal with CWD outliers
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'
# #for Fran
# wdir="C:/Users/fmoore/Google Drive/Treeconomics/Data/"

# 1. Historic species niche data
niche_csv <- paste0(wdir, 'clim_niche.csv')
niche_df <- read_csv(niche_csv) %>% 
  select(-X1) %>% 
  drop_na()

# 2. ITRDB data
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)
tree_db = as.data.frame(tbl(conn, 'trees'))
spp_db = as.data.frame(tbl(conn,'species'))
site_db = as.data.frame(tbl(conn, "sites"))
# obs_db = tbl(conn, 'observations_new')

# 3. Dendrochronologies
crn_csv = paste0(wdir, 'clean_crn.csv')
crn_df <- read.csv(crn_csv, sep=',') %>%
  select(-X)

# 4. Site-specific weather history
cwd_csv = paste0(wdir, 'essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) %>%
  select(-site)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site historic climate ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  group_by(site_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum(petm),
            temp.an = mean(tmean),
            ppt.an = sum(ppt))

# Calculate site-level historic climate
hist_clim_df <- clim_df %>%
  group_by(site_id) %>%
  filter(year<1980) %>%
  summarise(aet.ave = mean(aet.an),
            cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an),
            temp.ave = mean(temp.an),
            ppt.ave = mean(ppt.an))

# Merge chronology and site-level climate data
crn_df <- crn_df %>%
  inner_join(clim_df, by = c("site_id", "year")) #note that we loose a few sites because we are missing CWD data - probably becaue they are on the coast and more sites because we don't have cwd data before 1900
crn_df <- crn_df %>%
  left_join(hist_clim_df, by = c("site_id"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species historic climate -------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niche_smry <- niche_df %>% 
  group_by(sp_code) %>% 
  summarize(sp_cwd_mean = mean(cwd),
            sp_cwd_sd = sd(cwd),
            sp_aet_mean = mean(aet),
            sp_aet_sd = sd(aet),
            sp_pet_mean = mean(pet),
            sp_pet_sd = sd(pet)) %>% 
  ungroup()

# Merge chronology and species climate niche data
crn_df <- crn_df %>%
  inner_join(niche_smry, by = c("species_id" = "sp_code")) #note: currently losing a bunch of observations because we don't yet have range maps

# Calculate species-niche standardized climate
crn_df <- crn_df %>% 
  mutate(aet.spstd = (aet.ave - sp_aet_mean) / sp_aet_sd,
         pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep site-level data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Log transform crn
crn_df <- crn_df %>%
  mutate(lcrn = log(CRNstd))

# Add number of trees
n_trees <- tree_db %>%
  group_by(site_id, species_id) %>%
  summarise(n_trees = n_distinct(tree_id))

# Create error weight based on sample depth of chronology
crn_df <- crn_df %>%
  left_join(n_trees, by = c("site_id", "species_id")) %>%
  mutate(errorweights = samp.depth / sum(samp.depth))

# Remove outlier cwd values - NEED TO THINK THROUGH THIS
mod_df <- crn_df %>% 
  filter(aet.an>0)

# Run first stage regression
mod <- felm(lcrn ~ cwd.an + pet.an + temp.an + ppt.an |site_id|0|site_id+year, weights = crn_df$errorweights, data = crn_df)


mod <- felm(lcrn ~ pet.an + cwd.an |site_id|0|site_id+year, weights = mod_df$errorweights, data = mod_df)
mod <- felm(lcrn ~ aet.an + cwd.an + temp.an + ppt.an |site_id|0|site_id+year, weights = mod_df$errorweights, data = mod_df)
summary(mod)

mod <- felm(CRNstd ~ cwd.an + cwd.an : cwd.spstd + cwd.spstd |site_id|0|site_id+year, data = crn_df)
summary(mod)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Deprecated code --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# site_lm <- crn_df %>%
#   group_by(species_id, site_id) %>%
#   nest() %>%
#   mutate(mod = map(data, fit_mod),
#          cwd.pet.cov = map(mod, getcov),
#          mod = map(mod, tidy)) %>%
#   unnest(mod,cwd.pet.cov) %>%
#   filter(term %in% c('cwd.an', 'pet.an'))
# 
# siteCoef <- site_lm %>%
#   pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value")) %>%
#   select(-data)
# 
# s_id <- "bb"
# sp_id <- "abal"
# dat <- crn_df %>%
#   filter(species_id == sp_id, site_id == s_id)
# mod <- dat %>% fit_mod()
# summary(mod)
# 
# 
# # Create list of species
# sp_list <- spp_db %>%
#   pull(species_id)
# sp_id <- "psme"
# 
# df_list = list()
# for (i in 1:length(sp_list)) {
#   sp_id = sp_list[i]
# 
#   # ID trees from selected species
#   tree_ids = tree_db %>%
#     filter(species_id == sp_id) %>%
#     select('tree_id') %>%
#     collect()
#   
#   # Pull observations of identified trees
#   obs = obs_db %>%
#     filter(tree_id %in% local(tree_ids$tree_id)) %>%
#     arrange(tree_id, desc(year)) %>%
#     collect()
#   
#   # Add climate data
#   obs = obs %>%
#     inner_join(clim_df, by = c("site_id", "year")) #note that we loose a few sites because we are missing CWD data - probably becaue they are on the coast and more sites because we don't have cwd data before 1900
#   obs = obs %>%
#     left_join(hist_clim_df, by = c("site_id"))
#   obs <- obs %>%
#     inner_join(tree_db, by = c("tree_id", "site_id"))
#   
#   # Add tree age
#   obs$age=obs$year-obs$first_year
#   
#   ###### Set up for first stage #####
#   obs <- obs %>%
#     mutate(tree_id = paste0(site_id, tree_id))
#   
#   # Remove trees with negative ages
#   df_trim <- obs %>%
#     filter(age>0)
#   print(paste0(sp_id, " - dropping observations due to invalid age: ",(df_trim$age <=0) %>% sum()))
#   
#   
#   # Remove trees with a very short record (<5 years)
#   treeobs=df_trim %>%
#     group_by(tree_id) %>%
#     summarize(nyears=n())
#   df_trim=left_join(df_trim,treeobs, by = 'tree_id')
#   print(paste0(sp_id, " - dropping trees due to short record: ",(treeobs$nyears <=5) %>% sum()))
#   df_trim = df_trim %>%
#     filter(nyears>5)
#   
#   # Remove sites with few trees
#   ntrees=df_trim%>%
#     group_by(site_id,species_id)%>%
#     summarize(ntrees=length(unique(tree_id)))
#   df_trim=left_join(df_trim,ntrees, by = c("site_id", "species_id"))
#   print(paste0(sp_id, " - dropping sites due to few trees: ",(ntrees$ntrees <=5) %>% sum()))
#   df_trim=df_trim%>%
#     filter(ntrees>5)
#   
#   
#   # Remove NA and move on to next interation in loop if no observations remain
#   complete_df <- df_trim %>%
#     drop_na(c("cwd.an","aet.an","pet.an", "ring_width"))
#   if (dim(complete_df)[1]==0) {
#     print(paste0("No valid sites for: ", sp_id))
#     next
#   }
#   
#   # Count number of valid tree-year observations per site
#   nobs <- complete_df %>%
#     group_by(site_id, species_id) %>%
#     summarize(nobs = n())
#   
#   # Run first stage
#   site_lm <- complete_df %>% 
#     group_by(site_id, species_id) %>%
#     nest() %>%
#     mutate(mod = map(data, fit_mod),
#            cwd.pet.cov = map(mod, getcov),
#            mod = map(mod, tidy)) %>%
#     unnest(mod,cwd.pet.cov) %>%
#     filter(term %in% c('cwd.an', 'pet.an'))
#   
#   siteCoef <- site_lm %>%
#     pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value")) %>%
#     select(-data)
#   
#   # Attach climate, nobs and ntrees for each site
#   site_df <- obs %>%
#     select(c("site_id", "species_id")) %>%
#     distinct() %>%
#     left_join(hist_clim_df, by = c("site_id")) %>%
#     left_join(siteCoef, by = c("site_id", "species_id")) %>%
#     left_join(nobs, by = c("site_id", "species_id")) %>%
#     left_join(ntrees, by = c("site_id", "species_id"))
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