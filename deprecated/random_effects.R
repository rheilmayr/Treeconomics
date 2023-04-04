#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/25/20
# Purpose: Replace two-stage approach with random effects model? 
#
# Input files:
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(MASS)
library(tidyverse)
library(lfe)
library(broom)
library(purrr)
library(patchwork)
library(ggpubr)
library(ggiraphExtra)
library(patchwork)
library(hexbin)
library(dbplyr)
library(RSQLite)
library(modi)
library(margins)
library(lme4)

select <- dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "rwi_data\\")
dendro_sites <- read.csv(paste0(dendro_dir, "2_valid_sites.csv")) %>% 
  select(-X) %>% 
  mutate(file_name = paste0('sid-', site_id, '_spid-', species_id, '.csv'))

# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site))

cwd_sites <- cwd_df %>% 
  select(site) %>% 
  distinct()
cwd_dendro_sites <- dendro_sites %>% 
  inner_join(cwd_sites, by = c("site_id" = "site"))


# Connect to database
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tree_db = as.data.frame(tbl(conn, 'trees'))
sp_site_index <- tree_db %>% 
  select(species_id, site_id) %>% 
  distinct()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site and species historic climate ------------------------------
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

# Calculate species niche based on ITRDB sites
sp_clim_df <- clim_df %>% 
  left_join(sp_site_index, by = c("site_id"))

sp_clim_df <- sp_clim_df %>% 
  group_by(species_id) %>% 
  filter(year<1980) %>% 
  summarise(pet.sp.ave = mean(pet.an),
            cwd.sp.ave = mean(cwd.an),
            pet.sp.sd = sd(pet.an),
            cwd.sp.sd = sd(cwd.an))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create full tree by year dataset ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_obs <- function(s_id, sp_id){
  dendro_file <- paste0('sid-', s_id, '_spid-', sp_id, '.csv')
  dendro_dat <- read.csv(paste0(dendro_dir, dendro_file)) %>% 
    mutate(site_id = as.character(s_id)) %>% 
    select(-X, -site_id) %>% 
    drop_na()
  return(dendro_dat)
}

select_sites = cwd_dendro_sites[1:500,]

tree_df <- select_sites %>% 
  group_by(site_id, species_id) %>% 
  nest() %>% 
  mutate(data = map2(site_id, species_id, load_obs)) %>% 
  unnest(data)


tree_df <- tree_df %>% 
  left_join(clim_df, by = c("site_id", "year")) %>% 
  drop_na()

# Merge historic climate
tree_df <- tree_df %>% 
  left_join(hist_clim_df, by = c("site_id")) %>% 
  left_join(sp_clim_df, by = c("species_id"))


tree_df <- tree_df %>% 
  mutate(pet.spstd = (pet.ave - pet.sp.ave) / pet.sp.sd,
         cwd.spstd = (cwd.ave - cwd.sp.ave) / cwd.sp.sd)

tree_df <- tree_df %>% 
  mutate(ln_rwi = log(rwi)) %>% 
  drop_na()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run hierarchical random effects model ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod <- lmer(rwi ~ cwd.an + cwd.an:cwd.spstd + (1 | site_id) + (1 | site_id:tree_id), data = tree_df)
summary(mod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run bayesian hierarchical models ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bmod <- brm(rwi ~ cwd.an + pet.an + (1 | site_id) + (1 | site_id:tree_id), data = tree_df)
summary(bmod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run tree-level fixed effects models ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod <- felm(rwi ~ cwd.an:cwd.spstd + cwd.an | tree_id, data = tree_df)
summary(mod)



