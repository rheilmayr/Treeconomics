#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Run plot-level regressions of RWI sensitivity to annual weather variability
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using species_niche.R
# - tree_ring_data_V2.db: Compiled database of ITRDB observations 
# - essentialcwd_data.csv: File detailing plot-level weather history
# - dendro_dir: Directory containing processed RWI data from "1. Dendro preprocess.R"
#
# ToDo:
# - fix joins to prevent duplicate species_id
# - think through how to deal with CWD outliers
# - track down lost observations - currently dropping a lot due to NAN or failed RWI generation
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
library(fixest)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'out\\climate\\essentialcwd_data.csv')
cwd_df <- read_csv(cwd_csv)
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site))

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out\\dendro\\")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)
# old_dendro_df <- read_csv(paste0(dendro_dir, "rwi_long_old.csv"))
# old_dendro_df <- old_dendro_df %>% 
#   filter(year>1900) %>% 
#   select(-core_id)
# dendro_sites <- read.csv(paste0(dendro_dir, "2_valid_sites.csv")) %>% 
#   select(-X) %>% 
#   mutate(file_name = paste0('sid-', site_id, '_spid-', species_id, '.csv'))



# cwd_sites <- cwd_df %>% 
#   select(site) %>% 
#   distinct()
# cwd_dendro_sites <- dendro_sites %>% 
#   inner_join(cwd_sites, by = c("site_id" = "site"))

# 4. Site information
site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)

dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')

# 5. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))
psme_list <- site_smry %>% 
  filter(species_id == "psme") %>% 
  pull(collection_id)


# 1. Historic species niche data
niche_csv <- paste0(wdir, 'out/climate/clim_niche.csv')
niche_df <- read_csv(niche_csv) %>% 
  select(-X1)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site historic climate ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  rename(collection_id = site_id) %>% 
  group_by(collection_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum((aet+cwd)))
# temp.an = mean(tmean),
# ppt.an = sum(ppt))

dendro_df <- dendro_df %>% 
  left_join(clim_df, by = c("collection_id", "year"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore variance by site ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level historic climate
hist_clim_df <- clim_df %>%
  group_by(collection_id) %>%
  filter(year<1980) %>%
  summarise(aet.ave = mean(aet.an),
            aet.std = sd(aet.an),
            cwd.ave = mean(cwd.an),
            cwd.std = sd(cwd.an),
            pet.ave = mean(pet.an),
            pet.std = sd(pet.an),
            cwd.min = min(cwd.an)) %>% 
  left_join(site_smry, by = "collection_id") %>% 
  select(-genus, -gymno_angio, -species_id)

dendro_df <- dendro_df %>% 
  left_join(hist_clim_df, by = "collection_id")


niche_df <- niche_df %>% 
  select(sp_code, sp_pet_mean = pet_mean, sp_pet_sd = pet_sd, sp_cwd_mean = cwd_mean, sp_cwd_sd = cwd_sd)

# Merge chronology and species climate niche data
dendro_df <- dendro_df %>%
  inner_join(niche_df, by = c("species_id" = "sp_code")) #note: currently losing ~7k trees (10%) because we don't have range maps

# Calculate species-niche standardized climate
dendro_df <- dendro_df %>%
  mutate(pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)

dendro_df <- dendro_df %>% 
  mutate(pet.an.spstd = (pet.an - sp_pet_mean) / sp_pet_sd,
         cwd.an.spstd = (cwd.an - sp_cwd_mean) / sp_cwd_sd)


# Calculate site-standardized climate
dendro_df <- dendro_df %>% 
  mutate(pet.an.site.std = (pet.an - pet.ave) / pet.std,
         cwd.an.site.std = (cwd.an - cwd.ave) / cwd.std)





test_lm <- lm(rwi ~ cwd.an.spstd * cwd.spstd + cwd.an.spstd * pet.spstd + pet.an.spstd * cwd.spstd + pet.an.spstd * pet.spstd + I(cwd.an.spstd ** 2), dendro_df)
summary(test_lm)

test_lm <- lm(rwi ~ cwd.an.spstd * cwd.spstd + I(cwd.an.spstd ** 2), dendro_df)
summary(test_lm)


test_lm <- feols(rwi ~ cwd.an.spstd * cwd.spstd + cwd.an.spstd * pet.spstd + pet.an.spstd * cwd.spstd + pet.an.spstd * pet.spstd | collection_id + species_id, dendro_df)
print(test_lm, se = "twoway")

rwl_lm <- feols(rwl ~ cwd.an.spstd * cwd.spstd + cwd.an.spstd * pet.spstd + pet.an.spstd * cwd.spstd + pet.an.spstd * pet.spstd | collection_id + species_id, dendro_df)
print(rwl_lm, se = "twoway")




test_lm <- feols(rwi ~ cwd.an.site.std * cwd.spstd + cwd.an.site.std * pet.spstd + pet.an.site.std * cwd.spstd + pet.an.site.std * pet.spstd | collection_id + species_id, dendro_df)
print(test_lm, se = "twoway")

test_lm <- feols(rwi ~ cwd.an.site.std * cwd.spstd | collection_id + species_id, dendro_df)
print(test_lm, se = "twoway")



# TODO: Drop outliers

test_lm <- feols(rwi ~ cwd.an.spstd * cwd.spstd | species_id + collection_id, dendro_df)
print(test_lm, se = "twoway")
test_lm <- feols(rwi ~ cwd.an.spstd * cwd.spstd + cwd.an.spstd * pet.spstd | species_id + collection_id, dendro_df)
print(test_lm, se = "twoway")

library(modelsummary)
modelsummary(test_lm)




predictions <- prediction(test_lm, at = list(cwd.spstd = seq(-3, 3, .5)), calculate_se = T) %>% 
  summary() %>% 
  rename(cwd.spstd = "at(cwd.spstd)")


margins_plot <- ggplot(predictions, aes(x = cwd.spstd)) + 
  geom_line(aes(y = Prediction), color = "#404788FF", size = 2) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.8, fill = "#404788FF") +
  geom_line(aes(y = upper), linetype = 3) +
  geom_line(aes(y = lower), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") +
  # xlab(element_blank()) + 
  ylab("Predicted sensitivity\nto CWD") + 
  theme_bw(base_size = 22) +
  scale_y_continuous(labels = scales::scientific)



test_lm <- feols(rwi ~ cwd.an * cwd.spstd + cwd.an * pet.spstd | species_id + collection_id, dendro_df)
print(test_lm, se = "twoway")
test_lm <- feols(rwi ~ cwd.an * cwd.spstd + cwd.an * pet.spstd + pet.an * cwd.spstd + pet.an * pet.spstd | species_id + collection_id, dendro_df)
print(test_lm, se = "twoway")



test_lm <- feols(rwi ~ cwd.an * cwd.ave + cwd.an * pet.ave + pet.an * cwd.ave + pet.an * pet.ave | family + collection_id, dendro_df)
print(test_lm, se = "twoway")

test_lm <- feols(rwi ~ cwd.an * cwd.ave + cwd.an * pet.ave | family + collection_id, dendro_df)
print(test_lm, se = "twoway")





test_lm <- feols(rwl ~ cwd.an.spstd * cwd.spstd + cwd.an.spstd * pet.spstd + pet.an.spstd * cwd.spstd + pet.an.spstd * pet.spstd | species_id, dendro_df)
print(test_lm)

test_lm <- lm(rwl ~ cwd.an.spstd * cwd.spstd + cwd.an.spstd * pet.spstd + pet.an.spstd * cwd.spstd + pet.an.spstd * pet.spstd, dendro_df)
summary(test_lm)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Early life drought models  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
young_trees <- dendro_df %>% 
  group_by(collection_id, tree) %>% 
  summarise(min_year = min(year)) %>% 
  mutate(young = min_year>1901)

# flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv')) %>%
#   select(collection_id,tree,young)
young_trees <- young_trees %>%
  select(-min_year) %>% 
  right_join(dendro_df, by = c("collection_id", "tree")) %>%
  filter(young==TRUE)%>%
  select(collection_id,tree,year,cwd.an)%>%
  group_by(collection_id,tree)%>%
  filter(year<(min(year)+10))%>%
  summarize(earlylifecwd=mean(cwd.an))

youngtrees_df <- dendro_df %>%
#  select(collection_id, tree, year, species_id, ln_rwi, cwd.an, pet.an) %>% 
  inner_join(young_trees, by = c("collection_id", "tree")) %>% 
  mutate(cwd.early.spstd = (earlylifecwd - sp_cwd_mean) / sp_cwd_sd)


earlylifemod = feols(rwi ~ cwd.an * cwd.early.spstd + cwd.an * cwd.spstd | collection_id + species_id , data = youngtrees_df)
print(earlylifemod)

earlylifelist=list()
#interactions effect
for(i in 1:length(relgenus)){
  gendat=youngtrees_df%>%
    filter(genus == relgenus[i])%>%
    filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
  gendat$id=interaction(gendat$collection_id,gendat$tree)
  intmod=felm(ln_rwi~cwd.an*I(earlylifecwd>300)+pet.an-I(earlylifecwd>300)|id|0|collection_id,data=gendat)
  earlylifelist[[i]]=intmod
  print(i)
}

youngtrees_df$id=interaction(youngtrees_df$collection_id,youngtrees_df$tree)
youngtrees_df=youngtrees_df%>%
  filter(earlylifecwd<quantile(earlylifecwd,p=0.99,na.rm=T))%>%
  filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
earlylifemod=felm(ln_rwi~cwd.an*I(earlylifecwd>300)+pet.an-I(earlylifecwd>300)|id|0|collection_id,data=youngtrees_df)
earlylifemod=felm(ln_rwi~cwd.an*earlylifecwd + pet.an - earlylifecwd|id|0|collection_id,data=youngtrees_df)


youngtrees_df <- youngtrees_df %>% 
  left_join(hist_clim_df, by = "collection_id")

subset <- youngtrees_df %>% 
  filter(sp_id == "psme")
earlylifemod=felm(ln_rwi~cwd.an*earlylifecwd + cwd.an*cwd.ave + pet.an |id|0|collection_id,data=subset)
summary(earlylifemod)


earlylifemod=felm(ln_rwi~cwd.an*cwd.ave + pet.an |id|0|collection_id,data=subset)
summary(earlylifemod)
