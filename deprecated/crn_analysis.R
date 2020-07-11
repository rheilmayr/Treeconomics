#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run full analysis using chronologies rather than tree-level data
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import libraries ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(lfe)
library(broom)
library(patchwork)
library(hexbin)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Load processed chronologies
crn_csv <- paste0(wdir, 'clean_crn.csv')
crn_df <- read_csv(crn_csv) %>% 
  select(-X1)

# 2. Load weather data
cwd_csv <- paste0(wdir, 'essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) 

# 3. Load species niche data
niche_csv <- paste0(wdir, 'clim_niche.csv')
niche_df <- read_csv(niche_csv) %>% 
  select(-X1) %>% 
  drop_na()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site weather ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual weather
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

# Merge to crn
crn_df <- crn_df %>% 
  inner_join(clim_df, by = c("site_id", "year"))

crn_df <- crn_df %>% 
  drop_na()


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
  inner_join(hist_clim_df, by = c("site_id")) %>% 
  inner_join(niche_smry, by = c("species_id" = "sp_code")) #note: currently losing a bunch of observations because we don't yet have range maps


# Calculate species-niche standardized climate
crn_df <- crn_df %>%
  mutate(aet.spstd = (aet.ave - sp_aet_mean) / sp_aet_sd,
         pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd,
         cwd.spstd.an = (cwd.an - sp_cwd_mean) / sp_cwd_sd,
         pet.spstd.an = (pet.an - sp_pet_mean) / sp_pet_sd)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Big regression using site's deviation from species mean niche ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crn_df <- crn_df %>% 
  mutate(sp_s_id = paste0(species_id, site_id))
mod <- felm(CRNstd ~ cwd.spstd.an + pet.spstd.an | sp_s_id | 0 | sp_s_id + year, data = crn_df)
summary(mod)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Site-level regressions ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crn_df <- crn_df %>% 
  mutate(ln_crn = log(CRNstd))

lm_sites <- crn_df %>% 
  group_by(site_id, species_id) %>% 
  nest()
lm_sites <- lm_sites %>% 
  mutate(mod = map(data, ~lm(CRNres ~ cwd.an + pet.an, data = .)))

lm_sites <- lm_sites %>% 
  mutate(mod = map(mod, tidy)) %>% 
  unnest(c(mod)) %>% 
  filter(term %in% c('cwd.an', 'pet.an'))

site_coef <- lm_sites %>%
  pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))



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
site_coef <- site_coef %>% 
  inner_join(hist_clim_df, by = c("site_id"))
site_coef <- site_coef  %>%
  inner_join(niche_smry, by = c("species_id" = "sp_code")) #note: currently losing a bunch of observations because we don't yet have range maps


# Calculate species-niche standardized climate
site_coef <- site_coef %>%
  mutate(aet.spstd = (aet.ave - sp_aet_mean) / sp_aet_sd,
         pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd,
         cwd.spstd.an = (cwd.an - sp_cwd_mean) / sp_cwd_sd,
         pet.spstd.an = (pet.an - sp_pet_mean) / sp_pet_sd)


trim_df <- site_coef %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T)) %>%
  ungroup()
trim_df <- trim_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh)



mod <- lm(estimate_cwd.an ~ cwd.spstd, data = trim_df)
mod <- felm(estimate_cwd.an ~ cwd.spstd|species_id, data = trim_df)
summary(mod)



crn_df <- crn_df %>% 
  inner_join(hist_clim_df, by = c("site_id")) %>% 
  inner_join(niche_smry, by = c("species_id" = "sp_code")) #note: currently losing a bunch of observations because we don't yet have range maps

crn_df <- crn_df %>%
  mutate(aet.spstd = (aet.ave - sp_aet_mean) / sp_aet_sd,
         pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)

crn_df <- crn_df %>% 
  drop_na() %>% 
  mutate(ln_crn = log(CRNstd))

crn_df <- crn_df %>% 
  mutate(s_sp_id = paste0(site_id, species_id))


mod <- felm(ln_crn ~ cwd.an + lag(cwd.an, 1) + lag(cwd.an, 2) + lead(cwd.an, 1) + lead(cwd.an, 2) |s_sp_id|0|s_sp_id + year, data = crn_df)
summary(mod)


mod <- felm(ln_crn ~ cwd.an + cwd.an:cwd.spstd | s_sp_id | 0 | s_sp_id + year, data = crn_df)

mod <- felm(ln_crn ~ cwd.an + cwd.an:pet.an + pet.an |s_sp_id|0|s_sp_id + year, data = crn_df)
summary(mod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot sensitivities -------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Boxplots of effects by cwd / pet bins
plot_dat <- site_coef %>%
  group_by("species_id") %>%
  mutate(ppt.q = as.factor(ntile(ppt.ave, 5)),
         temp.q = as.factor(ntile(temp.ave, 5)),
         ppt.high = as.factor(ntile(ppt.ave, 2)-1),
         temp.high = as.factor(ntile(temp.ave, 2)-1),
         climate = paste0('ppt', ppt.high, '_temp', temp.high)) %>%
  ungroup()

p1 <- plot_dat %>% ggplot(aes(x=ppt.q, y=estimate_ppt.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p2 <- plot_dat %>% ggplot(aes(x=ppt.q, y=estimate_temp.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p3 <- plot_dat %>% ggplot(aes(x=temp.q, y=estimate_ppt.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p4 <- plot_dat %>% ggplot(aes(x=temp.q, y=estimate_temp.an)) +
  geom_boxplot() +
  ylim(-0.005, 0.005)

(p1 | p2) / (p3 | p4)










plot_dat <- site_coef %>%
  group_by("species_id") %>%
  mutate(cwd.q = as.factor(ntile(cwd.ave, 10)),
         pet.q = as.factor(ntile(pet.ave, 10)),
         cwd.high = as.factor(ntile(cwd.ave, 2)-1),
         pet.high = as.factor(ntile(pet.ave, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high)) %>%
  ungroup()

p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p2 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p3 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p4 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.005, 0.005)

(p1 | p2) / (p3 | p4)
