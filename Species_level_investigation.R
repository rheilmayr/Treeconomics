library(tidyverse)
library(broom)
library(feols)

### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv")) 
trim_df <- flm_df %>% filter(outlier == 0)

# 2. Species niche
niche_df <- read_csv(paste0(wdir, "out//climate//clim_niche.csv")) %>%
  select(-...1, species_id = sp_code)



sp_reg <- function(sp_data){
  sp_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = cwd_errorweights, data = sp_data)
  sp_mod <- tidy(sp_mod)
  cwd_result <- sp_mod %>% filter(term == "cwd.spstd")
  return(cwd_result)
}

sp_nest <- trim_df %>% 
  group_by(species_id) %>% 
  mutate(species_n = n()) %>%
  filter(species_n > 50) %>% 
  nest()

sp_regs <- sp_nest %>% 
  mutate(mod = map(data, .f = sp_reg)) %>% 
  unnest(mod) %>% 
  select(-data, term)

sp_regs <- sp_regs %>% 
  left_join(niche_df, by = "species_id")

sp_regs %>% 
  filter(p.value <= 0.05)
sp_regs %>% 
  filter(p.value > 0.05)


mod <- lm(estimate ~ pet_mean + cwd_mean, data = sp_regs)
summary(mod)

mod <- lm(estimate ~ pet_sd + cwd_sd, data = sp_regs)
summary(mod)
