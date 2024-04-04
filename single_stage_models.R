#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run alternate specifications of second stage model
#
# Input files:
#
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(MASS)
library(car)
library(tidyverse)
library(broom)
library(purrr)
library(patchwork)
library(patchwork)
library(dbplyr)
library(RSQLite)
library(lmtest)
library(sandwich)
library(prediction)
library(Hmisc)
library(hexbin)
library(ggpubr)
library(ggiraphExtra)
library(modi)
library(margins)
library(tidylog)
library(fixest)
library(biglm)
library(gstat)
library(sf)
library(units)
library(dtplyr)
library(marginaleffects)
source("f_spec_chart_function.R")
library(lme4)
library(multcomp)
library(tictoc)

select <- dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
fs_spl <- read_csv(paste0(wdir, '2_output/first_stage/site_pet_cwd_std.csv'))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz"))

# 3. Site information
dendro_dir <- paste0(wdir, "1_input_processed/dendro/")
site_df <- read_csv(paste0(dendro_dir, 'site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# 4. Species information
sp_info <- read_csv(paste0(wdir, '1_input_processed/species_ranges/species_metadata.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_df <- site_df %>% 
  left_join(sp_info, by = "species_id")


# 5. Dendrochronologies - used for single-stage model comparison
dendro_df <- read_csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            .groups = "drop") %>% 
  as_tibble()


# 6. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "2_output/climate/site_an_clim.gz"))
dendro_df <- dendro_df %>% 
  left_join(an_site_clim, by = c("collection_id", "year"))

dendro_df <- dendro_df %>% 
  left_join(ave_site_clim, by = "collection_id")

dendro_df <- dendro_df %>% 
  filter(year > 1901)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run model --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formula = "rwi ~ cwd.an.spstd + pet.an.spstd | collection_id"

formula = "rwi ~ cwd.an.spstd:pet.spstd + cwd.an.spstd:I(pet.spstd**2) +
               cwd.an.spstd:cwd.spstd + cwd.an.spstd:I(cwd.spstd**2) +
               pet.an.spstd:pet.spstd + pet.an.spstd:I(pet.spstd**2) +
               pet.an.spstd:cwd.spstd + pet.an.spstd:I(cwd.spstd**2) +
                                    (1 + cwd.an.spstd + pet.an.spstd | collection_id)"

mod_df <- dendro_df
cwd_median <- mod_df %>% select(cwd.spstd) %>% drop_na() %>% pull(cwd.spstd) %>% median()

formula = "rwi ~ cwd.an.spstd + pet.an.spstd | collection_id"
mod <- feols(as.formula(formula), data = mod_df)
summary(mod)

formula = "rwi ~ ppt.an.spstd + pet.an.spstd | collection_id"
mod <- feols(as.formula(formula), data = mod_df)
summary(mod)

formula = "rwi ~ cwd.an.spstd.cru + pet.an.spstd.cru | collection_id"
mod <- feols(as.formula(formula), data = mod_df)
summary(mod)

formula = "rwi ~ cwd.an.spstd.tc + pet.an.spstd.tc | collection_id"
mod <- feols(as.formula(formula), data = mod_df)
summary(mod)

formula = "rwi ~ ppt.an.spstd.tc + pet.an.spstd.tc | collection_id"
mod <- feols(as.formula(formula), data = mod_df)
summary(mod)







# mod <- lmer(formula, data = mod_df, control = lmerControl(optimizer ="Nelder_Mead"))
test_str <- paste0("cwd.spstd + (2 * \`I(cwd.spstd^2)\` * ", as.character(cwd_median), ") = 0")
lincom <- glht(mod, linfct = c(test_str))
lincom <- summary(lincom)


mod <- lmer(formula, data = mod_df, control = lmerControl(optimizer ="Nelder_Mead"))
test_str <- paste0("cwd.spstd + (2 * \`I(cwd.spstd^2)\` * ", as.character(cwd_median), ") = 0")
lincom <- glht(mod, linfct = c(test_str))
lincom <- summary(lincom)




mod <- feols(rwi ~ pet.an.spstd + cwd.an.spstd | collection_id, data = dendro_df %>% filter(year > 1958))
mod %>% summary()
mod <- feols(rwi ~ pet.an.spstd + ppt.an.spstd | collection_id, data = dendro_df %>% filter(year > 1958))
mod %>% summary()
mod <- feols(rwi ~ pet.an.spstd.tc + cwd.an.spstd.tc | collection_id, data = dendro_df %>% filter(year > 1958))
mod %>% summary()
mod <- feols(rwi ~ pet.an.spstd.tc + ppt.an.spstd.tc | collection_id, data = dendro_df %>% filter(year > 1958))
mod %>% summary()


formula = as.formula("rwi ~ cwd.an.spstd + pet.an.spstd + 
               cwd.an.spstd:pet.spstd + cwd.an.spstd:I(pet.spstd**2) +
               cwd.an.spstd:cwd.spstd + cwd.an.spstd:I(cwd.spstd**2) +
               pet.an.spstd:pet.spstd + pet.an.spstd:I(pet.spstd**2) +
               pet.an.spstd:cwd.spstd + pet.an.spstd:I(cwd.spstd**2)")
mod <- lm(formula, data = dendro_df %>% filter(year > 1958))
mod %>% summary()


formula = as.formula("rwi ~ cwd.an.spstd.tc + pet.an.spstd.tc + 
               cwd.an.spstd.tc:pet.spstd.tc + cwd.an.spstd.tc:I(pet.spstd.tc**2) +
               cwd.an.spstd.tc:cwd.spstd.tc + cwd.an.spstd.tc:I(cwd.spstd.tc**2) +
               pet.an.spstd.tc:pet.spstd.tc + pet.an.spstd.tc:I(pet.spstd.tc**2) +
               pet.an.spstd.tc:cwd.spstd.tc + pet.an.spstd.tc:I(cwd.spstd.tc**2)")
mod <- lm(formula, data = dendro_df %>% filter(year > 1958))
mod %>% summary()






mod_df <- dendro_df %>% 
  left_join(site_df %>% select(collection_id, species_id), by = "collection_id") %>% 
  # filter(species_id == "pipo") %>%
  # mutate(dry_class = ifelse(ppt.spstd < -0.5, "dry", ifelse(ppt.spstd > 0.5, "wet", "medium"))) %>% 
  mutate(dry_class = cut(ppt.spstd, quantile(ppt.spstd, 0:4/4, na.rm = TRUE))) %>% 
  drop_na()
  
mod <- lm(rwi ~ poly(ppt.an.spstd,2)*dry_class + poly(pet.an.spstd,2)*dry_class, data = mod_df)
summary(mod)

marg_fx_df <- function(mod, mod_df){
  classes = mod_df$dry_class %>% unique()
  inc <- 0.1
  slope_df = tibble()
  for (c in classes) {
    print(c)
    ppt_range <- mod_df %>% filter(dry_class == c) %>% pull(ppt.an.spstd) %>% range()
    min <- ppt_range[1]
    max <- ppt_range[2]
    class_slopes <- slopes(mod, newdata = datagrid(dry_class = c, pet.an.spstd = 0, ppt.an.spstd = seq(min,max,inc))) %>% 
      mutate(dry_class = c)
    slope_df <- rbind(slope_df, class_slopes)
  }
  return(slope_df)
}

slope_df <- marg_fx_df(mod, mod_df)



slope_df %>% 
  filter(term == "ppt.an.spstd") %>%
  ggplot(aes(x = ppt.an.spstd, group = dry_class, color = dry_class)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)




mod_df <- dendro_df %>% 
  left_join(site_df %>% select(collection_id, species_id), by = "collection_id") %>% 
  # filter(species_id == "pipo") %>%
  # mutate(dry_class = ifelse(ppt.spstd < -0.5, "dry", ifelse(ppt.spstd > 0.5, "wet", "medium"))) %>% 
  mutate(dry_class = cut(ppt.spstd.tc, quantile(ppt.spstd, 0:4/4, na.rm = TRUE))) %>% 
  drop_na()

mod <- lm(rwi ~ poly(ppt.an.spstd.tc,2)*dry_class + poly(pet.an.spstd.tc,2)*dry_class, data = mod_df)
summary(mod)

marg_fx_df <- function(mod, mod_df){
  classes = mod_df$dry_class %>% unique()
  inc <- 0.1
  slope_df = tibble()
  for (c in classes) {
    print(c)
    ppt_range <- mod_df %>% filter(dry_class == c) %>% pull(ppt.an.spstd.tc) %>% range()
    min <- ppt_range[1]
    max <- ppt_range[2]
    class_slopes <- slopes(mod, newdata = datagrid(dry_class = c, pet.an.spstd.tc = 0, ppt.an.spstd.tc = seq(min,max,inc))) %>% 
      mutate(dry_class = c)
    slope_df <- rbind(slope_df, class_slopes)
  }
  return(slope_df)
}

slope_df <- marg_fx_df(mod, mod_df)

slope_df %>% 
  filter(term == "ppt.an.spstd.tc") %>%
  ggplot(aes(x = ppt.an.spstd.tc, group = dry_class, color = dry_class)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)







mod_df <- dendro_df %>% 
  left_join(site_df %>% select(collection_id, species_id), by = "collection_id") %>% 
  filter(species_id == "pipo") %>%
  mutate(dry_class = ifelse(cwd.spstd < -0.5, "wet", ifelse(cwd.spstd > 0.5, "dry", "medium"))) %>% 
  drop_na()

mod <- lm(rwi ~ poly(cwd.an.spstd.tc,2)*dry_class + poly(pet.an.spstd,2)*dry_class, data = mod_df)


marg_fx_df <- function(mod, mod_df){
  classes = mod_df$dry_class %>% unique()
  inc <- 0.1
  slope_df = tibble()
  for (c in classes) {
    print(c)
    cwd_range <- mod_df %>% filter(dry_class == c) %>% pull(cwd.an.spstd.tc) %>% range()
    min <- cwd_range[1]
    max <- cwd_range[2]
    class_slopes <- slopes(mod, newdata = datagrid(dry_class = c, pet.an.spstd = 0, cwd.an.spstd.tc = seq(min,max,inc))) %>% 
      mutate(dry_class = c)
    slope_df <- rbind(slope_df, class_slopes)
  }
  return(slope_df)
}

slope_df <- marg_fx_df(mod, mod_df)


slope_df %>% 
  filter(term == "cwd.an.spstd.tc") %>%
  ggplot(aes(x = cwd.an.spstd.tc, group = dry_class, color = dry_class)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)

