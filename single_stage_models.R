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
formula = "rwi ~ cwd.an.spstd:pet.spstd + cwd.an.spstd:I(pet.spstd**2) +
               cwd.an.spstd:cwd.spstd + cwd.an.spstd:I(cwd.spstd**2) +
               pet.an.spstd:pet.spstd + pet.an.spstd:I(pet.spstd**2) +
               pet.an.spstd:cwd.spstd + pet.an.spstd:I(cwd.spstd**2) +
                                    (1 + cwd.an.spstd + pet.an.spstd | collection_id)"

mod_df <- dendro_df
cwd_median <- mod_df %>% select(cwd.spstd) %>% drop_na() %>% pull(cwd.spstd) %>% median()


mod <- lmer(formula, data = mod_df, control = lmerControl(optimizer ="Nelder_Mead"))
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
