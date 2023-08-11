#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run regressions to explore impact of historical climate on weather sensitivity
#
# Input files:
#   site_ave_clim.gz:
#   site_pet_cwd_std.csv:
#   site_summary.csv:
# 
# Output files:
#   ss_bootstrap.rds
#
# ToDo:
# - Update / finalize genus analyses
# - Rebuild robustness tests based on final baseline model
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
library(tidylog)
library(extrafont)


set.seed(5597)

select <- dplyr::select



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions from ITRDB
flm_itrdb <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

# 2. Site-level regressions from FIA
flm_fia <- read_csv(paste0(wdir, "out/first_stage/fia_pet_cwd_std.csv"))

# 3. Average site climates
ave_site_clim_df <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))
flm_fia <- flm_fia %>% 
  left_join(ave_site_clim_df, by = c("collection_id"))

# 3. Site information
site_df <- read_csv(paste0(wdir, 'out/dendro/site_summary_fia.csv'))
site_df <- site_df %>% 
  select(collection_id, species_id, latitude, longitude)
site_df <- site_df %>% 
  mutate(species_id = str_to_lower(species_id))
flm_fia <- flm_fia %>% 
  left_join(site_df, by = "collection_id")




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep FIA data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add weighting based on inverse of first stage variance
flm_fia <- flm_fia %>% 
  mutate(cwd_errorweights = 1 / (std.error_cwd.an),
         pet_errorweights = 1 / (std.error_pet.an))

# Identify and trim extreme outliers
cwd_est_bounds = quantile(flm_fia$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds = quantile(flm_fia$estimate_pet.an, c(0.01, 0.99),na.rm=T)

flm_fia <- flm_fia %>%
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]))

test_fia <- flm_fia %>%
  filter(species_id %in% c("pipo", "psme", "pied", "pifl"),
         outlier == 0)
test_fia %>% select(estimate_cwd.an, std.error_cwd.an, ntrees) %>% summary()


test_itrdb <- flm_itrdb %>%
  filter(species_id %in% c("pipo", "psme", "pied", "pifl"),
         outlier == 0)

test_itrdb %>% select(estimate_cwd.an, std.error_cwd.an, ntrees) %>% summary()



cwd_range_fia <- test_fia %>% 
  pull(cwd.spstd) %>% 
  quantile(c(0.025, 0.975))
cwd_range_itrdb <- test_itrdb %>% 
  pull(cwd.spstd) %>% 
  quantile(c(0.025, 0.975))
pet_range_fia <- test_fia %>% 
  pull(pet.spstd) %>% 
  quantile(c(0.025, 0.975))
pet_range_itrdb <- test_itrdb %>% 
  pull(pet.spstd) %>% 
  quantile(c(0.025, 0.975))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spatial autocorrelation of datasets ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_points=st_as_sf(test_itrdb,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
vg.range_itrdb = vg.fit[2,3]


site_points=st_as_sf(test_fia,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
vg.range_fia = vg.fit[2,3]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Estimate models ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_formula = as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
cwd_mod_itrdb <- feols(cwd_formula, data = test_itrdb, weights = test_itrdb$cwd_errorweights,
                       vcov = conley(cutoff = vg.range_itrdb, distance = "spherical"))
cwd_mod_fia <- feols(cwd_formula, data = test_fia, weights = test_fia$cwd_errorweights,
                     vcov = conley(cutoff = vg.range_fia, distance = "spherical"))

summary(cwd_mod_itrdb)
summary(cwd_mod_fia)


pet_formula = as.formula("estimate_pet.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
pet_mod_itrdb <- feols(pet_formula, data = test_itrdb, weights = test_itrdb$pet_errorweights,
                       vcov = conley(cutoff = vg.range_itrdb, distance = "spherical"))
pet_mod_fia <- feols(pet_formula, data = test_fia, weights = test_fia$pet_errorweights,
                     vcov = conley(cutoff = vg.range_fia, distance = "spherical"))

summary(pet_mod_itrdb)
summary(pet_mod_fia)


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Export models ---------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# out_mods <- list("cwd_fia" = cwd_mod_fia,
#                  "cwd_itrdb" = cwd_mod_itrdb,
#                  "pet_fia" = pet_mod_fia,
#                  "pet_itrdb" = pet_mod_itrdb)
# 
# write_rds(out_mods, paste0(wdir, "out/second_stage/fia_itrdb_mods.rds"))
# read_rds(paste0(wdir, "out/second_stage/fia_itrdb_mods.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A: ME of CWD by cwd.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
at_pet <- 0
seq_inc <- 0.1

### Marginal effects of CWD - FIA
pred_fia <- predictions(cwd_mod_fia, newdata = datagrid(pet.spstd = at_pet, cwd.spstd = seq(cwd_range_fia[1],cwd_range_fia[2],seq_inc))) %>% 
  mutate(Dataset = "FIA")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(cwd_mod_itrdb, newdata = datagrid(pet.spstd = at_pet, cwd.spstd = seq(cwd_range_itrdb[1],cwd_range_itrdb[2],seq_inc))) %>% 
  mutate(Dataset = "ITRDB")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
cwd_cwdstd_plot <- pred_all %>% 
  ggplot(aes(x = cwd.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = Dataset, fill = Dataset)) + 
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme_bw() +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD")  +
  scale_color_manual(values = c("ITRDB" = "darkred",
                                "FIA"="darkblue")) +
  scale_fill_manual(values = c("ITRDB" = "darkred",
                               "FIA"="darkblue"))
cwd_cwdstd_plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# B: ME of CWD by pet.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Marginal effects of CWD - ITRDB
pred_fia <- predictions(cwd_mod_fia, newdata = datagrid(pet.spstd = seq(pet_range_fia[1],pet_range_fia[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(Dataset = "FIA")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(cwd_mod_itrdb, newdata = datagrid(pet.spstd = seq(pet_range_itrdb[1],pet_range_itrdb[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(Dataset = "ITRDB")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
cwd_petstd_plot <- pred_all %>% 
  ggplot(aes(x = pet.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = Dataset, fill = Dataset)) + 
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme_bw() +  
  xlab("Historic PET\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD")  +
  scale_color_manual(values = c("ITRDB" = "darkred",
                                "FIA"="darkblue")) +
  scale_fill_manual(values = c("ITRDB" = "darkred",
                                "FIA"="darkblue"))
cwd_petstd_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# C: ME of PET by cwd.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Marginal effects of CWD - ITRDB
pred_fia <- predictions(pet_mod_fia, newdata = datagrid(pet.spstd = 0, cwd.spstd = seq(cwd_range_fia[1],cwd_range_fia[2],seq_inc))) %>% 
  mutate(Dataset = "FIA")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(pet_mod_itrdb, newdata = datagrid(pet.spstd = at_pet, cwd.spstd = seq(cwd_range_itrdb[1],cwd_range_itrdb[2],seq_inc))) %>% 
  mutate(Dataset = "ITRDB")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
pet_cwdstd_plot <- pred_all %>% 
  ggplot(aes(x = cwd.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = Dataset, fill = Dataset)) + 
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme_bw()  +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to PET")   +
  scale_color_manual(values = c("ITRDB" = "darkred",
                                "FIA"="darkblue")) +
  scale_fill_manual(values = c("ITRDB" = "darkred",
                               "FIA"="darkblue"))
pet_cwdstd_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# D: ME of PET by pet.spstd ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Marginal effects of CWD - ITRDB
pred_fia <- predictions(pet_mod_fia, newdata = datagrid(pet.spstd = seq(pet_range_fia[1],pet_range_fia[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(Dataset = "FIA")

### Marginal effects of CWD - ITRDB
pred_itrdb <- predictions(pet_mod_itrdb, newdata = datagrid(pet.spstd = seq(pet_range_itrdb[1],pet_range_itrdb[2],seq_inc), cwd.spstd = 0)) %>% 
  mutate(Dataset = "ITRDB")

## Combine into single plot
pred_all <- rbind(pred_itrdb, pred_fia)
pet_petstd_plot <- pred_all %>% 
  ggplot(aes(x = pet.spstd, y = estimate, ymin = conf.low, ymax = conf.high, color = Dataset, fill = Dataset)) + 
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme_bw() +
  xlab("Historic PET\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to PET")  +
  scale_color_manual(values = c("ITRDB" = "darkred",
                                "FIA"="darkblue")) +
  scale_fill_manual(values = c("ITRDB" = "darkred",
                               "FIA"="darkblue"))
pet_petstd_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine plots ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadfonts(device = "win")

theme_set(
  theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

combined_plot <- (cwd_cwdstd_plot / cwd_petstd_plot) | (pet_cwdstd_plot / pet_petstd_plot)  

combined_plot <- combined_plot +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") & 
  plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=12))

combined_plot

ggsave(paste0(wdir, 'figures/FigS6_fia_compare.svg'), plot = combined_plot, bg= 'transparent', width = 10, height = 7)
