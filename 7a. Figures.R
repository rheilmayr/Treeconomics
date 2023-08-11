#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run regressions to explore impact of historical climate on weather sensitivity
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using species_niche.R
# - tree_ring_data_V2.db: Compiled database of ITRDB observations 
# - essentialcwd_data.csv: File detailing plot-level weather history
#
# ToDo:
# - Think through places where Monte carlo uncertainties can be added to figures
# - track down lost observations - currently dropping a lot due to NAN or failed RWI generation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
library(stringr)
library(raster)
library(rgdal)
library(viridis)
library(patchwork)
library(Hmisc)
library(latticeExtra)
library(prediction)
library(colorspace)
library(ggnewscale)
library(reshape2)
library(ggExtra)
library(ggsn)
library(maptools)
library(broom)
library(ggExtra)
library(extrafont)
library(marginaleffects)
library(tmap)
library(fixest)
library(forcats)
library(car)
librarian::shelf(ggplotify)

loadfonts(device = "win")
theme(family="Serif")

select <- dplyr::select
summarize <- dplyr::summarize

options(scipen=999)

theme_set(
  theme_bw(base_size = 25)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)
site_loc <- site_smry %>%
  select(collection_id, latitude, longitude)
flm_df <- flm_df %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))

trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

# DNLM results
dnlm_results <- read_rds(paste0(wdir, "out/first_stage/dnlm_lagged_effects"))
# dnlm_results <- read_rds(paste0(wdir, "out/first_stage/dnlm_orig")) # Something has changed - old version creates positive pet effect, but gone after july tweaks to data...

# 5. Prediction rasters
rwi_list <- list.files(paste0(wdir, "out/predictions/pred_10000/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))

## 6. Second stage model
mod_df <- read_rds(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))

# 7. Genus second stage model
genus_models <- readRDS(paste0(wdir, "out/second_stage/ss_conley_genus.rds"))

sp_models <- readRDS(paste0(wdir, "out/second_stage/ss_conley_species.rds"))

# 2. Species range maps
range_file <- paste0(wdir, 'in/species_ranges/merged_ranges.shp')
range_sf <- st_read(range_file)

# 8. Change in CWD
cmip_end <- load(paste0(wdir, 'in/CMIP5 CWD/cmip5_cwdaet_end.Rdat'))
cwd_cmip_end <- cwd_raster %>% mean()
rm(cwd_raster)
rm(aet_raster)

cmip_start <- load(paste0(wdir, 'in/CMIP5 CWD/cmip5_cwdaet_start.Rdat'))
cwd_cmip_start <- cwd_raster %>% mean()
rm(cwd_raster)
rm(aet_raster)
cwd_cmip_change <- cwd_cmip_end - cwd_cmip_start
cwd_cmip_pchange <- cwd_cmip_change / cwd_cmip_start

# cwdlist=list();aetlist=list()
# j=1
# for(i in c("start","end")){
#   load(paste0(wdir,"in/CMIP5 CWD/cmip5_cwdaet_",i,".Rdat"))
#   cwdlist[[j]]=cwd_raster;aetlist[[j]]=aet_raster
#   j=j+1
# }

# # 2. Species range maps
# range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
# range_sf <- st_read(range_file) %>% 
#   filter(sp_code %in% (trim_df %>% pull(species_id) %>% unique()))
# 
# # 1. Historic climate raster
# clim_file <- paste0(wdir, 'in//CRUData//historic_raster//HistoricCWD_AETGrids_Annual.Rdat')
# load(clim_file)
# cwd_historic <- cwd_historic %>% mean(na.rm = TRUE)


# # 4. Site information
# site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
# site_smry <- site_smry %>% 
#   select(collection_id, sp_id, latitude, longitude) %>% 
#   mutate(species_id = tolower(sp_id)) %>% 
#   select(-sp_id)
# 
# # 5. Species information
# sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
# sp_info <- sp_info %>% 
#   select(species_id, genus, gymno_angio, family)
# site_smry <- site_smry %>% 
#   left_join(sp_info, by = c("species_id"))
# 
# # 5. Prediction rasters
# sp_predictions <- readRDS(paste0(wdir, "out/predictions/sq_sp_predictions.rds"))
# 
# 

# 2. Historic climate raster
clim_file <- paste0(wdir, 'in/CRUData/historic_raster/HistoricCWD_AETGrids.Rdat')
load(clim_file)
cwd_historic <- sum(cwd_historic)
names(cwd_historic) = "cwd"


# Average site conditions
ave_site_clim_df <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define palettes ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_palette <- scale_colour_brewer(
  type = "seq",
  palette = 5,
  direction = 1,
  aesthetics = "colour"
)











#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main summary figure ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_predictions <- sp_predictions %>% 
  mutate(cwd_change = cwd_cmip_end_mean - cwd_cmip_start_mean,
         pet_change = pet_cmip_end_mean - pet_cmip_start_mean)

plot_dat <- sp_predictions %>%
  filter(((abs(cwd_hist)<2.5) & (abs(pet_hist<2.5)))) %>%
  drop_na()

cwd_seq_min <- -2.125
cwd_seq_max <- 2.375
cwd_seq_inc <- 0.25
cwd_sequence <- seq(cwd_seq_min, cwd_seq_max, seq_inc)


pet_seq_min <- -2.125
pet_seq_max <- 2.375
pet_seq_inc <- 0.25
pet_sequence <- seq(pet_seq_min, pet_seq_max, seq_inc)
convert_bin <- function(n){
  cwd_sequence[n] + 0.125
}
plot_dat <- plot_dat %>% 
  mutate(cwd.q = cut(cwd_hist, breaks = cwd_sequence, labels = FALSE),
         cwd.q = convert_bin(cwd.q),
         pet.q = cut(pet_hist, breaks = pet_sequence, labels = FALSE),
         pet.q = convert_bin(pet.q))


plot_dat <- plot_dat %>% 
  group_by(cwd.q, pet.q) %>% 
  summarize(rwi_pred = mean(rwi_pred_mean, na.rm = TRUE),
            # rwi_change = mean(rwi_pred_mean, na.rm = TRUE),
            # rwi_change_lb = mean(rwi_pred_025, na.rm = TRUE),
            # rwi_change_ub = mean(rwi_pred_975, na.rm = TRUE),
            # rwi_change_pclim = mean(rwi_pclim_mean, na.rm = TRUE),
            # rwi_change_pclim_lb = mean(rwi_pclim_025),
            # rwi_change_pclim_ub = mean(rwi_pclim_975),
            rwi_change = mean(rwi_pred_change_mean, na.rm = TRUE),
            rwi_change_lb = mean(rwi_pred_change_025, na.rm = TRUE),
            rwi_change_ub = mean(rwi_pred_change_975, na.rm = TRUE),
            rwi_change_pclim = mean(rwi_pclim_change_mean, na.rm = TRUE),
            rwi_change_pclim_lb = mean(rwi_pclim_change_025),
            rwi_change_pclim_ub = mean(rwi_pclim_change_975),
            cwd_sens = mean(cwd_sens, na.rm = TRUE),
            pet_sens = mean(pet_sens, na.rm = TRUE),
            cwd_change = mean(cwd_change, na.rm = TRUE),
            pet_change = mean(pet_change, na.rm = TRUE),
            rwi_dif = mean(rwi_pred_pclim_change_dif_mean, na.rm = TRUE),
            rwi_dif_ub = mean(rwi_pred_pclim_change_dif_025, na.rm = TRUE),
            rwi_dif_lb = mean(rwi_pred_pclim_change_dif_975, na.rm = TRUE),
            n = n()) %>%
  filter(n>5)


# theme_set(
#   theme_bw(base_size = 24)+
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# )

### CWD sensitivity
cwd_sens_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_sens)) +
  geom_tile() +
  xlim(c(cwd_seq_min, cwd_seq_max)) +
  ylim(c(pet_seq_min, pet_seq_max)) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  # scale_fill_distiller(type = "div") +
  # scale_fill_viridis_c(direction = -1, option = "viridis") +
  theme_bw(base_size = 35)+
  theme(legend.position = c(.21,.82),
        legend.text = element_text(size=18),
        legend.title = element_text(size=23),
        legend.background = element_blank())+
  labs(fill = "Predicted\nsensitivity\nto CWD") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()+
  xlab("")
cwd_sens_bin


### PET sensitivity
pet_sens_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = pet_sens)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  # scale_fill_viridis_c(direction = -1, option = "viridis") +
  theme_bw(base_size = 35)+
  theme(legend.position = c(.21,.82),
        legend.text = element_text(size=18),
        legend.title = element_text(size=23),
        legend.background = element_blank())+
  labs(fill = "Predicted\nsensitivity\nto PET") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()
pet_sens_bin


### CWD change
cwd_change_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_change)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_viridis_c(direction = 1, option = "magma") +
  theme_bw(base_size = 35)+
  theme(legend.position = c(.19,.83),
        legend.text = element_text(size=18),
        legend.title = element_text(size=23),
        legend.background = element_blank())+
  labs(fill = "Predicted change\nin CWD (std)") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()+
  ylab("")+
  xlab("")
cwd_change_bin

### PET change
pet_change_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = pet_change)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_viridis_c(direction = 1, option = "magma") +
  theme_bw(base_size = 35)+
  theme(legend.position = c(.19,.83),
        legend.text = element_text(size=18),
        legend.title = element_text(size=23),
        legend.background = element_blank())+
  labs(fill = "Predicted change\nin PET (std)") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()+
  ylab("")
pet_change_bin


### RWI change
rwi_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = rwi_change)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_viridis_c(direction = -1, option = "viridis") +
  # scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  theme_bw(base_size = 35)+
  theme(legend.position = c(.15,.83),
        legend.text = element_text(size=18),
        legend.title = element_text(size=23),
        legend.background = element_blank())+
  labs(fill = "Predicted \nchange in RWI") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()
rwi_bin

lgd_pos <- c(.15, .8)
cwd_sens_bin <- cwd_sens_bin +
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t=0, r=0, b=0, l=0, "cm"))

pet_sens_bin <- pet_sens_bin +
  theme(
    plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank())

sens_plot <- cwd_sens_bin / pet_sens_bin
#ggsave(paste0(wdir, "figures\\", "pred_full_a.svg"), sens_plot, width = 9, height = 15)


lgd_pos <- c(.23, .8)
cwd_change_bin <- cwd_change_bin +
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank())


pet_change_bin <- pet_change_bin + 
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank())


lgd_pos <- c(.15, .8)
rwi_bin <- rwi_bin +
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
rwi_bin


pred_full <- cwd_sens_bin/ pet_sens_bin | cwd_change_bin/pet_change_bin | rwi_bin
pred_full
ggsave(paste0(wdir, "figures\\", "4_pred_full.svg"), pred_full, width = 40, height = 35)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Transect plots ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# pred_dat <- plot_dat %>% 
#   select(cwd.q, pet.q, 
#          rwi_change = rwi_change, 
#          rwi_change_lb = rwi_change_lb, 
#          rwi_change_ub = rwi_change_ub) %>% 
#   mutate(scenario = "rwi_change")
# 
# pclim_dat <- plot_dat %>% 
#   select(cwd.q, pet.q, 
#          rwi_change = rwi_change_pclim, 
#          rwi_change_lb = rwi_change_pclim_lb, 
#          rwi_change_ub = rwi_change_pclim_ub) %>% 
#   mutate(scenario = "rwi_change_pclim")
# 
# 
# transect_dat <- pred_dat %>% 
#   rbind(pclim_dat) %>% 
#   mutate(scenario = fct_relevel(scenario, "rwi_change", "rwi_change_pclim"))
# 
# transect_1 <- transect_dat %>% 
#   filter(pet.q == 0) %>%
#   group_by(cwd.q, scenario) %>% 
#   summarise(rwi_change = mean(rwi_change),
#             rwi_change_lb = mean(rwi_change_lb),
#             rwi_change_ub = mean(rwi_change_ub)) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_ribbon(aes(ymin = rwi_change_lb,
#                   ymax = rwi_change_ub,
#                   fill = scenario),
#               alpha = 0.2) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-1, 0.2)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted")) +
#   scale_fill_manual(name = "Scenario",
#                     labels = c("Full model", 
#                                "Variable shift in climate,\nconstant sensitivity"), 
#                     values = c("dark blue", "dark red", "dark green")) +
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std above mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# transect_1
# 
# 
# transect_2 <- transect_dat %>% 
#   filter(pet.q == 2) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_ribbon(aes(ymin = rwi_change_lb,
#                   ymax = rwi_change_ub,
#                   fill = scenario),
#               alpha = 0.2) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-1.1, 0.2)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted"))+
#   scale_fill_manual(name = "Scenario",
#                     labels = c("Full model", 
#                                "Variable shift in climate,\nconstant sensitivity"), 
#                     values = c("dark blue", "dark red", "dark green")) +
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std below mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.25),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# transect_2
# 
#   locator <- rwi_bin + 
#   theme_bw(base_size = 20)+
#   theme(legend.position = c(.18,.83),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank())+
#   geom_hline(yintercept = 1, size = 1) + 
#   geom_hline(yintercept = -1, size = 1)
# 
# locator | transect_1 / transect_2




transect_0 <- plot_dat %>% 
  filter(pet.q == 0) %>%
  group_by(cwd.q) %>% 
  summarise(rwi_change = mean(rwi_dif),
            rwi_change_lb = mean(rwi_dif_lb),
            rwi_change_ub = mean(rwi_dif_ub)) %>% 
  ggplot(aes(x = cwd.q, y = rwi_change)) +
  geom_ribbon(aes(ymin = rwi_change_lb,
                  ymax = rwi_change_ub),
              alpha = 0.2) +
  geom_line(size = 2) +
  theme_bw(base_size = 20)+
  ylim(c(-0.3, 0.2)) +
  xlim(c(-2, 2)) +
  # ggtitle("Historic PET = historic species mean") +
  # ylab("Predicted difference in RWI change - neutral model vs ourse") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  ylab("Difference in predicted RWI changes by 2100\n(Heterogeneous sensitivity model - Neutral model)") +
  theme(legend.position = c(.18,.75),
        legend.text = element_text(size=13),
        legend.title = element_text(size=18),
        legend.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)
transect_0
ggsave(paste0(wdir, "figures\\", "a5_dif_pred.svg"), transect_0, width = 8, height = 8)


# ## Alternate version that plots difference in two models
# transect_n1 <- plot_dat %>% 
#   filter(pet.q == -1) %>%
#   group_by(cwd.q) %>% 
#   summarise(rwi_change = mean(rwi_dif),
#             rwi_change_lb = mean(rwi_dif_lb),
#             rwi_change_ub = mean(rwi_dif_ub)) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change)) +
#   geom_ribbon(aes(ymin = rwi_change_lb,
#                   ymax = rwi_change_ub),
#               alpha = 0.2) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-0.4, 0.3)) +
#   xlim(c(-2, 2)) +
#   ggtitle("Historic PET = 1 std below mean") +
#   # ylab("Predicted difference in RWI change - neutral model vs ourse") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# transect_n1
# 
# 
# 
# transect_1 <- plot_dat %>% 
#   filter(pet.q == 1) %>%
#   group_by(cwd.q) %>% 
#   summarise(rwi_change = mean(rwi_dif),
#             rwi_change_lb = mean(rwi_dif_lb),
#             rwi_change_ub = mean(rwi_dif_ub)) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change)) +
#   geom_ribbon(aes(ymin = rwi_change_lb,
#                   ymax = rwi_change_ub),
#               alpha = 0.2) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-0.4, 0.3)) +
#   xlim(c(-2, 2)) +
#   ggtitle("Historic PET = 1 std above mean") +
#   # ylab("Predicted difference in RWI change - neutral model vs ours") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# transect_1
# 
# plot <- transect_1 / transect_0 / transect_n1
# plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Changes in CWD and PET  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#overlay species ranges to only show changes for areas with species we are considering
ranges_dissolve <- st_read(paste0(wdir,"in/species_ranges/merged_ranges_dissolve.shp"))

#clip raster to species range
cwd_cmip_change=raster::mask(cwd_cmip_change, ranges_dissolve)

data(World)


cwd_cmip_df <- as.data.frame(cwd_cmip_change, xy = TRUE)
world <- map_data("world")

theme_set(
  theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

cwd_map <- ggplot() +
  geom_map(data = world, map = world,
           aes(map_id = region),
           color = "lightgray", fill = "lightgray",size = 0.1) +
  geom_tile(data = cwd_cmip_df %>% filter(!is.na(layer)), aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(name = "Change in CWD\n(mm/year)") +
  coord_fixed() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
cwd_map
#dim 7x5
  
  
  
  



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Species-level changes in CWD  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gen_dat <- read_csv(paste0(wdir, "out/species_gen_gr.csv")) %>% 
  rename(sp_code = "species_id")

sp_predictions_gen <- sp_predictions %>% 
  left_join(gen_dat)

sp_plot_dat <- sp_predictions_gen %>% 
  group_by(sp_code) %>% 
  mutate(beyond_max_cwd = cwd_cmip_end_mean > max(cwd_hist)) %>% 
  ungroup() %>% 
  mutate(cwd_end_bin = case_when(
    (beyond_max_cwd == TRUE) ~ "Beyond prior max",
    (cwd_cmip_end_mean > 1.5) ~ "1-2 s.d. above prior mean",
    (cwd_cmip_end_mean >= 0) ~ "0-1 s.d. above prior mean",
    (cwd_cmip_end_mean < 0) ~ "Below prior mean"),
    cwd_end_bin = fct_relevel(cwd_end_bin, 
                              "Beyond prior max", "1-2 s.d. above prior mean", 
                              "0-1 s.d. above prior mean", 
                              "Below prior mean"))

bin_shares <- sp_plot_dat %>% # Add grouping by genera? 
  group_by(genus, cwd_end_bin) %>% 
  summarise(n = n()) %>% 
  group_by(genus) %>% 
  mutate(prop_cwd_bin = prop.table(n)) %>% 
  filter(cwd_end_bin == "Beyond prior max") %>% 
  arrange(prop_cwd_bin) %>% 
  select(genus, prop_cwd_bin)

sp_plot_dat <- sp_plot_dat %>% 
  left_join(bin_shares, by = "genus") %>% 
  mutate(genus = factor(genus))

sp_plot_dat$genus <- fct_reorder(sp_plot_dat$genus, sp_plot_dat$prop_cwd_bin, median) ## Probably want to define order based on multiple categories to make smoother plot


theme_set(
  theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


genus_cwd_change <- sp_plot_dat %>%
  ggplot(aes(x= genus, fill=cwd_end_bin)) +
  geom_bar(position = "fill", alpha=.9) +
  scale_fill_viridis_d(direction = -1) + # Use diverging color scheme? 
  ylab("Proportion of range")+
  xlab("Genera")+
  guides(fill=guide_legend("Change in CWD"))+
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))

genus_cwd_change

cwd_change <- sp_plot_dat %>% 
  ggplot(aes(x = cwd_cmip_start_mean, y = cwd_cmip_end_mean)) +
  geom_hex(bins = 60) +
  scale_fill_viridis(limits = c(0,5000), option = "magma") +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  xlim(-5, 20) +
  ylim(-5, 20) +
  xlab("Mean CWD in 1970-2000\n(Deviation from species mean)") +
  ylab("Mean CWD in 2091-2100\n(Deviation from species mean)") +
  labs(fill = "Count of\ngrid cells") +
  theme_bw(base_size = 12)

pet_change <- sp_plot_dat %>% 
  ggplot(aes(x = pet_cmip_start_mean, y = pet_cmip_end_mean)) +
  geom_hex(bins = 60) +
  scale_fill_viridis(limits = c(0,5000), option = "magma") +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  xlim(-5, 18) +
  ylim(-5, 18) +
  xlab("Mean PET in 1970-2000\n(Deviation from species mean)") +
  ylab("Mean PET in 2091-2100\n(Deviation from species mean)") +
  labs(fill = "Count of\ngrid cells") +
  theme_bw(base_size = 12)
# +
#   theme(legend.position = c(.9,.2),
#         legend.background = element_blank())

change_plot <- (cwd_change | pet_change) + 
  plot_layout(guides = "collect")
change_plot

# cwd_change_fig <- cwd_map / change_plot / genus_cwd_change +
#   plot_layout(heights = c(1.4, 1.7, 1))

cwd_change_fig <- change_plot / genus_cwd_change 


#ggsave(paste0(wdir, "figures\\", "3_cwd_change.svg"), cwd_change_fig, width = 7.5, height = 11)

cwd_change_fig


rwi_change <- sp_plot_dat %>% 
  filter(pet_hist > 1 & pet_hist < 1.5) %>% 
  ggplot(aes(x = cwd_hist, y = rwi_pred_change_mean)) +
  geom_hex(bins = 30) +
  scale_fill_viridis() +
  xlab("Mean CWD in 1970-2000\n(Deviation from species mean)") +
  ylab("RWI") +
  xlim(c(-2,3)) +
  labs(fill = "Count of\ngrid cells") +
  theme_bw(base_size = 12)
rwi_change

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Species-level changes in RWI  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_plot_dat <- sp_plot_dat %>% 
  mutate(rwi_change_bin = case_when(
    (rwi_pred_change_mean < -.5) ~ "1. RWI decline > 0.5",
    (rwi_pred_change_mean < -0.25) ~ "2. RWI decline between 0.25 and .5",
    (rwi_pred_change_mean < 0) ~ "3. RWI decline between 0 and 0.25",
    (rwi_pred_change_mean >= 0) ~ "4. RWI increase",
  ))
sp_plot_dat %>% 
  ggplot(aes(x = sp_code, fill = rwi_change_bin)) +
  geom_bar(position = "fill")


# transect_dat <- plot_dat %>% 
#   select(cwd.q, pet.q, rwi_change, rwi_change_psens, rwi_change_pclim) %>% 
#   pivot_longer(c(rwi_change, rwi_change_psens, rwi_change_pclim), names_to = 'scenario', values_to = "rwi_change") %>% 
#   mutate(scenario = fct_relevel(scenario, "rwi_change", "rwi_change_psens", "rwi_change_pclim"))
# 
# transect_1 <- transect_dat %>% 
#   filter(pet.q == 1) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-0.6, 0.2)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted")) +
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Constant shift in climate,\nvariable sensitivity",
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std above mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# 
# 
# transect_2 <- transect_dat %>% 
#   filter(pet.q == -1) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-0.4, 0.4)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted"))+
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Constant shift in climate,\nvariable sensitivity",
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std below mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#       legend.text = element_text(size=13),
#       legend.title = element_text(size=18),
#       legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# 
# 
# locator <- rwi_bin + 
#   theme_bw(base_size = 20)+
#   theme(legend.position = c(.18,.83),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank())+
#   geom_hline(yintercept = 1, size = 1) + 
#   geom_hline(yintercept = -1, size = 1)
# 
# locator | transect_1 / transect_2




# 
# lgd_pos <- c(.2, .8)
# cwd_sens_bin <- cwd_sens_bin +
#   theme(
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"))
# ggsave(paste0(wdir, "figures\\", "pres_cwd_sens.svg"), cwd_sens_bin, width = 9, height = 9)
# 
# pet_sens_bin <- pet_sens_bin +
#   theme(
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank())
# ggsave(paste0(wdir, "figures\\", "pres_pet_sens.svg"), pet_sens_bin, width = 9, height = 9)
# 

# cwd_change_bin <- cwd_change_bin +
#   theme(
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"))
# cwd_change_bin
# ggsave(paste0(wdir, "figures\\", "pres_cwd_change.svg"), cwd_change_bin, width = 9, height = 9)

# pet_change_bin <- pet_change_bin + 
#   theme(
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank())
# pet_change_bin
# ggsave(paste0(wdir, "figures\\", "pres_pet_change.svg"), pet_change_bin, width = 9, height = 9)
# 
# change_plot <- cwd_change_bin / pet_change_bin
# ggsave(paste0(wdir, "figures\\", "pred_full_b.svg"), change_plot, width = 9, height = 15)
# 
# 


# test +
#   scale_fill_continuous(type="viridis", name="Change in CWD\n1970-2000 to 2091-2100\n(mm per year)")
# 
# test %>% tmap_grob()
# 
# tm_shape(land) +
#   tm_raster(data = cwd_cmip_change, palette = terrain.colors(10)) +
# 
# map <-ggplot() +
#   geom_map(data = world, map = world,
#            aes(x=long, y=lat, map_id = region),
#            color = "black", fill = "lightgray", size = 0.1,alpha=.7) +
#   theme_bw(base_size =20) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_raster(data = cwd_cmip_change)
# 
# 
# 
#   geom_point(data = site_loc, aes(y = latitude, x = longitude), color = "#443A83FF")+
#   #geom_sf(data = flm_df, color = "#443A83FF", alpha = .2) +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   scale_x_continuous("Longitude", breaks = xr, labels = xlabels) +
#   scale_y_continuous("Latitude", breaks = yr, labels = ylabels)
# 
# map
# 
# 
# 
# 
# 
# # cwdmeanchange_end=mean(cwdlist[[3]])-mean(cwdlist[[1]])
# 
# 
# # ranges=st_read(paste0(wdir,"in/species_ranges/merged_ranges.shp"))
# r <- raster(nrows=nrow(cwd_cmip_change),ncols=ncol(cwd_cmip_change));extent(r)=extent(cwd_cmip_change)
# ranges_raster=rasterize(ranges_dissolve,r,field=1)
# test=rasterToPolygons(ranges_raster,dissolve = TRUE)
# testsf=st_as_sf(test)
# 
# plot(ranges)
# 
# world <- ne_download(scale = "medium", type="land",returnclass = "sf",category="physical")
# world=st_crop(world,extent(cwdmeanchange_end))
# world=st_transform(world,projection(cwdmeanchange_end))
# 
# library(reshape2)
# temp=as.matrix(cwd_cmip_change)
# colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# temp=melt(temp)
# colnames(temp)=c("lat","long","cwd")
# 
# area_grid=raster::area(cwdmeanchange_end)
# 
# #find global spatial average of cwd and pet changes - calculate areas of grid cells
# # cwdchanges=cellStats((mask(cwdlist[[3]],test)-mask(cwdlist[[1]],test))*area_grid,stat="sum")/cellStats(area_grid,stat="sum")
# # min=which.min(cwdchanges);max=which.max(cwdchanges)
# # mincwd=cwdlist[[3]][[min]]-cwdlist[[1]][[min]];maxcwd=cwdlist[[3]][[max]]-cwdlist[[1]][[max]]
# 
# world_raster=rasterize(world,r,field=1)
# world_raster=data.frame(temp[,1:2],melt(as.matrix(world_raster))[,3])
# colnames(world_raster)=c("lat","long","land")
# world_raster$land[which(world_raster$land==1)]="land"
# 
# a=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# a=a+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
# a=a+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# a=a+scale_fill_continuous(type="viridis",na.value="transparent",name="Change in CWD\n1970-2000 to 2091-2100\n(mm per year)")
# a=a+geom_raster()
# a=a+annotate("text",x=-150,y=-45,size=10,label="a)")
# # 
# # 
# # mincwd=raster::mask(mincwd,test)
# # temp=as.matrix(mincwd)
# # colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# # temp=melt(temp)
# # colnames(temp)=c("lat","long","cwd")
# # 
# # b=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# # b=b+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
# # b=b+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# # b=b+scale_fill_continuous(type="viridis",na.value="transparent",limits=c(-200,2700),guide=NULL)
# # b=b+annotate("text",x=-142,y=-45,size=8,label="Min")+geom_raster()
# # 
# # maxcwd=raster::mask(maxcwd,test)
# # temp=as.matrix(maxcwd)
# # colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# # temp=melt(temp)
# # colnames(temp)=c("lat","long","cwd")
# # 
# # d=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# # d=d+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="lightgrey")
# # d=d+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# # d=d+scale_fill_continuous(type="viridis",na.value="transparent",limits=c(-200,2700),guide=NULL)
# # d=d+annotate("text",x=-142,y=-45,size=8,label="Max")+geom_raster()
# 
# #a/(b+d)
# 
# #get spatial averages in cwd and pet across species range for all three time points for all model runs
# 
# petlist=list()
# for(i in 1:3) petlist[[i]]=cwdlist[[i]]+aetlist[[i]]
# 
# cwdchange_dist=matrix(nrow=dim(cwdlist[[1]])[3],ncol=length(cwdlist))
# petchange_dist=matrix(nrow=dim(petlist[[1]])[3],ncol=length(petlist))
# 
# area_mask=mask(area_grid,test)
# 
# for(i in 1:length(cwdlist)){
#   cwd_temp=cellStats(cwdlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
#   pet_temp=cellStats(petlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
#   cwdchange_dist[,i]=cwd_temp;petchange_dist[,i]=pet_temp
#   print(i)
# }
# 
# colnames(cwdchange_dist)=c("Historic","Mid-Century","End-Century")
# colnames(petchange_dist)=c("Historic","Mid-Century","End-Century")
# cwdchange_dist=as.data.frame(cwdchange_dist);cwdchange_dist$Model=1:nrow(cwdchange_dist)
# petchange_dist=as.data.frame(petchange_dist);petchange_dist$Model=1:nrow(petchange_dist)
# 
# cwdchange_dist=cwdchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="CWD")
# petchange_dist=petchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="PET")
# 
# changedist=cbind(cwdchange_dist,petchange_dist[,3])
# changedist=changedist%>%pivot_longer(c("CWD","PET"),names_to="Variable",values_to="CWD_PET")
# changedist$Time=ordered(changedist$Time,c("Historic","Mid-Century","End-Century"))
# changedist$Variable=ordered(changedist$Variable,c("PET","CWD"))
# historic=changedist%>%group_by(Variable,Time)%>%dplyr::summarise(hist_mean=quantile(CWD_PET,0.5))
# changedist$CWD_PET[which(changedist$Time=="Historic")]=NA
# 
# e=ggplot(changedist,aes(x=Time,y=CWD_PET,col=Variable))+geom_boxplot(position="identity",width=0.15,inherit.aes=FALSE,aes(x=Time,y=CWD_PET,group=interaction(Variable,Time)),col="black",outlier.shape = NA,lwd=0.75)
# e=e+scale_y_continuous(limits=c(400,1500),name="CWD or PET (mm per yer)")
# e=e+geom_jitter(width=0.1)+theme_bw()
# e=e+geom_line(data=historic,aes(x=Time,y=hist_mean,group=Variable),col="black",lty=2)
# e=e+geom_point(data=historic%>%filter(Time=="Historic"),aes(y=hist_mean,col=Variable),x="Historic",size=3)
# e=e+labs(x="",color="")
# e=e+scale_color_manual(values=c("#972c25","#8e9169"))
# e=e+scale_x_discrete(labels=c("Historic"="Historic\n1970-2000","Mid-Century"="Mid-Century\n2045-2055","End-Century"="End-Century\n2090-2100"))
# e=e+annotate("text",x=0.7,y=470,size=10,label="b)")
# x11()
# a/e
# 
# ###Supplementary Figure Showing PET map with AET and PET changes
# 
# petmeanchange_end=mean(petlist[[3]])-mean(petlist[[1]])
# 
# petmeanchange_end=raster::mask(petmeanchange_end,test)
# temp=as.matrix(petmeanchange_end)
# colnames(temp)=xFromCol(petmeanchange_end);rownames(temp)=yFromRow(petmeanchange_end)
# temp=melt(temp)
# colnames(temp)=c("lat","long","cwd")
# 
# a=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# a=a+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
# a=a+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# a=a+scale_fill_continuous(type="viridis",na.value="transparent",name="Change in PET\n1970-2000 to 2091-2100\n(mm per year)")
# a=a+geom_raster()
# a=a+annotate("text",x=-150,y=-45,size=10,label="a)")
# 
# #get model averages of AET changes
# aetchange_dist=matrix(nrow=dim(aetlist[[1]])[3],ncol=3)
# for(i in 1:length(cwdlist)){
#   aet_temp=cellStats(aetlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
#   aetchange_dist[,i]=aet_temp
#   print(i)
# }
# 
# colnames(aetchange_dist)=c("Historic","Mid-Century","End-Century")
# aetchange_dist=as.data.frame(aetchange_dist);aetchange_dist$Model=1:nrow(aetchange_dist)
# aetchange_dist=aetchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="AET")
# 
# changedist=cbind(petchange_dist,aetchange_dist[,3])
# changedist=changedist%>%pivot_longer(c("PET","AET"),names_to="Variable",values_to="PET_AET")
# changedist$Time=ordered(changedist$Time,c("Historic","Mid-Century","End-Century"))
# historic=changedist%>%group_by(Variable,Time)%>%dplyr::summarise(hist_mean=quantile(PET_AET,0.5))
# changedist$PET_AET[which(changedist$Time=="Historic")]=NA
# 
# e=ggplot(changedist,aes(x=Time,y=PET_AET,col=Variable))+geom_boxplot(position="identity",width=0.15,inherit.aes=FALSE,aes(x=Time,y=PET_AET,group=interaction(Variable,Time)),col="black",outlier.shape = NA,lwd=0.75)
# e=e+scale_y_continuous(limits=c(400,1500),name="PET or AET (mm per yer)")
# e=e+geom_jitter(width=0.1)+theme_bw()
# e=e+geom_line(data=historic,aes(x=Time,y=hist_mean,group=Variable),col="black",lty=2)
# e=e+geom_point(data=historic%>%filter(Time=="Historic"),aes(y=hist_mean,col=Variable),x="Historic",size=3)
# e=e+labs(x="",color="")
# e=e+scale_color_manual(values=c("#579473","#972c25"))
# e=e+scale_x_discrete(labels=c("Historic"="Historic\n1970-2000","Mid-Century"="Mid-Century\n2045-2055","End-Century"="End-Century\n2090-2100"))
# e=e+annotate("text",x=0.7,y=1400,size=10,label="b)")
# 
# x11()
# a/e
