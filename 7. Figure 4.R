#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 8/7/23
# Purpose: Create figure S1 illustrating data context
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

select <- dplyr::select
summarize <- dplyr::summarize

options(scipen=999)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Prediction rasters
rwi_list <- list.files(paste0(wdir, "out/predictions/pred_10000/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))
gen_dat <- read_csv(paste0(wdir, "out/species_gen_gr.csv")) %>% 
  rename(sp_code = "species_id")
sp_predictions_gen <- sp_predictions %>% 
  left_join(gen_dat)


# 2. Change in CWD
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define style  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_palette <- scale_colour_brewer(
  type = "seq",
  palette = 5,
  direction = 1,
  aesthetics = "colour"
)

theme_set(
  theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Panel A: Map of CWD changes ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data(World)
cwd_cmip_df <- as.data.frame(cwd_cmip_change, xy = TRUE)
world <- map_data("world")

cwd_map <- ggplot() +
  geom_map(data = world, map = world,
           aes(map_id = region),
           color = "lightgray", fill = "lightgray",linewidth = 0.1) +
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
# Panels B-C: Map of CWD changes ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  labs(fill = "Count of\ngrid cells")

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
  labs(fill = "Count of\ngrid cells")
# +
#   theme(legend.position = c(.9,.2),
#         legend.background = element_blank())

change_plot <- (cwd_change | pet_change) + 
  plot_layout(guides = "collect")
change_plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Genus-level changes in CWD  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


genus_cwd_change <- sp_plot_dat %>%
  ggplot(aes(x= genus, fill=cwd_end_bin)) +
  geom_bar(position = "fill", alpha=.9) +
  scale_fill_viridis_d(direction = -1) + # Use diverging color scheme? 
  ylab("Proportion of range")+
  xlab("Genera")+
  guides(fill=guide_legend("Change in CWD"))+
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))

genus_cwd_change

# cwd_change_fig <- cwd_map / change_plot / genus_cwd_change +
#   plot_layout(heights = c(1.4, 1.7, 1))

cwd_change_fig <- cwd_map / change_plot / genus_cwd_change  +
  plot_layout(heights = c(2, 2.5, 3)) &
  plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=12))

cwd_change_fig

ggsave(paste0(wdir, "figures/", "Fig4_cwd_change.svg"), cwd_change_fig, width = 10, height = 11)
