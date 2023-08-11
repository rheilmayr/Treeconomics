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
theme(family="Serif")

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
    theme(text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main summary figure ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_predictions <- sp_predictions %>% 
  mutate(cwd_change = cwd_cmip_end_mean - cwd_cmip_start_mean,
         pet_change = pet_cmip_end_mean - pet_cmip_start_mean)

plot_dat <- sp_predictions %>%
  filter(((abs(cwd_hist)<2.5) & (abs(pet_hist<2.5)))) %>%
  drop_na()

seq_inc <- 0.25
cwd_seq_min <- -2.125
cwd_seq_max <- 2.375

cwd_sequence <- seq(cwd_seq_min, cwd_seq_max, seq_inc)


pet_seq_min <- -2.125
pet_seq_max <- 2.375
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


### CWD sensitivity
cwd_sens_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_sens)) +
  geom_tile() +
  xlim(c(cwd_seq_min, cwd_seq_max)) +
  ylim(c(pet_seq_min, pet_seq_max)) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  # scale_fill_distiller(type = "div") +
  # scale_fill_viridis_c(direction = -1, option = "viridis") +
  theme(legend.position = c(.21,.82),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
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
  theme(legend.position = c(.21,.82),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
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
  theme(legend.position = c(.19,.83),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
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
  theme(legend.position = c(.19,.83),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
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
  theme(legend.position = c(.15,.83),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        legend.background = element_blank())+
  labs(fill = "Predicted \nchange in RWI") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()
rwi_bin

lgd_pos <- c(.15, .78)
cwd_sens_bin <- cwd_sens_bin +
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    # plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
    legend.key.size = unit(10, "pt"))

pet_sens_bin <- pet_sens_bin +
  theme(
    # plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(10, "pt"))

sens_plot <- cwd_sens_bin / pet_sens_bin
#ggsave(paste0(wdir, "figures\\", "pred_full_a.svg"), sens_plot, width = 9, height = 15)


lgd_pos <- c(.25, .8)
cwd_change_bin <- cwd_change_bin +
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    # plot.margin = margin(t=0, r=0, b=0, l=-1, "cm"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(10, "pt"))


pet_change_bin <- pet_change_bin + 
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    # plot.margin = margin(t=0, r=0, b=0, l=-1, "cm"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(10, "pt"))


lgd_pos <- c(.1, .9)
rwi_bin <- rwi_bin +
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(10, "pt")
  )
rwi_bin


pred_full <- (cwd_sens_bin / pet_sens_bin) | plot_spacer() | (cwd_change_bin / pet_change_bin) | plot_spacer() |rwi_bin
pred_full <- pred_full + plot_layout(widths = c(3, -1.01 , 4, -1.01, 6)) & 
  plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=12))
ggsave(paste0(wdir, "figures\\", "Fig5_pred_full.svg"), pred_full, width = 14, height = 7)
