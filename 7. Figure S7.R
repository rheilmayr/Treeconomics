#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 8/7/23
# Purpose: Create figure S7 illustrating transect of different predictions
#
# Input files:
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
rwi_list <- list.files(paste0(wdir, "2_output/predictions/sp_rwi/"), pattern = ".gz", full.names = TRUE)
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
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep data  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_predictions <- sp_predictions %>% 
  mutate(cwd_change = cwd_cmip_end_mean - cwd_cmip_start_mean,
         pet_change = pet_cmip_end_mean - pet_cmip_start_mean)

plot_dat <- sp_predictions %>%
  filter(((abs(cwd_hist)<2.5) & (abs(pet_hist<2.5)))) %>%
  drop_na()

cwd_seq_min <- -2.125
cwd_seq_max <- 2.375
seq_inc <- 0.25
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Transect plots ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  geom_line(linewidth = 2) +
  ylim(c(-0.3, 0.2)) +
  xlim(c(-2, 2)) +
  xlab("Historic CWD\n(Deviation from species mean)") +
  ylab("Difference in predicted RWI changes by 2100\n(Heterogeneous sensitivity model - Neutral model)") +
  theme(legend.position = c(.18,.75),
        legend.text = element_text(size=13),
        legend.title = element_text(size=18),
        legend.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)
transect_0
ggsave(paste0(wdir, "3_results/figures/FigS7_dif_pred.svg"), transect_0, width = 8, height = 6)


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
