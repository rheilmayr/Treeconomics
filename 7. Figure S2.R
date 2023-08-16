#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 8/11/23
# Purpose: Create figure S2 comparing climate distribution of species ranges and itrdb
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

# 2. Prediction rasters
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
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ITRDB bias plot --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CWD
fullrange_cwd <- sp_predictions %>% 
  select(cwd_hist) 
fullrange_quantiles <- (fullrange_cwd$cwd_hist) %>% 
  quantile(probs = seq(0, 1, 0.01))

itrdb_cwd <- flm_df %>% 
  select(cwd.spstd)
itrdb_quantiles <- (itrdb_cwd$cwd.spstd) %>% 
  quantile(probs = seq(0, 1, 0.01))

itrdb_hist <- itrdb_cwd %>% 
  ggplot(aes(x = cwd.spstd)) +
  geom_histogram(bins = 50) +
  xlim(-3, 5) +
  theme_bw() +
  ggtitle("CWD frequency among ITRDB sites") +
  xlab("Historic CWD (Deviation from species mean)")

fullrange_hist <- fullrange_cwd %>% 
  ggplot(aes(x = cwd_hist)) +
  geom_histogram(bins = 50) +
  xlim(-3, 5) +
  theme_bw() +
  ggtitle("CWD frequency across species ranges") +
  xlab("Historic CWD (Deviation from species mean)")

quantile_df <- tibble(itrdb = itrdb_quantiles, fullrange = fullrange_quantiles)

qq_plot <- quantile_df %>%
  ggplot(aes(x = fullrange, y = itrdb)) +
  geom_point() +
  xlim(c(-3, 5)) +
  ylim(c(-3, 5)) +
  # xlim(c(-4, 15)) +
  # ylim(c(-4, 15)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  coord_fixed() +
  ggtitle("QQ-plot comparing\nCWD distributions") +
  xlab("CWD quantiles in species ranges") +
  ylab("CWD quantiles in ITRDB sites")

cwd_qqplot <- (itrdb_hist / fullrange_hist) | qq_plot


# PET
fullrange_pet <- sp_predictions %>% 
  select(pet_hist) 
fullrange_pquantiles <- (fullrange_pet$pet_hist) %>% 
  quantile(probs = seq(0, 1, 0.01))

itrdb_pet <- flm_df %>% 
  select(pet.spstd)
itrdb_pquantiles <- (itrdb_pet$pet.spstd) %>% 
  quantile(probs = seq(0, 1, 0.01))

itrdb_pet_hist <- itrdb_pet %>% 
  ggplot(aes(x = pet.spstd)) +
  geom_histogram(bins = 50) +
  xlim(-3, 5) +
  theme_bw() +
  ggtitle("PET frequency among ITRDB sites") +
  xlab("Historic PET (Deviation from species mean)")

fullrange_pet_hist <- fullrange_pet %>% 
  ggplot(aes(x = pet_hist)) +
  geom_histogram(bins = 50) +
  xlim(-3, 5) +
  theme_bw() +
  ggtitle("PET frequency across species ranges") +
  xlab("Historic PET (Deviation from species mean)")

pquantile_df <- tibble(itrdb = itrdb_pquantiles, fullrange = fullrange_pquantiles)

pet_qq_plot <- pquantile_df %>%
  ggplot(aes(x = fullrange, y = itrdb)) +
  geom_point() +
  xlim(c(-3, 5)) +
  ylim(c(-3, 5)) +
  # xlim(c(-3.5, 5)) +
  # ylim(c(-3.5, 5)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  coord_fixed() +
  ggtitle("QQ-plot comparing\nPET distributions") +
  xlab("PET quantiles in species ranges") +
  ylab("PET quantiles in ITRDB sites")

pet_qqplot <- (itrdb_pet_hist / fullrange_pet_hist) | pet_qq_plot

qqplot <- cwd_qqplot / pet_qqplot
qqplot <- qqplot &
  plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=12))
qqplot
ggsave(paste0(wdir, 'figures\\FigS2_qqplots.svg'), plot = qqplot, width = 11, height = 7, units = "in")
