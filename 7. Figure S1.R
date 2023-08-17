#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 8/7/23
# Purpose: Create figure S1 illustrating data context
#
# Input files:
#
# 
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
library(ggh4x)
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

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))
flm_df <- flm_df %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# 2. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)
site_loc <- site_smry %>%
  select(collection_id, latitude, longitude)


trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

# 3. Species range maps
range_file <- paste0(wdir, 'in/species_ranges/merged_ranges.shp')
range_sf <- st_read(range_file)

# 4. Average site conditions
ave_site_clim_df <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))



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
          text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ITRDB map ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
world <- map_data("world")

site_loc <- site_smry %>% 
  select(collection_id, latitude, longitude)

xr <- seq(-180, 180, 45)
xlabels <- parse(text=str_c(abs(xr), "^o"))
yr <- seq(-90, 90, 35)
ylabels <- parse(text=str_c(abs(yr), "^o"))

map <- ggplot() +
  geom_map(data = world, map = world,
           aes(map_id = region),
           color = "black", fill = "lightgray", linewidth = 0.1, alpha=.7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_blank())+
  geom_point(data = site_loc, aes(y = latitude, x = longitude), color = "#443A83FF", size = 0.5)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_x_continuous("Longitude", breaks = xr, labels = xlabels) +
  scale_y_continuous("Latitude", breaks = yr, labels = ylabels) +
# +
#   force_panelsizes(rows = unit(1.5, "in"),
#                    cols = unit(2.5, "in")) +
  coord_equal(ratio=1) # square plot to avoid the distortion

map



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Range maps ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Generate plot illustrating range maps against historic CWD
sp_codes <- list("pipo", "tsca", "tadi", "qust", "qual")
sp_codes <- list("pipo", "tsca", "tadi")
sp_range <- range_sf %>% 
  filter(sp_code %in% sp_codes)
sp_bbox <- st_bbox(sp_range)
lon_lims <- c(sp_bbox$xmin - 1, sp_bbox$xmax + 1)
lat_lims <- c(sp_bbox$ymin - 1, sp_bbox$ymax + 1)

range_map <- ggplot() +
  geom_map(data = world, map = world,
           aes(map_id = region),
           color = "black", fill = "lightgray",alpha=.3, size = 0.1) +
  geom_sf(data = sp_range, aes(color = sp_code, fill = sp_code), alpha = .4) +
  scale_fill_viridis_d(name = bquote('Species')) +
  scale_color_viridis_d(name = bquote('Species')) +
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.85, 0.27),
        legend.background=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_blank())
# +
#   force_panelsizes(rows = unit(1.5, "in"),
#                    cols = unit(2.5, "in"))



range_map

mapsplot <- map / range_map
mapsplot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Observation frequency plot --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_min <- flm_df$cwd.spstd %>% quantile(0.01) %>% print()
cwd_max <- flm_df$cwd.spstd %>% quantile(0.99) %>% print()
pet_min <- flm_df$pet.spstd %>% quantile(0.01) %>% print()
pet_max <- flm_df$pet.spstd %>% quantile(0.99) %>% print()
cwd_min_man <- -2
cwd_max_man <- 4
pet_min_man <- -4
pet_max_man <- 2

xmin <- -3
xmax <- 4
ymin <- -4
ymax = 3

### Summary plot of sample distribution
# hex <- flm_df %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs / 1000)) +
hex <- flm_df %>% ggplot(aes(x = cwd.spstd, y = pet.spstd)) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  geom_point(alpha=.0)+
  geom_hex(bins = 50) +
  xlim(xmin, xmax) +
  ylim(ymin, ymax) +
  labs(fill = "Number of sites") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed() +
  scale_fill_viridis_c() +
  theme(legend.position = c(.24,.84),
        legend.background = element_blank(),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# +
#   force_panelsizes(rows = unit(3, "in"),
#                    cols = unit(3, "in"))
hex


hex2 <- ggMarginal(hex, type="histogram", fill ="#404788FF", alpha=.5)
hex2 <- hex2 %>% as.ggplot()
hex2

### Add summary plot of sample distribution based on raw CWD/PET
hex_raw <- flm_df %>% ggplot(aes(x = cwd.ave, y = pet.ave)) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  geom_point(alpha=.0)+
  geom_hex(bins = 50) +
  xlim(-100, 2500) +
  ylim(-100, 2500) +
  labs(fill = "Number of sites") +
  ylab(bquote("Historic PET (mm"*H[2]*"O)")) +
  xlab(bquote("Historic CWD (mm"*H[2]*"O)")) +
  coord_fixed() +
  scale_fill_viridis_c() +
  # scale_color_viridis_c(name = bquote('Species')) +
  theme(legend.position = c(.24, .84),
        legend.background = element_blank(),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# +
#   force_panelsizes(rows = unit(3, "in"),
#                    cols = unit(3, "in"))
hex_raw


hex2_raw <- ggMarginal(hex_raw, type="histogram", fill ="#404788FF", alpha=.5)
hex2_raw <- hex2_raw %>% as.ggplot()
# +
#   force_panelsizes(rows = unit(3, "in"),
#                    cols = unit(3, "in"))
hex2_raw



## Combine figures
fig <- (map / range_map) | hex2_raw | hex2
fig

ggsave(paste0(wdir, 'figures/FigS1_data_summary.pdf'), plot = fig, bg= 'transparent', width = 15, height = 8, units = "in")


# fig
# 
# layout <- "
# AACCDD
# BBCCDD
# "
# 
# layout <- c(
#   area(t = 1, l = 1, b = 2, r = 2),
#   area(t = 2, l = 1, b = 3, r = 2),
#   area(t = 1, l = 2, b = 3, r = 4),
#   area(t = 1, l = 4, b = 3, r = 6)  
# )
# 
# 
# layout <- c(
#   area(t = 1, l = 1, b = 2, r = 2),
#   area(t = 3, l = 1, b = 4, r = 2),
#   area(t = 1, l = 3, b = 4, r = 7),
#   area(t = 1, l = 8, b = 4, r = 12)  
# )
# 
# 
# fig <- map + range_map + hex_raw + hex +
#   plot_layout(design = layout) &
#   plot_annotation(tag_levels="A") &
#   theme(plot.tag = element_text(face = 'bold', size=12))

# figs1 <- map/range_map+plot_layout(heights = c(1,1))
# figs1+
#   plot_annotation(tag_levels = 'A') & theme(
#     plot.tag = element_text(face = 'bold', size=12, family ="Helvetica"),
#     text=element_text(family ="Helvetica"))
# 
# figs_maps <- mapsplot | hex2_raw | hex2
# figs_maps
# 
# figs_maps + plot_layout(heights = c(1,2))+
#   plot_annotation(tag_levels = 'A') & theme(
#     plot.tag = element_text(face = 'bold', size=12, family ="Helvetica"),
#     text=element_text(family ="Helvetica"))
# 
