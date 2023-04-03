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
wdir <- 'NewRemote/'

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
  left_join(site_loc, by = "collection_id") %>%
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define palettes ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_palette <- scale_colour_brewer(
  type = "seq",
  palette = 5,
  direction = 1,
  aesthetics = "colour"
)


# sp_predictions %>% 
#   ggplot(aes(x = cwd.spstd, y = cwd_sens)) +
#   geom_point(color = "dark blue", alpha = 0.5) +
#   xlim(c(-2,2)) +
#   # ylim(c(-.5,.3)) +
#   theme_bw() +
#   stat_smooth(method = "gam", formula = y ~ s(x), size = 1, color = "red")


# sp_predictions %>% 
#   ggplot(aes(x = cwd.spstd, y = rwi)) +
#   geom_point(color = "dark blue", alpha = 0.5) +
#   xlim(c(-2,2)) +
#   ylim(c(-.3,.3)) +
#   theme_bw() +
#   stat_smooth(method = "gam", formula = y ~ s(x), size = 1, color = "red")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ITRDB map, histogram, and climate with species ranges ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# map <- ggplot() +
#   theme_bw(base_size = 15)+
#   geom_sf(data = world, color = "darkgrey", fill ="lightgrey", alpha=.9) +
#   #geom_sf(data = sp_range, fill = 'lightblue', alpha = .9, colour = NA) +
#   geom_sf(data = flm_df, color = "#443A83FF", alpha = .2) +
#   ylab("Latitude")+
#   xlab("Longitude")
# map

world <- map_data("world")

site_loc <- site_smry %>% 
  select(collection_id, latitude, longitude)

xr <- seq(-180, 180, 45)
xlabels <- parse(text=str_c(abs(xr), "^o"))
yr <- seq(-90, 90, 35)
ylabels <- parse(text=str_c(abs(yr), "^o"))

map <-ggplot() +
  geom_map(data = world, map = world,
           aes(x=long, y=lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.1,alpha=.7) +
  theme_bw(base_size =25)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_point(data = site_loc, aes(y = latitude, x = longitude), color = "#443A83FF")+
  #geom_sf(data = flm_df, color = "#443A83FF", alpha = .2) +
  ylab("Latitude")+
  xlab("Longitude")+
  scale_x_continuous("Longitude", breaks = xr, labels = xlabels) +
  scale_y_continuous("Latitude", breaks = yr, labels = ylabels)

map

##exported dim 4x7


# map <- ggplot() +
#   geom_map(data = world, map = world,
#            aes(x=long, y=lat, map_id = region),
#            color = "black", fill = "lightgray", size = 0.1,alpha=.7) +
#   theme_bw(base_size =35)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   geom_sf(data = flm_df, color = "#443A83FF", alpha = .2) +
#   ylab("Latitude")+
#   xlab("Longitude")
# map



#ggsave(paste0(wdir, 'figures\\1a_itrdb_map.svg'), plot = map, width = 9, height = 6, units = "in")
# 
# mean(flm_df$cwd.ave)
# range(flm_df$cwd.ave)
# mean(flm_df$pet.ave)
# range(flm_df$pet.ave)
# 
# ggplot(flm_df)
# 
# histogram_conceptual <- ggplot(flm_df, aes(x = cwd.spstd)) + 
#   geom_histogram(bins = 40, alpha=0.5, fill = "#404788FF") +
#   xlim(c(-2.5, 2.5)) +
#   theme_bw(base_size = 22) + 
#   ylab("Number of sites")+
#   theme(legend.position = c(.1,.5),legend.title = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"),
#         panel.border = element_blank())+
#   xlab("")+
#   ylab("")+
#   scale_y_continuous(position = "right")
# 
# histogram_conceptual
#ggsave(paste0(wdir, 'figures\\1c_hist_conceptual.svg'), plot = histogram_conceptual, width = 9, height = 6, units = "in")


# range_sf %>% ggplot() +
#   geom_sf() +
#   ylab("Latitude") +
#   xlab("Longitude")
# 
# cwd_clip <- cwd_historic %>% 
#   raster::mask(range_sf)
# 
# plot(cwd_clip)


### Generate plot illustrating range maps against historic CWD
sp_codes <- list("pipo", "tsca", "tadi", "qust", "qual")
sp_range <- range_sf %>% 
  filter(sp_code %in% sp_codes)
sp_bbox <- st_bbox(sp_range)
lon_lims <- c(sp_bbox$xmin - 1, sp_bbox$xmax + 1)
lat_lims <- c(sp_bbox$ymin - 1, sp_bbox$ymax + 1)

cwd_historic_df <- as.data.frame(cwd_historic, xy = TRUE)

cwd_historic_df <- cwd_historic_df %>% 
  filter(x >= lon_lims[1],
         x <= lon_lims[2],
         y >= lat_lims[1],
         y <= lat_lims[2])
#world <- ne_coastline(scale = "medium", returnclass = "sf")

range_map <- ggplot() +
  #geom_tile(data = cwd_historic_df, aes(x = x, y = y, fill = layer)) +
  #scale_fill_viridis_c(name = bquote('Historic CWD (mmH2O)')) +
  #geom_sf(data = world, aes(fill ="lightgrey"), alpha=.9) +
  #new_scale_fill() +
  #geom_sf(data = sp_range, aes(colour = sp_code), fill = NA) +
  geom_map(data = world, map = world,
           aes(x=long, y=lat, map_id = region),
           color = "black", fill = "lightgray",alpha=.3, size = 0.1) +
  geom_sf(data = sp_range, aes(color = sp_code, fill = sp_code), alpha = .4) +
  scale_fill_viridis_d(name = bquote('Species')) +
  scale_color_viridis_d(name = bquote('Species')) +
  #scale_colour_discrete(name = "Species") +
  #scale_fill_discrete(name = "Species") +
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)+
  theme_bw(base_size = 20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.90, 0.27),
        legend.background=element_blank())


range_map

map/range_map + plot_layout(heights=c(1.3,2))

mapsplot <- map/range_map + plot_layout(heights=c(1.3,2))

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

base_text_size = 25
### Summary plot of sample distribution
# hex <- flm_df %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs / 1000)) +
hex <- flm_df %>% ggplot(aes(x = cwd.spstd, y = pet.spstd)) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  geom_point(alpha=.0)+
  geom_hex() +
  xlim(xmin, xmax) +
  ylim(ymin, ymax) +
  labs(fill = "Number of sites") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed() +
  scale_fill_viridis_c() +
  # scale_color_viridis_c(name = bquote('Species')) +
  theme_bw(base_size = base_text_size)+
  theme(legend.position = c(.24,.83),
        legend.text = element_text(size=base_text_size - 6),
        legend.title = element_text(size=base_text_size - 4),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
hex

  
hex2 <- ggMarginal(hex, type="histogram", fill ="#404788FF", alpha=.5)
hex2 <- hex2 %>% as.ggplot()
hex2

figs1 <- map/range_map+plot_layout(heights = c(1,1))
figs1+
  plot_annotation(tag_levels = 'A') & theme(
    plot.tag = element_text(face = 'bold', size=12, family ="Helvetica"),
    text=element_text(family ="Helvetica"))

figs_maps <- mapsplot | hex2 
figs_maps

figs_maps + plot_layout(heights = c(1,2))+
  plot_annotation(tag_levels = 'A') & theme(
    plot.tag = element_text(face = 'bold', size=12, family ="Helvetica"),
    text=element_text(family ="Helvetica"))



##hypothesis figure
hyp <- as.data.frame(rbind(h1="H1: Range-edge", h0="H0: Consistent", h2="H2: Drought-naive"))
colnames(hyp)="hyp"

h1 <- desc(10:1)
h0 <- rep(-5.5,times=10)
h2 <- desc(1:10)


hypdat <- rbind(h1,h0, h2)

newdat <- cbind(hypdat, hyp) %>% 
  pivot_longer(-hyp) %>% 
  mutate(cwd=rep(1:10, times=3))

hypfig <- ggplot(newdat, aes(x=cwd,y=value, color=hyp))+
  #geom_smooth()+
  scale_color_viridis_d(name = bquote('Sensitivity')) +
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y=element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("segment", x = 0, xend = 9, y = -10, yend = -1, size=1.1, color="#404788FF",
           arrow = arrow(type="closed", length = unit(.3,"cm")))+
  annotate("segment", x = 0, xend = 9, y = -5.5, yend = -5.5,size=1.1, color="grey",
           arrow = arrow(type="closed", length = unit(.3,"cm")))+
  annotate("segment", x = 0, xend = 9, y = -1, yend = -10,size=1.1, color="#20A387FF",
           arrow = arrow(type="closed", length = unit(.3,"cm")))+
  annotate("text", x = 1.5, y = -8, angle =50, 
           size=7, label = "H2: Drought-naive", family ="Helvetica")+
  annotate("text", x = 1.7, y = -5.1, 
           size=7, label = "H0: Consistent", family ="Helvetica")+
  annotate("text", x = 2, y = -2.4, angle =-50, 
           size=7, label = "H1: Range-edge", family ="Helvetica")+
  #xlab("Climate water deficit (mm/H20)")+
  labs(y="Marginal effect of CWD on RWI",x=expression(Climate~water~deficit~(mm/H[2]*O)))

hypfig

f1 <- (map/range_map) | hex + hypfig
#ggsave(paste0(wdir, 'figures\\1_data_and_hypoth.svg'), plot = f1, width = 8, height = 5)


figs2 <- figs1|(hex + hypfig)

figs2+plot_layout(widths = c(1,2,.5))+
  plot_annotation(tag_levels = 'A') & theme(
    plot.tag = element_text(face = 'bold', size=12, family ="Helvetica"),
    text=element_text(family ="Helvetica"))
figs2

#ggsave(paste0(wdir, 'figures\\1b_obs_density.svg'), plot = hex, width = 12, height = 12)


# ## Plot to show relative variation
# plot_dat <- plot_dat %>% 
#   mutate(site_variation = rwl_sd / rwl_mean)
# group_dat <- plot_dat %>% 
#   group_by(cwd.q, pet.q) %>% 
#   dplyr::summarize(wvar = wtd.var(site_variation, na.rm = TRUE),
#                    wsd = sqrt(wvar),
#                    wmean = weighted.mean(site_variation, na.rm = TRUE),
#                    n = n(),
#                    error = qt(0.975, df = n-1)*wsd/sqrt(n),
#                    lower = wmean - error,
#                    upper = wmean + error) %>% 
#   filter(n>10)
# 
# binned_margins <- group_dat %>% 
#   ggplot(aes(x = cwd.q, y = pet.q, fill = wmean)) +
#   geom_tile() +
#   # xlim(c(-3, 4)) +
#   #ylim(c(-1.5, 1.5))+
#   # scale_fill_gradientn (colours = c("darkblue","lightblue")) +
#   scale_fill_viridis_c(direction = -1) +
#   theme_bw(base_size = 22)+
#   ylab("Deviation from mean PET")+
#   xlab("Deviation from mean CWD")+
#   theme(legend.position = "left") +
#   labs(fill = "SD / Mean RWL") +
#   scale_x_continuous(labels = cwd.quantiles[label_pattern], breaks = cwd.breaks[label_pattern]) +
#   scale_y_continuous(labels = pet.quantiles[label_pattern], breaks = pet.breaks[label_pattern]) +
#   # scale_x_continuous(labels = cwd.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   # scale_y_continuous(labels = pet.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   coord_fixed()
# binned_margins
# 
# 
# group_dat <- plot_dat %>% 
#   group_by(cwd.q, pet.q) %>% 
#   dplyr::summarize(sum_weight = sum(errorweights, na.rm = TRUE),
#                    mean_weight = mean(errorweights, na.rm = TRUE),
#                    n = n()) %>% 
#   filter(n>15)
# 
# binned_margins <- group_dat %>% 
#   ggplot(aes(x = cwd.q, y = pet.q, fill = mean_weight)) +
#   geom_tile() +
#   # xlim(c(-3, 4)) +
#   #ylim(c(-1.5, 1.5))+
#   # scale_fill_gradientn (colours = c("darkblue","lightblue")) +
#   scale_fill_viridis_c(direction = 1) +
#   theme_bw(base_size = 22)+
#   ylab("Deviation from mean PET")+
#   xlab("Deviation from mean CWD")+
#   theme(legend.position = "left") +
#   labs(fill = "Sum of weights") +
#   scale_x_continuous(labels = cwd.quantiles[label_pattern], breaks = cwd.breaks[label_pattern]) +
#   scale_y_continuous(labels = pet.quantiles[label_pattern], breaks = pet.breaks[label_pattern]) +
#   # scale_x_continuous(labels = cwd.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   # scale_y_continuous(labels = pet.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   coord_fixed()
# binned_margins
# 


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
  ggtitle("QQ-plot comparing CWD distributions") +
  xlab("Species ranges") +
  ylab("ITRDB sites")

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
  ggtitle("QQ-plot comparing PET distributions") +
  xlab("Species ranges") +
  ylab("ITRDB sites")

pet_qqplot <- (itrdb_pet_hist / fullrange_pet_hist) | pet_qq_plot

qqplot <- cwd_qqplot / pet_qqplot
qqplot

ggsave(paste0(wdir, 'figures\\a1_qqplots.svg'), plot = qqplot, width = 11, height = 7, units = "in")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Histogram of first stage coefficients --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_trim_df <- trim_df %>% 
  select(collection_id, estimate = estimate_cwd.an, p = p.value_cwd.an) %>% 
  mutate(variable = "cwd")

pet_trim_df <- trim_df %>% 
  select(collection_id, estimate = estimate_pet.an, p = p.value_pet.an) %>% 
  mutate(variable = "pet")

long_trim_df <- rbind(cwd_trim_df, pet_trim_df) %>% 
  mutate(`p<0.05` = p<0.05)

cwd_median <- long_trim_df %>% filter(variable == "cwd") %>% pull(estimate) %>% median()
pet_median <- long_trim_df %>% filter(variable == "pet") %>% pull(estimate) %>% median()

cwd_est_plot <- long_trim_df %>%
  filter(variable == "cwd") %>% 
  # filter(estimate > -0.01 & estimate<0.0001) %>% 
  ggplot(aes(x = estimate, fill = `p<0.05`)) +
  geom_histogram(bins = 200) +
  # scale_x_continuous(trans="log1p") +
  theme_bw() +
  xlim(-1.5, 1.5) +
  scale_fill_manual(values = c("steelblue2", "dodgerblue4")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  # geom_vline(xintercept = cwd_median, color = "tomato3", size = 1.5) +
  xlab("Coefficient estimate") +
  ylab("Frequency") +
  # ylim(c(0, 150)) +
  theme(legend.position = c(.9,.85),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  ggtitle("Site-level estimates of contemporaneous, marginal effect of CWD")

pet_est_plot <- long_trim_df %>%
  filter(variable == "pet") %>% 
  ggplot(aes(x = estimate, fill = `p<0.05`)) +
  geom_histogram(bins = 200) +
  # scale_x_continuous(trans="log1p") +
  theme_bw() +
  xlim(-1.5, 1.5) +
  scale_fill_manual(values = c("steelblue2", "dodgerblue4")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  # geom_vline(xintercept = pet_median, color = "tomato3", size = 1.5) +
  xlab("Coefficient estimate") +
  ylab("Frequency") +
  # ylim(c(0, 150)) +
  theme(legend.position = c(.9,.85),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  ggtitle("Site-level estimates of contemporaneous, marginal effect of PET")

cwd_est_plot
pet_est_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dynamic lag plot --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_median_effect <- dnlm_results %>%
  # filter(cwd_quantile==4) %>%
  group_by(lag) %>%
  summarise(med_effect = median(cwd_effect),
            upper = quantile(cwd_effect, 0.66, na.rm = T),
            lower = quantile(cwd_effect, 0.33, na.rm = T))
cwd_plot <- cwd_median_effect %>%
  ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
  geom_line(color = "dodgerblue4") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.2, fill = "dodgerblue4", alpha = 0.2) +
  theme_bw() +
  xlim(0, 10) +
  ylim(-0.3, 0.1)
cwd_dynamic = cwd_plot +
  xlab("Lag (years)") +
  ylab("Median effect of CWD=1 shock on RWI") +
  ggtitle("Dynamic, marginal effect of CWD on RWI")



pet_median_effect <- dnlm_results %>%
  # filter(pet_quantile==4) %>%
  group_by(lag) %>%
  summarise(med_effect = median(pet_effect),
            upper = quantile(pet_effect, 0.66, na.rm = T),
            lower = quantile(pet_effect, 0.33, na.rm = T))
pet_plot <- pet_median_effect %>%
  ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
  geom_line(color = "dodgerblue4") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.2, fill = "dodgerblue4", alpha = 0.2) +
  theme_bw() +
  xlim(0, 10) +
  ylim(-0.15, 0.3)

pet_dynamic = pet_plot +
  xlab("Lag (years)") +
  ylab("Median effect of PET=1 shock on RWI") +
  ggtitle("Dynamic, marginal effect of PET on RWI")


first_stage_effects <- (cwd_est_plot / pet_est_plot) | (cwd_dynamic / pet_dynamic)
first_stage_effects
ggsave(paste0(wdir, 'figures\\a2_first_stage_effects.svg'), plot = first_stage_effects,
width = 15, height = 8, units = "in")

# plot_dnlm <- function(crosspredictions){
#   nlags = 15
#   bylag = 0.1
#   
#   plot_df <- data.frame(est = crosspredictions$matfit[3,],
#                         ci_low = crosspredictions$matlow[3,],
#                         ci_high = crosspredictions$mathigh[3,]) %>% 
#     mutate(lag = seq(0,nlags, bylag))
#   
#   plot_df %>% 
#     ggplot(aes(x = lag)) +
#     geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "dodgerblue4", alpha = 0.2) +
#     geom_line(aes(y = est), color = "dodgerblue4") +
#     geom_hline(yintercept = 0) +
#     theme_bw() +
#     xlab("Time lag (years)") +
#     ylab("Effect on ring width index") +
#     scale_x_continuous(breaks = seq(0,10,2), limits = c(0,10))
#   }
# 
# cwd_dynamic <- plot_dnlm(dnlm_results$cwd) +
#   ggtitle("Dynamic effect of 1 s.d. increase in CWD")
#   
# pet_dynamic <- plot_dnlm(dnlm_results$pet) +
#   ggtitle("Dynamic effect of 1 s.d. increase in PET")





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main sensitivity plot --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# test_df <- mc_df %>% 
#   left_join(flm_df %>% select(collection_id, cwd.spstd, pet.spstd), by = "collection_id")
# 
# trim_df <- mc_df %>% 
#   rename(estimate_cwd.an = cwd_coef,
#          estimate_pet.an = pet_coef)

# Note: PET 1-99% quantiles vary from -1.9 to 3.5; PET from -2.9 to 1.8. -3 to 3.5 seems like a good block for plots
cwd_pred_min <- -2
cwd_pred_max <- 2
pet_pred_min <- -2
pet_pred_max <- 2

### Binned plot of cwd sensitivity
seq_inc <- 0.25
half_inc <- (seq_inc / 2)
cwd_seq_min <- cwd_pred_min - half_inc
cwd_seq_max <- cwd_pred_max + half_inc
cwd_sequence <- seq(cwd_seq_min, cwd_seq_max, seq_inc)

pet_seq_min <- pet_pred_min - half_inc
pet_seq_max <- pet_pred_max + half_inc
pet_sequence <- seq(pet_seq_min, pet_seq_max, seq_inc)

convert_bin <- function(n){
  cwd_sequence[n] + half_inc
}

plot_dat <- trim_df
# plot_dat <- trim_df %>%
#   filter(cwd.spstd > pred_min, 
#          cwd.spstd < pred_max,
#          pet.spstd > pred_min,
#          cwd.spstd < pred_max) %>% 
#   # filter(((abs(cwd.spstd)<3) & (abs(pet.spstd<3)))) %>%
#   drop_na()

plot_dat_a <- plot_dat %>%
  mutate(cwd.q = cut(cwd.spstd, breaks = cwd_sequence, labels = FALSE),
         cwd.q = convert_bin(cwd.q),
         pet.q = cut(pet.spstd, breaks = pet_sequence, labels = FALSE),
         pet.q = convert_bin(pet.q)) %>% 
  filter(pet.spstd > pet_pred_min,
         pet.spstd < pet_pred_max,
         cwd.spstd > cwd_pred_min,
         cwd.spstd < cwd_pred_max)

plot_dat_b <- plot_dat_a %>%
  group_by(cwd.q, pet.q) %>%
  summarize(cwd_sens = mean(estimate_cwd.an, na.rm = TRUE),
            # cwd_sens = weighted.mean(estimate_cwd.an, w = cwd_errorweights, na.rm = TRUE),
            pet_sens = mean(estimate_pet.an, na.rm = TRUE),
            n = n()) %>%
  filter(n>=5)



# sym_log <- function(x){
#   y = sign(x) * log10(1 + abs(x))
#   return(y)
# }
# plot_dat_b <- plot_dat_b %>% mutate(cwd_sens = sym_log(cwd_sens))

base_text_size = 18

# binned_margins <- plot_dat_b %>%
#   ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_sens)) +
#   geom_tile() +
#   # xlim(c(pred_min, pred_max)) +
#   # ylim(c(pred_min, pred_max)) +
#   # scale_fill_viridis_c(direction = -1) +
#   scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
#   ylab("Deviation from mean PET")+
#   xlab("Deviation from mean CWD")+
#   theme_bw(base_size = base_text_size)+
#   labs(fill = "Marginal effect\nof CWD") +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   theme(legend.position = c(.18,.83),
#         legend.key = element_blank(),
#         legend.background = element_blank(), 
#         legend.title=element_text(size=base_text_size - 4),
#         legend.text = element_text(size = base_text_size - 6))+
#   #panel.grid.major = element_blank(), 
#   #panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))+
#   coord_fixed() +
#   geom_hline(yintercept = 0, size = 1, linetype = 2) +
#   geom_vline(xintercept = 0, size = 1, linetype = 2) +
#   xlim(c(pred_min, pred_max)) +
#   ylim(c(pred_min, pred_max))
# # +
# #   theme(panel.grid.major = element_blank(), 
# #         panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
# 
# 

cwd_binned_margins <- plot_dat_b %>%
  ggplot(aes(x = cwd.q, y = pet.q, z = cwd_sens)) +
  stat_summary_hex(fun = function(x) mean(x), bins=12)+
  scale_fill_gradient2(low = "#401552", mid = "grey93", high = "#82bead", midpoint = .98, 
                       na.value = NA, name="Mean RWI")+
  
  # xlim(c(pred_min, pred_max)) +
  # ylim(c(pred_min, pred_max)) +
  #scale_fill_viridis_c(direction = -1) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0,
                                  limits = c(-.33, .05),
                                  oob = scales::squish) +
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme_bw(base_size = base_text_size)+
  labs(fill = "Marginal effect\nof CWD") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  theme(
        legend.key = element_blank(),
        legend.background = element_blank(), 
        legend.title=element_text(size=base_text_size - 4),
        legend.text = element_text(size = base_text_size - 6))+
  #panel.grid.major = element_blank(), 
  #panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))+
  coord_fixed() +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  xlim(c(cwd_pred_min, cwd_pred_max)) +
  ylim(c(pet_pred_min, pet_pred_max))


cwd_binned_margins

# ggsave(paste0(wdir, 'figures\\binned_margins.svg'), plot = binned_margins)


# ### Binned plot of pet sensitivity
# binned_margins <- plot_dat_b %>%
#   ggplot(aes(x = cwd.q, y = pet.q, fill = pet_sens)) +
#   geom_tile() +
#   # scale_fill_viridis_c(direction = -1) +
#   scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
#   theme_bw(base_size = 22)+
#   xlim(c(-2, 1.1))+
#   ylim(c(-2,1.1))+
#   ylab("Deviation from mean PET")+
#   xlab("Deviation from mean CWD")+
#   theme(
#     legend.position = c(.15,.85),
#     legend.key = element_blank(),
#     legend.background = element_blank()
#   ) +
#   labs(fill = "Marginal effect\nof PET") +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   coord_fixed() +
#   geom_hline(yintercept = 0, size = 1, linetype = 2) +
#   geom_vline(xintercept = 0, size = 1, linetype = 2)
# 
# 
# binned_margins



# nbins = 8
# label_gaps <- 1
# label_pattern <- seq(1, nbins + 1,label_gaps)
# plot_dat <- trim_df %>%
#   filter(((abs(cwd.spstd)<3) & (abs(pet.spstd<3)))) %>%
#   drop_na()
# 
# plot_dat <- plot_dat %>%
#   mutate(cwd.q = as.numeric(ntile(cwd.spstd, nbins)),
#          pet.q = as.numeric(ntile(pet.spstd, nbins)))
# 
# cwd.quantiles = quantile(plot_dat$cwd.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2) %>% lapply(round, digits = 1)
# pet.quantiles = quantile(plot_dat$pet.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2) %>% lapply(round, digits = 1)
# cwd.breaks = seq(0.5, nbins+0.5, 1)
# pet.breaks = seq(0.5, nbins+0.5, 1)
# 
# 
# group_dat <- plot_dat %>% 
#   group_by(cwd.q, pet.q) %>% 
#   dplyr::summarize(wvar = wtd.var(estimate_cwd.an, errorweights, na.rm = TRUE),
#                    wsd = sqrt(wvar),
#                    wmean = weighted.mean(estimate_cwd.an, errorweights, na.rm = TRUE),
#                    n = n(),
#                    error = qt(0.975, df = n-1)*wsd/sqrt(n),
#                    lower = wmean - error,
#                    upper = wmean + error) %>% 
#   filter(n>15)
# 
# 
# binned_margins <- group_dat %>% 
#   ggplot(aes(x = cwd.q, y = pet.q, fill = wmean)) +
#   geom_tile() +
#   # xlim(c(-3, 4)) +
#   #ylim(c(-1.5, 1.5))+
#   # scale_fill_gradientn (colours = c("darkblue","lightblue")) +
#   scale_fill_viridis_c(direction = -1) +
#   theme_bw(base_size = 22)+
#   ylab("Deviation from mean PET")+
#   xlab("Deviation from mean CWD")+
#   theme(legend.position = "left") +
#   labs(fill = "Marginal effect\nof CWD") +
#   scale_x_continuous(labels = cwd.quantiles[label_pattern], breaks = cwd.breaks[label_pattern]) +
#   scale_y_continuous(labels = pet.quantiles[label_pattern], breaks = pet.breaks[label_pattern]) +
#   # scale_x_continuous(labels = cwd.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   # scale_y_continuous(labels = pet.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   coord_fixed()
# 
# binned_margins
# # ggsave(paste0(wdir, 'figures\\binned_margins.svg'), plot = binned_margins)
# 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot marginal effects from ss model ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pull_marg_fx <- function(at_pet, at_cwd, mod_df){
  cwd_me_predictions <- mod_df$cwd_int + (at_cwd * mod_df$cwd_cwd) + (at_cwd^2 * mod_df$cwd_cwd2) + (at_pet * mod_df$cwd_pet) + (at_pet^2 * mod_df$cwd_pet2)
  cwd_ci_min <- cwd_me_predictions %>% quantile(0.025)
  cwd_ci_max <- cwd_me_predictions %>% quantile(0.975)
  cwd_mean <- cwd_me_predictions %>% mean()
  
  pet_me_predictions <- mod_df$pet_int + (at_cwd * mod_df$pet_cwd) + (at_cwd^2 * mod_df$pet_cwd2) + (at_pet * mod_df$pet_pet) + (at_pet^2 * mod_df$pet_pet2)
  pet_ci_min <- pet_me_predictions %>% quantile(0.025)
  pet_ci_max <- pet_me_predictions %>% quantile(0.975)
  pet_mean <- pet_me_predictions %>% mean()
  return(tibble(cwd_mean = cwd_mean, cwd_ci_min = cwd_ci_min, cwd_ci_max = cwd_ci_max,
                pet_mean = pet_mean, pet_ci_min = pet_ci_min, pet_ci_max = pet_ci_max))
}



### Marginal effects of CWD df
at_pet <- 0
init_pull_marg_fx = partial(.f = pull_marg_fx, mod_df = mod_df)
cwd_me_df <- tibble(at_cwd = seq(cwd_min, cwd_max, seq_inc))
cwd_me_df <- cwd_me_df %>%
  mutate(cwd_me = pmap(list(at_pet = at_pet,
                            at_cwd = cwd_me_df$at_cwd),
                       .f = init_pull_marg_fx)) %>% 
  unnest(cwd_me)


### Marginal effects of PET df
at_cwd <- 0
pet_me_df <- tibble(at_pet = seq(pet_min, pet_max, seq_inc))
pet_me_df <- pet_me_df %>%
  mutate(pet_me = pmap(list(at_cwd = at_cwd,
                            at_pet = pet_me_df$at_pet),
                       .f = init_pull_marg_fx)) %>% 
  unnest(pet_me)


### CWD marginal effects plots
cwd_cwd_margins_plot <- ggplot(cwd_me_df, aes(x = at_cwd)) + 
  geom_line(aes(y = cwd_mean), size = 2) +
  geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkred") +
  geom_line(aes(y = cwd_ci_max), linetype = 3) +
  geom_line(aes(y = cwd_ci_min), linetype = 3) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD") + 
  xlim(c(cwd_pred_min, cwd_pred_max)) +
  theme_bw(base_size = 18)

cwd_cwd_margins_plot


pet_cwd_margins_plot <- ggplot(pet_me_df, aes(x = at_pet)) + 
  geom_line(aes(y = cwd_mean), size = 2) +
  geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkred") +
  geom_line(aes(y = cwd_ci_max), linetype = 3) +
  geom_line(aes(y = cwd_ci_min), linetype = 3) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  xlab("Historic PET\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD") + 
  xlim(c(pet_pred_min, pet_pred_max)) +
  theme_bw(base_size = 18)

pet_cwd_margins_plot



theme_set(
  theme_bw(base_size = 45)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


cwd_full_margins <- cwd_binned_margins / cwd_cwd_margins_plot  / pet_cwd_margins_plot+ 
  # plot_layout(widths = c(1,1))+
  plot_annotation(tag_levels="A") & theme(plot.tag = element_text(face = 'bold', size=23))


cwd_full_margins
# ggsave(paste0(wdir, 'figures\\2_cwd_margins.svg'), plot = cwd_full_margins, width = 9, height = 14, units = "in")
#ggsave(paste0(wdir, 'figures\\2_cwd_margins_only.svg'), plot = margins_plot, width = 15, height = 9, units = "in")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Supplemental figure for PET ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Binned plot of pet sensitivity
pet_binned_margins <- plot_dat_b %>%
  ggplot(aes(x = cwd.q, y = pet.q, z = pet_sens)) +
  stat_summary_hex(fun = function(x) mean(x), bins=12)+
  scale_fill_gradient2(low = "#401552", mid = "grey93", high = "#82bead", midpoint = .98, 
                       na.value = NA, name="Mean RWI")+
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme_bw(base_size = base_text_size)+
  labs(fill = "Marginal effect\nof PET") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  theme(
    legend.key = element_blank(),
    legend.background = element_blank(), 
    legend.title=element_text(size=base_text_size - 4),
    legend.text = element_text(size = base_text_size - 6))+
  coord_fixed() +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  xlim(c(cwd_pred_min, cwd_pred_max)) +
  ylim(c(pet_pred_min, pet_pred_max))


pet_binned_margins


## PET marginal effects plots
pet_pet_margins_plot <- ggplot(pet_me_df, aes(x = at_pet)) + 
  geom_line(aes(y = pet_mean), size = 2) +
  geom_ribbon(aes(ymin=pet_ci_min, ymax=pet_ci_max), alpha=0.2, fill = "darkblue") +
  geom_line(aes(y = pet_ci_max), linetype = 3) +
  geom_line(aes(y = pet_ci_min), linetype = 3) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  xlab("Historic PET\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to PET") + 
  xlim(c(pet_pred_min, pet_pred_max)) +
  ylim(c(-0.4, 0.3)) +
  theme_bw(base_size = 18)

pet_pet_margins_plot

cwd_pet_margins_plot <- ggplot(cwd_me_df, aes(x = at_cwd)) + 
  geom_line(aes(y = pet_mean), size = 2) +
  geom_ribbon(aes(ymin=pet_ci_min, ymax=pet_ci_max), alpha=0.2, fill = "darkblue") +
  geom_line(aes(y = pet_ci_max), linetype = 3) +
  geom_line(aes(y = pet_ci_min), linetype = 3) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to PET") + 
  xlim(c(cwd_pred_min, cwd_pred_max)) +
  ylim(c(-0.6, 0.1)) +
  theme_bw(base_size = 18)

cwd_pet_margins_plot


pet_full_margins <- pet_binned_margins / cwd_pet_margins_plot / pet_pet_margins_plot + 
  plot_annotation(tag_levels="A") & theme(plot.tag = element_text(face = 'bold', size=23))


pet_full_margins
# ggsave(paste0(wdir, 'figures\\a3_pet_margins.svg'), plot = pet_full_margins, width = 9, height = 14, units = "in")
#ggsave(paste0(wdir, 'figures\\2_cwd_margins_only.svg'), plot = margins_plot, width = 15, height = 9, units = "in")

margins_plot <- cwd_full_margins | pet_full_margins
margins_plot +
  plot_annotation(tag_levels="A") & theme(plot.tag = element_text(face = 'bold', size=23)) &
  plot_layout(guides = "collect")

ggsave(paste0(wdir, 'figures\\2_all_margins.svg'), plot = margins_plot, width = 15, height = 14, units = "in")




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Marginal effects by genera ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Binned plot of cwd sensitivity
gen_marg_fx_df <- function(gen_mod, gen_data){
  # cwd_min = 
  cwd_inc <- 0.1
  at_pet <- 0
  cwd_min <- gen_data %>% pull(cwd.spstd) %>% min()
  cwd_max <- gen_data %>% pull(cwd.spstd) %>% max()
  gen_pred <- predictions(gen_mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = seq(cwd_min,cwd_max,cwd_inc)))
  return(gen_pred)
}


pull_coefs <- function(gen_mod, gen_dat){
  median_cwd <- gen_dat %>% pull(cwd.spstd) %>% median()
  lht <- linearHypothesis(gen_mod, c(paste0('cwd.spstd + ', as.character(median_cwd), ' * I(cwd.spstd^2) = 0')))
  pvalue <- lht$`Pr(>Chisq)`[2] %>% round(digits = 3)
  coefs = gen_mod$coefficients
  me <- (coefs['cwd.spstd'] + median_cwd * 2 * coefs['I(cwd.spstd^2)']) %>% round(digits = 3)
  
  # coef_table <- gen_mod %>% coeftable()
  # coef <- coef_table[2,1]  %>% round(digits = 3)
  # p <- coef_table[2, 4] %>% round(digits = 3)
  n <- gen_mod$nobs
  
  label <-  paste0("   n sites: ", n, ";\n   slope: ", me, ";\n   p value: ", pvalue, "\n")
  return(label)
}


# genus_info <- sp_info %>% 
#   group_by(genus) %>% 
#   summarise(gymno_angio = first(gymno_angio))

coef_labels <- genus_models %>%
  mutate(labels = map2(model_estimates, data, pull_coefs)) %>%
  select(genus, labels) %>%
  unnest(labels) %>%
  arrange(genus)


genus_predictions <- genus_models %>% 
  mutate(predictions = map2(model_estimates, data, gen_marg_fx_df))

genus_predictions <- genus_predictions %>% 
  unnest(predictions) %>% 
  select(-data, -model_estimates) %>% 
  arrange(genus)

gen_plot <- genus_predictions %>% 
  # filter(genus %in% genus_keep) %>%
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  # geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkblue") +
  theme_bw(base_size = 22) + 
  facet_wrap(~genus, scales = "free") +
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)) +
  geom_line(aes(y = conf.low), linetype = 3) +
  geom_line(aes(y = conf.high), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") +
  xlim(c(-3, 3))


gen_plot <- gen_plot +
  geom_text(data = coef_labels, aes(label = labels, x = -Inf, y = -Inf),
            hjust = 0, vjust = 0)


gen_plot
ggsave(paste0(wdir, 'figures\\a3_genus_margins_nonlinear.svg'), gen_plot, width = 22, height = 12, units = "in")


# margins_plot <- margins_plot +
#   geom_text(data = genus_coefs, aes(x = -1.8, y = 0, label = lab))
# 
# margins_plot
# 
# 
# 
# gen_pred %>% ggplot(aes(x = cwd.spstd, y = predicted))
# gen_pred %>% plot()
# plot_cme(gen_mod, effect = "cwd.spstd", condition = "cwd.spstd")
# genus_lims <- genus_predictions %>% 
#   group_by(genus) %>% 
#   summarise(cwd_min = first(min_cwd),
#             cwd_max = first(max_cwd),
#             gymno_angio = first(gymno_angio))
# genus_df <- genus_predictions %>%
#   select(-min_cwd, -max_cwd, -range_cwd) %>% 
#   group_by(genus) %>% 
#   nest() %>% 
#   left_join(genus_lims, by = "genus")
# 
# genus_df <- genus_df %>% 
#   mutate(cwd_me = pmap(list(mod_df = data, cwd_min = cwd_min, cwd_max = cwd_max),
#                        .f = gen_marg_fx_df))
# 
# genus_df <- genus_df %>% 
#   select(genus, cwd_me, gymno_angio) %>% 
#   unnest(cwd_me)
# 
# 
# 
# 
# 
# genus_keep <- genus_predictions %>% 
#   # filter(range_cwd>2) %>% 
#   filter(n_collections>25) %>%
#   pull(genus) %>% 
#   unique()
# 
# margins_plot <- genus_df %>% 
#   filter(genus %in% genus_keep) %>%
#   ggplot(aes(x = at_cwd)) + 
#   geom_line(aes(y = cwd_mean)) +
#   geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max, fill = gymno_angio), alpha=0.2) +
#   # geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkblue") +
#   theme_bw(base_size = 22) + 
#   facet_wrap(~genus, scales = "free") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA)) +
#   geom_line(aes(y = cwd_ci_min), linetype = 3) +
#   geom_line(aes(y = cwd_ci_max), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") +
#   xlim(c(-2, 2))
# 
# # margins_plot <- margins_plot +
# #   geom_text(data = genus_coefs, aes(x = -1.8, y = 0, label = lab))
# 
# margins_plot
# #ggsave(paste0(wdir, 'figures\\3_genus_margins.svg'), margins_plot, width = 10, height = 8, units = "in")


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
  ylim(c(-0.4, 0.3)) +
  xlim(c(-2, 2)) +
  # ggtitle("Historic PET = historic species mean") +
  # ylab("Predicted difference in RWI change - neutral model vs ourse") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  ylab("Difference in predicted\nRWI changes by 2100 ") +
  theme(legend.position = c(.18,.75),
        legend.text = element_text(size=13),
        legend.title = element_text(size=18),
        legend.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)
transect_0
ggsave(paste0(wdir, "figures\\", "a5_dif_pred.svg"), transect_0, width = 4, height = 4)


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
