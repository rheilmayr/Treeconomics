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
library(prediction)
library(colorspace)
library(ggnewscale)
library(reshape2)
library(ggExtra)
library(ggsn)
library(maptools)
library(broom)
library(ggExtra)


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
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv")) %>%
  select(-X1)

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
dnlm_results <- read_rds(paste0(wdir, "out/first_stage/dnlm."))

# 5. Prediction rasters
rwi_list <- list.files(paste0(wdir, "out/predictions/sp_rwi_pred_10000/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))
# sp_predictions <- readRDS(paste0(wdir, "out/predictions/sp_predictions.rds"))

## 6. Second stage model
mod_df <- read_rds(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))

# 7. Genus second stage model
genus_predictions <- readRDS(paste0(wdir, "out/second_stage/ss_bootstrap_genus.rds"))


# 2. Species range maps
range_file <- paste0(wdir, 'in/species_ranges/merged_ranges.shp')
range_sf <- st_read(range_file)

# 8. Climate model data
cwdlist=list();aetlist=list()
j=1
for(i in c("start","mid","end")){
  load(paste0(wdir,"in/CMIP5 CWD/cmip5_cwdaet_",i,".Rdat"))
  cwdlist[[j]]=cwd_raster;aetlist[[j]]=aet_raster
  j=j+1
}

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
  theme_bw(base_size =20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_point(data = site_loc, aes(y = latitude, x = longitude), color = "#443A83FF")+
  #geom_sf(data = flm_df, color = "#443A83FF", alpha = .2) +
  ylab("Latitude")+
  xlab("Longitude")+
  scale_x_continuous("Longitude", breaks = xr, labels = xlabels) +
  scale_y_continuous("Latitude", breaks = yr, labels = ylabels)

map


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

mean(flm_df$cwd.ave)
range(flm_df$cwd.ave)
mean(flm_df$pet.ave)
range(flm_df$pet.ave)

ggplot(flm_df)

histogram_conceptual <- ggplot(flm_df, aes(x = cwd.spstd)) + 
  geom_histogram(bins = 40, alpha=0.5, fill = "#404788FF") +
  xlim(c(-2.5, 2.5)) +
  theme_bw(base_size = 22) + 
  ylab("Number of sites")+
  theme(legend.position = c(.1,.5),legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"),
        panel.border = element_blank())+
  xlab("")+
  ylab("")+
  scale_y_continuous(position = "right")

histogram_conceptual
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
sp_codes <- list("pipo", "tsca", "tadi", "pisy", "qust", "qual")
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Observation frequency plot --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin <- -2.5
xmax <- 2.5

### Summary plot of sample distribution
hex <- flm_df %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs)) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  geom_point(alpha=.0)+
  geom_hex() +
  xlim(xmin, xmax) +
  ylim(xmin, xmax) +
  labs(fill = "Number of tree-year\nobservations") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed() +
  scale_fill_viridis_c() +
  # scale_color_viridis_c(name = bquote('Species')) +
  theme_bw(base_size = 31)+
  theme(legend.position = c(.15,.83),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
hex
# hex2 <- ggMarginal(hex, type="histogram", fill ="#404788FF", alpha=.5)
# hex2

figs1 <- map/range_map
figs1+
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
  labs(y="RWI sensitivity",x=expression(Climate~water~deficit~(mm/H[2]*O)))

hypfig


map/range_map|hex2+hypfig


figs1

figs1+plot_layout(widths = c(1,1,.5))+
  plot_annotation(tag_levels = 'A') & theme(
    plot.tag = element_text(face = 'bold', size=12, family ="Helvetica"),
    text=element_text(family ="Helvetica"))
figs1

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

itrdb_cwd <- trim_df %>% 
  select(cwd.spstd)
itrdb_quantiles <- (itrdb_cwd$cwd.spstd) %>% 
  quantile(probs = seq(0, 1, 0.01))

itrdb_hist <- itrdb_cwd %>% 
  ggplot(aes(x = cwd.spstd)) +
  geom_histogram(bins = 50) +
  xlim(-2.5, 5) +
  theme_bw() +
  ggtitle("CWD frequency among ITRDB sites")

fullrange_hist <- fullrange_cwd %>% 
  ggplot(aes(x = cwd_hist)) +
  geom_histogram(bins = 50) +
  xlim(-2.5, 5) +
  theme_bw() +
  ggtitle("CWD frequency across species ranges")

quantile_df <- tibble(itrdb = itrdb_quantiles, fullrange = fullrange_quantiles)

qq_plot <- quantile_df %>%
  ggplot(aes(x = fullrange, y = itrdb)) +
  geom_point() +
  # xlim(c(-2, 2)) +
  # ylim(c(-2, 2)) +
  xlim(c(-4, 15)) +
  ylim(c(-4, 15)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  coord_fixed() +
  ggtitle("QQ-plot comparing CWD distributions")

(itrdb_hist / fullrange_hist) | qq_plot


# PET
fullrange_pet <- sp_predictions %>% 
  select(pet_hist) 
fullrange_pquantiles <- (fullrange_pet$pet_hist) %>% 
  quantile(probs = seq(0, 1, 0.01))

itrdb_pet <- trim_df %>% 
  select(pet.spstd)
itrdb_pquantiles <- (itrdb_pet$pet.spstd) %>% 
  quantile(probs = seq(0, 1, 0.01))

itrdb_pet_hist <- itrdb_pet %>% 
  ggplot(aes(x = pet.spstd)) +
  geom_histogram(bins = 50) +
  xlim(-2.5, 5) +
  theme_bw() +
  ggtitle("PET frequency among ITRDB sites")


fullrange_pet_hist <- fullrange_pet %>% 
  ggplot(aes(x = pet_hist)) +
  geom_histogram(bins = 50) +
  xlim(-2.5, 5) +
  theme_bw() +
  ggtitle("PET frequency across species ranges")

pquantile_df <- tibble(itrdb = itrdb_pquantiles, fullrange = fullrange_pquantiles)

pet_qq_plot <- pquantile_df %>%
  ggplot(aes(x = fullrange, y = itrdb)) +
  geom_point() +
  # xlim(c(-2, 2)) +
  # ylim(c(-2, 2)) +
  xlim(c(-3.5, 5)) +
  ylim(c(-3.5, 5)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  coord_fixed() +
  ggtitle("QQ-plot comparing PET distributions")

(itrdb_pet_hist / fullrange_pet_hist) | pet_qq_plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Histogram of first stage coefficients --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_trim_df <- trim_df %>% 
  select(collection_id, genus, species_id, estimate = estimate_cwd.an, p = p.value_cwd.an) %>% 
  mutate(variable = "cwd")

pet_trim_df <- trim_df %>% 
  select(collection_id, genus, species_id, estimate = estimate_pet.an, p = p.value_pet.an) %>% 
  mutate(variable = "pet")

long_trim_df <- rbind(cwd_trim_df, pet_trim_df) %>% 
  mutate(`p<0.05` = p<0.05)

cwd_median <- long_trim_df %>% filter(variable == "cwd") %>% pull(estimate) %>% median()
pet_median <- long_trim_df %>% filter(variable == "pet") %>% pull(estimate) %>% median()

cwd_est_plot <- long_trim_df %>%
  filter(variable == "cwd") %>% 
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
  ylim(c(0, 150)) +
  theme(legend.position = c(.9,.85),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  ggtitle("Site-level estimates of marginal effect of CWD")

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
  ylim(c(0, 150)) +
  theme(legend.position = c(.9,.85),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  ggtitle("Site-level estimates of marginal effect of PET")





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dynamic lag plot --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_dnlm <- function(crosspredictions){
  nlags = 15
  bylag = 0.1
  
  plot_df <- data.frame(est = crosspredictions$matfit[3,],
                        ci_low = crosspredictions$matlow[3,],
                        ci_high = crosspredictions$mathigh[3,]) %>% 
    mutate(lag = seq(0,nlags, bylag))
  
  plot_df %>% 
    ggplot(aes(x = lag)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "dodgerblue4", alpha = 0.2) +
    geom_line(aes(y = est), color = "dodgerblue4") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    xlab("Time lag (years)") +
    ylab("Effect on ring width index") +
    scale_x_continuous(breaks = seq(0,10,2), limits = c(0,10))
  }

cwd_dynamic <- plot_dnlm(dnlm_results$cwd) +
  ggtitle("Dynamic effect of 1 s.d. increase in CWD")
  
pet_dynamic <- plot_dnlm(dnlm_results$pet) +
  ggtitle("Dynamic effect of 1 s.d. increase in PET")

(cwd_est_plot / pet_est_plot) | (cwd_dynamic / pet_dynamic)


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
cwd_min <- flm_df$pet.spstd %>% quantile(0.01)
cwd_max <- flm_df$pet.spstd %>% quantile(0.99)
pred_min <- -3
pred_max <- 3


### Binned plot of cwd sensitivity
seq_inc <- 0.25
half_inc <- (seq_inc / 2)
seq_min <- pred_min - half_inc
seq_max <- pred_max + half_inc
sequence <- seq(seq_min, seq_max, seq_inc)

convert_bin <- function(n){
  sequence[n] + half_inc
}

plot_dat <- trim_df %>%
  filter(cwd.spstd > pred_min, 
         cwd.spstd < pred_max,
         pet.spstd > pred_min,
         cwd.spstd < pred_max) %>% 
  # filter(((abs(cwd.spstd)<3) & (abs(pet.spstd<3)))) %>%
  drop_na()

plot_dat_a <- plot_dat %>%
  mutate(cwd.q = cut(cwd.spstd, breaks = sequence, labels = FALSE),
         cwd.q = convert_bin(cwd.q),
         pet.q = cut(pet.spstd, breaks = sequence, labels = FALSE),
         pet.q = convert_bin(pet.q))

plot_dat_b <- plot_dat_a %>%
  group_by(cwd.q, pet.q) %>%
  summarize(cwd_sens = mean(estimate_cwd.an, na.rm = TRUE),
            pet_sens = mean(estimate_pet.an, na.rm = TRUE),
            n = n()) %>%
  filter(n>=10)

binned_margins <- plot_dat_b %>%
  ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_sens)) +
  geom_tile() +
  # xlim(c(pred_min, pred_max)) +
  # ylim(c(pred_min, pred_max)) +
  # scale_fill_viridis_c(direction = -1) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme_bw(base_size = 22)+
  theme(legend.position = c(.18,.83),
        legend.key = element_blank(),
        legend.background = element_blank())+
  #panel.grid.major = element_blank(), 
  #panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))+
  labs(fill = "Marginal effect\nof CWD") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  coord_fixed() +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2)
# +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))


binned_margins

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
  cwd_me_predictions <- mod_df$cwd_int + (at_cwd * mod_df$cwd_cwd) + (at_pet * mod_df$cwd_pet)
  cwd_ci_min <- cwd_me_predictions %>% quantile(0.025)
  cwd_ci_max <- cwd_me_predictions %>% quantile(0.975)
  cwd_mean <- cwd_me_predictions %>% mean()
  
  pet_me_predictions <- mod_df$pet_int + (at_cwd * mod_df$pet_cwd) + (at_pet * mod_df$pet_pet)
  pet_ci_min <- pet_me_predictions %>% quantile(0.025)
  pet_ci_max <- pet_me_predictions %>% quantile(0.975)
  pet_mean <- pet_me_predictions %>% mean()
  return(tibble(cwd_mean = cwd_mean, cwd_ci_min = cwd_ci_min, cwd_ci_max = cwd_ci_max,
                pet_mean = pet_mean, pet_ci_min = pet_ci_min, pet_ci_max = pet_ci_max))
}



### Binned plot of cwd sensitivity
cwd_inc <- 0.1
at_pet <- 0

init_pull_marg_fx = partial(.f = pull_marg_fx, mod_df = mod_df)
cwd_me_df <- tibble(at_cwd = seq(cwd_min, cwd_max, seq_inc))
cwd_me_df <- cwd_me_df %>%
  mutate(cwd_me = pmap(list(at_pet = at_pet,
                            at_cwd = cwd_me_df$at_cwd),
                       .f = init_pull_marg_fx)) %>% 
  unnest(cwd_me)

## Compare to observed sensitivities
plot_dat_a %>%
  group_by(cwd.q) %>%
  summarize(cwd_sens = mean(estimate_cwd.an, na.rm = TRUE),
            pet_sens = mean(estimate_pet.an, na.rm = TRUE),
            n = n())
# %>% 
#   ggplot(aes(x = cwd.q, y = pet_sens)) +
#   geom_point()



# cwd_me_predictions <- prediction(cwd_mod, at = list(cwd.spstd = seq(cwd_min, cwd_max, .1)), 
#                              vcov = cwd_vcov, calculate_se = T, data = mod_df) %>% 
#   summary() %>% 
#   rename(cwd.spstd = "at(cwd.spstd)")

margins_plot <- ggplot(cwd_me_df, aes(x = at_cwd)) + 
  # stat_smooth(data = trim_df, aes(x = cwd.spstd, y = estimate_cwd.an)) +
  geom_line(aes(y = cwd_mean), size = 2) +
  geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkblue") +
  geom_line(aes(y = cwd_ci_max), linetype = 3) +
  geom_line(aes(y = cwd_ci_min), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD") + 
  xlim(c(pred_min, pred_max)) +
  theme_bw(base_size = 23)
# +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))

margins_plot


histogram <- ggplot(trim_df, aes(x = cwd.spstd)) + 
  geom_histogram(bins = 40, alpha = 0.8, fill = "#404788FF", color="white") +
  xlim(c(pred_min, pred_max)) +
  theme_bw(base_size = 23) + 
  ylab("# sites") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  theme(aspect.ratio = 0.3,
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
histogram


out_fig <- binned_margins / margins_plot + 
  # plot_layout(widths = c(1,1))+
  plot_annotation(tag_levels="A") & theme(plot.tag = element_text(face = 'bold', size=23))


out_fig
#ggsave(paste0(wdir, 'figures\\2_cwd_margins.svg'), plot = out_fig, width = 20, height = 12, units = "in")
#ggsave(paste0(wdir, 'figures\\2_cwd_margins.png'), plot = out_fig, width = 20, height = 12, units = "in")
#ggsave(paste0(wdir, 'figures\\2_cwd_margins_only.svg'), plot = margins_plot, width = 15, height = 9, units = "in")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Marginal effects by genera ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Binned plot of cwd sensitivity
gen_marg_fx_df <- function(mod_df, cwd_min, cwd_max){
  cwd_inc <- 0.1
  at_pet <- 0
  
  init_pull_marg_fx <- partial(.f = pull_marg_fx, 
                               mod_df = mod_df)
  cwd_me_df <- tibble(at_cwd = seq(cwd_min, cwd_max, seq_inc))
  cwd_me_df <- cwd_me_df %>%
    mutate(cwd_me = pmap(list(at_pet = at_pet,
                              at_cwd = cwd_me_df$at_cwd),
                         .f = init_pull_marg_fx)) %>% 
    unnest(cwd_me)
  return(cwd_me_df)
}

genus_lims <- genus_predictions %>% 
  group_by(genus) %>% 
  summarise(cwd_min = first(min_cwd),
            cwd_max = first(max_cwd),
            gymno_angio = first(gymno_angio))
genus_df <- genus_predictions %>%
  select(-min_cwd, -max_cwd, -range_cwd) %>% 
  group_by(genus) %>% 
  nest() %>% 
  left_join(genus_lims, by = "genus")

genus_df <- genus_df %>% 
  mutate(cwd_me = pmap(list(mod_df = data, cwd_min = cwd_min, cwd_max = cwd_max),
                       .f = gen_marg_fx_df))

genus_df <- genus_df %>% 
  select(genus, cwd_me, gymno_angio) %>% 
  unnest(cwd_me)





genus_keep <- genus_predictions %>% 
  # filter(range_cwd>2) %>% 
  filter(n_collections>25) %>%
  pull(genus) %>% 
  unique()

margins_plot <- genus_df %>% 
  filter(genus %in% genus_keep) %>%
  ggplot(aes(x = at_cwd)) + 
  geom_line(aes(y = cwd_mean)) +
  geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max, fill = gymno_angio), alpha=0.2) +
  # geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkblue") +
  theme_bw(base_size = 22) + 
  facet_wrap(~genus, scales = "free") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  geom_line(aes(y = cwd_ci_min), linetype = 3) +
  geom_line(aes(y = cwd_ci_max), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") +
  xlim(c(-2, 2))

# margins_plot <- margins_plot +
#   geom_text(data = genus_coefs, aes(x = -1.8, y = 0, label = lab))

margins_plot
#ggsave(paste0(wdir, 'figures\\3_genus_margins.svg'), margins_plot, width = 10, height = 8, units = "in")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main summary figure ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## ToDo - rwi_null should probably be calculated as part of the MC analysis in prediction script?
sp_predictions <- sp_predictions %>% 
  mutate(rwi_null = cwd_hist * cwd_sens + pet_hist * pet_sens + int_sens,
         rwi_change_mean = rwi_pred_mean - rwi_null,
         rwi_change_975 = rwi_pred_975 - rwi_null,
         rwi_change_025 = rwi_pred_025 - rwi_null,
         rwi_change_pclim_mean = rwi_pclim_mean - rwi_null,
         rwi_change_pclim_975 = rwi_pclim_975 - rwi_null,
         rwi_change_pclim_025 = rwi_pclim_025 - rwi_null,
         cwd_change = cwd_fut - cwd_hist,
         pet_change = pet_fut - pet_hist)

plot_dat <- sp_predictions %>%
  filter(((abs(cwd_hist)<2.5) & (abs(pet_hist<2.5)))) %>%
  drop_na()

seq_min <- -2.625
seq_max <- 2.625
seq_inc <- 0.25
sequence <- seq(seq_min, seq_max, seq_inc)
convert_bin <- function(n){
  sequence[n] + 0.125
}
plot_dat <- plot_dat %>% 
  mutate(cwd.q = cut(cwd_hist, breaks = sequence, labels = FALSE),
         cwd.q = convert_bin(cwd.q),
         pet.q = cut(pet_hist, breaks = sequence, labels = FALSE),
         pet.q = convert_bin(pet.q))


plot_dat <- plot_dat %>% 
  group_by(cwd.q, pet.q) %>% 
  summarize(rwi_pred = mean(rwi_pred_mean, na.rm = TRUE),
            rwi_change = mean(rwi_change_mean, na.rm = TRUE),
            rwi_change_lb = mean(rwi_change_025, na.rm = TRUE),
            rwi_change_ub = mean(rwi_change_975, na.rm = TRUE),
            rwi_change_pclim = mean(rwi_change_pclim_mean, na.rm = TRUE),
            rwi_change_pclim_lb = mean(rwi_change_pclim_025),
            rwi_change_pclim_ub = mean(rwi_change_pclim_975),
            cwd_sens = mean(cwd_sens, na.rm = TRUE),
            pet_sens = mean(pet_sens, na.rm = TRUE),
            cwd_change = mean(cwd_change, na.rm = TRUE),
            pet_change = mean(pet_change, na.rm = TRUE),
            n = n()) %>% 
  filter(n>10)


### CWD sensitivity
cwd_sens_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_sens)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  # scale_fill_distiller(type = "div") +
  # scale_fill_viridis_c(direction = -1, option = "viridis") +
  theme_bw(base_size = 55)+
  theme(legend.position = c(.13,.83),
        legend.text = element_text(size=25),
        legend.title = element_text(size=29),
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
  theme_bw(base_size = 55)+
  theme(legend.position = c(.13,.83),
        legend.text = element_text(size=25),
        legend.title = element_text(size=29),
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
  theme_bw(base_size = 55)+
  theme(legend.position = c(.19,.83),
        legend.text = element_text(size=25),
        legend.title = element_text(size=29),
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
  theme_bw(base_size = 55)+
  theme(legend.position = c(.19,.83),
        legend.text = element_text(size=25),
        legend.title = element_text(size=29),
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


lgd_pos <- c(.15, .8)
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


lgd_pos <- c(.23, .8)
rwi_bin <- rwi_bin +
  theme(
    legend.position = c(lgd_pos),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
rwi_bin
#ggsave(paste0(wdir, "figures\\", "pred_full_c.svg"), rwi_bin, width = 9, height = 9)


cwd_change_bin/pet_change_bin | cwd_sens_bin/ pet_sens_bin | rwi_bin


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Transect plots ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pred_dat <- plot_dat %>% 
  select(cwd.q, pet.q, 
         rwi_change = rwi_change, 
         rwi_change_lb = rwi_change_lb, 
         rwi_change_ub = rwi_change_ub) %>% 
  mutate(scenario = "rwi_change")

pclim_dat <- plot_dat %>% 
  select(cwd.q, pet.q, 
         rwi_change = rwi_change_pclim, 
         rwi_change_lb = rwi_change_pclim_lb, 
         rwi_change_ub = rwi_change_pclim_ub) %>% 
  mutate(scenario = "rwi_change_pclim")


transect_dat <- pred_dat %>% 
  rbind(pclim_dat) %>% 
  mutate(scenario = fct_relevel(scenario, "rwi_change", "rwi_change_pclim"))

transect_1 <- transect_dat %>% 
  filter(pet.q == 1) %>% 
  ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
  geom_ribbon(aes(ymin = rwi_change_lb,
                  ymax = rwi_change_ub,
                  fill = scenario),
              alpha = 0.2) +
  geom_line(size = 2) +
  theme_bw(base_size = 20)+
  ylim(c(-2, 0.4)) +
  xlim(c(-2, 2)) +
  scale_linetype_manual(values=c("solid", "dotted", "dotted")) +
  scale_fill_manual(name = "Scenario",
                    labels = c("Full model", 
                               "Variable shift in climate,\nconstant sensitivity"), 
                    values = c("dark blue", "dark red", "dark green")) +
  scale_color_manual(name = "Scenario",
                     labels = c("Full model", 
                                "Variable shift in climate,\nconstant sensitivity"), 
                     values = c("dark blue", "dark red", "dark green")) +
  ggtitle("Historic PET = 1 std above mean") +
  ylab("Predicted change in RWI") +
  xlab("Historic CWD (Deviation from species mean)") +
  theme(legend.position = c(.18,.75),
        legend.text = element_text(size=13),
        legend.title = element_text(size=18),
        legend.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)
transect_1


transect_2 <- transect_dat %>% 
  filter(pet.q == -1) %>% 
  ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
  geom_ribbon(aes(ymin = rwi_change_lb,
                  ymax = rwi_change_ub,
                  fill = scenario),
              alpha = 0.2) +
  geom_line(size = 2) +
  theme_bw(base_size = 20)+
  ylim(c(-2, 0.2)) +
  xlim(c(-2, 2)) +
  scale_linetype_manual(values=c("solid", "dotted", "dotted"))+
  scale_fill_manual(name = "Scenario",
                    labels = c("Full model", 
                               "Variable shift in climate,\nconstant sensitivity"), 
                    values = c("dark blue", "dark red", "dark green")) +
  scale_color_manual(name = "Scenario",
                     labels = c("Full model", 
                                "Variable shift in climate,\nconstant sensitivity"), 
                     values = c("dark blue", "dark red", "dark green")) +
  ggtitle("Historic PET = 1 std below mean") +
  ylab("Predicted change in RWI") +
  xlab("Historic CWD (Deviation from species mean)") +
  theme(legend.position = c(.18,.25),
        legend.text = element_text(size=13),
        legend.title = element_text(size=18),
        legend.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)


locator <- rwi_bin + 
  theme_bw(base_size = 20)+
  theme(legend.position = c(.18,.83),
        legend.text = element_text(size=13),
        legend.title = element_text(size=18),
        legend.background = element_blank())+
  geom_hline(yintercept = 1, size = 1) + 
  geom_hline(yintercept = -1, size = 1)

locator | transect_1 / transect_2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Changes in CWD and PET  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwdmeanchange_end=mean(cwdlist[[3]])-mean(cwdlist[[1]])

#overlay species ranges to only show changes for areas with species we are considering
range_dissolve <- st_read(paste0(wdir,"in/species_ranges/merged_ranges_dissolve.shp"))

ranges=st_read(paste0(wdir,"in/species_ranges/merged_ranges.shp"))
r=raster(nrows=nrow(cwdlist[[1]]),ncols=ncol(cwdlist[[1]]));extent(r)=extent(cwdlist[[1]])
ranges_raster=rasterize(ranges,r,field=1)
test=rasterToPolygons(ranges_raster,dissolve = TRUE)
testsf=st_as_sf(test)

plot(ranges)

world <- ne_download(scale = "medium", type="land",returnclass = "sf",category="physical")
world=st_crop(world,extent(cwdmeanchange_end))
world=st_transform(world,projection(cwdmeanchange_end))

library(reshape2)
#clip raster to species range
cwdmeanchange_end=raster::mask(cwdmeanchange_end,test)
temp=as.matrix(cwdmeanchange_end)
colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
temp=melt(temp)
colnames(temp)=c("lat","long","cwd")

area_grid=raster::area(cwdmeanchange_end)

#find global spatial average of cwd and pet changes - calculate areas of grid cells
# cwdchanges=cellStats((mask(cwdlist[[3]],test)-mask(cwdlist[[1]],test))*area_grid,stat="sum")/cellStats(area_grid,stat="sum")
# min=which.min(cwdchanges);max=which.max(cwdchanges)
# mincwd=cwdlist[[3]][[min]]-cwdlist[[1]][[min]];maxcwd=cwdlist[[3]][[max]]-cwdlist[[1]][[max]]

world_raster=rasterize(world,r,field=1)
world_raster=data.frame(temp[,1:2],melt(as.matrix(world_raster))[,3])
colnames(world_raster)=c("lat","long","land")
world_raster$land[which(world_raster$land==1)]="land"

a=ggplot(temp,aes(x=long,y=lat,fill=cwd))
a=a+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
a=a+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
a=a+scale_fill_continuous(type="viridis",na.value="transparent",name="Change in CWD\n1970-2000 to 2091-2100\n(mm per year)")
a=a+geom_raster()
a=a+annotate("text",x=-150,y=-45,size=10,label="a)")
# 
# 
# mincwd=raster::mask(mincwd,test)
# temp=as.matrix(mincwd)
# colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# temp=melt(temp)
# colnames(temp)=c("lat","long","cwd")
# 
# b=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# b=b+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
# b=b+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# b=b+scale_fill_continuous(type="viridis",na.value="transparent",limits=c(-200,2700),guide=NULL)
# b=b+annotate("text",x=-142,y=-45,size=8,label="Min")+geom_raster()
# 
# maxcwd=raster::mask(maxcwd,test)
# temp=as.matrix(maxcwd)
# colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# temp=melt(temp)
# colnames(temp)=c("lat","long","cwd")
# 
# d=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# d=d+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="lightgrey")
# d=d+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# d=d+scale_fill_continuous(type="viridis",na.value="transparent",limits=c(-200,2700),guide=NULL)
# d=d+annotate("text",x=-142,y=-45,size=8,label="Max")+geom_raster()

#a/(b+d)

#get spatial averages in cwd and pet across species range for all three time points for all model runs

petlist=list()
for(i in 1:3) petlist[[i]]=cwdlist[[i]]+aetlist[[i]]

cwdchange_dist=matrix(nrow=dim(cwdlist[[1]])[3],ncol=length(cwdlist))
petchange_dist=matrix(nrow=dim(petlist[[1]])[3],ncol=length(petlist))

area_mask=mask(area_grid,test)

for(i in 1:length(cwdlist)){
  cwd_temp=cellStats(cwdlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
  pet_temp=cellStats(petlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
  cwdchange_dist[,i]=cwd_temp;petchange_dist[,i]=pet_temp
  print(i)
}

colnames(cwdchange_dist)=c("Historic","Mid-Century","End-Century")
colnames(petchange_dist)=c("Historic","Mid-Century","End-Century")
cwdchange_dist=as.data.frame(cwdchange_dist);cwdchange_dist$Model=1:nrow(cwdchange_dist)
petchange_dist=as.data.frame(petchange_dist);petchange_dist$Model=1:nrow(petchange_dist)

cwdchange_dist=cwdchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="CWD")
petchange_dist=petchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="PET")

changedist=cbind(cwdchange_dist,petchange_dist[,3])
changedist=changedist%>%pivot_longer(c("CWD","PET"),names_to="Variable",values_to="CWD_PET")
changedist$Time=ordered(changedist$Time,c("Historic","Mid-Century","End-Century"))
changedist$Variable=ordered(changedist$Variable,c("PET","CWD"))
historic=changedist%>%group_by(Variable,Time)%>%dplyr::summarise(hist_mean=quantile(CWD_PET,0.5))
changedist$CWD_PET[which(changedist$Time=="Historic")]=NA

e=ggplot(changedist,aes(x=Time,y=CWD_PET,col=Variable))+geom_boxplot(position="identity",width=0.15,inherit.aes=FALSE,aes(x=Time,y=CWD_PET,group=interaction(Variable,Time)),col="black",outlier.shape = NA,lwd=0.75)
e=e+scale_y_continuous(limits=c(400,1500),name="CWD or PET (mm per yer)")
e=e+geom_jitter(width=0.1)+theme_bw()
e=e+geom_line(data=historic,aes(x=Time,y=hist_mean,group=Variable),col="black",lty=2)
e=e+geom_point(data=historic%>%filter(Time=="Historic"),aes(y=hist_mean,col=Variable),x="Historic",size=3)
e=e+labs(x="",color="")
e=e+scale_color_manual(values=c("#972c25","#8e9169"))
e=e+scale_x_discrete(labels=c("Historic"="Historic\n1970-2000","Mid-Century"="Mid-Century\n2045-2055","End-Century"="End-Century\n2090-2100"))
e=e+annotate("text",x=0.7,y=470,size=10,label="b)")
x11()
a/e

###Supplementary Figure Showing PET map with AET and PET changes

petmeanchange_end=mean(petlist[[3]])-mean(petlist[[1]])

petmeanchange_end=raster::mask(petmeanchange_end,test)
temp=as.matrix(petmeanchange_end)
colnames(temp)=xFromCol(petmeanchange_end);rownames(temp)=yFromRow(petmeanchange_end)
temp=melt(temp)
colnames(temp)=c("lat","long","cwd")

a=ggplot(temp,aes(x=long,y=lat,fill=cwd))
a=a+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
a=a+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
a=a+scale_fill_continuous(type="viridis",na.value="transparent",name="Change in PET\n1970-2000 to 2091-2100\n(mm per year)")
a=a+geom_raster()
a=a+annotate("text",x=-150,y=-45,size=10,label="a)")

#get model averages of AET changes
aetchange_dist=matrix(nrow=dim(aetlist[[1]])[3],ncol=3)
for(i in 1:length(cwdlist)){
  aet_temp=cellStats(aetlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
  aetchange_dist[,i]=aet_temp
  print(i)
}

colnames(aetchange_dist)=c("Historic","Mid-Century","End-Century")
aetchange_dist=as.data.frame(aetchange_dist);aetchange_dist$Model=1:nrow(aetchange_dist)
aetchange_dist=aetchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="AET")

changedist=cbind(petchange_dist,aetchange_dist[,3])
changedist=changedist%>%pivot_longer(c("PET","AET"),names_to="Variable",values_to="PET_AET")
changedist$Time=ordered(changedist$Time,c("Historic","Mid-Century","End-Century"))
historic=changedist%>%group_by(Variable,Time)%>%dplyr::summarise(hist_mean=quantile(PET_AET,0.5))
changedist$PET_AET[which(changedist$Time=="Historic")]=NA

e=ggplot(changedist,aes(x=Time,y=PET_AET,col=Variable))+geom_boxplot(position="identity",width=0.15,inherit.aes=FALSE,aes(x=Time,y=PET_AET,group=interaction(Variable,Time)),col="black",outlier.shape = NA,lwd=0.75)
e=e+scale_y_continuous(limits=c(400,1500),name="PET or AET (mm per yer)")
e=e+geom_jitter(width=0.1)+theme_bw()
e=e+geom_line(data=historic,aes(x=Time,y=hist_mean,group=Variable),col="black",lty=2)
e=e+geom_point(data=historic%>%filter(Time=="Historic"),aes(y=hist_mean,col=Variable),x="Historic",size=3)
e=e+labs(x="",color="")
e=e+scale_color_manual(values=c("#579473","#972c25"))
e=e+scale_x_discrete(labels=c("Historic"="Historic\n1970-2000","Mid-Century"="Mid-Century\n2045-2055","End-Century"="End-Century\n2090-2100"))
e=e+annotate("text",x=0.7,y=1400,size=10,label="b)")

x11()
a/e

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