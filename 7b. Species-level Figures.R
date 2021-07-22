#===============================================================================
# Species deep dive plots
# 
# 
#===============================================================================

#===============================================================================
# 1) Pkg imports ---------
#===============================================================================
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
library(effects)
select <- dplyr::select



#===============================================================================
# 2) Data imports  ---------
#===============================================================================
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv")) %>%
  select(-X1)

# 2. Species range maps
range_file <- paste0(wdir, 'in/species_ranges/merged_ranges.shp')
range_sf <- st_read(range_file)

# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)
site_loc <- site_smry %>% 
  select(collection_id, latitude, longitude)
flm_df <- flm_df %>% 
  left_join(site_loc, by = "collection_id")

# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))

# 5. Prediction rasters
sp_predictions <- readRDS(paste0(wdir, "out/predictions/sq_sp_predictions.rds"))

# 6. Dendro examples - note: currently just exporting two pipo sites in first stage script
dendro_ex <- read.csv(paste0(wdir, "out/dendro/example_sites.csv"))


#===============================================================================
# 3) Prep climate / prediction data  ---------
#===============================================================================
# Define species
flm_df %>% group_by(species_id) %>% tally() %>% arrange(desc(n))

spp_code <- 'pipo'


trim_df <- flm_df %>% 
  filter(species_id == spp_code,
         outlier != 1)

spp_predictions <- sp_predictions %>% 
  filter(sp_code == spp_code)

sp_fut <- (spp_predictions %>% 
             pull(clim_future_sp))[[1]]
names(sp_fut) <- c("cwd.fut", "pet.fut")

sp_hist <- (spp_predictions %>% 
              pull(clim_historic_sp))[[1]]
sp_sens  <- (spp_predictions %>% 
               pull(sensitivity))[[1]]
sp_rwi  <- (spp_predictions %>% 
              pull(rwi_predictions))[[1]]
names(sp_rwi) <- "rwi"


clim_compare <- brick(c(sp_fut, sp_hist, sp_sens, sp_rwi))
clim_compare <- clim_compare %>% 
  as.data.frame(xy = TRUE)


#===============================================================================
# 3) Data summary plots  ---------
#===============================================================================
### Map of ITRDB sites and species range
# Pull relevant ITRDB sites
trim_df <- trim_df %>%
  drop_na() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Pull relevant range map
sp_range <- range_sf %>% 
  filter(sp_code == spp_code)
sp_bbox <- st_bbox(sp_range)

lon_lims <- c(sp_bbox$xmin - 1, sp_bbox$xmax + 1)
lat_lims <- c(sp_bbox$ymin - 1, sp_bbox$ymax + 1)

# Plot species ranges
world <- ne_coastline(scale = "medium", returnclass = "sf")
map <- ggplot() +
  geom_sf(data = world) +
  geom_sf(data = sp_range, fill = '#21908CFF', alpha = .5, colour = NA) +
  geom_sf(data = trim_df, color = '#440154FF', fill = 'red', alpha = .8) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
map



### Plot of climatic range with ITRDB sites
xmin <- -2
xmax <- 2
hex <- clim_compare %>% ggplot(aes(x = cwd.spstd, y = pet.spstd)) +
  geom_density_2d(colour = "black") +
  geom_density_2d_filled(alpha = 0.5) +
  xlim(xmin, xmax) +
  ylim(xmin, xmax) +
  # labs(fill = "Number of cells") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  theme_bw(base_size = 22) +
  coord_fixed() +
  geom_point(data = trim_df, aes(x = cwd.spstd, y = pet.spstd), colour = "black", alpha = 0.5)
hex


#===============================================================================
# Plot of first stage for several sites  ---------
#===============================================================================
high_sens = "CO559"
low_sens = "CA585"

high_coords <- trim_df %>% 
  filter(collection_id == high_sens) %>% 
  pull(geometry)
low_coords <- trim_df %>% 
  filter(collection_id == low_sens) %>% 
  pull(geometry)

high_val <- trim_df %>% 
  filter(collection_id == high_sens) %>% 
  pull(estimate_cwd.an) %>% 
  round(digits = 3)

low_val <- trim_df %>% 
  filter(collection_id == low_sens) %>% 
  pull(estimate_cwd.an) %>% 
  round(digits = 3)

high_lab <- paste0("Sensitivity = ", as.character(high_val))
low_lab <- paste0("Sensitivity = ", as.character(low_val))

# ex_plots <- trim_df %>% 
#   filter(collection_id == high_sens | collection_id == low_sens)
map_ex <- ggplot() +
  geom_sf(data = world) +
  geom_sf(data = sp_range, fill = 'lightgrey', alpha = 1, colour = NA) +
  geom_sf(data = trim_df, color = 'darkblue', alpha = 1) +
  annotate("text", x = high_coords[[1]][1], y = high_coords[[1]][2], label = "A", color = "darkred", size = 12) +
  annotate("text", x = low_coords[[1]][1], y = low_coords[[1]][2], label = "B", color = "darkred", size = 12) +
  # geom_sf(data = ex_plots, color = 'red', fill = 'red', alpha = .8) +
  theme_bw(base_size = 20)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  # ylab("Latitude")+
  # xlab("Longitude")+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)


high_ex <- dendro_ex %>% 
  filter(collection_id == high_sens) %>% 
  mutate(label="high")
low_ex <- dendro_ex %>% 
  filter(collection_id == low_sens) %>% 
  mutate(label="low")

both_ex <- high_ex %>% 
  full_join(low_ex)

both_fig <- both_ex %>%
  ggplot(aes(x = cwd.an.spstd, y = rwi, color=label, fill=label)) +
  geom_point( alpha=.1) +
  scale_fill_manual(values = c("#404788", "#efca2a")) +
  scale_color_manual(values = c("#404788", "#efca2a")) +
  geom_smooth(method=lm) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ylab("Ring width index")+
  xlab("Climatic water deficit") +
  ylim(c(0,3)) +
  xlim(c(-1.25, 0))+
  guides(fill=F, color=F)
#ggtitle("Site A") +
#annotate("text", x = -0.85, y = 1.25, label = high_lab, size = 7)
both_fig

# low_ex <- dendro_ex %>% 
#   filter(collection_id == low_sens) %>% 
#   ggplot(aes(x = cwd.an.spstd, y = rwi)) +
#   geom_point(alpha=.1, color="#efca2a") +
#   geom_smooth(method=lm, color = "#efca2a", fill = "#efca2a") +
#   theme_bw(base_size = 20) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5))+
#   ylab("Ring width index")+
#   xlab("Climatic water deficit") +
#   ylim(c(0,3)) +
#   xlim(c(-1.25, 0)) 
#   #ggtitle("Site B") 
#   #annotate("text", x = -1.1, y = 1.3, label = low_lab, size = 7)
# low_ex
# 
# (map_ex | high_ex / low_ex) +
#   plot_layout(widths = c(2,1))

#===============================================================================
# 4) Observed sensitivity  ---------
#===============================================================================
### Map of CWD sensitivity
sens_map <- ggplot() +
  geom_sf(data = world) +
  geom_sf(data = sp_range, fill = 'gray', alpha = .9, colour = NA) +
  geom_sf(data = trim_df, aes(color = estimate_cwd.an)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_color_viridis(direction = -1)
sens_map


### Plot of observed CWD sensitivity on climatic range 
lgd_pos <- c(.18, .75)
sens_niche <- clim_compare %>% 
  filter(cwd_sens <= 1&cwd.spstd<=1&pet.spstd<=1) %>% 
  ggplot(aes(x = cwd.spstd, y = pet.spstd)) +
  geom_density_2d(colour = "black", alpha=.5) +
  # geom_density_2d_filled(alpha = 0.5) +
  #xlim(xmin, xmax) +
  #ylim(xmin, xmax) +
  # labs(fill = "Number of cells") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  theme_bw(base_size = 20) +
  xlim(-1.7,1)+
  ylim(-1.7,1)+
  coord_fixed() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = filter(trim_df, estimate_cwd.an <= 1&cwd.spstd<=1&pet.spstd<=1), aes(x = cwd.spstd, y = pet.spstd, 
                                                                                         color = estimate_cwd.an), size = 3) +
  scale_color_viridis(direction = -1)+
  labs(color = "Sensitivity\nto CWD") +
  theme(legend.position = c(lgd_pos),
        legend.key = element_blank(),
        legend.background = element_blank())
sens_niche



### Plot of cwd sensitivity against cwd.spstd
mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd , data = trim_df)
eff = effect("cwd.spstd", mod, partial.residuals=T)
closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, cwd = eff$x$cwd.spstd)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$cwd)] + eff$residuals)

partial_plot <- ggplot(filter(x,cwd <=1&fit<=1), aes(x = cwd, y = fit)) +
  theme_bw(base_size = 20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_line(size = 1,color="#21908CFF", fill="#21908CFF") +
  #scale_color_gradient()+
  #geom_smooth(method="lm", color="#21908CFF", fill="#21908CFF")+
  geom_point(data = filter(xy, x<=1&y<=1), aes(x = x, y = y, color=y), alpha=.5) +
  scale_colour_viridis(direction = -1)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#21908CFF") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  ylab("Partial effect of historical CWD\non CWD sensitivity")+
  guides(scale = "none", color=F)
# geom_smooth(data = xy, aes(x = trans(x), y = y), 
#             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)
partial_plot


combined_plot <- (map_ex | both_fig) / (sens_niche | partial_plot)+ plot_layout(widths = c(1,4))
combined_plot
combined_plot <- (map_ex | high_ex / low_ex) / (sens_niche | partial_plot)
ggsave(paste0(wdir, "figures\\", "sp_example.svg"), combined_plot, width = 15, height = 12)

#===============================================================================
# 5) Predictions  ---------
#===============================================================================
### Plot of cwd.spstd vs cwd.fut
clim_compare %>% 
  ggplot(aes(x = cwd.spstd, y = cwd.fut)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  xlim(c(-4,4)) +
  ylim(c(-4,4)) +
  theme_bw()

### Plot of cwd.spstd vs cwd.sens
clim_compare %>% 
  ggplot(aes(x = cwd.spstd, y = cwd_sens)) +
  geom_point() +
  xlim(c(-2,2)) +
  ylim(c(-.3,.3)) +
  theme_bw()

### Plot of cwd.spstd vs rwi
clim_compare %>% 
  ggplot(aes(x = cwd.spstd, y = rwi)) +
  geom_point(color = "dark blue", alpha = 0.5) +
  xlim(c(-2,2)) +
  ylim(c(-.5,.3)) +
  theme_bw() +
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, color = "red")


### Map of historic CWD
cwd_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = clim_compare %>% drop_na(), aes(x = x, y = y, fill = cwd.spstd)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
cwd_map

### Map of historic PET
pet_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = clim_compare %>% drop_na(), aes(x = x, y = y, fill = pet.spstd)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
pet_map

### Map of CWD change
clim_compare <- clim_compare %>% 
  mutate(cwd_change = cwd.fut - cwd.spstd,
         pet_change = pet.fut - pet.spstd)

cwd_change_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = clim_compare %>% drop_na(), aes(x = x, y = y, fill = cwd_change)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
cwd_change_map


### Map of PET change
pet_change_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = clim_compare %>% drop_na(), aes(x = x, y = y, fill = pet_change)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
pet_change_map


### Map of CWD sensitivity
cwd_sens_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = clim_compare %>% drop_na(), aes(x = x, y = y, fill = cwd_sens)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
cwd_sens_map


### Map of PET sensitivity
pet_sens_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = clim_compare %>% drop_na(), aes(x = x, y = y, fill = pet_sens)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
pet_sens_map


### Map of predicted RWI
rwi_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = clim_compare %>% drop_na(), aes(x = x, y = y, fill = rwi)) +
  theme_bw(base_size = 22)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)
rwi_map


