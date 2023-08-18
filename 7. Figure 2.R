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
library(dplR)
select <- dplyr::select


base_text_size = 12
theme_set(
  theme_bw(base_size = base_text_size)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt

#===============================================================================
# 2) Data imports  ---------
#===============================================================================
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv")) 

# 2. Species range maps
range_file <- paste0(wdir, '1_input_processed/species_ranges/merged_ranges_dissolve.shp')
range_sf <- st_read(range_file)

# 3. Site information
site_smry <- read_csv(paste0(wdir, '1_input_processed/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)
# site_loc <- site_smry %>% 
#   select(collection_id, latitude, longitude)
# flm_df <- flm_df %>% 
#   left_join(site_loc, by = "collection_id")

# # 4. Species information
# sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
# sp_info <- sp_info %>% 
#   select(species_id, genus, gymno_angio, family)
# site_smry <- site_smry %>% 
#   left_join(sp_info, by = c("species_id"))

# 5. Prediction rasters
rwi_list <- list.files(paste0(wdir, "2_output/predictions/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))

# 6. Dendro examples - note: exporting two pipo sites in first stage script
dendro_ex <- read_csv(paste0(wdir, "1_input_processed/dendro/example_sites.csv"))

# 7. Raw dendro file for one site
rwl_path <- paste0(wdir, "1_input_processed/dendro/ca585.rwl")
rwl_dat <- read.tucson(paste0(rwl_path))



#===============================================================================
# Prep climate / prediction data  ---------
#===============================================================================
# Define species
flm_df %>% group_by(species_id) %>% tally() %>% arrange(desc(n))

spp_code <- 'pipo'


trim_df <- flm_df %>% 
  filter(species_id == spp_code)

spp_predictions <- sp_predictions %>% 
  filter(sp_code == spp_code)

# sp_fut <- (spp_predictions %>% 
#              pull(clim_future_sp))[[1]]
# names(sp_fut) <- c("cwd.fut", "pet.fut")
# 
# sp_hist <- (spp_predictions %>% 
#               pull(clim_historic_sp))[[1]]
# sp_sens  <- (spp_predictions %>% 
#                pull(sensitivity))[[1]]
# sp_rwi  <- (spp_predictions %>% 
#               pull(rwi_predictions))[[1]]
# names(sp_rwi) <- "rwi"
# 
# 
# clim_compare <- brick(c(sp_fut, sp_hist, sp_sens, sp_rwi))
# clim_compare <- clim_compare %>% 
#   as.data.frame(xy = TRUE)


#===============================================================================
# Define example sites  ---------
#===============================================================================
# Pull relevant ITRDB sites
trim_df <- trim_df %>%
  drop_na() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

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


high_fs <- trim_df %>% filter(collection_id == high_sens)
low_fs <- trim_df %>% filter(collection_id == low_sens)

high_lab <- paste0("sensitivity = ", as.character(high_val))
low_lab <- paste0("sensitivity = ", as.character(low_val))

high_color <- "#404788"
low_color <- "#efca2a"
low_color <- "#1b9e77"

#440154FF

#===============================================================================
# Step 1: data and detrending  ---------
#===============================================================================
### Map of ITRDB sites and species range

# Pull relevant range map
sp_range <- range_sf %>% 
  filter(sp_code == spp_code)
sp_bbox <- st_bbox(sp_range)

lon_lims <- c(sp_bbox$xmin - 1, sp_bbox$xmax + 1)
lat_lims <- c(sp_bbox$ymin - 1, sp_bbox$ymax + 1)

# Plot species ranges
world <- ne_coastline(scale = "medium", returnclass = "sf")
map <- ggplot(trim_df, aes(x = Longitude, y = Latitude))
  
  
map <- ggplot() +
  geom_sf(data = world) +
  geom_sf(data = sp_range, fill = 'black', alpha = .2, colour = NA) +
  geom_sf(data = trim_df, color = '#B0357B', fill = '#B0357B', alpha = .4) +
  ylab("Latitude")+
  xlab("Longitude")+
  # theme(axis.text.x=element_text(size=base_text_size - 2),
  #       axis.text.y=element_text(size = base_text_size - 2))+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_x_continuous(breaks=seq(-120,100,10)) +
  geom_point(aes(x = high_coords[[1]][1], y = high_coords[[1]][2]), color = high_color, size = 3) +
  geom_point(aes(x = low_coords[[1]][1], y = low_coords[[1]][2]), color = low_color, size = 3) +
  geom_segment(
    data = trim_df,
    x = -105, y = 50,
    xend = high_coords[[1]][1], yend = high_coords[[1]][2] + 1.3,
    lineend = "round",
    linejoin = "round",
    linewidth = 1,
    arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
    colour = high_color
  ) +
  geom_label(aes(x = -105, y = 50, label = "Site A"), fill = high_color, color = "white", size = 5, label.size = NA) +
  geom_segment(
    data = trim_df,
    # x = -low_coords[[1]][1] - 20, y = low_coords[[1]][2],
    x = -120, y = 25,
    # x = -120, y = 25,
    xend = low_coords[[1]][1], yend = low_coords[[1]][2] - 1.3,
    lineend = "round",
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
    colour = low_color
  ) +
  geom_label(aes(x = -119, y = 25, label = "Site B"), fill = low_color, color = "white", size = 5, label.size = NA)

map
# ggsave(paste0(wdir, 'figures/Methods figure/TransparentFigs/map_ex.png'), plot = map, bg= 'transparent', width = 2.25, height = 2.9)



### Illustrate detrending
rwl <- rwl_dat %>% select(PPP02B)
rwi <- rwl %>% 
  select(PPP02B) %>% 
  detrend(method = "Spline", make.plot = TRUE)
detrend_df <- tibble("RWL" = rwl$PPP02B, "RWI" = rwi$PPP02B) %>% 
  mutate(Spline = RWL / RWI) %>% 
  drop_na() %>% 
  mutate(Age = row_number(),
         Detrend_spline = 1)

rwl_plot <- ggplot(data = detrend_df) +
  geom_line(aes(x = Age, y = RWL), size = 0.5) +
  geom_line(aes(x = Age, y = Spline),linetype="dashed", linewidth = 0.25) +
  ylab("RWL (mm)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


rwi_plot <- ggplot(data = detrend_df) +
  geom_line(aes(x = Age, y = RWI), size = 0.5) +
  geom_line(aes(x = Age, y = Detrend_spline),linetype="dashed", size = 0.25) +
  ylab("RWI") +
  xlab("Age (years)")

despline_plot <- rwl_plot / rwi_plot
despline_plot

# ggsave(paste0(wdir, 'figures/Methods figure/TransparentFigs/despline.png'), plot = despline_plot, bg= 'transparent', width = 2, height = 2.9)


panel_a <- map | despline_plot
panel_a
ggsave(paste0(wdir, '3_results/figures/methods_panels/panel_A.png'), plot = panel_a, bg= 'transparent', width = 4.25, height = 2.9)

#===============================================================================
# Step 2: Climatic range  ---------
#===============================================================================
### Plot of climatic range with ITRDB sites
xmin <- -1.5
xmax <- 2
hex <- spp_predictions %>% 
  ggplot(aes(x = cwd_hist, y = pet_hist)) +
  geom_density_2d(color = 'grey') +
  # geom_density_2d(colour = "black") +
  #geom_density_2d_filled(alpha = .6) +
  xlim(xmin, xmax) +
  ylim(xmin, xmax) +
  # labs(fill = "Number of cells") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed() +
  guides(fill=F, colour=F)+
  geom_point(data = trim_df, aes(x = cwd.spstd, y = pet.spstd), colour = '#B0357B', alpha = 0.4, size = 0.75) +
  geom_point(aes(x = high_fs$cwd.spstd, y = high_fs$pet.spstd), color = high_color, size = 3) +
  geom_point(aes(x = low_fs$cwd.spstd, y = low_fs$pet.spstd), color = low_color, size = 3) +
  geom_segment(
    x = 0.5, y = -1.3,
    xend = high_fs$cwd.spstd + 0.2, yend = high_fs$pet.spstd - 0.04,
    lineend = "round",
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
    colour = high_color
  ) +
  geom_label(aes(x = 0.5, y = -1.3, label = "Site A"), fill = high_color, color = "white", size = 5, label.size = NA) +
  geom_segment(
    x = 1, y = -0.8,
    xend = low_fs$cwd.spstd + 0.2, yend = low_fs$pet.spstd - 0.04,
    lineend = "round",
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
    colour = low_color
  ) +
  geom_label(aes(x = 1, y = -0.8, label = "Site B"), fill = low_color, color = "white", size = 5, label.size = NA)
hex

ggsave(paste0(wdir, '3_results/figures/methods_panels/hex.png'), plot = hex, bg= 'transparent', width = 2.9, height = 2.9)



#===============================================================================
# Step 3: First stage for example sites  ---------
#===============================================================================
high_ex <- dendro_ex %>% 
  filter(collection_id == high_sens) %>% 
  mutate(label="high")
low_ex <- dendro_ex %>% 
  filter(collection_id == low_sens) %>% 
  mutate(label="low")

both_ex <- high_ex %>% 
  rbind(low_ex)

both_fig <- both_ex %>%
  ggplot(aes(x = cwd.an.spstd, y = rwi, color=label, fill=label)) +
  geom_point( alpha=.1) +
  scale_fill_manual(values = c(high_color, low_color)) +
  scale_color_manual(values = c(high_color, low_color)) +
  geom_smooth(method="lm") +
  #theme_bw(base_size = 25) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #plot.title = element_text(hjust = 0.5))+
  ylab("Ring width index")+
  xlab("Climatic water deficit") +
  ylim(c(0,3)) +
  xlim(c(-1.25, 0))+
  guides(fill=F, color=F)+
  geom_label(aes(x = -1.05, y = 2.5, label = "Site A"), fill = high_color, color = "white", size = 5, label.size = NA) +
  geom_label(aes(x = -0.4, y = 2.5, label = "Site B"), fill = low_color, color = "white", size = 5, label.size = NA)
# 
#   annotate("text", x = -1.05, y = 2.5, label = paste("Site A"), color = high_color, size = 12/ pt_size)+
#   annotate("text", x = -0.4, y = 2.5, label = paste("Site B"), color = low_color, size = 12/ pt_size)

both_fig
ggsave(paste0(wdir, '3_results/figures/methods_panels/both_fig.png'), plot = both_fig, bg= 'transparent', width = 3.5, height = 2.9)


#===============================================================================
# Step 4: Second stage variation in sensitivity  ---------
#===============================================================================
### Plot of cwd sensitivity against cwd.spstd
mod <- lm(estimate_cwd.an ~ cwd.spstd + I(cwd.spstd^2) + pet.spstd + I(pet.spstd^2), data = trim_df)
eff = effect("cwd.spstd", mod, partial.residuals=T)
closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, cwd = eff$x$cwd.spstd)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$cwd)] + eff$residuals)

partial_plot <- ggplot(filter(x,cwd <=1&fit<=1), aes(x = cwd, y = fit, color="black")) +
  #theme_bw(base_size = 20)+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #plot.title = element_text(hjust = 0.5))+
  geom_line(size = 1,color="black", fill="black") +
  #scale_color_gradient()+
  #geom_smooth(method="lm", color="#21908CFF", fill="#21908CFF")+
  geom_point(data = filter(xy, x<=1&y<=1), aes(x = x, y = y, color= "black"), alpha=.5) +
  scale_colour_manual(values = "black")+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "black") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  ylab("Partial effect of CWD\non CWD sensitivity")+
  guides(scale = "none", color=F) +
  geom_point(aes(x = high_fs$cwd.spstd, y = high_fs$estimate_cwd.an), color = high_color, size = 3) +
  geom_point(aes(x = low_fs$cwd.spstd, y = low_fs$estimate_cwd.an), color = low_color, size = 3) +
  geom_segment(
    x = -0.2, y = -2.4,
    xend = high_fs$cwd.spstd + 0.1, yend = high_fs$estimate_cwd.an - 0.05,
    lineend = "round",
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
    colour = high_color
  ) +
  geom_label(aes(x = -0.2, y = -2.4, label = "Site A"), fill = high_color, color = "white", size = 5, label.size = NA) +
  geom_segment(
    x = -1, y = 0.8,
    xend = low_fs$cwd.spstd - 0.1, yend = low_fs$estimate_cwd.an + 0.1,
    lineend = "round",
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
    colour = low_color
  ) +
  geom_label(aes(x = -1, y = 0.7, label = "Site B"), fill = low_color, color = "white", size = 5, label.size = NA)
  # geom_smooth(data = xy, aes(x = trans(x), y = y), 
#             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)
partial_plot

ggsave(paste0(wdir, '3_results/figures/methods_panels/partial_plot.png'), plot = partial_plot, bg= 'transparent', width = 3.75, height = 3.05)

# combined_plot <- (map_ex | both_fig) / (sens_niche | partial_plot)+ plot_layout(widths = c(1,4))
# combined_plot
# combined_plot <- (map_ex | high_ex / low_ex) / (sens_niche | partial_plot)
# #ggsave(paste0(wdir, "figures\\", "5_sp_example.svg"), combined_plot, width = 15, height = 12)


#===============================================================================
# Step 5: Prediction of sensitivity  ---------
#===============================================================================
spp_predictions <- spp_predictions %>% filter(abs(cwd_hist) < 2) ## TODO - FIgure out correct cut-off for predictions

### Map of CWD sensitivity
cwd_sens_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = cwd_sens)) +
  #theme_bw(base_size = 22)+
  theme(legend.position = c(.15,.18))+
  ylab("Latitude")+
  xlab("Longitude")+
  guides(fill=guide_legend("Sens."))+
  #scale_fill_viridis_c(direction = -1) +
  scale_fill_viridis(option="mako", direction = -1)+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_x_continuous(breaks=seq(-120,100,10)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x=element_text(size=base_text_size - 6),
        # axis.text.y=element_text(size = base_text_size - 6),
        legend.key.size = unit(8, "pt"),
        legend.title=element_text(size=base_text_size - 2), 
        legend.text=element_text(size=base_text_size - 4))+
  theme()
cwd_sens_map

ggsave(paste0(wdir, '3_results/figures/methods_panels/cwd_sens_map.png'), plot = cwd_sens_map, bg= 'transparent', width = 2.25, height = 2.9)


#===============================================================================
# Step 5: Prediction of RWI change  ---------
#===============================================================================
### Map of CWD change
spp_predictions <- spp_predictions %>% 
  mutate(cwd_change = cwd_cmip_end_mean - cwd_cmip_start_mean,
         pet_change = pet_cmip_end_mean - pet_cmip_start_mean)

cwd_change_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = cwd_change)) +
  #theme_bw(base_size = 22)+
  guides(fill=guide_legend("Δ CWD"))+
  theme(legend.position = c(.18,.15))+
  ylab("Latitude")+
  xlab("Longitude")+
  #scale_fill_viridis_c(direction = -1) +
  scale_fill_viridis(option="magma")+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_x_continuous(breaks=seq(-120,100,10)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x=element_text(size=base_text_size - 6),
        # axis.text.y=element_text(size = base_text_size - 6),
        legend.key.size = unit(8, "pt"),
        legend.title=element_text(size=base_text_size - 2), 
        legend.text=element_text(size=base_text_size - 4))
cwd_change_map

ggsave(paste0(wdir, '3_results/figures/methods_panels/cwd_change_map.png'), plot = cwd_change_map, bg= 'transparent', width = 2.25, height = 2.9)


### Map of predicted RWI
rwi_map <- ggplot() +
  geom_sf(data = world) +
  geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = rwi_pred_change_mean)) +
  # theme_bw(base_size = 12)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_fill_viridis_c(direction = -1) +
  #scale_fill_viridis(option="mako")+
  guides(fill=guide_legend(title="Δ RWI"))+
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_x_continuous(breaks=seq(-120,100,10)) +
  theme(legend.position = c(.18,.15),
        # axis.text.x=element_text(size=base_text_size - 6),
        # axis.text.y=element_text(size = base_text_size - 6),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(8, "pt"),
        legend.title=element_text(size=base_text_size - 2), 
        legend.text=element_text(size=base_text_size - 4))
rwi_map

ggsave(paste0(wdir, '3_results/figures/methods_panels/rwi_map.png'), plot = rwi_map, bg= 'transparent', width = 2.25, height = 2.9)

