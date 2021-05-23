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
select <- dplyr::select



#===============================================================================
# 2) Data imports  ---------
#===============================================================================
### Define path
wdir <- 'remote\\'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv")) %>%
  select(-X1)

# 2. Species range maps
range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
range_sf <- st_read(range_file)

# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
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
  geom_sf(data = sp_range, fill = 'lightblue', alpha = .9, colour = NA) +
  geom_sf(data = trim_df, color = 'darkblue', fill = 'red', alpha = .8) +
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
sens_niche <- clim_compare %>% ggplot(aes(x = cwd.spstd, y = pet.spstd)) +
  geom_density_2d(colour = "black") +
  # geom_density_2d_filled(alpha = 0.5) +
  xlim(xmin, xmax) +
  ylim(xmin, xmax) +
  # labs(fill = "Number of cells") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  theme_bw(base_size = 22) +
  coord_fixed() +
  geom_point(data = trim_df, aes(x = cwd.spstd, y = pet.spstd, color = estimate_cwd.an)) +
  scale_color_viridis(direction = -1)
sens_niche


### Plot of cwd sensitivity against cwd.spstd
trim_df <- trim_df %>% 
  mutate(significant =   p.value_cwd.an<0.05)
trim_df %>% 
  ggplot(aes(x = cwd.spstd, y = estimate_cwd.an, color = std.error_cwd.an)) +
  scale_color_viridis() +
  geom_point() +
  ylim(c(-.5, .1)) +
  xlim(c(-2,2)) +
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, color = "black") +
  theme_bw()


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


