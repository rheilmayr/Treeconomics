#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Creat predictions of growth impacts from climate change
#
# Input files:
# - ss_mod: R model object saved from Second stage
# - 
#
# ToDo:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(fixest)
library(raster)
library(sp)
library(sf)

library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Second stage model
mod <- readRDS(paste0(wdir, "out\\second_stage\\psme_mod.rds"))
mod <- readRDS(paste0(wdir, "out\\second_stage\\ss_mod.rds"))

# 2. Historic climate raster
clim_file <- paste0(wdir, 'CRU\\HistoricCWD_AETGrids.Rdat')
load(clim_file)
cwd_historic <- sum(cwd_historic)
aet_historic <- sum(aet_historic)
pet_historic <- aet_historic + cwd_historic

# 3. Climate projections from CMIP5
cmip <- load(paste0(wdir, 'CMIP5 CWD\\cmip5_cwdaet_end.Rdat'))
pet_raster <- aet_raster + cwd_raster
pet_future <- pet_raster
cwd_future <- cwd_raster
rm(pet_raster)
rm(cwd_raster)
rm(aet_raster)

# 4. Species range maps
range_file <- paste0(wdir, 'range_maps//merged_ranges.shp')
range_sf <- st_read(range_file)

# 5. Species climate niche
niche <- read.csv(paste0(wdir, "out//clim_niche.csv")) %>% 
  select(-X)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate deviation from species' historic range ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define species
# spp_code <- "pipo"
spp_code <- "psme"

sp_range <- range_sf %>%
  filter(sp_code == spp_code)

sp_niche <- niche %>%
  drop_na() %>% 
  filter(sp_code == spp_code) %>% 
  select(-sp_code) 
# %>% 
#   summarise_all(list(mean = mean, sd = sd))


sp_range <- range_sf %>% 
  filter(sp_code == spp_code)
pet_historic_sp <- pet_historic %>%  
  mask(sp_range)
pet_historic_spstd <- (pet_historic_sp - sp_niche$pet_mean) / sp_niche$pet_sd
names(pet_historic_spstd) = "pet.spstd"

cwd_historic_sp <- cwd_historic %>% 
  mask(sp_range)
cwd_historic_spstd <- (cwd_historic_sp - sp_niche$cwd_mean) / sp_niche$cwd_sd
names(cwd_historic_spstd) = "cwd.spstd"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create climate sensitivity rasters ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clim_historic_sp <- brick(list(cwd_historic_spstd, pet_historic_spstd))
cwd_sens <- raster::predict(clim_historic_sp, mod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot climate and sensitivity ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlims <- c(-132, -100)
ylims <- c(30, 56)
cwd_sensitivity_df <- cwd_sens %>% 
  as.data.frame(xy = TRUE) %>% 
  drop_na() %>% 
  filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
  rename(sensitivity = layer)

cwd_historic_df <- cwd_historic_spstd %>% 
  as.data.frame(xy = TRUE) %>% 
  filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>%
  drop_na()

pet_historic_df <- pet_historic_spstd %>% 
  as.data.frame(xy = TRUE) %>% 
  filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
  drop_na()



world <- ne_coastline(scale = "medium", returnclass = "sf")
sens_map <- ggplot() +
  geom_raster(data = cwd_sensitivity_df, aes(x = x, y = y, fill = sensitivity)) +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  scale_fill_viridis_c(direction = -1,
                       guide = guide_legend(
    label.theme = element_text(angle = 90),
    label.position = "bottom"
  )) +
  theme_bw(base_size = 22)+
  theme(legend.position = "bottom") +
  ylab("Latitude")+
  xlab("Longitude")+
  geom_sf(data = world, color = 'black') +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US

pet_map <- ggplot() +
  geom_raster(data = pet_historic_df, aes(x = x, y = y, fill = pet.spstd)) +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  scale_fill_viridis_c(guide = guide_legend(
    label.theme = element_text(angle = 90),
    label.position = "bottom"
  )) +
  theme_bw(base_size = 22)+
  theme(legend.position = "bottom") +
  ylab("Latitude")+
  xlab("Longitude")+
  geom_sf(data = world, color = 'black') +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US

cwd_map <- ggplot() +
  geom_raster(data = cwd_historic_df, aes(x = x, y = y, fill = cwd.spstd)) +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  scale_fill_viridis_c(guide = guide_legend(
    label.theme = element_text(angle = 90),
    label.position = "bottom"
  )) +
  theme_bw(base_size = 22) +
  theme(legend.position = "bottom") +
  ylab("Latitude")+
  xlab("Longitude")+
  geom_sf(data = world, color = 'black') +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US

sens_map | cwd_map



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Predict growth deviation from future climate ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pull_cmip_model <- function(mod_n){
  cwd_lyr <- cwd_future[[mod_n]]
  names(cwd_lyr) = "cwd"
  pet_lyr <- pet_future[[mod_n]]
  names(pet_lyr) = "pet"
  future_clim <- brick(list(cwd_lyr, pet_lyr))
  return(future_clim)
}

sp_clip <- function(rast){
  rast_clip <- mask(rast, sp_range)
  return(rast_clip)
}

calc_rwi <- function(cmip_rast){
  rwi_rast <- cmip_rast$cwd * cwd_sens + cmip_rast$pet * pet_sens
  return(rwi_rast)
}

n_cmip_mods <- dim(pet_future)[3]
projections <- tibble(mod_n = 1:n_cmip_mods)
projections <- projections %>% 
  mutate(cmip_rast = map(mod_n, pull_cmip_model),
         cmip_rast = map(cmip_rast, sp_clip),
         rwi_rast = map(cmip_rast, calc_rwi))



rwi_rast <- projections %>% 
  filter(mod_n == 5) %>% 
  pull(rwi_rast)

rwi_df <- rwi_rast[[1]] %>% 
  as.data.frame(xy = TRUE) %>% 
  drop_na() %>% 
  filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
  rename(log_rwi = layer)
rwi_map <- ggplot() +
  geom_raster(data = rwi_df, aes(x = x, y = y, fill = log_rwi)) +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  scale_fill_viridis_c(direction = -1,
                       guide = guide_legend(
                         label.theme = element_text(angle = 90),
                         label.position = "bottom"
                       )) +
  theme_bw(base_size = 22)+
  theme(legend.position = "bottom") +
  ylab("Latitude")+
  xlab("Longitude")+
  geom_sf(data = world, color = 'black') +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US

rwi_map


cwd_prj <- projections %>% 
  filter(mod_n==1) %>% 
  pull(cmip_rast)

cwd_change <- cwd_prj[[1]]$cwd - cwd_historic

change_df <- cwd_change %>% 
  as.data.frame(xy = TRUE) %>% 
  drop_na() %>% 
  filter(x>xlims[1], x<xlims[2], y>ylims[1], y<ylims[2]) %>% 
  rename(change = layer)

change_map <- ggplot() +
  geom_raster(data = change_df, aes(x = x, y = y, fill = change)) +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9) +
  scale_fill_viridis_c(guide = guide_legend(
                         label.theme = element_text(angle = 90),
                         label.position = "bottom"
                       )) +
  theme_bw(base_size = 22)+
  theme(legend.position = "bottom") +
  ylab("Latitude")+
  xlab("Longitude")+
  geom_sf(data = world, color = 'black') +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) ## Western US

change_map
