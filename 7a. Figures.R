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
# - add nobs
# - fix joins to prevent duplicate species_id
# - think through how to deal with CWD outliers
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



select <- dplyr::select
summarize <- dplyr::summarize

options(scipen=999)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

# 5. Prediction rasters
sp_predictions <- readRDS(paste0(wdir, "out/predictions/sp_predictions.rds"))


# 6. Second stage model
cwd_mod <- readRDS(paste0(wdir, "out\\second_stage\\cwd_mod.rds"))
cwd_vcov <- readRDS(paste0(wdir, "out\\second_stage\\cwd_mod_vcov.rds"))
mod_df <- trim_df

pet_mod <- readRDS(paste0(wdir, "out\\second_stage\\pet_mod.rds"))
int_mod <- readRDS(paste0(wdir, "out\\second_stage\\int_mod.rds"))
  
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define palettes ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_palette <- scale_colour_brewer(
  type = "seq",
  palette = 5,
  direction = 1,
  aesthetics = "colour"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stack rasters ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
create_prediction_df <- function(spp_predictions){
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
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  return(clim_compare)
}

sp_predictions <- sp_predictions %>% 
  group_by(sp_code) %>% 
  nest() %>% 
  mutate(pred_df = map(data, create_prediction_df)) %>% 
  select(-data) %>% 
  unnest(cols = pred_df) %>% 
  mutate(cwd_change = cwd.fut - cwd.spstd,
         pet_change = pet.fut - pet.spstd)

sp_predictions <- sp_predictions %>% 
  mutate(rwi_null = cwd.spstd * cwd_sens + pet.spstd * pet_sens + intercept,
         rwi_change = rwi - rwi_null)

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
world <- ne_countries(scale = "medium", returnclass = "sf")
map <- ggplot() +
  theme_bw(base_size = 22)+
  geom_sf(data = world, color = "lightgrey", fill = "lightgrey") +
  # geom_sf(data = sp_range, fill = 'lightblue', alpha = .9, colour = NA) +
  geom_sf(data = flm_df, color = 'darkblue', alpha = .5) +
  ylab("Latitude")+
  xlab("Longitude")
map
ggsave(paste0(wdir, 'figures\\1a_itrdb_map.svg'), plot = map, width = 9, height = 6, units = "in")



histogram_conceptual <- ggplot(flm_df, aes(x = cwd.spstd)) + 
  geom_histogram(bins = 40, alpha=0.8, fill = "#404788FF") +
  xlim(c(-2.5, 2.5)) +
  theme_bw(base_size = 22) + 
  ylab("Number of sites")
histogram_conceptual
ggsave(paste0(wdir, 'figures\\1c_hist_conceptual.svg'), plot = histogram_conceptual, width = 9, height = 6, units = "in")


# range_sf %>% ggplot() +
#   geom_sf() +
#   ylab("Latitude") +
#   xlab("Longitude")
# 
# cwd_clip <- cwd_historic %>% 
#   raster::mask(range_sf)
# 
# plot(cwd_clip)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Observation frequency plot --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin <- -2.5
xmax <- 2.5

### Summary plot of sample distribution
hex <- flm_df %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs)) +
  geom_hex() +
  xlim(xmin, xmax) +
  ylim(xmin, xmax) +
  labs(fill = "Number of tree-year\nobservations") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  theme_bw(base_size = 22) +
  coord_fixed() +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  theme(
    legend.position = c(.15,.9),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
hex
ggsave(paste0(wdir, 'figures\\1b_obs_density.svg'), plot = hex, width = 12, height = 12)


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
# Main sensitivity plot --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Binned plot of cwd sensitivity
seq_min <- -2.625
seq_max <- 2.625
seq_inc <- 0.25
sequence <- seq(seq_min, seq_max, seq_inc)
convert_bin <- function(n){
  sequence[n] + 0.125
}

plot_dat <- trim_df %>%
  filter(((abs(cwd.spstd)<3) & (abs(pet.spstd<3)))) %>%
  drop_na()
plot_dat <- plot_dat %>%
  mutate(cwd.q = cut(cwd.spstd, breaks = sequence, labels = FALSE),
         cwd.q = convert_bin(cwd.q),
         pet.q = cut(pet.spstd, breaks = sequence, labels = FALSE),
         pet.q = convert_bin(pet.q))


plot_dat <- plot_dat %>%
  group_by(cwd.q, pet.q) %>%
  summarize(cwd_sens = mean(estimate_cwd.an, na.rm = TRUE),
            pet_sens = mean(estimate_pet.an, na.rm = TRUE),
            n = n()) %>%
  filter(n>10)


binned_margins <- plot_dat %>%
  ggplot(aes(x = cwd.q, y = pet.q, fill = cwd_sens)) +
  geom_tile() +
  xlim(c(-2, 1.1))+
  ylim(c(-2,1.1))+
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme_bw(base_size = 22)+
  theme(
    legend.position = c(.15,.85),
    legend.key = element_blank(),
    legend.background = element_blank()
  ) +
  scale_fill_viridis_c(direction = -1) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  # theme_bw(base_size = 30)+
  # theme(
  #   legend.position = c(.25,.85),
  #   legend.key = element_blank(),
  #   legend.background = element_blank()
  # ) +
  # scale_fill_continuous_diverging(rev = TRUE, mid = 0,
  #                                 limits = c(-.5, .2),
  #                                 breaks = c(.2, 0, -.2, -.4),
  #                                 labels = c(.2, 0, -.2, -.4)) +
  labs(fill = "Marginal effect\nof CWD") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  coord_fixed() +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2)
  

binned_margins
ggsave(paste0(wdir, 'figures\\binned_margins.svg'), plot = binned_margins)


### Binned plot of pet sensitivity
binned_margins <- plot_dat %>%
  ggplot(aes(x = cwd.q, y = pet.q, fill = pet_sens)) +
  geom_tile() +
  # scale_fill_viridis_c(direction = -1) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  theme_bw(base_size = 22)+
  xlim(c(-2, 1.1))+
  ylim(c(-2,1.1))+
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme(
    legend.position = c(.15,.85),
    legend.key = element_blank(),
    legend.background = element_blank()
  ) +
  labs(fill = "Marginal effect\nof PET") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  coord_fixed() +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2)


binned_margins



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
pred_min <- mod_df$cwd.spstd %>% quantile(0.01)
pred_max <- mod_df$cwd.spstd %>% quantile(0.99)
sq_predictions <- prediction(cwd_mod, at = list(cwd.spstd = seq(pred_min, pred_max, .1)), vcov = cwd_vcov, calculate_se = T, data = mod_df) %>% 
  summary() %>% 
  rename(cwd.spstd = "at(cwd.spstd)")
margins_plot <- ggplot(sq_predictions, aes(x = cwd.spstd)) + 
  # stat_smooth(data = trim_df, aes(x = cwd.spstd, y = estimate_cwd.an)) +
  geom_line(aes(y = Prediction), size = 2) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
  geom_line(aes(y = upper), linetype = 3) +
  geom_line(aes(y = lower), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity\nto CWD") + 
  xlim(c(pred_min, pred_max)) +
  theme_bw(base_size = 22)
margins_plot


histogram <- ggplot(mod_df, aes(x = cwd.spstd)) + 
  geom_histogram(bins = 40, alpha=0.8, fill = "#404788FF") +
  xlim(c(pred_min, pred_max)) +
  theme_bw(base_size = 22) + 
  ylab("# sites") +
  theme(aspect.ratio = 0.3,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
histogram


out_fig <- binned_margins | (histogram / margins_plot)
out_fig
ggsave(paste0(wdir, 'figures\\2_cwd_margins.svg'), plot = out_fig, width = 20, height = 12, units = "in")
ggsave(paste0(wdir, 'figures\\2_cwd_margins.png'), plot = out_fig, width = 20, height = 12, units = "in")



ggsave(paste0(wdir, 'figures\\2_cwd_margins_only.svg'), plot = margins_plot, width = 15, height = 9, units = "in")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main summary figure ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_dat <- sp_predictions %>%
  filter(((abs(cwd.spstd)<2.5) & (abs(pet.spstd<2.5)))) %>%
  drop_na()

seq_min <- -2.625
seq_max <- 2.625
seq_inc <- 0.25
sequence <- seq(seq_min, seq_max, seq_inc)
convert_bin <- function(n){
  sequence[n] + 0.125
}
plot_dat <- plot_dat %>% 
  mutate(cwd.q = cut(cwd.spstd, breaks = sequence, labels = FALSE),
         cwd.q = convert_bin(cwd.q),
         pet.q = cut(pet.spstd, breaks = sequence, labels = FALSE),
         pet.q = convert_bin(pet.q))


plot_dat <- plot_dat %>% 
  group_by(cwd.q, pet.q) %>% 
  summarize(rwi = mean(rwi, na.rm = TRUE),
            rwi_change = mean(rwi_change, na.rm = TRUE),
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
  theme_bw(base_size = 20)+
  labs(fill = "Predicted\nsensitivity\nto CWD") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()
cwd_sens_bin


### PET sensitivity
pet_sens_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = pet_sens)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  # scale_fill_viridis_c(direction = -1, option = "viridis") +
  theme_bw(base_size = 20)+
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
  theme_bw(base_size = 20)+
  labs(fill = "Predicted change\nin CWD (std)") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()
cwd_change_bin

### PET change
pet_change_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = pet_change)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_viridis_c(direction = 1, option = "magma") +
  theme_bw(base_size = 20)+
  labs(fill = "Predicted change\nin PET (std)") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()
pet_change_bin


### RWI change
rwi_bin <- plot_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = rwi_change)) +
  geom_tile() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  scale_fill_viridis_c(direction = -1, option = "viridis") +
  # scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  theme_bw(base_size = 25)+
  labs(fill = "Predicted change\nin RWI") +
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
ggsave(paste0(wdir, "figures\\", "pred_full_a.svg"), sens_plot, width = 9, height = 15)


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
ggsave(paste0(wdir, "figures\\", "pred_full_c.svg"), rwi_bin, width = 9, height = 9)



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