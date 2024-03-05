#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 8/7/23
# Purpose: Create figure S3 illustrating first-stage estimates
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
flm_df <- read_csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv"))

# 2. Site information
site_smry <- read_csv(paste0(wdir, '1_input_processed/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)
site_loc <- site_smry %>%
  select(collection_id, latitude, longitude)
flm_df <- flm_df %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

# 3. DNLM results
dnlm_results <- read_rds(paste0(wdir, "2_output/first_stage/dnlm_lagged_effects.rmd"))


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
  xlim(-1.25, 1) +
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
  xlim(-1.25, 1) +
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
  geom_ribbon(alpha = 0.2, fill = "dodgerblue4") +
  xlim(0, 10) +
  ylim(-0.2, 0.1)
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
  geom_ribbon(alpha = 0.2, fill = "dodgerblue4") +
  xlim(0, 10) +
  ylim(-0.2, 0.25)

pet_dynamic = pet_plot +
  xlab("Lag (years)") +
  ylab("Median effect of PET=1 shock on RWI") +
  ggtitle("Dynamic, marginal effect of PET on RWI")


first_stage_effects <- (cwd_est_plot / pet_est_plot) | (cwd_dynamic / pet_dynamic)
first_stage_effects <- first_stage_effects  + 
  plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=12))


first_stage_effects
ggsave(paste0(wdir, '3_results/figures/FigS3_first_stage_effects.svg'), plot = first_stage_effects,
       width = 15, height = 8, units = "in")




