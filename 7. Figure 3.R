#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 8/7/23
# Purpose: Create figure S1 illustrating data context
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'
# wdir <- 'G:/.shortcut-targets-by-id/10TtqG9P3BY70rcYp-WACmO38J5zBeflA/Treeconomics/Data/replication - original/'


# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "2_output/first_stage/site_pet_cwd_std_augmented.csv"))
trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

## 2. Second stage model
mod_df <- read_rds(paste0(wdir, "2_output/second_stage/ss_bootstrap.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define style  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_palette <- scale_colour_brewer(
  type = "seq",
  palette = 5,
  direction = 1,
  aesthetics = "colour"
)

base_text_size = 12

theme_set(
  theme_bw(base_size = base_text_size)+
    theme(text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main sensitivity plot --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_pred_min <- -2
cwd_pred_max <- 2
pet_pred_min <- -2
pet_pred_max <- 2

cwd_pred_min <- -6
cwd_pred_max <- 6
pet_pred_min <- -6
pet_pred_max <- 6

cwd_min <- flm_df$cwd.spstd %>% quantile(0.01, na.rm = TRUE) %>% print()
cwd_max <- flm_df$cwd.spstd %>% quantile(0.99, na.rm = TRUE) %>% print()
pet_min <- flm_df$pet.spstd %>% quantile(0.01, na.rm = TRUE) %>% print()
pet_max <- flm_df$pet.spstd %>% quantile(0.99, na.rm = TRUE) %>% print()



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


cwd_binned_margins <- plot_dat_b %>%
  ggplot(aes(x = cwd.q, y = pet.q, z = cwd_sens)) +
  stat_summary_hex(fun = function(z) mean(z), bins=12)+
  # scale_fill_continuous_diverging(rev = TRUE, mid = 0,
  #                                 limits = c(-.33, .05),
  #                                 oob = scales::squish) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0,
                                  oob = scales::squish) +
  labs(fill = "Marginal effect\nof CWD") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  theme(legend.position = c(0.18, 0.78),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.title=element_text(size= 8),
    legend.text = element_text(size = 8))+
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  coord_fixed()
#   xlim(c(cwd_pred_min, cwd_pred_max)) +
#   ylim(c(pet_pred_min, pet_pred_max))


cwd_binned_margins <- plot_dat_a %>%
  ggplot(aes(x = cwd.spstd, y = pet.spstd, z = estimate_cwd.an)) +
  stat_summary_hex(fun = function(z) mean(z), bins=10)+
  # scale_fill_continuous_diverging(rev = TRUE, mid = 0,
  #                                 limits = c(-.33, .05),
  #                                 oob = scales::squish) +
  scale_fill_continuous_diverging(rev = TRUE, mid = 0,
                                  oob = scales::squish) +
  labs(fill = "Marginal effect\nof CWD") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  theme(legend.position = c(0.18, 0.78),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(10, "pt"),
        legend.title=element_text(size= 8),
        legend.text = element_text(size = 8))+
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  geom_vline(xintercept = 0, size = 1, linetype = 2) +
  coord_fixed()
#   xlim(c(cwd_pred_min, cwd_pred_max)) +
#   ylim(c(pet_pred_min, pet_pred_max))


cwd_binned_margins

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
  xlim(c(cwd_pred_min, cwd_pred_max))

cwd_cwd_margins_plot


pet_cwd_margins_plot <- ggplot(pet_me_df, aes(x = at_pet)) + 
  geom_line(aes(y = cwd_mean), size = 2) +
  geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkred") +
  geom_line(aes(y = cwd_ci_max), linetype = 3) +
  geom_line(aes(y = cwd_ci_min), linetype = 3) +
  geom_hline(yintercept = 0, size = 1, linetype = 2) +
  xlab("Historic PET\n(Deviation from species mean)") + 
  ylab("Pred. sensitivity to CWD") + 
  xlim(c(pet_pred_min, pet_pred_max))

pet_cwd_margins_plot


cwd_full_margins <- cwd_binned_margins / cwd_cwd_margins_plot  / pet_cwd_margins_plot


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
  # scale_fill_gradient2(low = "#401552", mid = "grey93", high = "#82bead", midpoint = .98, 
  #                      na.value = NA, name="Mean RWI")+
  scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
  labs(fill = "Marginal effect\nof PET") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  coord_fixed() +
  theme(legend.position = c(0.18, 0.78),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.title=element_text(size= 8),
    legend.text = element_text(size = 8))+
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
  ylim(c(-0.4, 0.3))

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
  ylim(c(-0.6, 0.1))

cwd_pet_margins_plot


pet_full_margins <- pet_binned_margins / cwd_pet_margins_plot / pet_pet_margins_plot


pet_full_margins
# ggsave(paste0(wdir, 'figures\\a3_pet_margins.svg'), plot = pet_full_margins, width = 9, height = 14, units = "in")
#ggsave(paste0(wdir, 'figures\\2_cwd_margins_only.svg'), plot = margins_plot, width = 15, height = 9, units = "in")

margins_plot <- cwd_full_margins | pet_full_margins
margins_plot <- margins_plot +
  plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=12))
margins_plot

ggsave(paste0(wdir, '3_results/figures/Fig3_all_margins.svg'), plot = margins_plot, width = 7.5, height = 10, units = "in")

