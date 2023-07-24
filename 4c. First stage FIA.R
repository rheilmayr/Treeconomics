#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/17/23
# Purpose: Run plot-level regressions of RWI sensitivity to annual weather variability using FIA data
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using "3b. Species_niche.R"
# - site_an_clim.gz: File detailing site-level weather history. Generated using "3b. Species_niche.R"
# - rwi_long.csv: Directory containing processed RWI data from "1b. Parse ITRDB.R"
# - species_gen_gr.csv: Annotated data about species.
# - site_summary.csv: Generated using "1b. Parse ITRDB.R"
#
# ToDo:
# - fix joins to prevent duplicate species_id
# - think through how to deal with CWD outliers
# - track down lost observations - currently dropping a lot due to NAN or failed RWI generation
# - Incorporate code from tree-level analysis script to generate DLNM plots? 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyr)
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(lfe)
library(broom)
library(purrr)
library(fixest)
library(dtplyr)
library(rstanarm)
library(broom.mixed)
library(furrr)
library(tidylog)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'


# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out/dendro/")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long_fia.csv"))
dendro_df <- dendro_df %>% 
  mutate(species_id = tolower(species_id)) 

# 2. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "out/climate/site_an_clim.gz"))
dendro_df <- dendro_df %>% 
  left_join(an_site_clim, by = c("collection_id", "year"))






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run site-level regressions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_mod <- function(site_data, outcome = "rwi", energy_var = "pet.an"){
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree_id) %>%  n_distinct()
  no_cwd_var <- (site_data %>% select(cwd.an) %>% n_distinct() == 1)
  no_pet_var <- (site_data %>% select(energy_var) %>% n_distinct() == 1)
  
  if (no_cwd_var | no_pet_var) {
    message(paste0("Site has no variation in cwd.an or ", energy_var))
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / pet data 
    tryCatch(
      expr = {
        formula <- as.formula(paste0(outcome, " ~ ", energy_var, " + cwd.an"))
        mod <- lm(formula, data = site_data)
        
        mod_sum <- summary(mod)
        mod_vcov <- vcov(mod)
        # cov <- list(int_cwd = mod_vcov[1, 2], 
        #             int_pet = mod_vcov[1, 3], 
        #             pet_cwd = mod_vcov[2, 3])
        nobs <- nobs(mod)
        mod <- tidy(mod) %>%
          mutate(term = term %>% str_replace("\\(Intercept\\)", "intercept")) %>% 
          filter(term %in% c('intercept', 'cwd.an', energy_var)) %>% 
          pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
        # mod <- mod %>% 
        #   rename_all(funs(stringr::str_replace_all(., energy_var, 'energy.an')))
        mod$cov_int_cwd = mod_vcov[c("(Intercept)"), c("cwd.an")]
        cov_var_name <- paste0("cov_int_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("(Intercept)"), c(energy_var)]
        cov_var_name <- paste0("cov_cwd_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("cwd.an"), c(energy_var)]
        mod$r2 = mod_sum$r.squared
      },
      error = function(e){ 
        message("Returned regression error")
        reg_error <<- e[1]
        failed <<- T
      }
    )    
  }
  if (failed){
    return(NA)
  }
  return(tibble(mod = list(mod), nobs = nobs, ntrees = ntrees, error = reg_error))
}

site_df <- dendro_df %>% 
  drop_na() %>% 
  rename(cwd.an = cwd.an.spstd,
         pet.an = pet.an.spstd,
         temp.an = temp.an.spstd) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  filter(nobs>10) %>% 
  nest()

fs_mod_bl <- partial(fs_mod, outcome = "rwi", energy_var = "pet.an")

site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = fs_mod_bl))

data_df <- site_df %>% 
  select(collection_id,data)

fs_df <- site_df %>% 
  select(collection_id, fs_result) %>% 
  unnest(fs_result)

fs_df <- fs_df[which(!(fs_df %>% pull(mod) %>% is.na())),]
fs_df <- fs_df %>% 
  unnest(mod)

fs_df <- fs_df %>% 
  select(-error)

fs_df %>% write_csv(paste0(wdir, 'out\\first_stage\\fia_pet_cwd_std.csv'))



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
seq_inc <- 0.1
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
            n = n())

base_text_size = 18


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
