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
# - Mechanisms: Think through this model in more detail
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(MASS)
library(tidyverse)
library(broom)
library(purrr)
library(patchwork)
library(patchwork)
library(dbplyr)
library(RSQLite)
library(lmtest)
library(sandwich)
library(prediction)
library(Hmisc)
library(hexbin)
library(ggpubr)
library(ggiraphExtra)
library(modi)
library(margins)
library(tidylog)
library(fixest)
library(biglm)

set.seed(5597)


select <- dplyr::select

mc_n <- 10000


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote\\'

# 1. Historic species niche data
niche_csv <- paste0(wdir, 'out/climate/clim_niche.csv')
niche_df <- read_csv(niche_csv) %>% 
  select(-X1)

# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'out/climate/essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) %>%
  select(-site)

# 3a. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\site_pet_cwd_std.csv')) %>%
  select(-X1)

# # 3b. Tree-level regressions
# tree_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv')) %>%
#   select(-X1)

# ggplot(flm_df %>% filter(abs(estimate_cwd.an)<0.01, cwd.min>5), aes(x = estimate_cwd.an)) + 
#   geom_histogram(bins = 100, alpha=0.4)
# dim(flm_df %>% filter(abs(estimate_cwd.an)<100))[1] / dim(flm_df)[1] 
# dim(tree_df %>% filter(abs(estimate_cwd.an)<100))[1] / dim(tree_df)[1] 


# old_flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd_old.csv')) %>%
#   select(-X1) %>%
#   select(-c(year, aet.an, cwd.an, pet.an)) %>%
#   unique()
# old_sites <- old_flm_df %>% select(collection_id) %>% unique() %>% pull()

# 4. Site information
site_df <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id)
flm_df <- flm_df %>% 
  left_join(site_df, by = "collection_id")
flm_df <- flm_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# 5. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
flm_df <- flm_df %>% 
  left_join(sp_info, by = "species_id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Characterize site-level climate in relation to species range -------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  group_by(site_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum(aet+cwd)) %>% 
  drop_na()

# Calculate site-level historic climate
hist_clim_df <- clim_df %>%
  group_by(site_id) %>%
  filter(year<1980) %>%
  summarise(aet.ave = mean(aet.an),
            cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an),
            cwd.min = min(cwd.an))

# Merge historic climate to flm data
flm_df <- flm_df %>% 
  inner_join(hist_clim_df, by = c("collection_id" = "site_id"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species historic climate -------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niche_df <- niche_df %>% 
  select(sp_code, sp_pet_mean = pet_mean, sp_pet_sd = pet_sd, sp_cwd_mean = cwd_mean, sp_cwd_sd = cwd_sd)

# Merge chronology and species climate niche data
flm_df <- flm_df %>%
  inner_join(niche_df, by = c("species_id" = "sp_code")) #note: currently losing ~7k trees (10%) because we don't have range maps

# Calculate species-niche standardized climate
flm_df <- flm_df %>%
  mutate(pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)

# # Calculate tree-species-niche standardized climate
# flm_df <- flm_df %>%
#   mutate(pet.trstd = (pet.trmean - sp_pet_mean) / sp_pet_sd,
#          cwd.trstd = (cwd.trmean - sp_cwd_mean) / sp_cwd_sd)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Trim data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add weighting based on inverse of first stage variance
flm_df <- flm_df %>% 
  mutate(errorweights = 1 / (std.error_cwd.an),
         errorweights2 = sqrt(ntrees),
         pet_errorweights = 1 / (std.error_pet.an),
         int_errorweights = 1 / (std.error_intercept))

# Identify and trim extreme outliers
cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)


flm_df <- flm_df %>% 
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (cwd.spstd<cwd_spstd_bounds[1]) |
           (cwd.spstd>cwd_spstd_bounds[2]) |
           (pet.spstd<pet_spstd_bounds[1]) |
           (pet.spstd>pet_spstd_bounds[2]))
           

flm_df %>% write.csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

# flm_df <- flm_df %>%
#   # group_by(species_id) %>%
#   mutate() %>%
#   ungroup()

trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()


  #        abs(cwd.spstd)<5,
  #        abs(pet.spstd)<5,
  #        abs(std.error_cwd.an) < 50 * abs(estimate_cwd.an),
  #        abs(std.error_cwd.an) < 50 * abs(estimate_cwd.an)) %>% 
  # drop_na()


# increment <- 0.25
# breaks <- seq(xmin, xmax, increment)
# 
# plot_dat <- flm_df %>%
#   filter(abs(estimate_cwd.an)<10) %>% 
#   drop_na()
# 
# plot_dat <- plot_dat %>%
#   mutate(cwd.q = as.numeric(cut(cwd.spstd, breaks = breaks)),
#          pet.q = as.numeric(cut(pet.spstd, breaks = breaks)))
# cwd.quantiles = breaks
# pet.quantiles = breaks
# cwd.breaks = breaks
# pet.breaks = breaks
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
#   filter(n>10)
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
#   # scale_x_continuous(labels = cwd.quantiles[label_pattern], breaks = cwd.breaks[label_pattern]) +
#   # scale_y_continuous(labels = pet.quantiles[label_pattern], breaks = pet.breaks[label_pattern]) +
#   # scale_x_continuous(labels = cwd.quantiles, breaks = breaks) +
#   # scale_y_continuous(labels = pet.quantiles, breaks = breaks) +
#   ylab("Historic PET\n(Deviation from species mean)") +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   coord_fixed()
# 
# binned_margins


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Random draws of coefs from first stage ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
draw_coefs <- function(n, cwd_est, pet_est, int_est, cwd_ste, pet_ste, int_ste, 
                       cwd_pet_cov, pet_int_cov, cwd_int_cov){
  mu <- c("cwd_coef" = cwd_est, 
          "pet_coef" = pet_est, 
          "int_coef" = int_est)
  vcov <- matrix(data = NA, nrow = 3, ncol = 3)
  vcov[1,1] <- cwd_ste^2
  vcov[2,2] <- pet_ste^2
  vcov[3,3] <- int_ste^2
  vcov[1,2] <- cwd_pet_cov
  vcov[2,1] <- cwd_pet_cov
  vcov[1,3] <- cwd_int_cov
  vcov[3,1] <- cwd_int_cov
  vcov[2,3] <- pet_int_cov
  vcov[3,2] <- pet_int_cov
  
  draw <- mvrnorm(n, mu, vcov)
  draw <- as_tibble(draw)
  draw$iter_idx <- seq(1,n)
  draw <- draw %>% select(iter_idx, cwd_coef, pet_coef, int_coef)
  return(draw)
}

run_ss <- function(data, outcome = "cwd_coef"){
  formula <- as.formula(paste(outcome, " ~ cwd.spstd + pet.spstd"))
  mod <- lm(formula, data=data)
  vcov <- vcovCL(mod, cluster = data$collection_id)
  coefs <- mod$coefficients
  draw <- mvrnorm(1, coefs, vcov) %>% 
    enframe(name = paste("var", outcome, sep = "_"), 
            value = outcome)
  return(draw)
}

## Create n random draws of first stage coefficients for each site
mc_df <- trim_df %>%
  mutate(coef_draws = pmap(list(n = mc_n, 
                                cwd_est = trim_df$estimate_cwd.an, 
                                pet_est = trim_df$estimate_pet.an,
                                int_est = trim_df$estimate_intercept, 
                                cwd_ste = trim_df$std.error_cwd.an,
                                pet_ste = trim_df$std.error_pet.an, 
                                int_ste = trim_df$std.error_intercept,
                                cwd_pet_cov = trim_df$cov_cwd_pet, 
                                cwd_int_cov = trim_df$cov_int_cwd,
                                pet_int_cov = trim_df$cov_int_pet), 
                           draw_coefs))


## Nest draws into n separate datasets
mc_df <- mc_df %>% 
  unnest(coef_draws) %>% 
  select(collection_id, iter_idx, cwd_coef, pet_coef, int_coef, cwd.spstd, 
         pet.spstd, errorweights) %>% 
  group_by(iter_idx) %>% 
  nest()


## Run second stage model for each of the n datasets
mc_df <- mc_df %>%
  mutate(ss_cwd_mod = data %>% map(run_ss, outcome = "cwd_coef"),
         ss_pet_mod = data %>% map(run_ss, outcome = "pet_coef"),
         ss_int_mod = data %>% map(run_ss, outcome = "int_coef")) %>% 
  unnest(c(ss_cwd_mod, ss_pet_mod, ss_int_mod)) %>% 
  select(iter_idx, parameter = var_cwd_coef, cwd_coef, pet_coef, int_coef)


## Save out coefficients that reflect uncertainty from both first and second stage models
saveRDS(mc_df, paste0(wdir, "out/second_stage/ss_mc_mods.rds"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Primary second stage model --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run regression and cluster s.e.

# trim_df <- trim_df %>% 
#   filter(gymno_angio=="gymno")

mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=flm_df)
cluster_vcov <- vcovCL(mod, cluster = flm_df$collection_id)
coeftest(mod, vcov = vcovCL, cluster = flm_df$collection_id)

mod_df = trim_df
mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=mod_df)
cluster_vcov <- vcovCL(mod, cluster = mod_df$species_id)
coeftest(mod, vcov = vcovCL, cluster = mod_df$species_id)
saveRDS(mod, paste0(wdir, "out\\second_stage\\cwd_mod.rds"))
saveRDS(cluster_vcov, paste0(wdir, "out\\second_stage\\cwd_mod_vcov.rds"))
# saveRDS(sq_pet_mod, paste0(wdir, "out\\second_stage\\ss_sq_pet_mod.rds"))

pet_mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd, weights = pet_errorweights, data=mod_df)
pet_cluster_vcov <- vcovCL(pet_mod, cluster = mod_df$collection_id)
coeftest(pet_mod, cluster = mod_df$collection_id)
saveRDS(pet_mod, paste0(wdir, "out\\second_stage\\pet_mod.rds"))
saveRDS(pet_cluster_vcov, paste0(wdir, "out\\second_stage\\pet_mod_vcov.rds"))


int_mod <- lm(estimate_intercept ~ cwd.spstd + pet.spstd, weights = int_errorweights, data = mod_df)
int_cluster_vcov <- vcovCL(int_mod, cluster = mod_df$collection_id)
coeftest(int_mod, cluster = mod_df$collection_id)
saveRDS(int_mod, paste0(wdir, "out\\second_stage\\int_mod.rds"))
saveRDS(int_cluster_vcov, paste0(wdir, "out\\second_stage\\int_mod_vcov.rds"))


# Squared terms
sq_cwd_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd^2) + I(pet.spstd^2), weights = errorweights, data=mod_df)
coeftest(sq_cwd_mod, vcov = vcovCL, cluster = trim_df$collection_id)
sq_cwd_cluster_vcov <- vcovCL(sq_cwd_mod, cluster = trim_df$collection_id)
saveRDS(sq_cwd_mod, paste0(wdir, "out\\second_stage\\sq_cwd_mod.rds"))
saveRDS(sq_cwd_cluster_vcov, paste0(wdir, "out\\second_stage\\sq_cwd_mod_vcov.rds"))

sq_pet_mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd + I(cwd.spstd^2) + I(pet.spstd^2), weights = pet_errorweights, data=mod_df)
sq_pet_cluster_vcov <- vcovCL(sq_pet_mod, cluster = trim_df$collection_id)
saveRDS(sq_pet_mod, paste0(wdir, "out\\second_stage\\sq_pet_mod.rds"))
saveRDS(sq_pet_cluster_vcov, paste0(wdir, "out\\second_stage\\sq_pet_mod_vcov.rds"))

sq_int_mod <- lm(estimate_intercept ~ cwd.spstd + pet.spstd + I(cwd.spstd^2) + I(pet.spstd^2), weights = int_errorweights, data=mod_df)
sq_int_cluster_vcov <- vcovCL(sq_int_mod, cluster = trim_df$collection_id)
saveRDS(sq_int_mod, paste0(wdir, "out\\second_stage\\sq_int_mod.rds"))
saveRDS(sq_int_cluster_vcov, paste0(wdir, "out\\second_stage\\sq_int_mod_vcov.rds"))


mod_df %>% 
  mutate(pred_rwi = estimate_intercept + (estimate_pet.an * pet.spstd) + (estimate_cwd.an * cwd.spstd)) %>% 
  select(pred_rwi) %>% 
  summary() #CAUTION: Shouldn't mean here be much closer to 1? Investigate in first stage script....







# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Margins plots --------------------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# trim_df <- trim_df %>% 
#   mutate(significant =   p.value_cwd.an<0.05)
# trim_df %>% 
#   ggplot(aes(x = cwd.spstd, y = estimate_cwd.an, color = std.error_cwd.an)) +
#   stat_smooth(method = "loess") +
#   theme_bw()
# 
# xmin = -2
# xmax = 2
# sq_predictions <- prediction(sq_mod, at = list(cwd.spstd = seq(xmin, xmax, .1)), vcov = sq_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(cwd.spstd = "at(cwd.spstd)")
# margins_plot <- ggplot(sq_predictions, aes(x = cwd.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto CWD") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot
# 
# 
# sq_predictions <- prediction(sq_mod, at = list(pet.spstd = seq(xmin, xmax, .1)), vcov = sq_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(pet.spstd = "at(pet.spstd)")
# margins_plot <- ggplot(sq_predictions, aes(x = pet.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic PET\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto CWD") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot
# 
# sq_pet_predictions <- prediction(sq_pet_mod, at = list(pet.spstd = seq(xmin, xmax, .1)), vcov = sq_pet_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(pet.spstd = "at(pet.spstd)")
# 
# margins_plot <- ggplot(sq_pet_predictions, aes(x = pet.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic PET\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto PET") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot
# 
# sq_pet_predictions <- prediction(sq_pet_mod, at = list(cwd.spstd = seq(xmin, xmax, .1)), vcov = sq_pet_cluster_vcov, calculate_se = T) %>%
#   summary() %>%
#   rename(cwd.spstd = "at(cwd.spstd)")
# 
# margins_plot <- ggplot(sq_pet_predictions, aes(x = cwd.spstd)) +
#   geom_line(aes(y = Prediction), size = 2) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") +
#   ylab("Predicted sensitivity\nto PET") +
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::scientific)
# margins_plot



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Robustness tests --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# No weights
flm_df <- flm_df %>% 
  mutate(equal_weights = 1)

# Trimming
flm_df <- flm_df %>% 
  mutate(trim_y = estimate_cwd.an>cwd_est_bounds[1] & estimate_cwd.an<cwd_est_bounds[2],
         trim_x = cwd.spstd>cwd_spstd_bounds[1] & cwd.spstd<cwd_spstd_bounds[2] & 
                  pet.spstd>pet_spstd_bounds[1] & pet.spstd<pet_spstd_bounds[2])

trim_y <- list(TRUE)
trim_x <- list(TRUE, FALSE)
# cluster_level <- list(list(TRUE, "species_id"), list(FALSE, "collection_id"))
cluster_level <- list(list(FALSE, "collection_id"))
weighting <- list(list(TRUE, FALSE, "errorweights"), list(FALSE, TRUE, "errorweights2"), list(FALSE, FALSE, "equalweights"))
# gymno_angio <- list(list(TRUE, FALSE, "gymno"), list(FALSE, TRUE, "angio"), list(TRUE, TRUE, c("gymno", "angio")))
gymno_angio <- list(list(TRUE, TRUE, c("gymno", "angio")))
equations <- list(list(FALSE, FALSE, "estimate_cwd.an ~ cwd.spstd + pet.spstd"), 
                  list(FALSE, TRUE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd**2) + I(pet.spstd**2)"), 
                  list(TRUE, FALSE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + factor(species_id)"),
                  list(TRUE, TRUE, "estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd**2) + I(pet.spstd**2) + factor(species_id)"))


specs <- data.frame(coef=NaN, 
                    se=NaN, 
                    trim_x = NaN, trim_y = NaN, 
                    weight_se = NaN, weight_n = NaN,
                    include_gymno = NaN, include_angio = NaN,
                    species_control = NaN, squared_term = NaN,
                    species_cluster = NaN)
i = 1
for (t_y in trim_y){
  for (t_x in trim_x){
    for(wght in weighting){
      for (g_a in gymno_angio){
        for (eqn in equations){
          for(c_lvl in cluster_level){
            data <- flm_df %>% filter(gymno_angio %in% g_a[[3]])
            if (t_x) {data <- data %>% filter(trim_x==TRUE)}
            if (t_y) {data <- data %>% filter(trim_y==TRUE)}
            mod <- lm(eqn[[3]], weights = data[[wght[[3]]]], data=data)
            cluster_vcov <- vcovCL(mod, cluster = data[[c_lvl[[2]]]])
            mod <- coeftest(mod, vcov = vcovCL, cluster = data[[c_lvl[[2]]]])
            specs[i,] <- data_frame(coef = mod[2,1],
                                    se = mod[2,2],
                                    trim_x = t_x,
                                    trim_y = t_y,
                                    weight_se = wght[[1]],
                                    weight_n = wght[[2]],
                                    include_gymno = g_a[[1]],
                                    include_angio = g_a[[2]],
                                    species_control = eqn[[1]],
                                    squared_term = eqn[[2]],
                                    species_cluster = c_lvl[[1]])
            i <- i+1
          }
        }
      }
    }
  }
}

specs <- specs %>% drop_na()
specs <- specs %>% 
  select(-include_gymno, -include_angio, -trim_y)

specs <- specs %>% 
  select(-species_cluster)

highlight_n <- which(specs$trim_x == TRUE  &
                     # specs$species_cluster == FALSE &
                     # specs$trim_y == TRUE &
                     specs$weight_se == TRUE &
                     specs$weight_n == FALSE &
                     specs$species_control == FALSE &
                     specs$squared_term == FALSE)


schart(specs, highlight=highlight_n, order = "increasing", ci = c(.95, .99))


labels <- list("Trimming:" = c("Trim outliers"),
               "Weighting:" = c("Inverse s.e.", "Square root of n"),
               "Controls" = c("Squared CWD and PET", "Species fixed effects"))

schart(specs, labels, highlight=highlight_n, order = "increasing", ci = c(.95, .99), adj=c(0,0), offset=c(5.5,5), ylab = "Slope of line relating\nsites' historic CWD to\nCWD's impact on growth")







# 
# sq_all_predictions = sq_pet_predictions %>% rename(cwd.spstd = pet.spstd) + sq_predictions
# 
# flm_df %>% filter(species_id == "tsme") %>%  ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) + geom_point() + ylim(c(-.5,.5))





# Drop PET
nopet_mod <- lm(estimate_cwd.an ~ cwd.spstd, weights = errorweights, data=trim_df)
nopet_cluster_vcov <- vcovCL(nopet_mod, cluster = trim_df$collection_id)
coeftest(nopet_mod, vcov = nopet_cluster_vcov)

# Don't trim dataset
allobs_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=flm_df)
allobs_cluster_vcov <- vcovCL(allobs_mod, cluster = flm_df$collection_id)
coeftest(allobs_mod, vcov = allobs_cluster_vcov)

# Limit to gymnosperms
subset_df <- trim_df %>% filter(gymno_angio=="gymno")
# subset_df <- trim_df %>% filter(gymno_angio=="gymno") %>% filter(collection_id %in% old_sites)
subset_df <- trim_df %>% filter(genus=="Pinus")
# subset_df <- trim_df %>% filter(genus=="Pinus") %>% filter(collection_id %in% old_sites)
subset_df %>% dim()
subset_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=subset_df)
subset_cluster_vcov <- vcovCL(subset_mod, cluster = subset_df$collection_id)
coeftest(subset_mod, vcov = subset_cluster_vcov)


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Investigate mechanism - tree or site --------------------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Focus on trees from genera with effect, and for which we have full weather history
# young_df <- trim_df %>%
#   filter(young==TRUE) %>%
#   drop_na()
# 
# mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trstd + pet.spstd + pet.trstd, weights = errorweights, data=young_df)
# cluster_vcov <- vcovCL(mod, cluster = young_df$collection_id)
# coeftest(mod, vcov = vcovCL, cluster = young_df$collection_id)
# 
# # Could include all trees from sites with any variation in tree-level drought history (including trees born <1901)
# young_sites <- trim_df %>%
#   filter(young == T) %>%
#   pull(collection_id) %>%
#   unique()
# youngsite_df <- trim_df %>%
#   filter(collection_id %in% young_sites) %>%
#   drop_na()
# mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trstd + pet.spstd + pet.trstd, weights = errorweights, data=youngsite_df)
# coeftest(mod, vcov = vcovCL, cluster = youngsite_df$collection_id)
# 
# # Include collection fixed effects
# mod <- feols(estimate_cwd.an ~ cwd.trstd + pet.trstd | collection_id, weights = trim_df$errorweights, data=trim_df)
# summary(mod)
# 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate variation by genus  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# predict_cwd_effect <- function(trim_df){
#   trim_df <- trim_df %>%
#     filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
#            abs(cwd.spstd)<5,
#            abs(pet.spstd)<5) %>%
#     drop_na()
#   gen_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=trim_df)
#   cluster_vcov <- vcovCL(gen_mod, cluster = trim_df$collection_id)
#   print(coeftest(gen_mod, vcov = vcovCL, cluster = trim_df$collection_id))
#   predictions <- prediction(gen_mod, at = list(cwd.spstd = seq(-3, 3, .5)), vcov = cluster_vcov, calculate_se = T) %>% 
#     summary() %>% 
#     rename(cwd.spstd = "at(cwd.spstd)") 
#   return(predictions)
# }
# # ## Plot marginal effects by gymno / angio
# grp_freq <- flm_df %>% 
#   group_by(gymno_angio) %>% 
#   summarise(n_collections = n_distinct(collection_id)) %>% 
#   arrange(desc(n_collections))
# 
# grp_keep <- grp_freq %>% 
#   filter(n_collections>50) %>% 
#   pull(gymno_angio)
# 
# grp_df <- flm_df %>%
#   filter(gymno_angio %in% grp_keep) %>% 
#   group_by(gymno_angio) %>% 
#   nest() %>%
#   mutate(prediction = map(data, predict_cwd_effect))
# 
# grp_df <- grp_df %>% 
#   unnest(prediction) %>% 
#   select(-data)
# 
# grp_df$gymno_angio <- grp_df$gymno_angio %>% recode("angio" = "Angiosperm", "gymno" = "Gymnosperm")
# 
# margins_plot <- ggplot(grp_df, aes(x = cwd.spstd)) + 
#   geom_line(aes(y = Prediction)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   facet_wrap(~gymno_angio, scales = "free", ncol = 1) +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") + 
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::comma)
# margins_plot
# ggsave(paste0(wdir, 'figures\\gymno_angio.svg'), plot = margins_plot)
# 
# ## Plot marginal effects by gymno / angio
# grp_freq <- flm_df %>% 
#   group_by(gymno_angio) %>% 
#   summarise(n_collections = n_distinct(collection_id)) %>% 
#   arrange(desc(n_collections))
# 
# grp_keep <- grp_freq %>% 
#   filter(n_collections>50) %>% 
#   pull(gymno_angio)
# 
# grp_df <- flm_df %>%
#   filter(gymno_angio %in% grp_keep) %>% 
#   group_by(gymno_angio) %>% 
#   nest() %>%
#   mutate(prediction = map(data, predict_cwd_effect))
# 
# grp_df <- grp_df %>% 
#   unnest(prediction) %>% 
#   select(-data)
# 
# grp_df$gymno_angio <- grp_df$gymno_angio %>% recode("angio" = "Angiosperm", "gymno" = "Gymnosperm")
# 
# margins_plot <- ggplot(grp_df, aes(x = cwd.spstd)) + 
#   geom_line(aes(y = Prediction)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   facet_wrap(~gymno_angio, scales = "free", ncol = 1) +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") + 
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::comma)
# margins_plot
# ggsave(paste0(wdir, 'figures\\gymno_angio.svg'), plot = margins_plot)


predict_cwd_effect <- function(trim_df){
  trim_df <- trim_df %>%
    filter(abs(cwd.spstd)<3,
           abs(pet.spstd)<3) %>%
    drop_na()
  gen_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=trim_df)
  cluster_vcov <- vcovCL(gen_mod, cluster = trim_df$collection_id)
  mod_cl <- tidy(coeftest(gen_mod, vcov = vcovCL, cluster = trim_df$collection_id))
  pred_min <- trim_df$cwd.spstd %>% quantile(0.025)
  pred_max <- trim_df$cwd.spstd %>% quantile(0.975)
  predictions <- prediction(gen_mod, at = list(cwd.spstd = seq(pred_min, pred_max, 0.1)), vcov = cluster_vcov, calculate_se = T) %>% 
    summary() %>% 
    rename(cwd.spstd = "at(cwd.spstd)")
  out <- tibble(predictions = list(predictions), model = list(mod_cl))
  return(out)
}



## Plot marginal effects by genus
genus_freq <- mod_df %>% 
  group_by(genus) %>% 
  summarise(n_collections = n_distinct(collection_id),
            min_cwd = min(cwd.spstd),
            max_cwd = max(cwd.spstd),
            range_cwd = max_cwd - min_cwd) %>% 
  arrange(desc(n_collections))

genus_keep <- genus_freq %>% 
  filter(n_collections>2) %>%
  pull(genus)


genus_key <- sp_info %>% 
  select(genus, gymno_angio) %>% 
  unique()

genus_df <- mod_df %>% 
  filter(genus %in% genus_keep) %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(prediction = map(data, predict_cwd_effect)) %>% 
  unnest(prediction) %>% 
  left_join(genus_key, by = "genus")


genus_coefs <- genus_df %>%
  unnest(model) %>% 
  filter(term == "cwd.spstd") %>% 
  select(genus, estimate, p.value)

genus_coefs <- genus_coefs %>% 
  mutate(lab = paste0("Slope: ", as.character(format(estimate, digits = 3, scientific = F)), 
                      "\nP value: ", as.character(format(p.value, digits = 3, scientific = T))))

genus_predictions <- genus_df %>% 
  unnest(predictions) %>% 
  select(-data, -model) %>% 
  left_join(genus_freq, by = "genus") %>% 
  left_join(genus_coefs, by = "genus")

saveRDS(genus_predictions, paste0(wdir, "out/second_stage/genus_mods.rds"))

# histogram <- ggplot(flm_df %>% filter(genus %in% genus_keep), aes(x = cwd.spstd)) + 
#   facet_wrap(~genus) +
#   geom_histogram(aes(fill = gymno_angio), bins = 40, alpha=0.4) +
#   xlim(c(-3, 3)) +
#   theme_bw(base_size = 22) + 
#   ylab("# sites") +
#   xlab("Historic CWD\n(Deviation from species mean)")
# histogram


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Deprecated code below this  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
margins(mod, at  = list(pet.spstd = -2:2))
cplot(mod, 'pet.spstd', what = "prediction")
persp(mod, 'cwd.spstd', 'pet.spstd')
img <- image(mod, 'cwd.spstd', 'pet.spstd')




group_dat <- plot_dat %>% 
  group_by(cwd.q) %>% 
  summarize(wvar = weighted.var(estimate_cwd.an, errorweights),
            wsd = sqrt(wvar),
            wmean = weighted.mean(estimate_cwd.an, errorweights),
            n = n(),
            error = qt(0.975, df = n-1)*wsd/sqrt(n),
            lower = wmean - error,
            upper = wmean + error)


ggplot(group_dat, aes(x=cwd.q, y=wmean)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) 


pet_dat <- plot_dat %>% 
  group_by(pet.q) %>% 
  summarize(wvar.cwd = weighted.var(estimate_cwd.an, errorweights),
            wsd.cwd = sqrt(wvar),
            wmean = weighted.mean(estimate_cwd.an, errorweights),
            n = n(),
            error = qt(0.975, df = n-1)*wsd/sqrt(n),
            lower = wmean - error,
            upper = wmean + error)

ggplot(group_dat, aes(x=pet.q, y=wmean)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) 


p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p1



## Plot estimates by historic cwd / pet - would be awesome to combine these and make prettier
plot_dat <- flm_df
coef_plot1 <- 
  plot_dat %>% 
  ggplot(aes(x = cwd.spstd, y = pet.spstd, 
             z = estimate_pet.an, weight = errorweights)) +
  stat_summary_hex() +
  xlim(c(-3, 4)) +
  #ylim(c(-1.5, 1.5))+
  scale_fill_gradientn (colours = c("darkblue","lightblue", "lightgrey")) +
  theme_bw()+
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme(legend.title = element_blank(),
        legend.position = "bottom")

coef_plot1 

coef_plot2 <- plot_dat %>% ggplot(aes(x = cwd.spstd, y = estimate_cwd.an, weight = errorweights)) +
  stat_summary_bin(bins = 20)+
  theme_bw()+
  xlim(-3, 4)+
  ylab("Estimate CWD")+
  xlab("Deviation from mean CWD")


coef_plot2

coef_plot3 <- plot_dat %>% ggplot(aes(x = pet.spstd, y = estimate_cwd.an, weight = nobs)) +
  stat_summary_bin(bins = 20)+
  xlim(-3,4)+
  xlab("Deviation from mean PET")+
  ylab("Estimate CWD")+
  theme_bw()+
  rotate()

coef_plot3
coef_plot1/coef_plot2


ggarrange(coef_plot2, NULL,coef_plot1, coef_plot3, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = T)




## Plot prediction along one axis
cwd.x <- seq(-2, 2, length.out=1000)
pred <- predict(cwd.mod, data.frame(cwd.spstd=cwd.x, pet.spstd = 0),
                interval="prediction",level=0.90)
pred <- cbind(pred, cwd.x) %>%
  as_tibble()
plot(cwd.x, pred$fit, type="l", xlab="Historic CWD",ylab="Estimated effect of CWD",las=1,lwd=2,col="#0dcfca",
     ylim = c(pred$lwr %>% min(), pred$upr %>% max()))
polygon(c(cwd.x,rev(cwd.x)),c(pred$lwr,rev(pred$upr)),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)


## Plot prediction along two axes
cwd.x <- seq(-2, 2, length.out=50)
pet.x <- seq(-2, 2, length.out=50)
grid.x <- expand.grid(cwd.x, pet.x)

pred <- predict(cwd.mod,data.frame(cwd.spstd=grid.x$Var1, pet.spstd = grid.x$Var2),
                interval="prediction",level=0.90)
pred <- cbind(pred,grid.x) %>% 
  as_tibble()

p <- pred %>% ggplot(aes(x = Var1, y = Var2, fill = fit)) +
  geom_tile()
p



## Bootstrap distributions from first stage
N <- 100

draw_est_smpl <- function(dat){
  var1 <- as.numeric(dat[which(names(dat)=="std.error_cwd.an")])^2
  var2 <- as.numeric(dat[which(names(dat)=="std.error_pet.an")])^2
  cov <- as.numeric(dat[which(names(dat)=="pet_cwd_cov")])^2
  vcov <- matrix(c(var1, cov, cov, var2), 2)
  mu <- c(data$estimate_cwd.an, data$estimate_pet.an)
  smpl <- mvrnorm(N, mu = mu, Sigma = vcov )
  colnames(smpl)=c("cwd.est.smpl","pet.est.smpl")
  smpl <- smpl %>% as_tibble()
  return(smpl)
}

draw_df <- flm_df %>% 
  rowwise %>% 
  do(draw = draw_est_smpl(.))
draw_df <- cbind(flm_df, draw_df) %>% 
  unnest(draw)
  

## Boxplots of effects by cwd / pet bins


sig_df <- flm_df %>% 
  filter(p.value_cwd.an<0.05)

plot_dat <- flm_df %>%
  # filter(genus %in% c("pi", "ps", "ju")) %>% 
  group_by("species_id") %>%
  mutate(cwd.q = as.factor(ntile(cwd.ave, 5)),
         pet.q = as.factor(ntile(pet.ave, 5)),
         cwd.high = as.factor(ntile(cwd.ave, 2)-1),
         pet.high = as.factor(ntile(pet.ave, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high)) %>%
  ungroup()



plot_dat %>% 
  mutate(inv_std = 1/ std.error_cwd.an) %>% 
  group_by(cwd.q) %>% 
  summarize(wmean = weighted.mean(estimate_cwd.an, inv_std, na.rm = TRUE))

plot_dat %>% 
  group_by(cwd.q) %>% 
  summarize(mean = mean(estimate_cwd.an, na.rm = TRUE))


p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p2 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p3 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p4 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)

(p1 | p2) / (p3 | p4)





sp <- "pipo"
sp_data <- flm_df %>% filter(species_id == sp)
sp_data %>% filter(cwd.spstd<1) %>% 
  ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) + 
  geom_point() + 
  geom_point(data = sp_data %>% filter(collection_id=="CO559"), color = 'red', alpha = 1, size = 5) + 
  geom_smooth(method = lm) +
  theme_bw(base_size = 25) +
  ylab("Sensitivity to CWD (\U03B2)")+
  xlab("Historic CWD\n(Deviation from species mean)") +
  geom_hline(yintercept = 0)

sp_data %>% arrange(estimate_cwd.an) %>% select(collection_id, estimate_cwd.an, p.value_cwd.an, cwd.ave, cwd.spstd)
(sp_data %>% arrange(desc(estimate_cwd.an)) %>% select(collection_id, estimate_cwd.an, p.value_cwd.an, cwd.ave, cwd.spstd))[10:20,]
