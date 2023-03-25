#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 9/21/20
# Purpose: Assess tree-level evolution of drought sensitivity
#
# Input files:
# - 
#
# ToDo:
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
library(dlnm)
library(tidyverse)
library(data.table)
library(tidylog)
library(dplyr)
library(data.table)
library(stringr)
library(dtplyr)
library(ggplotify)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out\\dendro\\")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            rwl = mean(rwl),
            .groups = "drop") %>% 
  mutate(tree_id = paste0(collection_id, "_", tree)) %>% 
  as_tibble()


# 2. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_an_clim.gz"))
# dendro_df <- dendro_df %>% 
#   left_join(an_site_clim, by = c("collection_id", "year"))


# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)


# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))
dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')


# 5. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out\\climate\\site_ave_clim.gz"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create lagged cwd  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlags=15
shock = 1

cwd_lag_names <- str_c("cwd_L", str_pad(0:nlags, width = 2, pad = 0))
pet_lag_names <- str_c("pet_L", str_pad(0:nlags, width = 2, pad = 0))
clim_lagged <- an_site_clim %>% 
  group_by(collection_id) %>% 
  arrange(collection_id, year) %>% 
  do(data.frame(., setNames(shift(.$cwd.an.spstd, 0:nlags), cwd_lag_names))) %>%
  do(data.frame(., setNames(shift(.$pet.an.spstd, 0:nlags), pet_lag_names)))

# Drop NAs from missing lags
clim_lagged <- clim_lagged %>% 
  ungroup() %>% 
  filter(complete.cases(clim_lagged[,grep("cwd_L", colnames(clim_lagged))]),
         complete.cases(clim_lagged[,grep("pet_L", colnames(clim_lagged))]))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create cross-basis matrices  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_lagged = data.frame(clim_lagged[,grep("cwd_L",colnames(clim_lagged))])
pet_lagged = data.frame(clim_lagged[,grep("pet_L",colnames(clim_lagged))])


cb_knots=logknots(nlags,3)
cwd_cb = crossbasis(an_site_clim$cwd.an.spstd,
                    lag=c(0,nlags),
                    group = an_site_clim$collection_id,
                    argvar=list(fun = "lin"),
                    arglag=list(fun = "ns",
                                knots=cb_knots,
                                Boundary.knots=c(0,nlags)))
pet_cb = crossbasis(an_site_clim$pet.an.spstd,
                    lag=c(0,nlags),
                    group = an_site_clim$collection_id,
                    argvar=list(fun = "lin"),
                    arglag=list(fun = "ns",
                                knots=cb_knots,
                                Boundary.knots=c(0,nlags)))


# cwd_cb = crossbasis(cwd_lagged,
#                     lag=c(0,nlags),
#                     argvar=list(fun = "ns",
#                                 Boundary.knots=c(-3,3)),
#                     arglag=list(fun = "ns",
#                                 knots=cb_knots,
#                                 Boundary.knots=c(0,10)))
# pet_cb = crossbasis(pet_lagged,
#                     lag=c(0,nlags),
#                     argvar=list(fun = "ns",
#                                 Boundary.knots=c(-3,3)),
#                     arglag=list(fun = "ns",
#                                 knots=cb_knots,
#                                 Boundary.knots=c(0,10)))


cwd_cb_dat <- cwd_cb %>% 
  as_tibble() %>% 
  rename_with(~ stringr::str_replace(.x, 
                                     pattern = "v1", 
                                     replacement = "cwd_cbv1"), 
              matches("v1")) 
pet_cb_dat <- pet_cb %>% 
  as_tibble()  %>% 
  rename_with(~ stringr::str_replace(.x, 
                                     pattern = "v1", 
                                     replacement = "pet_cbv1"), 
              matches("v1")) 

cb_df <- an_site_clim %>% 
  cbind(cwd_cb_dat, pet_cb_dat)


# cb_df <- clim_lagged %>% 
#   cbind(cwd_cb_dat, pet_cb_dat)
cb_df <- dendro_df %>% 
  inner_join(cb_df, by = c("collection_id", "year")) %>% 
  left_join(ave_site_clim, by = "collection_id") 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Filter outliers ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_spstd_bounds = quantile(ave_site_clim$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(ave_site_clim$pet.spstd, c(0.01, 0.99), na.rm = T)
rwi_bounds = quantile(dendro_df$rwi, c(0.01, 0.99), na.rm = T)

clim_sample <- ave_site_clim %>% 
  mutate(outlier = 
           (cwd.spstd<cwd_spstd_bounds[1]) |
           (cwd.spstd>cwd_spstd_bounds[2]) |
           (pet.spstd<pet_spstd_bounds[1]) |
           (pet.spstd>pet_spstd_bounds[2])) %>% 
  filter(outlier==0) %>% 
  pull(collection_id)

dendro_sample <- dendro_df %>% 
  mutate(outlier =
           (rwi<rwi_bounds[1]) |
           (rwi>rwi_bounds[2])) %>% 
  filter(outlier==0) %>% 
  pull(collection_id)

mod_df <- cb_df %>%
  filter()
  # filter(collection_id %in% clim_sample) %>%
  # filter(collection_id %in% dendro_sample) %>% 
  # filter(pet.spstd > 0.5 & cwd.spstd > 0.5)
  # filter(pet.spstd < -0 & cwd.spstd < -0)
  # filter((collection_id %in% old_incompletes))


mod_df <- cb_df %>% 
  left_join(trim_df %>% select(collection_id, cwd_errorweights), by = "collection_id")


# flm_df %>% arrange(desc(estimate_pet.an))
# 
# mod_df <- cb_df %>%
# # #   filter(cwd.spstd < 0)
# #   # filter(collection_id %in% site_sample) %>%
# #   # filter(collection_id == "CO559")
#   filter(collection_id == "RUSS094")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prep dnlm regression ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_vars <- dim(cwd_cb_dat)[2]
bylag = 0.1

petnam <- paste("pet_cbv1.l", 1:n_vars, sep="")
cwdnam <- paste("cwd_cbv1.l", 1:n_vars, sep="")
xnam <- append(cwdnam, petnam)
# fmla <- as.formula(paste("rwi ~ ", paste(xnam, collapse= "+")) %>% paste0(" | collection_id"))
# vg.range <- 363. # km, calculated in 5b. Second stage

# lagmod=feols(fmla, vcov = conley(cutoff = vg.range, distance = "spherical"), data=mod_df)
# # lagmod=feols(fmla, vcov = conley(cutoff = vg.range, distance = "spherical"), data=mod_df, weights = mod_df$cwd_errorweights)
# # fmla <- as.formula(paste("rwi ~ ", paste(xnam, collapse= "+")))
# # lagmod=feols(fmla, data=mod_df)
# 
# cwd_cp = crosspred(cwd_cb, lagmod, cen = 0, at = -1:1, bylag = bylag, cumul = TRUE)
# pet_cp = crosspred(pet_cb, lagmod, cen = 0, at = -1:1, bylag = bylag, cumul = TRUE)
# 
# lag_effects <- tibble(lag = seq(0,nlags,bylag), 
#                       cwd_effect = cwd_cp$matfit[3,], 
#                       pet_effect = pet_cp$matfit[3,],
#                       cwd_cum =  cwd_cp$cumfit[3,nlags + 1],
#                       pet_cum =  pet_cp$cumfit[3,nlags + 1])
# 
# dynamic_cwd <- plot(cwd_cp,
#      var=shock,
#      cumul=FALSE,
#      xlab=paste0("Lagged Effect of CWD=", as.integer(shock)),
#      ylab="Ring Width Growth",main=paste(""),
#      xlim=c(0,15))
# 
# 
# dynamic_pet <- plot(pet_cp,
#      var=shock,
#      cumul = FALSE,
#      xlab=paste0("Lagged Effect of PET=", as.integer(shock)),
#      ylab="Ring Width Growth",main=paste(""),
#      xlim=c(0,15))






# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Site-level lagged model  ------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# lag_df <- dendro_df %>% 
#   inner_join(clim_lagged, by = c("collection_id", "year")) %>% 
#   left_join(ave_site_clim, by = "collection_id") 
# 
# n_vars = 15
# cwd_lags <- c("cwd_L00", "cwd_L01", "cwd_L02", "cwd_L03", "cwd_L04", "cwd_L05", "cwd_L06", "cwd_L07", "cwd_L08")
# pet_lags <- c("pet_L00", "pet_L01", "pet_L02", "pet_L03", "pet_L04", "pet_L05", "pet_L06", "pet_L07", "pet_L08")
# all_lags <- append(cwd_lags, pet_lags)
# lags_fmla <- as.formula(paste("rwi ~ ", paste(all_lags, collapse= "+")))
# 
# run_lags <- function(mod_df){
#   failed <- F
#   tryCatch(
#     expr = {
#       lagmod=feols(lags_fmla, data=mod_df)
#     },
#     error = function(e) {
#       message("Returned dlnm error")
#       reg_error <<- e[1]
#       failed <<- T
#     }
#   )
#   if (failed){
#     return(NA)
#   }
#   return(lagmod %>% tidy())
# }
# 
# lag_df <- lag_df %>% 
#   group_by(collection_id) %>% 
#   nest()
# 
# lag_df <- lag_df %>% 
#   mutate(dnlm = map(data, run_lags))
# 
# lag_df <- lag_df %>% 
#   mutate(mod = map(data, run_lags)) 
# 
# test <- lag_df %>% 
#   select(-data) %>% 
#   unnest(mod) %>% 
#   filter(term %>% str_detect("cwd")) %>% 
#   mutate(lag = str_sub(term, -2, -1) %>% as.integer()) %>% 
#   group_by(lag) %>% 
#   summarize(median = median(estimate),
#             low = quantile(estimate, 0.33),
#             high = quantile(estimate, 0.66))
# test %>% 
#   ggplot() +
#   geom_line(aes(x = lag, y = median)) +
#   geom_ribbon(aes(x = lag, ymax = high, ymin = low), alpha = 0.2)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Site-level DNLM  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_dnlm <- function(mod_df){
  failed <- F
  tryCatch(
    expr = {
      fmla <- as.formula(paste("rwi ~ ", paste(xnam, collapse= "+")))
      lagmod=feols(fmla, data=mod_df)
      
      cwd_cp = crosspred(cwd_cb, lagmod, cen = 0, at = -1:1, bylag = bylag, cumul = TRUE)
      pet_cp = crosspred(pet_cb, lagmod, cen = 0, at = -1:1, bylag = bylag, cumul = TRUE)
      
      lag_effects <- tibble(lag = seq(0,nlags,bylag),
                            cwd_effect = cwd_cp$matfit[3,],
                            pet_effect = pet_cp$matfit[3,],
                            cwd_cum =  cwd_cp$cumfit[3,nlags + 1],
                            pet_cum =  pet_cp$cumfit[3,nlags + 1],
                            cwd_cum_se = cwd_cp$cumse[3,nlags + 1])
    },
    error = function(e) {
      message("Returned dlnm error")
      reg_error <<- e[1]
      failed <<- T
    }
  )
  if (failed){
    return(NA)
  }
  return(lag_effects)
}

mod_df <- cb_df %>% 
  group_by(collection_id) %>% 
  nest()

mod_df <- mod_df %>% 
  mutate(dnlm = map(data, run_dnlm))


lagged_effects <- mod_df %>%
  select(collection_id, dnlm) %>%
  unnest(dnlm) %>%
  left_join(ave_site_clim, by = "collection_id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export final cumulative effects (for robustness table run)  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cum_effects <- lagged_effects %>%
  select(collection_id, lag, cwd_cum, pet_cum, cwd_cum_se, cwd.spstd, pet.spstd) %>%
  distinct()  %>% 
  filter(lag == 15)

y_lims <- cum_effects$cwd_cum %>% quantile(c(0.025, 0.975))

cum_mod_df <- cum_effects %>% 
  filter(cwd_cum > y_lims[1],
         cwd_cum < y_lims[2]) %>%
  mutate(errorweights = 1 / cwd_cum_se)

cum_mod_df %>% 
  write_rds(paste0(wdir, "out/first_stage/dnlm_cum_effects.rds"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export lagged effects (for appendix sub-plot)  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lagged_effects <- lagged_effects %>%
  # filter(cwd.spstd>0) %>% 
  select(-cwd_cum, -pet_cum) %>%
  drop_na()

lagged_effects$cwd_quantile <- ntile(lagged_effects$cwd.spstd, 3) %>% as.factor()
lagged_effects$pet_quantile = ntile(lagged_effects$pet.spstd, 3) %>% as.factor()

lagged_effects %>% 
  write_rds(paste0(wdir, "out/first_stage/dnlm_lagged_effects.rmd"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Illustrate appendix figure ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_median_effect <- lagged_effects %>%
  # filter(cwd_quantile==4) %>%
  group_by(lag) %>%
  summarise(med_effect = median(cwd_effect),
            upper = quantile(cwd_effect, 0.66, na.rm = T),
            lower = quantile(cwd_effect, 0.33, na.rm = T))
cwd_plot <- cwd_median_effect %>%
  ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.2) +
  theme_bw() +
  xlim(0, 10) +
  ylim(-0.3, 0.1)
cwd_plot +
  xlab("Lag (years)") +
  ylab("Median effect of CWD=1 shock on RWI")


pet_median_effect <- lagged_effects %>%
  # filter(pet_quantile==4) %>%
  group_by(lag) %>%
  summarise(med_effect = median(pet_effect),
            upper = quantile(pet_effect, 0.66, na.rm = T),
            lower = quantile(pet_effect, 0.33, na.rm = T))
pet_plot <- pet_median_effect %>%
  ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.2) +
  theme_bw() +
  xlim(0, 10) +
  ylim(-0.15, 0.3)

pet_plot +
  xlab("Lag (years)") +
  ylab("Median effect of PET=1 shock on RWI")




















# # median_effect <- lagged_effects %>%
# #   group_by(lag) %>%
# #   summarise(med_effect = median(cwd_effect),
# #             upper = quantile(cwd_effect, 0.66, na.rm = T),
# #             lower = quantile(cwd_effect, 0.33, na.rm = T),
# #             upper2 = med_effect + sd(cwd_effect),
# #             lower2 = med_effect - sd(cwd_effect))
# # 
# # median_effect %>%
# #   ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
# #   geom_line() +
# #   geom_ribbon(alpha = 0.2) +
# #   theme_bw()
# # 
# # 
# # cwd_p <- ggplot() +
# #   geom_line(aes(x = lag, y = med_effect), data = median_effect, color = "cornflowerblue") +
# #   geom_ribbon(aes(x = lag, ymax = upper, ymin = lower), data = median_effect,
# #               fill = "cornflowerblue", alpha = 0.3) +
# #   theme_bw()
# # cwd_p
# # 
# # 
# # p <- ggplot() +
# #   geom_line(aes(x = lag, y = cwd_effect, group = collection_id), data = lagged_effects,
# #             alpha = 0.05) +
# #   geom_line(aes(x = lag, y = med_effect), data = median_effect, color = "cornflowerblue") +
# #   geom_ribbon(aes(x = lag, ymax = upper, ymin = lower), data = median_effect,
# #               fill = "cornflowerblue", alpha = 0.3) +
# #   ylim(c(-0.6, 0.25)) +
# #   theme_bw()
# # p
# # 
# # 
# # 
# # 
# # Look at variation across terciles. Particularly interesting for pet
# lagged_effects <- lagged_effects %>% 
#   group_by(collection_id) %>% 
#   mutate(cwd_tercile = ntile(x = cwd.spstd, 3),
#          pet_tercile = ntile(x = pet.spstd, 3))
# 
# median_effect <- lagged_effects %>%
#   group_by(cwd_quantile, lag) %>%
#   summarise(med_effect = median(cwd_effect),
#             upper = quantile(cwd_effect, 0.66, na.rm = T),
#             lower = quantile(cwd_effect, 0.33, na.rm = T),
#             upper2 = med_effect + sd(cwd_effect),
#             lower2 = med_effect - sd(cwd_effect))
# median_effect %>%
#   ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
#   facet_grid(rows = median_effect$cwd_quantile) +
#   geom_line() +
#   geom_ribbon(alpha = 0.2) +
#   # geom_hline(yintercept = 0) +
#   theme_bw()
# 
# 
# 
# median_effect <- lagged_effects %>%
#   group_by(lag, pet_quantile) %>%
#   summarise(med_effect = median(pet_effect),
#             upper = quantile(pet_effect, 0.66, na.rm = T),
#             lower = quantile(pet_effect, 0.33, na.rm = T))
# median_effect %>%
#   ggplot(aes(x = lag, y = med_effect, ymax = upper, ymin = lower)) +
#   geom_line() +
#   # geom_hline(yintercept = 0) +
#   geom_ribbon(alpha = 0.2) +
#   facet_grid(rows = median_effect$pet_quantile) +
#   theme_bw()
# 
# 
# 
# cwd_lagged
# 
# 
# # 
# # 
# # thresh_low <- cum_effects$cwd_cum %>% quantile(0.02, na.rm = T)
# # thresh_high <- cum_effects$cwd_cum %>% quantile(0.98, na.rm = T)
# # mod_df <- cum_effects %>% filter(cwd_cum>thresh_low, cwd_cum<thresh_high)
# # mod <- lm(cwd_cum ~ cwd.ave + pet.ave, data = mod_df)
# # summary(mod)
# # 
# # thresh_low <- cum_effects$pet_cum %>% quantile(0.02, na.rm = T)
# # thresh_high <- cum_effects$pet_cum %>% quantile(0.98, na.rm = T)
# # mod_df <- cum_effects %>% filter(pet_cum>thresh_low, pet_cum<thresh_high)
# # mod <- lm(pet_cum ~ cwd.ave + pet.ave, data = mod_df)
# # summary(mod)
# # 
# # 
# # cum_effects %>% 
# #   filter(cwd_cum>thresh_low, cwd_cum<thresh_high) %>% 
# #   ggplot(aes(y = cwd_cum)) +
# #   geom_boxplot()
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # # Aggregated DLNM  ------------------------------
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # gendat <- dendro_lagged
# # 
# # gendat <- gendat %>% 
# #   group_by(collection_id) %>% 
# #   nest()
# # 
# # keep_sites = site_smry %>% 
# #   filter(species_id == "pcgl") %>% 
# #   pull(collection_id) %>% 
# #   unique()
# # 
# # cwd_thresh <- hist_clim_df$cwd.ave %>% 
# #   quantile(0.66)
# # pet_thresh <- hist_clim_df$pet.ave %>% 
# #   quantile(0.66)
# # 
# # 
# # keep_sites = hist_clim_df %>% 
# #   filter(pet.ave > pet_thresh,
# #          cwd.ave > cwd_thresh) %>%
# #   # sample_frac(0.25) %>%
# #   pull(collection_id) %>% 
# #   unique()
# # 
# # gendat = dendro_lagged %>%
# #   filter(collection_id %in% keep_sites)
# # 
# # mod <- feols(rwi ~ pet.an + cwd.an + pet_L01 + cwd_L01 + pet_L02 + cwd_L02 + pet_L03 + cwd_L03 +
# #                pet_L05 + cwd_L05 + pet_L05 + cwd_L05| id, data = gendat)
# # mod <- lm(rwi ~ pet.an + cwd.an * cwd.ave + pet.an * pet.ave, data = gendat)
# # summary(mod)
# # 
# # at_vals = seq((shock - 1), (shock + 1), 0.01)
# # 
# # cwd_lagged = data.frame(gendat$cwd.an,gendat[,grep("cwd_L",colnames(gendat))])
# # pet_lagged = data.frame(gendat$pet.an,gendat[,grep("pet_L",colnames(gendat))])
# # 
# # cwd_cb = crossbasis(cwd_lagged,lag=c(0,nlags),argvar=list(fun = "lin"),arglag=list(knots=logknots(30, 2)))
# # pet_cb = crossbasis(pet_lagged, lag=c(0,nlags), argvar=list(fun = "lin"), arglag=list(knots=logknots(30, 2)))
# # 
# # 
# # lagmod=lm(rwi ~ cwd_cb + pet_cb, data=gendat)
# # summary(lagmod)
# # cwd_cp = crosspred(cwd_cb, lagmod, cen = 0, at = -1:1, bylag = 0.1, cumul = TRUE)
# # plot(cwd_cp,
# #      var=shock,
# #      xlab=paste0("Lagged Effect of CWD=", as.integer(shock)),
# #      ylab="Ring Width Growth",main=paste("Site"),cumul=FALSE)
# # 
# # pet_cp = crosspred(pet_cb, lagmod, cen = 0, at = -1:1, bylag = 0.1, cumul = TRUE)
# # plot(pet_cp,
# #      var=shock,
# #      xlab=paste0("Lagged Effect of PET=", as.integer(shock)),
# #      ylab="Ring Width Growth",main=paste("Site"),cumul=FALSE)
# # 
# # 
# # 
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # # Genus-specific distributed lag models  ------------------------------
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # shock <- 300
# # at_vals = 0:1000
# # relgenus=c("Juniperus", "Pinus", "Pseudotsuga", "Abies", "Fagus", "Picea", "Larix", "Quercus", "Tsuga")
# # genlist=list()
# # 
# # i = 1
# # for(i in 1:length(relgenus)){
# #   gendat=dendro_lagged%>%
# #     filter(genus==relgenus[i])%>%
# #     filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
# #   # shock = gendat$cwd.an %>% quantile(.9, na.rm = TRUE)
# #   # at_vals = (shock - 50):(shock + 50)
# #   lagged=data.frame(gendat$cwd.an,gendat[,grep("cwd_L",colnames(gendat))])
# #   cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(30,4)))
# #   gendat$id=interaction(gendat$collection_id,gendat$tree)
# #   lagmod=lm(gendat$width~cblagged+gendat$pet.an,data=gendat)
# #   # crosspredict = crosspred(cblagged, lagmod, cen=0, at=0:cblim[i]*1, cumul = FALSE)
# #   genlist[[i]]=list(cblagged,lagmod)
# #   genlist[[i]][[3]]=crosspred(cblagged,lagmod,cen=0,at=at_vals,cumul=FALSE)
# #   print(i)
# # }
# # 
# # +x11()
# # par(mfrow=c(3,3))
# # for(i in 1:length(relgenus)){
# #   plot(genlist[[i]][[3]],
# #        var=300,
# #        xlab=paste0("Lagged Effect of CWD=", as.integer(shock)),
# #        ylab="Log Ring Width Growth",main=paste("Genus=",relgenus[i]),cumul=FALSE)
# # }
# # 
# # intlist=list()
# # 
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # # Site-level distributed lag model  ------------------------------
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # sites <- dendro_lagged %>% 
# #   pull(collection_id) %>% 
# #   unique()
# # n_site = indices[i]
# # indices = c(1703, 302, 1252, 1367, 2028)
# # genlist=list()
# # 
# # for(i in 1:length(indices)){
# #   n_site <- indices[i]
# #   gendat=dendro_lagged%>%
# #     filter(collection_id==sites[n_site])%>%
# #     filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
# #   lagged=data.frame(gendat$cwd.an,gendat[,grep("cwd_L",colnames(gendat))])
# #   cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(30,4)))
# #   gendat$id=interaction(gendat$collection_id,gendat$tree)
# #   lagmod=lm(gendat$width~cblagged+gendat$pet.an,data=gendat)
# #   genlist[[i]]=list(cblagged,lagmod)
# #   genlist[[i]][[3]]=crosspred(cblagged,lagmod,cen=0,at=0:cblim[i]*1,cumul=TRUE)
# #   print(i)
# # }
# # 
# # 
# # x11()
# # par(mfrow=c(3,3))
# # for(i in 1:length(indices)){
# #   plot(genlist[[i]][[3]],var=300,xlab="Lagged Effect of CWD=800",ylab="Log Ring Width Growth",main=paste("Genus=",relgenus[i]),cumul=FALSE)
# # }
# # # subset_dat <- dendro_lagged%>%
# # #   filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T),
# # #          genus %in% c("Juniperus", "Pinus", "Pseudotsuga", "Abies", "Fagus", "Picea", "Larix"))
# # # lagged=data.frame(subset_dat$cwd.an,subset_dat[,grep("cwd_L",colnames(subset_dat))])
# # # cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(nlags,4)))
# # # subset_dat$id=interaction(subset_dat$collection_id,subset_dat$tree)
# # # lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
# # # prediction <- crosspred(cblagged,lagmod,cen=0,at=0:3000*1,cumul=TRUE)
# # # plot(prediction, var = 100, xlab = "Lagged Effect of CWD=800", ylab = "log Ring Width Growth")
# # 
# # 
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # # Lag interaction models  ------------------------------
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # for(i in 1:length(relgenus)){
# #   gendat=dendro_lagged%>%
# #     filter(genus==relgenus[i])%>%
# #     filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))%>%
# #     mutate(lag_5=rowMeans(select(.,cwd_L1,cwd_L2,cwd_L3,cwd_L4,cwd_L5)))%>%
# #     mutate(lag_6_10=rowMeans(select(.,cwd_L6,cwd_L7,cwd_L8,cwd_L9,cwd_L10)))
# #   gendat$id=interaction(gendat$collection_id,gendat$tree)
# #   intmod=felm(ln_rwi~cwd.an*lag_5+cwd.an*lag_6_10+pet.an|id|0|collection_id,data=gendat)
# #   intlist[[i]]=intmod
# #   print(i)
# # }
# # 
# # #try interactions model for whole dataset
# # dendro_lagged <- dendro_lagged %>%
# #   mutate(lag_5 = rowMeans(select(.,cwd_L01,cwd_L02,cwd_L03,cwd_L04,cwd_L05)),
# #          lag_2_5 = rowMeans(select(.,cwd_L02,cwd_L03,cwd_L04,cwd_L05)),
# #          lag_6_10 = rowMeans(select(.,cwd_L06,cwd_L07,cwd_L08,cwd_L09,cwd_L10)))
# # dendro_lagged$id=interaction(dendro_lagged$collection_id,dendro_lagged$tree)
# # intmod=felm(rwi~cwd.an*lag_2_5+cwd.an*lag_6_10+pet.an|id|0|collection_id,data=dendro_lagged)
# # summary(intmod)
# # intmod=felm(rwi~cwd.an + lag_5 + lag_6_10 + pet.an|id|0|collection_id,data=dendro_lagged)
# # summary(intmod)
# # 
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # # Early life drought models  ------------------------------
# # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # young_trees <- dendro_lagged %>% 
# #   group_by(collection_id, tree) %>% 
# #   summarise(min_year = min(year)) %>% 
# #   mutate(young = min_year>1901)
# # 
# # # flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv')) %>%
# # #   select(collection_id,tree,young)
# # young_trees <- young_trees %>%
# #   select(-min_year) %>% 
# #   right_join(dendro_lagged, by = c("collection_id", "tree")) %>%
# #   filter(young==TRUE)%>%
# #   select(collection_id,tree,year,cwd.an)%>%
# #   group_by(collection_id,tree)%>%
# #   filter(year<(min(year)+10))%>%
# #   summarize(earlylifecwd=mean(cwd.an))
# # 
# # youngtrees_df <- dendro_lagged %>%
# #   select(collection_id, tree, year, sp_id, ln_rwi, cwd.an, pet.an) %>% 
# #   inner_join(young_trees, by = c("collection_id", "tree"))
# # 
# # earlylifelist=list()
# # #interactions effect
# # for(i in 1:length(relgenus)){
# #   gendat=youngtrees_df%>%
# #     filter(genus == relgenus[i])%>%
# #     filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
# #   gendat$id=interaction(gendat$collection_id,gendat$tree)
# #   intmod=felm(ln_rwi~cwd.an*I(earlylifecwd>300)+pet.an-I(earlylifecwd>300)|id|0|collection_id,data=gendat)
# #   earlylifelist[[i]]=intmod
# #   print(i)
# # }
# # 
# # youngtrees_df$id=interaction(youngtrees_df$collection_id,youngtrees_df$tree)
# # youngtrees_df=youngtrees_df%>%
# #   filter(earlylifecwd<quantile(earlylifecwd,p=0.99,na.rm=T))%>%
# #   filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
# # earlylifemod=felm(ln_rwi~cwd.an*I(earlylifecwd>300)+pet.an-I(earlylifecwd>300)|id|0|collection_id,data=youngtrees_df)
# # earlylifemod=felm(ln_rwi~cwd.an*earlylifecwd + pet.an - earlylifecwd|id|0|collection_id,data=youngtrees_df)
# # 
# # 
# # youngtrees_df <- youngtrees_df %>% 
# #   left_join(hist_clim_df, by = "collection_id")
# # 
# # subset <- youngtrees_df %>% 
# #   filter(sp_id == "psme")
# # earlylifemod=felm(ln_rwi~cwd.an*earlylifecwd + cwd.an*cwd.ave + pet.an |id|0|collection_id,data=subset)
# # summary(earlylifemod)
# # 
# # 
# # earlylifemod=felm(ln_rwi~cwd.an*cwd.ave + pet.an |id|0|collection_id,data=subset)
# # summary(earlylifemod)
# # 
# # 
# # hist_clim_temp <- hist_clim_df %>% 
# #   left_join(site_df, by = "collection_id") 
# # hist_clim_temp <- hist_clim_temp %>% 
# #   group_by(sp_id) %>% 
# #   summarise(sp.cwd.median = median(cwd.ave, na.rm = TRUE)) %>% 
# #   right_join(hist_clim_temp, by = "sp_id") %>% 
# #   mutate(high_cwd_site = cwd.ave > sp.cwd.median)
# # 
# # cblim <- 1000
# # subset_dat <- dendro_lagged %>%
# #   left_join(hist_clim_temp %>% select(collection_id, high_cwd_site), by = "collection_id") %>% 
# #   filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T),
# #          sp_id =="pipo")
# # lagged=data.frame(subset_dat$cwd.an,subset_dat[,grep("cwd_L",colnames(subset_dat))])
# # cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(nlags,4)))
# # subset_dat$id=interaction(subset_dat$collection_id,subset_dat$tree)
# # lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
# # prediction <- crosspred(cblagged,lagmod,cen=0,at=0:3000*1,cumul=TRUE)
# # plot(prediction, var = 100, xlab = "Lagged Effect of CWD=800", ylab = "log Ring Width Growth")
# # 
# # 
# # subset_dat <- dendro_lagged %>%
# #   left_join(hist_clim_temp %>% select(collection_id, high_cwd_site), by = "collection_id") %>% 
# #   filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T),
# #          sp_id =="psme",
# #          high_cwd_site==0)
# # lagged=data.frame(subset_dat$cwd.an,subset_dat[,grep("cwd_L",colnames(subset_dat))])
# # cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(nlags,4)))
# # subset_dat$id=interaction(subset_dat$collection_id,subset_dat$tree)
# # lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
# # prediction <- crosspred(cblagged,lagmod,cen=0,at=0:3000*1,cumul=TRUE)
# # plot(prediction, var = 100, xlab = "Lagged Effect of CWD=800", ylab = "log Ring Width Growth")
# # 
# # 
# # 
# # 
# # 
# # calc_germ_clim <- function(site, germ_year){
# #   pre_germ_clim <- clim_df %>% 
# #     filter(collection_id == site,
# #            year < germ_year,
# #            year < 1970) %>% 
# #     summarise(pre_cwd = mean(cwd.an),
# #               pre_pet = mean(cwd.an))
# #   post_germ_clim <- clim_df %>% 
# #     filter(collection_id == site,
# #            year>germ_year,
# #            year<1970) %>% 
# #     summarise(post_cwd = mean(cwd.an),
# #               post_pet = mean(pet.an))
# #   germ_clim <- pre_germ_clim %>% 
# #     left_join(post_germ_clim, by = "collection_id")
# #   return(germ_clim)
# # }
# # 
# # young_trees <- young_trees %>% 
# #   filter(young==TRUE) %>% 
# #   mutate(germ_clim = map2(collection_id, min_year, calc_germ_clim))
# # 
# # youngtrees_df <- dendro_lagged %>%
# #   select(collection_id, tree, year, sp_id, ln_rwi, cwd.an, pet.an) %>% 
# #   inner_join(young_trees, by = c("collection_id", "tree"))
# # 
# # lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
# # 
# # 
# # 
# # 
# # 
# # 
# # # ### Binned plot of cwd sensitivity - could adapt to dlnm?
# # # seq_min <- -2.625
# # # seq_max <- 2.625
# # # seq_inc <- 0.25
# # # sequence <- seq(seq_min, seq_max, seq_inc)
# # # 
# # # convert_bin <- function(n){
# # #   sequence[n] + 0.125
# # # }
# # # 
# # # plot_dat <- cum_effects %>%
# # #   filter(cwd_cum>thresh_low, cwd_cum<thresh_high) %>%
# # #   filter(((abs(cwd.ave)<3) & (abs(pet.ave<3)))) %>%
# # #   drop_na()
# # # 
# # # plot_dat_a <- plot_dat %>%
# # #   mutate(cwd.q = cut(cwd.ave, breaks = sequence, labels = FALSE),
# # #          cwd.q = convert_bin(cwd.q),
# # #          pet.q = cut(pet.ave, breaks = sequence, labels = FALSE),
# # #          pet.q = convert_bin(pet.q))
# # # 
# # # 
# # # plot_dat_b <- plot_dat_a %>%
# # #   group_by(cwd.q, pet.q) %>%
# # #   dplyr::summarize(cwd_sens = mean(cwd_cum, na.rm = TRUE),
# # #                    pet_sens = mean(pet_cum, na.rm = TRUE),
# # #                    n = n()) %>%
# # #   filter(n>10)
# # # 
# # # 
# # # binned_margins <- plot_dat_b %>%
# # #   ggplot(aes(x = cwd.q, y = pet.q, fill = pet_sens)) +
# # #   geom_tile() +
# # #   xlim(c(-2, 1.1))+
# # #   ylim(c(-2,1.1))+
# # #   # scale_fill_viridis_c(direction = -1) +
# # #   scale_fill_continuous_diverging(rev = TRUE, mid = 0) +
# # #   ylab("Deviation from mean PET")+
# # #   xlab("Deviation from mean CWD")+
# # #   theme_bw(base_size = 22)+
# # #   theme(legend.position = c(.18,.83),
# # #         legend.key = element_blank(),
# # #         legend.background = element_blank())+
# # #   #panel.grid.major = element_blank(),
# # #   #panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))+
# # #   labs(fill = "Marginal effect\nof CWD") +
# # #   ylab("Historic PET\n(Deviation from species mean)") +
# # #   xlab("Historic CWD\n(Deviation from species mean)") +
# # #   coord_fixed() +
# # #   geom_hline(yintercept = 0, size = 1, linetype = 2) +
# # #   geom_vline(xintercept = 0, size = 1, linetype = 2)
# # # 
# # # 
# # # binned_margins
