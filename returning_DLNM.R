# Plan for initial re-analysis:
# 1. Start by focusing on a single (dominant) species to avoid having to standardize climate
# 2. Start with FIA data so that we aren't dealing with any sampling biases?
# 3. Use raw terraclim data to avoid challenges about input data
# 4. Start with precip / temp rather than CWD/PET?
# 5. Focus on differentiating non-linear effects (precip has declining returns) from adaptations (drier locations are better adapted to dry conditions)
# 6. Pay attention to relative vs absolute dry conditions


library(fixest)
library(tidyverse)
library(marginaleffects)
library(tidylog)
library(arrow)
library(splines)
library(patchwork)
library(dlnm)
library(data.table)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'
dendro_dir <- paste0(wdir, '1_input_processed/dendro/')

dendro_df <- read_csv(paste0(dendro_dir, "rwi_long_fia.csv"))
site_df <- read_csv(paste0(dendro_dir, "site_summary_fia.csv"))


# # 2. Historic site-level climate
# an_site_clim <- read_rds(paste0(wdir, "2_output/climate/site_an_clim_nospstd.gz"))
# dendro_df <- dendro_df %>% 
#   mutate(plot_cn = as.character(plot_cn)) %>%
#   left_join(an_site_clim, by = c("plot_cn" = "location_id", "year"))


# # 3. Historic site-level climate
# ave_site_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz"))
# dendro_df <- dendro_df %>% 
#   left_join(ave_site_clim, by = "collection_id")

terraclim_df <- read_parquet(paste0(wdir, "0_raw/TerraClimate/site_climate_full.parquet"))
terraclim_pts <- read_parquet(paste0(wdir, "0_raw/TerraClimate/site_pixel_map.parquet"))

terraclim_df %>% pull(site_id) %>% unique()


site_ids <- dendro_df %>% pull(plot_cn) %>% unique()
terraclim_df <- terraclim_df %>%
  filter(water_year > 1958, # drop first year since it doesn't include full water year
         site_id %in% site_ids) %>% # drop non-fia sites
  pivot_wider(names_from = variable, values_from = value)

tc_an <- terraclim_df %>%
  select(-year) %>%
  rename(year = water_year) %>%
  mutate(tmean = (tmmx + tmmn) / 2) %>%
  group_by(site_id, year) %>%
  summarise(temp.an = mean(tmean),
            ppt.an = sum(pr),
            cwd.an = sum(def),
            pet.an = sum(pet))

tc_ave <- tc_an %>%
  filter(year < 1980) %>%
  group_by(site_id) %>%
  summarise(temp.ave = mean(temp.an),
            ppt.ave = mean(ppt.an),
            cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an))

dendro_df <- dendro_df %>% 
  mutate(plot_cn = as.character(plot_cn)) %>%
  inner_join(tc_an, by = c("plot_cn" = "site_id", "year"))

dendro_df <- dendro_df %>% 
  inner_join(tc_ave, by = c("plot_cn" = "site_id"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Filter to single species --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dendro_df %>%
  group_by(species_id) %>%
  summarise(n_sites = n_distinct(collection_id))

mod_df <- dendro_df %>% filter(species_id == "PSME")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pooled distributed-lag model (DLNM) of CWD / PET dynamics ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Purpose: Estimate how a CWD (and PET) shock affects RWI and how long that
#   effect persists, pooling all PSME sites. Exposure-response is linear in
#   CWD/PET; the lag dimension is a flexible natural spline (non-parametric).
#   Shocks are defined relative to the site-mean climate (a +1 / +2 SD dry-year
#   anomaly) and the model absorbs cross-site differences via collection_id
#   fixed effects.
# Note: A later extension re-runs this per site (re-introduce the group_by /
#   nest / map pattern from "4b. DLNM first stage.R").

## 1. Tree-year RWI panel -----------------------------------------------------
## Average cores to a single series per tree x year (cf. "4b" lines 57-64)
tree_df <- mod_df %>%
  group_by(collection_id, tree_id, year) %>%
  summarise(rwi = mean(rwi), .groups = "drop")


## 2. Site-year climate series for lagging ------------------------------------
## Climate is constant within a site-year across trees, so collapse to one row
## per collection_id x year and order within site so lags are well defined.
nlags <- 10  # max lag (years); FIA series are shorter than ITRDB, so < 15

site_clim <- mod_df %>%
  distinct(collection_id, year, cwd.an, pet.an) %>%
  arrange(collection_id, year) %>%
  drop_na(cwd.an, pet.an)

## Keep only sites with enough years to support the lag basis
keep_sites <- site_clim %>%
  group_by(collection_id) %>%
  tally() %>%
  filter(n > (nlags + 5)) %>%
  pull(collection_id)

site_clim <- site_clim %>%
  filter(collection_id %in% keep_sites)


## 3. Crossbasis matrices (cf. "4b" lines 124-155) ----------------------------
## Linear in the exposure (CWD / PET), natural-spline across lags.
cb_knots <- logknots(nlags, 3)

cwd_cb <- crossbasis(site_clim$cwd.an,
                     lag = c(0, nlags),
                     group = site_clim$collection_id,  # no lag bleed across sites
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "ns",
                                   knots = cb_knots,
                                   Boundary.knots = c(0, nlags)))
pet_cb <- crossbasis(site_clim$pet.an,
                     lag = c(0, nlags),
                     group = site_clim$collection_id,
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "ns",
                                   knots = cb_knots,
                                   Boundary.knots = c(0, nlags)))

## Rename basis columns so each regressor is uniquely identified (cf. "4b" 141-152)
cwd_cb_dat <- cwd_cb %>%
  as_tibble() %>%
  rename_with(~ str_replace(.x, "v1", "cwd_cbv1"), matches("v1"))
pet_cb_dat <- pet_cb %>%
  as_tibble() %>%
  rename_with(~ str_replace(.x, "v1", "pet_cbv1"), matches("v1"))

## Attach crossbasis columns to the site-year frame, then join to tree-year RWI
cb_df <- site_clim %>%
  cbind(cwd_cb_dat, pet_cb_dat)

cb_df <- tree_df %>%
  inner_join(cb_df, by = c("collection_id", "year"))


## 4. Pooled regression with site fixed effects -------------------------------
## One pooled model (no nesting by site); collection_id FE absorb baseline growth.
cwdnam <- str_subset(colnames(cwd_cb_dat), "cwd_cbv1")
petnam <- str_subset(colnames(pet_cb_dat), "pet_cbv1")
xnam <- c(cwdnam, petnam)

fmla <- as.formula(paste0("rwi ~ ", paste(xnam, collapse = " + "),
                          " | collection_id"))
lagmod <- feols(fmla, data = cb_df, cluster = ~collection_id)
summary(lagmod)


## 5. Predict lag-response & persistence (cf. "4b" lines 185-193) --------------
## Shock defined relative to the site-mean climate. With a linear exposure,
## centering at the mean and evaluating at mean + k*SD gives the distributed
## effect of a k-SD dry-year anomaly. crosspred can't read feols objects
## directly, so pass coef / vcov explicitly.
cwd_mean <- mean(site_clim$cwd.an); cwd_sd <- sd(site_clim$cwd.an)
pet_mean <- mean(site_clim$pet.an); pet_sd <- sd(site_clim$pet.an)
bylag <- 0.1

cwd_cp <- crosspred(cwd_cb,
                    coef = coef(lagmod)[cwdnam],
                    vcov = vcov(lagmod)[cwdnam, cwdnam],
                    model.link = "identity",
                    cen = cwd_mean,
                    at = cwd_mean + c(0, 1, 2) * cwd_sd,
                    bylag = bylag, cumul = TRUE)
pet_cp <- crosspred(pet_cb,
                    coef = coef(lagmod)[petnam],
                    vcov = vcov(lagmod)[petnam, petnam],
                    model.link = "identity",
                    cen = pet_mean,
                    at = pet_mean + c(0, 1, 2) * pet_sd,
                    bylag = bylag, cumul = TRUE)

## Tidy tables for the +1 SD shock (row 2 of the `at` grid). dlnm returns the
## lag-specific effect (matfit) on the fine `bylag` grid, but the cumulative
## effect (cumfit) only at integer lags -- so keep them in separate tables and
## read the lag axis off the matrix column names ("lag0", "lag0.1", ...).
shock_row <- 2  # mean + 1*SD

lag_label <- function(m) as.numeric(sub("lag", "", colnames(m)))

## Lag-specific effect (fine grid): how a +1 SD shock plays out year by year
lag_response <- bind_rows(
  tibble(var = "cwd", lag = lag_label(cwd_cp$matfit),
         effect = cwd_cp$matfit[shock_row, ],
         lower  = cwd_cp$matlow[shock_row, ],
         upper  = cwd_cp$mathigh[shock_row, ]),
  tibble(var = "pet", lag = lag_label(pet_cp$matfit),
         effect = pet_cp$matfit[shock_row, ],
         lower  = pet_cp$matlow[shock_row, ],
         upper  = pet_cp$mathigh[shock_row, ]))

## Cumulative effect (integer lags): total accumulated impact through each lag
lag_cumulative <- bind_rows(
  tibble(var = "cwd", lag = lag_label(cwd_cp$cumfit),
         cum   = cwd_cp$cumfit[shock_row, ],
         lower = cwd_cp$cumlow[shock_row, ],
         upper = cwd_cp$cumhigh[shock_row, ]),
  tibble(var = "pet", lag = lag_label(pet_cp$cumfit),
         cum   = pet_cp$cumfit[shock_row, ],
         lower = pet_cp$cumlow[shock_row, ],
         upper = pet_cp$cumhigh[shock_row, ]))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualize lag dynamics ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot A: per-lag effect of a +1 SD CWD shock (decay / persistence)
pA <- lag_response %>%
  filter(var == "cwd") %>%
  ggplot(aes(x = lag, y = effect, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = "Lagged effect of a +1 SD CWD shock on RWI",
       subtitle = "Pooled PSME, linear exposure x natural-spline lag, site FE",
       x = "Lag (years)", y = "Effect on RWI") +
  theme_bw()

## Plot B: cumulative effect accumulating over lags (total persistent impact)
pB <- lag_cumulative %>%
  filter(var == "cwd") %>%
  ggplot(aes(x = lag, y = cum, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = "Cumulative effect of a +1 SD CWD shock",
       subtitle = "Flattening indicates the persistence horizon",
       x = "Lag (years)", y = "Cumulative effect on RWI") +
  theme_bw()

pA + pB +
  plot_annotation(title = "CWD shock dynamics (pooled DLNM)")

## PET companion panel
pC <- lag_response %>%
  filter(var == "pet") %>%
  ggplot(aes(x = lag, y = effect, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = "Lagged effect of a +1 SD PET shock on RWI",
       x = "Lag (years)", y = "Effect on RWI") +
  theme_bw()
pC


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export lag-effect tables --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list(lag_response = lag_response, lag_cumulative = lag_cumulative) %>%
  write_rds(paste0(wdir, "2_output/first_stage/dnlm_pooled_lag_effects.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CWD lag dynamics by historic-CWD tercile (dry vs wet plots) ---------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Purpose: Compare the dynamic lag response to a CWD shock across plots that
#   differ in historic baseline aridity (cwd.ave). Sites are split into
#   terciles of baseline CWD and a separate DLNM is fit to each. Each tercile's
#   response is predicted under two shock definitions:
#     - a COMMON shock (same mm anomaly for all groups; vertical gaps between
#       curves reflect differences in sensitivity / persistence), and
#     - a GROUP-SPECIFIC +1 SD shock (the variability each group actually
#       experiences; mixes shock magnitude with sensitivity).
#   Because the exposure is linear, the fitted lag *shape* per tercile is the
#   same under both shocks -- only the vertical scaling differs.

## 1. Assign site-level terciles of baseline CWD -----------------------------
## Terciles are assigned per site (not per observation) so data-rich sites do
## not skew the cut points. Higher CWD = drier, so tercile 3 is the dry group.
site_terciles <- mod_df %>%
  distinct(collection_id, cwd.ave) %>%
  filter(collection_id %in% keep_sites) %>%
  mutate(cwd_tercile = factor(ntile(cwd.ave, 3), levels = 1:3,
                              labels = c("Low CWD (wet)", "Mid CWD", "High CWD (dry)")))


## 2. Helper: fit one stratified DLNM for a set of sites ----------------------
## Mirrors the pooled crossbasis -> rename -> join -> feols steps above, but on
## a site subset. Returns the CWD crossbasis plus the pieces crosspred needs.
fit_dlnm <- function(site_ids) {
  sc <- site_clim %>% filter(collection_id %in% site_ids)

  cwd_cb <- crossbasis(sc$cwd.an, lag = c(0, nlags), group = sc$collection_id,
                       argvar = list(fun = "lin"),
                       arglag = list(fun = "ns", knots = cb_knots,
                                     Boundary.knots = c(0, nlags)))
  pet_cb <- crossbasis(sc$pet.an, lag = c(0, nlags), group = sc$collection_id,
                       argvar = list(fun = "lin"),
                       arglag = list(fun = "ns", knots = cb_knots,
                                     Boundary.knots = c(0, nlags)))

  cwd_cb_dat <- cwd_cb %>% as_tibble() %>%
    rename_with(~ str_replace(.x, "v1", "cwd_cbv1"), matches("v1"))
  pet_cb_dat <- pet_cb %>% as_tibble() %>%
    rename_with(~ str_replace(.x, "v1", "pet_cbv1"), matches("v1"))

  cb_df <- sc %>%
    cbind(cwd_cb_dat, pet_cb_dat) %>%
    inner_join(tree_df, ., by = c("collection_id", "year"))

  cwdnam <- str_subset(colnames(cwd_cb_dat), "cwd_cbv1")
  petnam <- str_subset(colnames(pet_cb_dat), "pet_cbv1")
  fmla <- as.formula(paste0("rwi ~ ", paste(c(cwdnam, petnam), collapse = " + "),
                            " | collection_id"))
  mod <- feols(fmla, data = cb_df, cluster = ~collection_id)

  list(cwd_cb = cwd_cb,
       coef   = coef(mod)[cwdnam],
       vcov   = vcov(mod)[cwdnam, cwdnam],
       mean   = mean(sc$cwd.an),
       sd     = sd(sc$cwd.an),
       nsites = n_distinct(sc$collection_id))
}


## 3. Helper: crosspred a given CWD shock into tidy rows ----------------------
predict_lag <- function(fit, delta, shock_type, tercile) {
  cp <- crosspred(fit$cwd_cb, coef = fit$coef, vcov = fit$vcov,
                  model.link = "identity",
                  cen = fit$mean, at = fit$mean + c(0, delta),
                  bylag = bylag, cumul = TRUE)
  resp <- tibble(tercile = tercile, shock_type = shock_type,
                 lag = lag_label(cp$matfit),
                 est = cp$matfit[2, ], lower = cp$matlow[2, ], upper = cp$mathigh[2, ])
  cum  <- tibble(tercile = tercile, shock_type = shock_type,
                 lag = lag_label(cp$cumfit),
                 est = cp$cumfit[2, ], lower = cp$cumlow[2, ], upper = cp$cumhigh[2, ])
  list(response = resp, cumulative = cum)
}


## 4. Fit each tercile once, then predict under both shock definitions --------
tercile_sites <- split(site_terciles$collection_id, site_terciles$cwd_tercile)
fits <- map(tercile_sites, fit_dlnm)

delta_common <- cwd_sd  # pooled SD (line ~183): same mm anomaly for all groups

preds <- imap(fits, function(fit, terc) {
  list(
    predict_lag(fit, delta_common, "Common shock (pooled SD)", terc),
    predict_lag(fit, fit$sd,       "Group-specific +1 SD",     terc)
  )
}) %>% flatten()

lag_response_terc   <- map_dfr(preds, "response")
lag_cumulative_terc <- map_dfr(preds, "cumulative")

## Lock tercile ordering (wet -> dry) for plotting
terc_levels <- levels(site_terciles$cwd_tercile)
lag_response_terc   <- lag_response_terc   %>% mutate(tercile = factor(tercile, terc_levels))
lag_cumulative_terc <- lag_cumulative_terc %>% mutate(tercile = factor(tercile, terc_levels))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Comparison figure: relative vs absolute shocks by tercile -----------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
terc_cols <- c("Low CWD (wet)" = "#0072B2",
               "Mid CWD"       = "#999999",
               "High CWD (dry)" = "#D55E00")

## Builds a combined lag-response + cumulative panel for one shock definition
plot_dynamics <- function(shock, show_x = TRUE) {
  resp <- lag_response_terc   %>% filter(shock_type == shock)
  cum  <- lag_cumulative_terc %>% filter(shock_type == shock)
  xlab <- if (show_x) "Lag (years)" else NULL

  pA <- ggplot(resp, aes(lag, est, color = tercile, fill = tercile,
                         ymin = lower, ymax = upper)) +
    geom_ribbon(alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = terc_cols) +
    scale_fill_manual(values = terc_cols) +
    labs(subtitle = paste0(shock, ": lag-specific"),
         x = xlab, y = "Effect on RWI", color = "Baseline CWD", fill = "Baseline CWD") +
    theme_bw()

  pB <- ggplot(cum, aes(lag, est, color = tercile, fill = tercile,
                        ymin = lower, ymax = upper)) +
    geom_ribbon(alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = terc_cols) +
    scale_fill_manual(values = terc_cols) +
    labs(subtitle = paste0(shock, ": cumulative"),
         x = xlab, y = "Cumulative effect on RWI", color = "Baseline CWD", fill = "Baseline CWD") +
    theme_bw()

  pA + pB
}

## Caption documenting the mm magnitude of each shock
terc_sd <- map_dbl(fits, "sd")  # named by tercile label
shock_caption <- paste0(
  "CWD shock magnitudes — common: ", round(delta_common), " mm (pooled SD of annual CWD); ",
  "group-specific +1 SD: ",
  paste(sprintf("%s %.0f mm", terc_levels, terc_sd[terc_levels]), collapse = ", "), ".")

## Top row = common shock, bottom row = group-specific shock; shared legend
tercile_fig <- (plot_dynamics("Common shock (pooled SD)", show_x = FALSE) /
                plot_dynamics("Group-specific +1 SD", show_x = TRUE)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "CWD shock dynamics by historic-CWD tercile",
    subtitle = "Top: common shock (gaps = sensitivity). Bottom: group-specific +1 SD (gaps mix shock size + sensitivity).",
    caption = shock_caption,
    theme = theme(plot.caption = element_text(hjust = 0, size = 9)))
tercile_fig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export tercile lag-effect tables ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list(lag_response = lag_response_terc, lag_cumulative = lag_cumulative_terc) %>%
  write_rds(paste0(wdir, "2_output/first_stage/dnlm_tercile_lag_effects.rds"))