#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 6/11/26
# Purpose: Characterize the non-parametric response of absolute ring width (RWL)
#   to annual CWD, and test whether historically dry sites respond differently
#   to dry conditions than historically wet sites (climate adaptation).
#
#   Approach borrows from the climate-economics literature on non-parametric
#   dose-response and adaptation. With one annual CWD value per site-year, this
#   is a flexible STEP FUNCTION of the annual CWD level (coarse fixed-width
#   bins) with site/tree fixed effects -- the annual-panel form of Burke,
#   Hsiang & Miguel (2015) / Dell, Jones & Olken (2012), NOT the sub-annual
#   "days-in-bin" estimator of Deschenes & Greenstone (2011) / Schlenker &
#   Roberts (2009). Adaptation is identified by interacting the CWD bins with
#   long-run baseline climate (cwd.ave); per Merel & Gammans (2021), the
#   nonlinearity lets cross-sectional (climate) variation enter identification.
#   See also Carleton et al. (2022 QJE) and Heutel, Miller & Molitor (2021).
#
# Input files (data_source = "itrdb"):
# - rwi_long.csv: ITRDB ring width data incl. absolute RWL. From "1b. Parse ITRDB.R".
# - site_summary.csv: Site metadata incl. species (sp_id). From "1b. Parse ITRDB.R".
# - site_ave_clim.gz: Site historic baseline climate + collection_id/location_id
#     crosswalk. From "3b. Species niche.R".
# - site_an_clim_nospstd.gz: Raw (non-standardized) annual climate. From "3b. Species niche.R".
#
# Output files:
# - binned_cwd_response.rds: Tidy bin-response tables (pooled + by tercile) and
#     the adaptation interaction test.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(fixest)
library(marginaleffects)
library(splines)
library(patchwork)
library(broom)
library(data.table)
library(arrow)

select <- dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Configuration ----------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wdir <- "remote/"
data_source <- "itrdb" # "itrdb" | "fia"  <- single switch to swap datasets
species <- "PSME"
binwidth <- 50 # mm; width of the narrow annual-CWD bins in the core
core_q <- c(0.01, 0.99) # trim quantiles: narrow `binwidth` bins span this
# central range of the selected species' cwd.an; the
# sparse outer tails (typically single-tercile extremes)
# are lumped into catch-all bins. core_range is derived
# dynamically from the data below.
ns_df <- 4 # spline df for PET and ring-age controls


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data loading (abstracted across data sources) -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Returns a standardized tibble:
#   collection_id, tree_id, core_id, year, growth, ring_age,
#   cwd.an, pet.an, cwd.ave, pet.ave
load_growth_data <- function(data_source, species) {
  sp_lower <- tolower(species)

  if (data_source == "itrdb") {
    dendro_dir <- paste0(wdir, "1_input_processed/dendro/")

    # Species lookup (rwi_long.csv has no species column)
    site_smry <- read_csv(
      paste0(dendro_dir, "site_summary.csv"),
      show_col_types = FALSE
    ) %>%
      mutate(species_id = tolower(sp_id)) %>%
      select(collection_id, species_id)
    sp_collections <- site_smry %>%
      filter(species_id == sp_lower) %>%
      pull(collection_id)

    # Ring widths (1.1 GB file -> read selected cols with fread, then filter early)
    dendro <- fread(
      paste0(dendro_dir, "rwi_long.csv"),
      select = c("collection_id", "core_id", "tree", "year", "rwl")
    ) %>%
      as_tibble() %>%
      filter(collection_id %in% sp_collections) %>%
      filter(!is.na(rwl), rwl > 0) %>%
      mutate(tree_id = paste0(collection_id, "_", tree)) %>%
      group_by(core_id) %>%
      mutate(ring_age = year - min(year)) %>% # cambial-age proxy (assumes pith)
      ungroup() %>%
      rename(growth = rwl)

    # Climate: baseline (carries the collection_id <-> location_id crosswalk)
    ave_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz")) %>%
      select(collection_id, location_id, cwd.ave, pet.ave)
    an_clim <- read_rds(paste0(
      wdir,
      "2_output/climate/site_an_clim_nospstd.gz"
    )) %>%
      select(location_id, year, cwd.an, pet.an)

    df <- dendro %>%
      inner_join(ave_clim, by = "collection_id") %>%
      inner_join(an_clim, by = c("location_id", "year")) %>%
      select(
        collection_id,
        tree_id,
        core_id,
        year,
        growth,
        ring_age,
        cwd.an,
        pet.an,
        cwd.ave,
        pet.ave
      )
  } else if (data_source == "fia") {
    ## FIA placeholder. Mirrors the TerraClimate load in returning_DLNM.R, but
    ## requires an FIA ring-width file that includes absolute RWL. The current
    ## repo file (rwi_long_fia.csv) has only RWI -- drop in the RWL dataset here.
    dendro_dir <- paste0(wdir, "1_input_processed/dendro/")
    fia_rwl_path <- paste0(dendro_dir, "rwl_long_fia.csv") # TODO: user-supplied FIA RWL file
    if (!file.exists(fia_rwl_path)) {
      stop(
        "FIA RWL file not found at ",
        fia_rwl_path,
        ". Provide an FIA ring-width file with an `rwl` column to use data_source = 'fia'."
      )
    }
    dendro <- read_csv(fia_rwl_path, show_col_types = FALSE) %>%
      filter(species_id == species, !is.na(rwl), rwl > 0) %>%
      group_by(core_cn) %>%
      mutate(ring_age = year - min(year)) %>%
      ungroup() %>%
      rename(growth = rwl, core_id = core_cn)

    terraclim_df <- read_parquet(paste0(
      wdir,
      "0_raw/TerraClimate/site_climate_full.parquet"
    ))
    site_ids <- dendro %>% pull(plot_cn) %>% unique()
    terraclim_df <- terraclim_df %>%
      filter(water_year > 1958, site_id %in% site_ids) %>%
      pivot_wider(names_from = variable, values_from = value)
    tc_an <- terraclim_df %>%
      select(-year) %>%
      rename(year = water_year) %>%
      group_by(site_id, year) %>%
      summarise(cwd.an = sum(def), pet.an = sum(pet), .groups = "drop")
    tc_ave <- tc_an %>%
      filter(year < 1980) %>%
      group_by(site_id) %>%
      summarise(
        cwd.ave = mean(cwd.an),
        pet.ave = mean(pet.an),
        .groups = "drop"
      )

    df <- dendro %>%
      mutate(plot_cn = as.character(plot_cn)) %>%
      inner_join(tc_an, by = c("plot_cn" = "site_id", "year")) %>%
      inner_join(tc_ave, by = c("plot_cn" = "site_id")) %>%
      select(
        collection_id,
        tree_id,
        core_id,
        year,
        growth,
        ring_age,
        cwd.an,
        pet.an,
        cwd.ave,
        pet.ave
      )
  } else {
    stop("Unknown data_source: ", data_source)
  }

  df %>% drop_na(growth, cwd.an, pet.an, cwd.ave, pet.ave)
}

df <- load_growth_data(data_source, species)
cat(
  "Loaded",
  nrow(df),
  "core-year obs from",
  n_distinct(df$collection_id),
  "sites for",
  species,
  "(",
  data_source,
  ")\n"
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define absolute annual-CWD bins (BHM-style step function) --------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Narrow (binwidth) bins across the dense, common-support core; the sparse tails
# below/above core_range are lumped into single catch-all bins so the model isn't
# fit on near-empty extreme bins. Each obs falls in one bin. Bins are relabelled
# to clean tokens (b01..bN, ordered by CWD) so the bracket/comma labels from
# cut() don't break fixest term names or string parsing.

# Derive the core range dynamically from the selected species' cwd.an distribution
# (trim the sparse outer tails), snapped to the bin width.
core_range <- c(
  floor(quantile(df$cwd.an, core_q[1]) / binwidth) * binwidth,
  ceiling(quantile(df$cwd.an, core_q[2]) / binwidth) * binwidth
) %>%
  unname()
cat(
  "\nDynamic core range (",
  core_q[1] * 100,
  "-",
  core_q[2] * 100,
  "pctile of cwd.an): ",
  core_range[1],
  "-",
  core_range[2],
  " mm\n",
  sep = ""
)

cwd_breaks <- unique(c(
  -Inf,
  seq(core_range[1], core_range[2], by = binwidth),
  Inf
))
df <- df %>%
  mutate(cwd_bin_raw = cut(cwd.an, breaks = cwd_breaks, dig.lab = 5))

raw_levels <- levels(df$cwd_bin_raw) # already in CWD order
bin_tokens <- setNames(sprintf("b%02d", seq_along(raw_levels)), raw_levels)
df <- df %>%
  mutate(
    cwd_bin = factor(bin_tokens[as.character(cwd_bin_raw)], levels = bin_tokens)
  )

# Reference (omitted) bin = the token of the bin containing the pooled median CWD
ref_raw <- as.character(cut(
  median(df$cwd.an),
  breaks = cwd_breaks,
  include.lowest = TRUE,
  dig.lab = 5
))
ref_bin <- unname(bin_tokens[ref_raw])
df <- df %>% mutate(cwd_bin = relevel(cwd_bin, ref = ref_bin))

# Bin lookup (token -> mean CWD x-position + readable label) and obs-per-bin table
bin_lookup <- df %>%
  group_by(cwd_bin, cwd_bin_raw) %>%
  summarise(bin_mid = mean(cwd.an), n = n(), .groups = "drop") %>%
  arrange(bin_mid)
cat(
  "\n=== Observations per CWD bin (reference bin =",
  ref_bin,
  "=",
  ref_raw,
  "mm) ===\n"
)
print(bin_lookup)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper: fit binned-response model and tidy the bin coefficients --------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Age control only matters for absolute RWL (RWI is already detrended).
age_term <- if (data_source %in% c("itrdb", "fia")) {
  paste0(" + ns(ring_age, ", ns_df, ")")
} else {
  ""
}

fit_binned <- function(data) {
  fmla <- as.formula(paste0(
    "growth ~ i(cwd_bin, ref = '",
    ref_bin,
    "') + ns(pet.an, ",
    ns_df,
    ")",
    age_term,
    " | tree_id"
  ))
  feols(fmla, data = data, cluster = ~collection_id)
}

# Extract i(cwd_bin) coefficients as a response curve relative to the ref bin
tidy_bins <- function(mod, label) {
  tidy(mod) %>%
    filter(str_detect(term, "^cwd_bin::"), !str_detect(term, ":cwd.ave")) %>%
    mutate(cwd_bin = str_remove(term, "^cwd_bin::")) %>%
    bind_rows(tibble(cwd_bin = ref_bin, estimate = 0, std.error = 0)) %>% # ref bin = 0
    left_join(bin_lookup, by = "cwd_bin") %>%
    mutate(
      group = label,
      lower = estimate - 1.96 * std.error,
      upper = estimate + 1.96 * std.error
    ) %>%
    arrange(bin_mid)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 1: pooled binned response ---------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod_pooled <- fit_binned(df)
summary(mod_pooled)
pooled_resp <- tidy_bins(mod_pooled, "Pooled")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 2: stratified response by historic-CWD tercile ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Terciles assigned per site (distinct collection_id) on baseline cwd.ave.
site_terciles <- df %>%
  distinct(collection_id, cwd.ave) %>%
  mutate(
    cwd_tercile = factor(
      ntile(cwd.ave, 3),
      levels = 1:3,
      labels = c("Low CWD (wet)", "Mid CWD", "High CWD (dry)")
    )
  )
df <- df %>% left_join(site_terciles, by = c("collection_id", "cwd.ave"))

terc_levels <- levels(site_terciles$cwd_tercile)
tercile_resp <- map_dfr(terc_levels, function(g) {
  tidy_bins(fit_binned(df %>% filter(cwd_tercile == g)), g)
}) %>%
  mutate(group = factor(group, terc_levels))

# Common-support: which bins are populated across all three terciles?
support <- df %>%
  count(cwd_tercile, cwd_bin) %>%
  pivot_wider(names_from = cwd_tercile, values_from = n, values_fill = 0)
cat("\n=== Common support: obs per bin x tercile ===\n")
print(support)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 3: continuous bin x baseline-CWD interaction (adaptation test) ---
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main effect of cwd.ave is absorbed by tree_id FE; its interaction with the
# (time-varying) CWD bin is identified. A positive coef on the highest-CWD bin
# x cwd.ave => drier-baseline sites suffer a SMALLER penalty at high CWD (adaptation).
fmla_int <- as.formula(paste0(
  "growth ~ i(cwd_bin, cwd.ave, ref = '",
  ref_bin,
  "') + i(cwd_bin, ref = '",
  ref_bin,
  "')",
  " + ns(pet.an, ",
  ns_df,
  ")",
  age_term,
  " | tree_id"
))
mod_int <- feols(fmla_int, data = df, cluster = ~collection_id)

# Adaptation test: interaction coefficients (bin x cwd.ave)
adapt_test <- tidy(mod_int, conf.int = TRUE) %>%
  filter(str_detect(term, "cwd_bin") & str_detect(term, "cwd.ave")) %>%
  mutate(cwd_bin = str_remove(str_remove(term, "^cwd_bin::"), ":cwd.ave$")) %>%
  left_join(bin_lookup, by = "cwd_bin") %>%
  arrange(bin_mid) %>%
  select(cwd_bin, bin_mid, estimate, std.error, conf.low, conf.high, p.value)
cat("\n=== Adaptation test: i(cwd_bin) x cwd.ave interaction ===\n")
print(adapt_test)

# Carleton-style predicted response functions at low/mid/high baseline aridity.
# feols + FE => marginaleffects can't propagate FE uncertainty, so vcov = FALSE
# (point predictions only). The formal uncertainty is in adapt_test above.
ave_q <- quantile(site_terciles$cwd.ave, c(0.1, 0.9))
int_pred_raw <- predictions(
  mod_int,
  vcov = FALSE,
  newdata = datagrid(
    cwd_bin = levels(df$cwd_bin),
    cwd.ave = ave_q,
    pet.an = mean(df$pet.an),
    ring_age = mean(df$ring_age)
  )
)
# Rebuild as a plain tibble: marginaleffects attaches the model/jacobian to its
# columns as attributes, so coerce each column to strip them (else the 42-row
# object balloons to ~750 MB).
int_pred <- tibble(
  cwd_bin = as.character(int_pred_raw$cwd_bin),
  cwd.ave = as.numeric(int_pred_raw$cwd.ave),
  estimate = as.numeric(int_pred_raw$estimate)
) %>%
  left_join(bin_lookup, by = "cwd_bin") %>%
  mutate(
    baseline = factor(
      cwd.ave,
      levels = ave_q,
      labels = c("P10 (wet)", "P90 (dry)")
    )
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure -----------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
terc_cols <- c(
  "Low CWD (wet)" = "#0072B2",
  "Mid CWD" = "#999999",
  "High CWD (dry)" = "#D55E00"
)
baseline_cols <- c("P10 (wet)" = "#0072B2", "P90 (dry)" = "#D55E00")

## Panel A: pooled binned response
pA <- ggplot(pooled_resp, aes(bin_mid, estimate, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.15) +
  geom_line() +
  geom_point() +
  labs(
    subtitle = "A. Pooled response",
    x = "Annual CWD (mm)",
    y = "RWL effect vs. reference bin (mm)"
  ) +
  coord_cartesian(xlim = core_range) + # focus on dense, common-support core
  theme_bw()

## Panel B: stratified tercile response curves (common CWD axis)
pB <- ggplot(
  tercile_resp,
  aes(
    bin_mid,
    estimate,
    color = group,
    fill = group,
    ymin = lower,
    ymax = upper
  )
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.12, color = NA) +
  geom_line() +
  geom_point(size = 1) +
  scale_color_manual(values = terc_cols) +
  scale_fill_manual(values = terc_cols) +
  labs(
    subtitle = "B. By historic-CWD tercile",
    x = "Annual CWD (mm)",
    y = "RWL effect vs. reference bin (mm)",
    color = "Baseline CWD",
    fill = "Baseline CWD"
  ) +
  coord_cartesian(xlim = core_range) + # focus on dense, common-support core
  theme_bw()

## Panel C: interaction-model predicted curves at P10/P90 baseline CWD
## (point predictions; see adapt_test for the formal interaction uncertainty)
pC <- ggplot(int_pred, aes(bin_mid, estimate, color = baseline)) +
  geom_line() +
  geom_point(size = 1) +
  scale_color_manual(values = baseline_cols) +
  labs(
    subtitle = "C. Interaction model: predicted at baseline P10/P90",
    x = "Annual CWD (mm)",
    y = "Predicted RWL (mm)",
    color = "Baseline CWD"
  ) +
  coord_cartesian(xlim = core_range) + # focus on dense, common-support core
  theme_bw()

binned_fig <- (pA + pB) /
  pC +
  plot_annotation(
    title = paste0(
      "Non-parametric RWL response to annual CWD (",
      species,
      ", ",
      data_source,
      ")"
    ),
    subtitle = "Step function of annual CWD, controlling for PET + ring age, with tree fixed effects",
    caption = paste0(
      "Reference bin: ",
      ref_bin,
      " mm. Outcome: absolute ring width (RWL, mm). ",
      "Bin width: ",
      binwidth,
      " mm."
    )
  )
binned_fig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export -----------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list(
  pooled_response = pooled_resp,
  tercile_response = tercile_resp,
  adaptation_test = adapt_test,
  interaction_pred = int_pred,
  bin_lookup = bin_lookup,
  common_support = support
) %>%
  write_rds(paste0(wdir, "2_output/first_stage/binned_cwd_response.rds"))
