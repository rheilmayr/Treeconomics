#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 6/11/26
# Purpose: Pooled, MULTI-SPECIES version of the binned CWD response + adaptation
#   analysis (cf. binned_cwd_response.R, which runs one species at a time). Tests
#   whether the non-parametric CWD -> growth response and the climate-adaptation
#   pattern are GENERAL across species, while only ever comparing sites within the
#   same species against one another.
#
#   Design:
#   - Outcome is log(RWL): a proportional (%) growth response, comparable across
#     species of very different absolute ring size, and additive with tree FE.
#   - All sites/species/years pooled into single models with tree fixed effects
#     (tree c site c species), so identification is within-species/within-tree.
#   - Historic-CWD terciles are defined WITHIN EACH SPECIES (a site's baseline
#     aridity ranked among its conspecific sites), so "dry-adapted" is comparable
#     across species. This also restores common support: a given absolute-CWD bin
#     contains both "dry-for-its-species" and "wet-for-its-species" sites.
#   - CWD exposure bins stay ABSOLUTE (mm): water deficit is the physical stressor;
#     adaptation is relative to what the population is used to (climate-econ logic;
#     Burke-Hsiang-Miguel 2015, Merel & Gammans 2021, Carleton et al. 2022).
#
# Input files (data_source = "itrdb"):
# - rwi_long.csv: ITRDB ring width data incl. absolute RWL. From "1b. Parse ITRDB.R".
# - site_summary.csv: Site metadata incl. species (sp_id). From "1b. Parse ITRDB.R".
# - site_ave_clim.gz: Site historic baseline climate + collection_id/location_id crosswalk.
# - site_an_clim_nospstd.gz: Raw (non-standardized) annual climate.
#
# Output files:
# - binned_cwd_response_pooled.rds: Tidy response tables (pooled + by species-
#     relative tercile) and the adaptation interaction tests.
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
wdir        <- "remote/"
data_source <- "itrdb"   # "itrdb" | "fia"  <- single switch to swap datasets
min_sites   <- 10        # keep species with >= this many sites (for within-sp terciles)
binwidth    <- 50          # mm; width of the narrow annual-CWD bins in the core
core_q      <- c(0.01, 0.99) # trim quantiles defining the narrow-bin core (pooled cwd.an)
ns_df       <- 4           # spline df for PET and ring-age controls


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data loading (all qualifying species pooled) --------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Returns a standardized tibble:
#   species_id, collection_id, tree_id, core_id, year, growth, ring_age,
#   cwd.an, pet.an, cwd.ave, pet.ave
load_growth_data_pooled <- function(data_source, min_sites) {
  if (data_source == "itrdb") {
    dendro_dir <- paste0(wdir, "1_input_processed/dendro/")

    # Species per collection; keep species with enough sites for within-sp terciles
    site_smry <- read_csv(paste0(dendro_dir, "site_summary.csv"), show_col_types = FALSE) %>%
      mutate(species_id = tolower(sp_id)) %>%
      select(collection_id, species_id)
    keep_species <- site_smry %>% count(species_id, name = "n_sites") %>%
      filter(n_sites >= min_sites) %>% pull(species_id)
    keep_collections <- site_smry %>% filter(species_id %in% keep_species) %>% pull(collection_id)

    # Ring widths (1.1 GB -> read selected cols with fread, filter to kept collections)
    dendro <- fread(paste0(dendro_dir, "rwi_long.csv"),
                    select = c("collection_id", "core_id", "tree", "year", "rwl")) %>%
      as_tibble() %>%
      filter(collection_id %in% keep_collections, !is.na(rwl), rwl > 0) %>%
      left_join(site_smry, by = "collection_id") %>%
      mutate(tree_id = paste0(collection_id, "_", tree)) %>%
      group_by(core_id) %>%
      mutate(ring_age = year - min(year)) %>%   # cambial-age proxy (assumes pith)
      ungroup() %>%
      rename(growth = rwl)

    # Climate: baseline (carries the collection_id <-> location_id crosswalk)
    ave_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz")) %>%
      select(collection_id, location_id, cwd.ave, pet.ave)
    an_clim <- read_rds(paste0(wdir, "2_output/climate/site_an_clim_nospstd.gz")) %>%
      select(location_id, year, cwd.an, pet.an)

    df <- dendro %>%
      inner_join(ave_clim, by = "collection_id") %>%
      inner_join(an_clim, by = c("location_id", "year")) %>%
      select(species_id, collection_id, tree_id, core_id, year, growth, ring_age,
             cwd.an, pet.an, cwd.ave, pet.ave)

  } else if (data_source == "fia") {
    ## FIA placeholder. Requires an FIA ring-width file that includes absolute RWL
    ## with a species_id column (the repo file rwi_long_fia.csv has only RWI).
    dendro_dir <- paste0(wdir, "1_input_processed/dendro/")
    fia_rwl_path <- paste0(dendro_dir, "rwl_long_fia.csv")  # TODO: user-supplied FIA RWL file
    if (!file.exists(fia_rwl_path)) {
      stop("FIA RWL file not found at ", fia_rwl_path,
           ". Provide an FIA ring-width file with `rwl` + `species_id` to use data_source = 'fia'.")
    }
    dendro <- read_csv(fia_rwl_path, show_col_types = FALSE) %>%
      filter(!is.na(rwl), rwl > 0)
    keep_species <- dendro %>% distinct(species_id, collection_id) %>%
      count(species_id, name = "n_sites") %>% filter(n_sites >= min_sites) %>% pull(species_id)
    dendro <- dendro %>%
      filter(species_id %in% keep_species) %>%
      group_by(core_cn) %>% mutate(ring_age = year - min(year)) %>% ungroup() %>%
      rename(growth = rwl, core_id = core_cn)

    terraclim_df <- read_parquet(paste0(wdir, "0_raw/TerraClimate/site_climate_full.parquet"))
    site_ids <- dendro %>% pull(plot_cn) %>% unique()
    terraclim_df <- terraclim_df %>%
      filter(water_year > 1958, site_id %in% site_ids) %>%
      pivot_wider(names_from = variable, values_from = value)
    tc_an <- terraclim_df %>% select(-year) %>% rename(year = water_year) %>%
      group_by(site_id, year) %>%
      summarise(cwd.an = sum(def), pet.an = sum(pet), .groups = "drop")
    tc_ave <- tc_an %>% filter(year < 1980) %>%
      group_by(site_id) %>%
      summarise(cwd.ave = mean(cwd.an), pet.ave = mean(pet.an), .groups = "drop")

    df <- dendro %>%
      mutate(plot_cn = as.character(plot_cn)) %>%
      inner_join(tc_an, by = c("plot_cn" = "site_id", "year")) %>%
      inner_join(tc_ave, by = c("plot_cn" = "site_id")) %>%
      select(species_id, collection_id, tree_id, core_id, year, growth, ring_age,
             cwd.an, pet.an, cwd.ave, pet.ave)
  } else {
    stop("Unknown data_source: ", data_source)
  }

  df %>% drop_na(growth, cwd.an, pet.an, cwd.ave, pet.ave)
}

df <- load_growth_data_pooled(data_source, min_sites)
cat("Loaded", nrow(df), "core-year obs from", n_distinct(df$collection_id),
    "sites across", n_distinct(df$species_id), "species (", data_source, ")\n")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Species-relative historic-CWD terciles --------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Rank each site's baseline CWD WITHIN its own species, so "dry"/"wet" mean
# dry/wet *for that species*. cwd_ave_rank (0-1) is the continuous analog.
site_terc <- df %>%
  distinct(species_id, collection_id, cwd.ave) %>%
  group_by(species_id) %>%
  mutate(cwd_tercile = factor(ntile(cwd.ave, 3), levels = 1:3,
                              labels = c("Wet for sp.", "Mid for sp.", "Dry for sp.")),
         cwd_ave_rank = percent_rank(cwd.ave)) %>%
  ungroup()
df <- df %>% left_join(site_terc, by = c("species_id", "collection_id", "cwd.ave"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Absolute annual-CWD bins (pooled across species) ----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Narrow bins across the dense core (pooled cwd.an), sparse tails lumped into
# catch-all bins. Relabel to clean tokens b01..bN to keep fixest term names parseable.
core_range <- c(floor(quantile(df$cwd.an, core_q[1]) / binwidth) * binwidth,
                ceiling(quantile(df$cwd.an, core_q[2]) / binwidth) * binwidth) %>%
  unname()
cat("\nDynamic core range (", core_q[1] * 100, "-", core_q[2] * 100, "pctile of cwd.an): ",
    core_range[1], "-", core_range[2], " mm\n", sep = "")

cwd_breaks <- unique(c(-Inf, seq(core_range[1], core_range[2], by = binwidth), Inf))
df <- df %>% mutate(cwd_bin_raw = cut(cwd.an, breaks = cwd_breaks, dig.lab = 5))

raw_levels <- levels(df$cwd_bin_raw)
bin_tokens <- setNames(sprintf("b%02d", seq_along(raw_levels)), raw_levels)
df <- df %>% mutate(cwd_bin = factor(bin_tokens[as.character(cwd_bin_raw)], levels = bin_tokens))

ref_raw <- as.character(cut(median(df$cwd.an), breaks = cwd_breaks, include.lowest = TRUE, dig.lab = 5))
ref_bin <- unname(bin_tokens[ref_raw])
df <- df %>% mutate(cwd_bin = relevel(cwd_bin, ref = ref_bin))

bin_lookup <- df %>%
  group_by(cwd_bin, cwd_bin_raw) %>%
  summarise(bin_mid = mean(cwd.an), n = n(), .groups = "drop") %>%
  arrange(bin_mid)
cat("\n=== Observations per CWD bin (reference bin =", ref_bin, "=", ref_raw, "mm) ===\n")
print(bin_lookup)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper: fit pooled binned-response model and tidy bin coefficients -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Outcome is log(growth) -> coefficients are proportional growth responses.
# Ring-age trend is SPECIES-SPECIFIC: each species gets its own age curve. We do
# this via fixest VARYING SLOPES on the FE side -- species_id[rage1..rageD] absorbs
# species-specific slopes on the natural-spline age basis in the demeaning step.
# (A dense ns(ring_age):species_id interaction would materialize ~n_species*ns_df
# columns over ~11M rows and blow past memory; varying slopes avoid that.)
# PET control stays pooled.
rage_basis <- ns(df$ring_age, ns_df)
colnames(rage_basis) <- paste0("rage", seq_len(ns_df))
df <- bind_cols(df, as_tibble(rage_basis))
age_fe <- paste0("species_id[", paste0("rage", seq_len(ns_df), collapse = ", "), "]")

fit_binned <- function(data) {
  fmla <- as.formula(paste0(
    "log(growth) ~ i(cwd_bin, ref = '", ref_bin, "') + ns(pet.an, ", ns_df, ")",
    " | tree_id + ", age_fe))
  feols(fmla, data = data, cluster = ~collection_id)
}

# Extract i(cwd_bin) coefficients as a response curve relative to the ref bin
tidy_bins <- function(mod, label) {
  tidy(mod) %>%
    filter(str_detect(term, "^cwd_bin::"), !str_detect(term, ":cwd_ave_rank")) %>%
    mutate(cwd_bin = str_remove(term, "^cwd_bin::")) %>%
    bind_rows(tibble(cwd_bin = ref_bin, estimate = 0, std.error = 0)) %>%  # ref bin = 0
    left_join(bin_lookup, by = "cwd_bin") %>%
    mutate(group = label,
           lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error) %>%
    arrange(bin_mid)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model A: pooled general response (all species, all terciles) -----------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod_pooled <- fit_binned(df)
summary(mod_pooled)
pooled_resp <- tidy_bins(mod_pooled, "All species")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model B: response by species-relative historic-CWD tercile ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# One pooled fit per species-relative tercile (each pools all species). Bins are
# shared/absolute, so the three curves sit on a common CWD axis and overlap.
terc_levels <- levels(site_terc$cwd_tercile)
tercile_resp <- map_dfr(terc_levels, function(g) {
  tidy_bins(fit_binned(df %>% filter(cwd_tercile == g)), g)
}) %>% mutate(group = factor(group, terc_levels))

# Common support: obs per bin x species-relative tercile
support <- df %>%
  count(cwd_tercile, cwd_bin) %>%
  pivot_wider(names_from = cwd_tercile, values_from = n, values_fill = 0) %>%
  left_join(bin_lookup %>% select(cwd_bin, bin_mid), by = "cwd_bin") %>%
  arrange(bin_mid)
cat("\n=== Common support: obs per bin x species-relative tercile ===\n")
print(support)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model C: continuous bin x within-species baseline-CWD rank (test) ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Single pooled model. cwd_ave_rank (within-species percentile of baseline CWD)
# is site-level => its main effect is absorbed by tree FE; the bin x rank
# interaction is identified. Positive coef on the highest-CWD bin x rank =>
# species-relatively-drier sites are LESS harmed at high CWD (adaptation).
fmla_int <- as.formula(paste0(
  "log(growth) ~ i(cwd_bin, cwd_ave_rank, ref = '", ref_bin, "') + i(cwd_bin, ref = '", ref_bin, "')",
  " + ns(pet.an, ", ns_df, ") | tree_id + ", age_fe))
mod_int <- feols(fmla_int, data = df, cluster = ~collection_id)

adapt_test <- tidy(mod_int, conf.int = TRUE) %>%
  filter(str_detect(term, "cwd_bin") & str_detect(term, "cwd_ave_rank")) %>%
  mutate(cwd_bin = str_remove(str_remove(term, "^cwd_bin::"), ":cwd_ave_rank$")) %>%
  left_join(bin_lookup, by = "cwd_bin") %>%
  arrange(bin_mid) %>%
  select(cwd_bin, bin_mid, estimate, std.error, conf.low, conf.high, p.value)
cat("\n=== Adaptation test: i(cwd_bin) x within-species baseline-CWD rank ===\n")
print(adapt_test)

# Predicted response curves at within-species P10 (wet) vs P90 (dry) baseline,
# reconstructed directly from coefficients: effect_b(rank) = beta_bin_b + gamma_b * rank
# (relative to the reference bin, so FE / controls drop out). The formal uncertainty
# is in adapt_test above.
int_coefs <- tidy(mod_int)
beta_bin <- int_coefs %>%
  filter(str_detect(term, "^cwd_bin::"), !str_detect(term, ":cwd_ave_rank")) %>%
  transmute(cwd_bin = str_remove(term, "^cwd_bin::"), beta = estimate)
gamma_bin <- int_coefs %>%
  filter(str_detect(term, ":cwd_ave_rank")) %>%
  transmute(cwd_bin = str_remove(str_remove(term, "^cwd_bin::"), ":cwd_ave_rank$"), gamma = estimate)
int_pred <- expand_grid(cwd_bin = bin_lookup$cwd_bin, rank = c(0.1, 0.9)) %>%
  mutate(cwd_bin = as.character(cwd_bin)) %>%
  left_join(beta_bin, by = "cwd_bin") %>%
  left_join(gamma_bin, by = "cwd_bin") %>%
  mutate(beta = replace_na(beta, 0), gamma = replace_na(gamma, 0),  # ref bin -> 0
         estimate = beta + gamma * rank) %>%
  left_join(bin_lookup, by = "cwd_bin") %>%
  mutate(baseline = factor(rank, levels = c(0.1, 0.9),
                           labels = c("P10 (wet for sp.)", "P90 (dry for sp.)"))) %>%
  arrange(bin_mid)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure -----------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
terc_cols <- c("Wet for sp." = "#0072B2", "Mid for sp." = "#999999", "Dry for sp." = "#D55E00")
baseline_cols <- c("P10 (wet for sp.)" = "#0072B2", "P90 (dry for sp.)" = "#D55E00")

## Panel A: pooled general response
pA <- ggplot(pooled_resp, aes(bin_mid, estimate, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.15) +
  geom_line() + geom_point() +
  labs(subtitle = "A. General response (all species pooled)", x = "Annual CWD (mm)",
       y = "log(RWL) effect vs. reference bin") +
  coord_cartesian(xlim = core_range) +
  theme_bw()

## Panel B: response by species-relative tercile (common CWD axis)
pB <- ggplot(tercile_resp, aes(bin_mid, estimate, color = group, fill = group,
                               ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = 0.12, color = NA) +
  geom_line() + geom_point(size = 1) +
  scale_color_manual(values = terc_cols) + scale_fill_manual(values = terc_cols) +
  labs(subtitle = "B. By species-relative historic-CWD tercile", x = "Annual CWD (mm)",
       y = "log(RWL) effect vs. reference bin", color = "Baseline CWD\n(within species)",
       fill = "Baseline CWD\n(within species)") +
  coord_cartesian(xlim = core_range) +
  theme_bw()

## Panel C: continuous interaction, predicted at within-species P10 vs P90 baseline
pC <- ggplot(int_pred, aes(bin_mid, estimate, color = baseline)) +
  geom_line() + geom_point(size = 1) +
  scale_color_manual(values = baseline_cols) +
  labs(subtitle = "C. Interaction model: predicted at within-species P10/P90 baseline",
       x = "Annual CWD (mm)", y = "Predicted log(RWL)", color = "Baseline CWD\n(within species)") +
  coord_cartesian(xlim = core_range) +
  theme_bw()

binned_fig <- (pA + pB) / pC +
  plot_annotation(
    title = paste0("Pooled multi-species RWL response to annual CWD (", data_source, ")"),
    subtitle = "log(RWL) step function of annual CWD; PET + ring age controlled; tree FE; within-species terciles",
    caption = paste0("Outcome: log(ring width). Reference bin: ", ref_raw, " mm. Bin width: ",
                     binwidth, " mm. ", n_distinct(df$species_id), " species, ",
                     n_distinct(df$collection_id), " sites."))
binned_fig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export -----------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list(pooled_response  = pooled_resp,
     tercile_response = tercile_resp,
     adaptation_test  = adapt_test,
     interaction_pred = int_pred,
     bin_lookup       = bin_lookup,
     common_support   = support) %>%
  write_rds(paste0(wdir, "2_output/first_stage/binned_cwd_response_pooled.rds"))
