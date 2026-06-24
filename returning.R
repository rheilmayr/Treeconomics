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
# mod_df <- mod_df %>%
#   mutate(ppt.an = cwd.an,
#          ppt.ave = cwd.ave)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run model --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formula = "rwi ~ cwd.an + pet.an | collection_id"

# formula = "rwi ~ cwd.an.spstd:pet.spstd + cwd.an.spstd:I(pet.spstd**2) +
#                cwd.an.spstd:cwd.spstd + cwd.an.spstd:I(cwd.spstd**2) +
#                pet.an.spstd:pet.spstd + pet.an.spstd:I(pet.spstd**2) +
#                pet.an.spstd:cwd.spstd + pet.an.spstd:I(cwd.spstd**2) +
#                                     (1 + cwd.an.spstd + pet.an.spstd | collection_id)"

# cwd_median <- mod_df %>% select(cwd.spstd) %>% drop_na() %>% pull(cwd.spstd) %>% median()
# formula = "rwi ~ temp.an + poly(ppt.an, 2) | collection_id"
mod <- feols(as.formula(formula), data = mod_df)
summary(mod)



formula = "rwi ~ cwd.an + pet.an | collection_id"
mod <- feols(as.formula(formula), data = mod_df)
summary(mod)

formula = "rwi ~ ns(pet.an, df = 4) + ns(cwd.an, df = 4) | collection_id"
mod_ns <- feols(as.formula(formula), data = mod_df)
summary(mod_ns)
p1 <- plot_predictions(mod_ns,
                       condition = "cwd.an",
                       vcov = FALSE) +
  labs(title = "Predicted RWI vs CWD",
       subtitle = "Natural splines (df = 4) with site fixed effects",
       x = "CWD (mm)",
       y = "Predicted RWI") +
  theme_minimal()
p1

p2 <- plot_predictions(mod_ns,
                       condition = "pet.an",
                       vcov = FALSE) +
  labs(title = "Predicted RWI vs PET",
       subtitle = "Natural splines (df = 4) with site fixed effects",
       x = "PET (mm)",
       y = "Predicted RWI") +
  theme_minimal()
p2


formula = "rwi ~ ns(temp.an, df = 4) + ns(ppt.an, df = 4) | collection_id"
mod_ns <- feols(as.formula(formula), data = mod_df)
summary(mod_ns)
p3 <- plot_predictions(mod_ns,
                       condition = "ppt.an",
                       vcov = FALSE) +
  labs(title = "Predicted RWI vs PPT",
       subtitle = "Natural splines (df = 4) with site fixed effects",
       x = "PPT (mm)",
       y = "Predicted RWI") +
  theme_minimal()
p3

p4 <- plot_predictions(mod_ns,
                       condition = "temp.an",
                       vcov = FALSE) +
  labs(title = "Predicted RWI vs temp",
       subtitle = "Natural splines (df = 4) with site fixed effects",
       x = "Temp (degrees)",
       y = "Predicted RWI") +
  theme_minimal()
p4

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Natural splines for non-parametric CWD effects ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Estimate model with natural splines (4 df for flexibility)
formula = "rwi ~ pet.an + ns(cwd.an, df = 4) | collection_id"
mod_ns <- feols(as.formula(formula), data = mod_df)
summary(mod_ns)

# 1. Prediction plot: Shows the overall shape of the relationship
p1 <- plot_predictions(mod_ns,
                       condition = "cwd.an",
                       vcov = FALSE) +
  labs(title = "Predicted RWI vs Annual CWD",
       subtitle = "Natural splines (df = 4) with site fixed effects",
       x = "Annual CWD (mm)",
       y = "Predicted RWI") +
  theme_minimal()

# 2. Marginal effects plot: Shows how the SLOPE changes with CWD
# This reveals where CWD has the strongest/weakest effects
p2 <- plot_slopes(mod_ns,
                  variables = "cwd.an",
                  condition = "cwd.an") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Marginal Effect of CWD on RWI",
       subtitle = "How does the growth-CWD slope vary? (with 95% CI)",
       x = "Annual CWD (mm)",
       y = "Marginal Effect (∂RWI/∂CWD)") +
  theme_minimal()

# Display combined plot
p1 + p2 + plot_annotation(
  title = "Non-linear CWD Effects on Tree Growth",
  subtitle = "Douglas-fir (PSME) from FIA data"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Test Local Adaptation: Do drought-adapted sites show greater resilience?
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hypothesis: Sites with higher historic CWD (cwd.ave) should show
# less severe growth declines during dry years (local adaptation to drought)

# 1. Create categorical variable for historic CWD
# Split into dry vs wet sites based on median
mod_df <- mod_df %>%
  mutate(cwd_ave_cat = ifelse(cwd.ave >= median(cwd.ave, na.rm = TRUE),
                               "High CWD Historic Climate",
                               "Low CWD Historic Climate"),
         cwd_ave_cat = factor(cwd_ave_cat, levels = c("Low CWD Historic Climate", "High CWD Historic Climate")))

# Check the split
cat("\n=== Historic Climate Categories ===\n")
table(mod_df$cwd_ave_cat) %>% print()
cat(sprintf("\nLow CWD sites: cwd.ave < %.0f mm\n", median(mod_df$cwd.ave, na.rm = TRUE)))
cat(sprintf("High CWD sites: cwd.ave >= %.0f mm\n", median(mod_df$cwd.ave, na.rm = TRUE)))

# Estimate interaction model with categorical moderator
formula = "rwi ~ pet.an + ns(cwd.an, df = 4) * cwd_ave_cat | collection_id"
mod_interact <- feols(as.formula(formula), data = mod_df)
summary(mod_interact)

# Check data ranges
cat("\n=== Data Summary ===\n")
cat(sprintf("Annual CWD range: %.0f - %.0f mm\n",
            min(mod_df$cwd.an, na.rm = TRUE), max(mod_df$cwd.an, na.rm = TRUE)))
cat(sprintf("Historic avg CWD range: %.0f - %.0f mm\n",
            min(mod_df$cwd.ave, na.rm = TRUE), max(mod_df$cwd.ave, na.rm = TRUE)))
cat(sprintf("Correlation between cwd.an and cwd.ave: %.3f\n",
            cor(mod_df$cwd.an, mod_df$cwd.ave, use = "complete.obs")))


# 2. Create group-specific prediction grids (avoid extrapolation)
# Only predict within the range each group has actually experienced
low_cwd_range <- mod_df %>%
  filter(cwd_ave_cat == "Low CWD Historic Climate") %>%
  summarise(min_cwd = min(cwd.an, na.rm = TRUE),
            max_cwd = max(cwd.an, na.rm = TRUE))

high_cwd_range <- mod_df %>%
  filter(cwd_ave_cat == "High CWD Historic Climate") %>%
  summarise(min_cwd = min(cwd.an, na.rm = TRUE),
            max_cwd = max(cwd.an, na.rm = TRUE))

cat("\n=== Group-specific CWD ranges ===\n")
cat(sprintf("Low CWD sites: %.0f - %.0f mm\n", low_cwd_range$min_cwd, low_cwd_range$max_cwd))
cat(sprintf("High CWD sites: %.0f - %.0f mm\n", high_cwd_range$min_cwd, high_cwd_range$max_cwd))

# Create separate grids for each group
low_cwd_grid <- datagrid(
  model = mod_interact,
  cwd.an = seq(low_cwd_range$min_cwd, low_cwd_range$max_cwd, length.out = 100),
  cwd_ave_cat = "Low CWD Historic Climate"
)

high_cwd_grid <- datagrid(
  model = mod_interact,
  cwd.an = seq(high_cwd_range$min_cwd, high_cwd_range$max_cwd, length.out = 100),
  cwd_ave_cat = "High CWD Historic Climate"
)

combined_grid <- bind_rows(low_cwd_grid, high_cwd_grid)


# 3. Visualization: Predicted growth curves (no extrapolation)
# Get predictions manually, then plot
preds_df <- predictions(mod_interact,
                        newdata = combined_grid,
                        vcov = FALSE)

p_adapt1 <- ggplot(preds_df, aes(x = cwd.an, y = estimate, color = cwd_ave_cat, fill = cwd_ave_cat)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Low CWD Historic Climate" = "#0072B2",
                                 "High CWD Historic Climate" = "#D55E00")) +
  scale_fill_manual(values = c("Low CWD Historic Climate" = "#0072B2",
                                "High CWD Historic Climate" = "#D55E00")) +
  labs(title = "Growth Response to CWD Varies by Historic Climate",
       subtitle = "Only showing predictions within observed range for each group",
       x = "Annual CWD (mm)",
       y = "Predicted RWI",
       color = "Site Type",
       fill = "Site Type") +
  theme_minimal() +
  theme(legend.position = "right")


# 4. KEY PLOT: Marginal effects (no extrapolation)
# Get slopes manually, then plot
slopes_df <- slopes(mod_interact,
                    variables = "cwd.an",
                    newdata = combined_grid,
                    vcov = FALSE)

p_adapt2 <- ggplot(slopes_df, aes(x = cwd.an, y = estimate, color = cwd_ave_cat, fill = cwd_ave_cat)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Low CWD Historic Climate" = "#0072B2",
                                 "High CWD Historic Climate" = "#D55E00")) +
  scale_fill_manual(values = c("Low CWD Historic Climate" = "#0072B2",
                                "High CWD Historic Climate" = "#D55E00")) +
  labs(title = "Marginal Effect of CWD: Low vs High Historic CWD Sites",
       subtitle = "Hypothesis: High-CWD-adapted sites show less growth decline at high CWD",
       x = "Annual CWD (mm)",
       y = "Marginal Effect (∂RWI/∂CWD)",
       color = "Site Type",
       fill = "Site Type") +
  annotate("text", x = min(slopes_df$cwd.an, na.rm = TRUE),
           y = Inf, vjust = 1.5, hjust = 0,
           label = "Look for: Orange line (high CWD sites) ABOVE blue line at high CWD",
           size = 3, fontface = "italic", color = "gray30") +
  theme_minimal() +
  theme(legend.position = "right")


# 5. Display plots
cat("\n=== Generating Visualizations ===\n")
p_adapt1
p_adapt2

# Combined view
(p_adapt1 / p_adapt2) +
  plot_annotation(
    title = "Local Adaptation to Historic Climate in Douglas-fir (PSME)",
    subtitle = "Testing if drought-adapted populations show greater drought resilience"
  )


# 6. Statistical interpretation helper
cat("\n=== Interpretation Guide ===\n")
cat("If LOCAL ADAPTATION exists:\n")
cat("  - Plot 1: Orange curve (high CWD sites) should be ABOVE blue curve at HIGH annual CWD\n")
cat("  - Plot 2: Orange line should have SHALLOWER (less negative) slope at high cwd.an\n")
cat("  - Meaning: Drought-adapted trees maintain better growth under high water deficit\n\n")

cat("If NO ADAPTATION (universal constraint):\n")
cat("  - Plot 1: Curves should be parallel (no divergence)\n")
cat("  - Plot 2: Orange and blue lines should overlap\n")
cat("  - Meaning: All populations respond identically to drought\n\n")

cat("Color coding:\n")
cat("  - Orange = High historic CWD (above median cwd.ave)\n")
cat("  - Blue = Low historic CWD (below median cwd.ave)\n\n")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Two-stage analysis: site-level sensitivities ~ historic climate ---------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Stage 1: fit rwi ~ cwd.an + pet.an for each site -------------------------
site_sens <- mod_df %>%
  group_by(collection_id) %>%
  filter(n() >= 10) %>%  # require at least 10 years per site
  nest() %>%
  mutate(
    mod    = map(data, ~lm(rwi ~ cwd.an + pet.an, data = .x)),
    tidied = map(mod, broom::tidy),
    n_obs  = map_int(data, nrow),
    r2     = map_dbl(mod, ~summary(.x)$r.squared)
  ) %>%
  unnest(tidied) %>%
  filter(term %in% c("cwd.an", "pet.an")) %>%
  select(collection_id, term, estimate, std.error, n_obs, r2)

# Join historic climate from tc_ave (keyed on plot_cn)
site_key <- mod_df %>%
  distinct(collection_id, plot_cn)

site_sens <- site_sens %>%
  left_join(site_key, by = "collection_id") %>%
  left_join(tc_ave, by = c("plot_cn" = "site_id"))


# Stage 1 diagnostics -------------------------------------------------------

# 1a. Distribution of observations per site
ggplot(site_sens %>% distinct(collection_id, n_obs),
       aes(x = n_obs)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Years of data per site",
       x = "Number of observations", y = "Count") +
  theme_minimal()

# 1b. Distribution of within-site R²
ggplot(site_sens %>% distinct(collection_id, r2),
       aes(x = r2)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Within-site model R² (rwi ~ cwd.an + pet.an)",
       x = "R²", y = "Count") +
  theme_minimal()

# 1c. Distribution of site-level coefficients
ggplot(site_sens, aes(x = estimate)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~term, scales = "free") +
  labs(title = "Distribution of site-level sensitivities",
       x = "Coefficient estimate", y = "Count") +
  theme_minimal()

# 1d. Estimate vs. std.error (flag unreliable sites)
ggplot(site_sens, aes(x = std.error, y = estimate)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~term, scales = "free") +
  labs(title = "Coefficient estimate vs. uncertainty",
       subtitle = "Points far right have high uncertainty",
       x = "Std. error", y = "Estimate") +
  theme_minimal()


# Stage 2: sensitivity ~ historic climate -----------------------------------
sens_cwd <- site_sens %>% filter(term == "cwd.an")
sens_pet <- site_sens %>% filter(term == "pet.an")

mod2_cwd <- lm(estimate ~ cwd.ave + pet.ave, data = sens_cwd,
               weights = 1 / std.error^2)
mod2_pet <- lm(estimate ~ cwd.ave + pet.ave, data = sens_pet,
               weights = 1 / std.error^2)

cat("\n=== Stage 2: CWD sensitivity ~ historic climate ===\n")
summary(mod2_cwd)

cat("\n=== Stage 2: PET sensitivity ~ historic climate ===\n")
summary(mod2_pet)


# Stage 2 diagnostics -------------------------------------------------------

# 2a. β_cwd ~ historic CWD
p2a <- ggplot(sens_cwd, aes(x = cwd.ave, y = estimate)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#D55E00") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "CWD sensitivity vs. historic CWD",
       x = "Historic CWD (mm)", y = "β_cwd") +
  theme_minimal()

# 2b. β_cwd ~ historic PET
p2b <- ggplot(sens_cwd, aes(x = pet.ave, y = estimate)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#D55E00") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "CWD sensitivity vs. historic PET",
       x = "Historic PET (mm)", y = "β_cwd") +
  theme_minimal()

# 2c. β_pet ~ historic CWD
p2c <- ggplot(sens_pet, aes(x = cwd.ave, y = estimate)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "PET sensitivity vs. historic CWD",
       x = "Historic CWD (mm)", y = "β_pet") +
  theme_minimal()

# 2d. β_pet ~ historic PET
p2d <- ggplot(sens_pet, aes(x = pet.ave, y = estimate)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "PET sensitivity vs. historic PET",
       x = "Historic PET (mm)", y = "β_pet") +
  theme_minimal()

(p2a + p2b) / (p2c + p2d) +
  plot_annotation(
    title = "Stage 2: How do site-level sensitivities vary with historic climate?",
    subtitle = "Douglas-fir (PSME)"
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Segmented regression: CWD threshold analysis ----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(segmented)

# Stage 1: fit segmented regression per site --------------------------------

# Diagnostic: run on one site to inspect errors before full run
test_site <- mod_df %>%
  group_by(collection_id) %>%
  filter(n() >= 20) %>%
  ungroup() %>%
  filter(collection_id == first(collection_id))

base_test <- lm(rwi ~ pet.an + cwd.an, data = test_site)
seg_test  <- segmented(base_test, seg.Z = ~cwd.an,
                       psi = median(test_site$cwd.an),
                       control = seg.control(it.max = 50, n.boot = 20))
summary(seg_test)

# Inspect actual structure to confirm column names before running at scale
cat("=== seg_test$psi structure ===\n")
print(seg_test$psi)
cat("\nColumn names of psi:", colnames(seg_test$psi), "\n")

cat("\n=== slope(seg_test)$cwd.an structure ===\n")
slopes_test <- slope(seg_test)$cwd.an
print(slopes_test)
cat("\nColumn names of slope output:", colnames(slopes_test), "\n")

fit_segmented <- function(df) {
  tryCatch({
    # Require meaningful CWD variation; skip near-constant sites
    if (sd(df$cwd.an, na.rm = TRUE) < 10) return(list(converged = FALSE, err = "low variance"))
    base <- lm(rwi ~ pet.an + cwd.an, data = df)
    seg  <- segmented(base, seg.Z = ~cwd.an,
                      psi = median(df$cwd.an),               # explicit start at median
                      control = seg.control(it.max = 50, n.boot = 20))
    slopes_out <- slope(seg)[["cwd.an"]]  # [[ ]] is safer than $ for names with dots
    bp    <- seg$psi[, "Est."]
    bp_se <- seg$psi[, "St.Err"]   # note: psi uses "St.Err" (no dot), slope() uses "St.Err."
    list(
      bp        = bp,
      bp_se     = bp_se,
      beta_lo   = slopes_out[1, "Est."],
      beta_hi   = slopes_out[2, "Est."],
      delta     = slopes_out[2, "Est."] - slopes_out[1, "Est."],
      delta_se  = sqrt(slopes_out[1, "St.Err."]^2 + slopes_out[2, "St.Err."]^2),
      converged = TRUE,
      err       = NA_character_
    )
  }, error = function(e) list(converged = FALSE, err = conditionMessage(e)))
}

seg_results <- mod_df %>%
  group_by(collection_id) %>%
  filter(n() >= 20) %>%
  nest() %>%
  mutate(result = map(data, fit_segmented)) %>%
  mutate(
    converged = map_lgl(result, "converged"),
    err_msg   = map_chr(result, ~.x$err %||% NA_character_),
    n_obs     = map_int(data, nrow)
  )

# Inspect failure reasons before proceeding
cat("\n=== Convergence diagnostics ===\n")
seg_results %>%
  filter(!converged) %>%
  count(err_msg, sort = TRUE) %>%
  print()

cat(sprintf("\nSegmented regression: %d sites attempted, %d converged\n",
            nrow(seg_results),
            sum(seg_results$converged)))

seg_results <- seg_results %>%
  filter(converged) %>%
  mutate(
    bp       = map_dbl(result, "bp"),
    bp_se    = map_dbl(result, "bp_se"),
    beta_lo  = map_dbl(result, "beta_lo"),
    beta_hi  = map_dbl(result, "beta_hi"),
    delta    = map_dbl(result, "delta"),
    delta_se = map_dbl(result, "delta_se")
  ) %>%
  select(-data, -result)

# Flag sites matching expected pattern (slope steepens above threshold) -----
seg_results <- seg_results %>%
  mutate(
    expected_pattern = delta < 0,
    delta_z          = delta / delta_se,
    sig_pattern      = expected_pattern & abs(delta_z) > 1.96
  )

cat("\n=== Pattern summary (H1 check) ===\n")
cat(sprintf("Sites with steeper slope above threshold (expected): %d / %d (%.0f%%)\n",
            sum(seg_results$expected_pattern), nrow(seg_results),
            100 * mean(seg_results$expected_pattern)))
cat(sprintf("Sites with significant steepening (p < 0.05):        %d / %d (%.0f%%)\n",
            sum(seg_results$sig_pattern), nrow(seg_results),
            100 * mean(seg_results$sig_pattern)))

# Join historic climate from tc_ave
site_key <- mod_df %>% distinct(collection_id, plot_cn)

seg_results <- seg_results %>%
  left_join(site_key, by = "collection_id") %>%
  left_join(tc_ave, by = c("plot_cn" = "site_id"))


# Stage 1 diagnostics -------------------------------------------------------

# 1a. Histogram of breakpoints
ggplot(seg_results, aes(x = bp)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Distribution of CWD breakpoints (ψ) across sites",
       x = "Breakpoint CWD (mm)", y = "Count") +
  theme_minimal()

# 1b. β_lo vs β_hi — steepening confirmed if points fall below y=x line
ggplot(seg_results, aes(x = beta_lo, y = beta_hi, color = expected_pattern)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                     labels = c("TRUE" = "Expected (steeper above)", "FALSE" = "Opposite")) +
  labs(title = "Slope below vs. above breakpoint",
       subtitle = "Points below diagonal: slope steepens above threshold (expected pattern)",
       x = "β below threshold", y = "β above threshold", color = NULL) +
  theme_minimal()

# 1c. Histogram of δ (change in slope)
ggplot(seg_results, aes(x = delta)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Distribution of slope change at threshold (δ = β_hi − β_lo)",
       subtitle = "Negative values = slope becomes more negative above threshold (H1)",
       x = "δ (change in slope)", y = "Count") +
  theme_minimal()


# Stage 2: breakpoint ~ historic climate (sig_pattern sites only) -----------
seg_stage2 <- seg_results %>% filter(sig_pattern)

cat(sprintf("\nStage 2 sample: %d sites with significant expected pattern\n",
            nrow(seg_stage2)))

mod2_bp <- lm(bp ~ cwd.ave + pet.ave, data = seg_stage2,
              weights = 1 / bp_se^2)

cat("\n=== Stage 2: CWD threshold (ψ) ~ historic climate ===\n")
summary(mod2_bp)


# Stage 2 diagnostics -------------------------------------------------------

# 2a. ψ ~ historic CWD (key test of H2)
p_bp1 <- ggplot(seg_stage2, aes(x = cwd.ave, y = bp)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#D55E00") +
  labs(title = "CWD threshold vs. historic CWD",
       subtitle = "Positive slope → H2 supported: drier sites tolerate more CWD before rapid decline",
       x = "Historic CWD (mm)", y = "Threshold ψ (mm)") +
  theme_minimal()

# 2b. ψ ~ historic PET
p_bp2 <- ggplot(seg_stage2, aes(x = pet.ave, y = bp)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  labs(title = "CWD threshold vs. historic PET",
       x = "Historic PET (mm)", y = "Threshold ψ (mm)") +
  theme_minimal()

p_bp1 + p_bp2 +
  plot_annotation(
    title = "Stage 2: Does historic climate predict the CWD tolerance threshold?",
    subtitle = sprintf("n = %d sites with significant threshold (expected pattern)", nrow(seg_stage2))
  )
