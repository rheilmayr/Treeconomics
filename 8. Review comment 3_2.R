#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/11/2023
# Review comment 3.3: 
#   The error distributions in equation 1 and 2 need to be writen out. It does 
#   not appear that they factor in spatial or temporal dependence in the errors 
#   in equation 1 nor spatial dependence in the errors of equation 2. Given the 
#   maps and time series presented, I assume there must be substantial residual 
#   dependence. This can be incorporated through usual spatio-temporal statistical methods.   
# 
#
# Approach:
# - 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(tidylog)
library(tseries)
library(dtplyr)
library(furrr)
library(nlme)
library(gstat)
library(sf)
library(fixest)


n_cores <- availableCores() - 6
future::plan(multisession, workers = n_cores)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import and integrate data (same as first stage) ------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'


# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out/dendro/")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            rwl = mean(rwl),
            rwi_ar = mean(rwi_ar),
            rwi_nb = mean(rwi_nb),
            .groups = "drop") %>% 
  as_tibble()


# 2. Historic site-level climate
an_site_clim <- read_rds(paste0(wdir, "out/climate/site_an_clim.gz"))
dendro_df <- dendro_df %>% 
  left_join(an_site_clim, by = c("collection_id", "year"))


# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)

dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore temporal autocorrelation in first stage data ------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


test <- dendro_df %>% 
  select(rwl) %>% 
  drop_na() %>% 
  acf()

$acf

dendro_df %>% 
  select(rwi) %>% 
  drop_na() %>% 
  acf()

dendro_df %>% 
  select(rwi_ar) %>% 
  drop_na() %>% 
  acf()



test_df <- dendro_df[1:1000,]

ts_df <- dendro_df %>% 
  group_by(collection_id, tree) %>% 
  nest() %>% 
  mutate(acf = future_map(data, calc_acf))

ts_df <- ts_df %>% 
  select(-data) %>% 
  unnest(acf)

ts_df %>% 
  ggplot(aes(group = lag, y = rwl_corr)) +
  geom_boxplot() +
  ylim(-0.5, 1)

ts_df %>% 
  ggplot(aes(group = lag, y = rwi_corr)) +
  geom_boxplot() +
  ylim(-0.5, 1)

ts_df %>% 
  ggplot(aes(group = lag, y = rwiar_corr)) +
  geom_boxplot() +
  ylim(-0.5, 1)




fs_mod <- function(site_data, outcome = "rwi", energy_var = "pet.an"){
  rwl_mean <- site_data$rwl %>% mean(na.rm = TRUE)
  rwl_sd <- site_data$rwl %>% sd(na.rm = TRUE)
  
  failed <- F
  reg_error <- NA
  femod <- NA
  pet_cwd_cov <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
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
  return(tibble(mod = list(mod), nobs = nobs, ntrees = ntrees, rwl_mean = rwl_mean, rwl_sd = rwl_sd, error = reg_error))
}



install.packages("nlme")
library(nlme)
library(lme4)
library(lmtest)
n_idx <- 900
test_data <- dendro_df %>% 
  drop_na() %>% 
  group_by(collection_id) %>% 
  nest()
test_data <- test_data[n_idx,2][[1]][[1]]

formula <- as.formula(paste0("rwl ~ pet.an.spstd + cwd.an.spstd"))
mod <- lm(formula, data = test_data)
acf(residuals(mod))
dwtest(formula = mod,  alternative = "two.sided")

formula <- as.formula(paste0("rwi ~ pet.an.spstd + cwd.an.spstd"))
mod <- lm(formula, data = test_data)
acf(residuals(mod))
dwtest(formula = mod,  alternative = "two.sided")

formula <- as.formula(paste0("rwi_ar ~ pet.an.spstd + cwd.an.spstd"))
mod <- lm(formula, data = test_data)
acf(residuals(mod))
dwtest(formula = mod,  alternative = "two.sided")

formula <- as.formula(paste0("rwi ~ pet.an.spstd + cwd.an.spstd"))
mod <- lm(formula, data=test_data)
summary(mod)

mod <- nlme::lme(formula,
                data=test_data, method="REML",
                random = ~ 1 | tree)
summary(mod)
plot(nlme::ACF(mod))
mod <- nlme::lme(formula,
                 data=test_data, method="REML",
                 random = ~ 1 | tree,
                 correlation = nlme::corAR1(form=~year|tree))
summary(mod)
plot(nlme::ACF(mod, resType = "normalized"))
dwtest(formula = mod,  alternative = "two.sided")


formula <- as.formula(paste0("rwi_ar ~ pet.an.spstd + cwd.an.spstd"))
mod <- nlme::lme(formula,
                 data=test_data, method="REML",
                 random = ~ 1 + year | tree)
summary(mod)
plot(nlme::ACF(mod, resType = "normalized"))


library(DHARMa)

library(lmtest)
lmtest::dwtest(model)




## Plot autocorrelation figures across plots
calc_acf <- function(site_data, nlags = 20){
  if (dim(data)[1] < 30){
    return(NULL)
  } else{
    mod_rwl <- lm("rwl ~ pet.an.spstd + cwd.an.spstd", data = site_data)
    acf_rwl <- acf(residuals(mod_rwl), lag.max = nlags)
    
    mod_rwi <- lm("rwi ~ pet.an.spstd + cwd.an.spstd", data = site_data)
    acf_rwi <- acf(residuals(mod_rwi), lag.max = nlags)
    
    mod_rwiar <- lm("rwi_ar ~ pet.an.spstd + cwd.an.spstd", data = site_data)
    acf_rwiar <- acf(residuals(mod_rwiar), lag.max = nlags)

    acf_df <- tibble(lag = seq(0, nlags), rwi_corr = acf_rwi$acf, rwl_corr = acf_rwl$acf, rwiar_corr = acf_rwiar$acf)
    return(acf_df)    
  }
}

test_df <- dendro_df %>% 
  drop_na() %>% 
  nest_by(collection_id) %>%
  ungroup() %>% 
  slice_sample(n = 300)

site_data <- test_df[1,2][[1]][[1]]

test_df <- test_df %>% 
  mutate(acf = map(data, calc_acf)) %>% 
  unnest(acf)

test_df %>% 
  ggplot(aes(group = lag, y = rwl_corr)) +
  geom_boxplot() +
  theme_bw() +
  ylim(-.4, 1)

test_df %>% 
  ggplot(aes(group = lag, y = rwi_corr)) +
  geom_boxplot() +
  theme_bw() +
  ylim(-.4, 1)


test_df %>% 
  ggplot(aes(group = lag, y = rwiar_corr)) +
  geom_boxplot() +
  theme_bw() +
  ylim(-.4, 1)



## Compare spl vs ar first stage results
temp <- fs_spl %>% 
  left_join(fs_ar %>% select(collection_id, estimate_cwd.an_ar = estimate_cwd.an, std.error_cwd.an_ar = std.error_cwd.an), by = "collection_id")
temp %>% 
  ggplot(aes(x = estimate_cwd.an, y = estimate_cwd.an_ar)) +
  geom_point() +
  xlim(-2, 1) +
  ylim(-2, 1) +  
  geom_smooth(method='lm', formula= y~x, linetype="dashed") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed()

temp %>% 
  ggplot(aes(x = std.error_cwd.an, y = std.error_cwd.an_ar)) +
  geom_point() +
  xlim(0, 1) +
  ylim(0, 1) +  
  geom_smooth(method='lm', formula= y~x, linetype="dashed") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed()



library(fixest)
formula <- as.formula("estimate_cwd.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
data <- fs_spl %>% filter(outlier == 0)
mod <- feols(formula, weights = data$cwd_errorweights, data = data,
             vcov = conley(cutoff = 370, distance = "spherical"))

data <- fs_ar %>% filter(outlier == 0)
mod <- feols(formula, weights = data$cwd_errorweights, data = data,
             vcov = conley(cutoff = 370, distance = "spherical"))
summary(mod)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore spatial autocorrelation in second stage data ------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

# Trim outliers
trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()


# Spatial autocorrelation of trim_df outcome variable
site_points=st_as_sf(trim_df,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
plot(vg, vg.fit)
vg.range_km = vg.fit[2,3]
vg.range_deg <- vg.range_km / 100 # Convert from km to degrees. Approximate since this changes across latitudes. But don't want to project


# Conley SE model
formula <- as.formula(estimate_cwd.an ~ cwd.spstd + cwd.spstd**2 + pet.spstd + pet.spstd**2)
bl_mod <- feols(formula,data=trim_df, weights = ~I(cwd_errorweights))
summary(bl_mod)

con_mod <- feols(formula,data=trim_df, weights = ~I(cwd_errorweights), vcov = conley(cutoff = vg.range_km, distance = "spherical"))
summary(con_mod)


# LME spatial autocorrelation model
formula <- as.formula(estimate_cwd.an ~ cwd.spstd + I(cwd.spstd**2) + pet.spstd + I(pet.spstd**2))
spat_re_mod <- nlme::lme(formula, data = mod_df,
                         weights = ~I(1/cwd_errorweights),
                         random = ~ 1 | dummy)
spat_re_mod <- spat_re_mod %>% 
  update(correlation = nlme::corExp(value = vg.range_deg, form = ~ longitude + latitude), method = "ML")
summary(spat_re_mod)


models <- list(
  "No correction" = bl_mod,
  "Conley standard errors"     = con_mod,
  "Spatially correlated errors" = spat_re_mod
)

modelsummary(models, gof_omit = 'DF|Deviance|R2|AIC|BIC|ICC|Log.Lik.|aicc|Std. errors')




