#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/11/2023
# Comment: 
  # My major concern is around potential statistical artifacts that drive 
  # the results. Specifically, I’m concerned that areas that experience wetter 
  # conditions will inherently have experienced more variation in moisture 
  # conditions and thus more variation in ring width index. Trees in areas that 
  # experience more climate variability will appear to be more sensitive to climate 
  # just by virtue of having experienced more variability. Standardizing by the 
  # z-score makes sense in terms of facilitating comparisons across sites, but it 
  # also obscures differences in the magnitude of variability in water availability 
  # across sites, potentially causing some statistical artifacts. The manuscript 
  # could test for these statistical artifacts in a number of ways; first, by 
  # testing for associations between the site mean climate variables and the 
  # standard deviation of those variables, and second, by controlling for SD of 
  # the climate variables in the second-stage regression. 
#
# Approach:
# - Re-run first stage model using RE, Bayesian RE and FE structures. Compare coefficients.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Explore concerns --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Relationship between cwd and sd(cwd)
flm_df %>% 
  ggplot(aes(x = cwd.spstd, y = cwd.sd)) +
  geom_point(alpha = 0.2) +
  geom_smooth() +
  xlim(-5, 10) +
  ylim(0, 30) +
  ylab("Standard deviation of CWD") +
  xlab("Historic CWD\n(deviation from species mean)") +
  theme_bw()


### Relationship between beta and sd(cwd)
flm_df %>% 
  filter(outlier == FALSE) %>% 
  ggplot(aes(x = cwd.sd, y = estimate_cwd.an)) +
  geom_point(alpha = 0.2) +
  geom_smooth() +
  xlab("Standard deviation of CWD") +
  ylab("Sensitivity to CWD") +
  # ylim(-2,2) +
  xlim(0,5) +
  theme_bw()


trim_df <- flm_df %>% filter(outlier == FALSE)
mod <- feols(estimate_cwd.an ~ cwd.sd + cwd.spstd + (cwd.spstd^2) + pet.sd + pet.spstd + (pet.spstd^2), 
           data = trim_df, weights = trim_df$cwd_errorweights,
           vcov = conley(cutoff = 460, distance = "spherical"))
summary(mod)



mod1 <- lm(cwd.sd ~ cwd.spstd, data = trim_df)
summary(mod1)


mod2 <- feols(estimate_cwd.an ~ cwd.sd, trim_df)
summary(mod2)

library(modelsummary)
models <- list(
  "sd(CWD) ~ ave(CWD)" = mod1,
  "beta ~ sd(CWD)" = mod2
)
modelsummary(models, coef_rename = c("cwd.spstd" = "Average CWD", "cwd.sd" = "sd(CWD)"), stars = TRUE, gof_omit = 'F|DF|Deviance|R2|AIC|BIC|ICC|Log.Lik.|aicc|Std. errors')

# #### REVIEWER COMMENTS - DRIVEN BY VARIABILITY?
# sd_df <- dendro_df %>% 
#   group_by(collection_id) %>% 
#   summarise(cwd_sd = sd(cwd.an.spstd, na.rm = TRUE),
#             cwd_mean = mean(cwd.an.spstd, na.rm = TRUE))
# 
# mod <- lm(cwd_mean ~ cwd_sd, sd_df)
# summary(mod)
# sd_df %>% ggplot(aes(x = cwd_mean, y = cwd_sd)) +
#   geom_point()
# 
# 
# #### REVIEWER COMMENTS - DRIVEN BY VARIABILITY?