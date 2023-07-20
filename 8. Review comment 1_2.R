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