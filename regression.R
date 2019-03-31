install.packages("plm")
install.packages("lfe")

library(tidyverse)
library(lfe)
library(broom)

library(ggiraphExtra)
### Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
setwd(wdir)

### Read data
spp <- "psme"
df <- read.csv(paste0(wdir,spp,'_data.csv'))
df <- df %>%
  mutate(tree_id = paste0(site_id, tree_id))

### Calculate deviation in growth rate                     ## Alternative is polynomial (logitic) age control interacted with site in single stage regression
find_deviation <- function(df) {
  fit <- lm(ring_width ~ age, data = df)                   ## Replace with logistic growth function
  df$predict <- fit %>% predict(type = "response")
  df$deviation <- df$ring_width - df$predict 
  return(df$deviation)
  }

df$deviation <- find_deviation(df)

df %>%                              ## Rewrite to run function by site
  group_by(site_id) %>%
  summarise(deviation = find_deviation())
  do(site_dev = find_deviation(.))

### Calculate site average aet and cwd
site_clim <- df %>%
  filter(year < 1980) %>%                                                                      ## Switch to raw climate data rather than using tree dataframe
  group_by(site_id, year) %>%
  summarize(cwd.an = max(cwd.an, na.rm = TRUE), aet.an = max(aet.an, na.rm = TRUE)) %>%
  filter((cwd.an > -Inf) & (aet.an > -Inf)) %>%
  group_by(site_id) %>%
  summarize(cwd.ave = mean(cwd.an, na.rm = TRUE), aet.ave = mean(aet.an, na.rm = TRUE))

df <- merge(x = df, y = site_clim, by = "site_id", all.x = TRUE)

### CWD data binned                                            ## Probably not necessary since our CWD/AET variables are already integrals
df <- cut(x, breaks = quantile(x, probs = seq(0, 1, 0.1)), 
          labels = 1:10, include.lowest = TRUE)

# Climate regressions
mod_interact <- felm(deviation ~ cwd.ave * cwd.an + cwd.ave * I(cwd.an^2)-cwd.ave|site_id|0|site_id, data = df )
summary(mod_interact)

mod_interact <- felm(deviation ~ aet.ave * aet.an + aet.ave * I(aet.an^2)-aet.ave|site_id|0|site_id, data = df )
summary(mod_interact)

complete_df <- df[which(complete.cases(df[,c(7,10)])),]

site_lm <- complete_df %>% 
  group_by(site_id) %>%
  do(fit_site = lm(deviation ~ cwd.an, data = . ))
siteCoef = tidy(site_lm, fit_site) %>%
  filter(term == "cwd.an")
siteCoef=merge(siteCoef,unique(complete_df[,c(1,11)]),all.x=T,all.y=F)
grandmodel=lm(estimate~cwd.ave+I(cwd.ave^2),data=siteCoef)
grandmodel=lm(estimate~cwd.ave*I(cwd.ave>250)-I(cwd.ave>250),data=siteCoef)

# Explore regression results
x = 0:1000
plot(x,x*mod_interact$coefficients[1]+x^2*mod_interact$coefficients[2])
waldtest(mod_interact,c(F,F,T,T))