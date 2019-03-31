install.packages("plm")

library(tidyverse)
library(plm)

### Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
setwd(wdir)

### Read data
spp <- "psme"
df <- read.csv(paste0(wdir,spp,'_data.csv'))
df <- df %>%
  mutate(tree_id = paste0(site_id, tree_id))

### Calculate deviation in growth rate
find_deviation <- function(df) {
  fit <- lm(ring_width ~ age, data = df)                   ## Replace with logistic growth function
  df$predict <- fit %>% predict(type = "response")
  df$deviation <- df$ring_width - df$predict 
  return(df$deviation)
  }

df$deviation <- find_deviation(df)

df$deviation <- df %>%                              ## Ask Fran about this - How to run function by site?
  group_by(site_id) %>%
  do(find_deviation)

### Calculate site average aet and cwd
site_clim <- df %>%
  group_by(site_id, year) %>%
  summarize(cwd.an = max(cwd.an, na.rm = TRUE), aet.an = max(aet.an, na.rm = TRUE)) %>%
  filter((cwd.an > -Inf) & (aet.an > -Inf)) %>%
  group_by(site_id) %>%
  summarize(cwd.ave = mean(cwd.an, na.rm = TRUE), aet.ave = mean(aet.an, na.rm = TRUE))

df <- merge(x = df, y = site_clim, by = "site_id", all.x = TRUE)

# Climate regressions
mod <- lm(deviation ~ cwd.an + aet.an, data = df)
summary(mod)

mod_squared <- lm(deviation ~ cwd.an + I(cwd.an^2) + aet.an + I(aet.an^2), data = df)
summary(mod_squared)


mod_interact <- lm(deviation ~ cwd.ave * cwd.an + cwd.ave * I(cwd.an^2), data = df)
summary(mod_interact)

mod_interact <- lm(deviation ~ cwd.ave * cwd.an, data = df)
summary(mod_interact)