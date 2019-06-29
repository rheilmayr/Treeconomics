library(tidyverse)
library(lfe)
library(broom)
library(ggiraphExtra)

### Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
#for Fran
wdir="C:/Users/fmoor/Google Drive/Treeconomics/Data/"
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
  summarise(deviation = find_deviation(.))
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
mod_interact <- felm(ring_width ~ cwd.ave * cwd.an + cwd.ave * I(cwd.an^2)+age*site_id+I(age^2)*site_id+I(age^3)*site_id-site_id-cwd.ave|site_id|0|site_id, data = df )
summary(mod_interact)

mod_interact <- felm(deviation ~ aet.ave * aet.an + aet.ave * I(aet.an^2)-aet.ave|site_id|0|site_id, data = df )
summary(mod_interact)

complete_df <- df[which(complete.cases(df[,c(7,10)])),]

ntrees=complete_df%>%
  group_by(site_id)%>%
  summarize(ntrees=length(unique(tree_id)))
complete_df=merge(complete_df,ntrees)
complete_df=complete_df%>%
  filter(ntrees>5)

site_lm <- complete_df %>% 
  group_by(site_id) %>%
  do(fit_site = lm(ring_width ~ cwd.an+age+I(age^2)+I(age^3)+I(age^4)+tree_id, data = . ))
SiteCoef=tidy(site_lm[[2]][[1]]) %>%
  filter(term=="cwd.an")
for(i in 2:length(site_lm[[2]])){
  SiteCoef=rbind(SiteCoef,tidy(site_lm[[2]][[i]])%>%filter(term=="cwd.an"))
}
SiteCoef$site_id=site_lm$site_id
siteCoef=merge(SiteCoef,unique(complete_df[,c(1,10,11)]),all.x=T,all.y=F)
siteCoef$weights=1/siteCoef$std.error/(sum(1/siteCoef$std.error))
grandmodel=lm(estimate~cwd.ave+I(cwd.ave^2),weights=weights,data=siteCoef)

# Explore regression results
x = 0:1000
y=predict(grandmodel,data.frame(cwd.ave=x),interval="prediction",level=0.95)
x11()
par(mar=c(5.1,5.1,4.1,2.1),cex=1.5)
plot(x,y[,1],type="l",xlab="Average Climatic Water Deficit",ylab="",las=1,lwd=2,col="#0dcfca",ylim=c(-11e-4,0))
polygon(c(x,rev(x)),c(y[,2],rev(y[,3])),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)
title(ylab="Effect of Water Deficit on Tree Growth", line=4)
abline(h=0,lty=2,lwd=2)



SiteCoef <- merge(x = SiteCoef, y = site_clim, by = "site_id", all.x = TRUE)
x11()
par(mar=c(5.1,5.1,4.1,2.1),cex=1.5)
plot(SiteCoef$cwd.ave, SiteCoef$estimate,xlab="Average Climatic Water Deficit",ylab="Growth sensitivity to CWD",ylim=c(-11e-3,11e-3), xlim=c(1,1e3))

ggplot(site_clim, aes(x=cwd.ave)) + geom_histogram()
