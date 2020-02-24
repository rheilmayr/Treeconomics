#install.packages("plm")
#install.packages("lfe")

library(MASS)
library(tidyverse)
library(lfe)
library(broom)
library(purrr)
library(ggiraphExtra)
select <- dplyr::select


### Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
#for Fran
# wdir="C:/Users/fmoor/Google Drive/Treeconomics/Data/"
setwd(wdir)

### Plotting parameters
pal=colorRampPalette(c('white','blue','yellow','red','darkred'))

### Read and prep data  ####
# spp=c("psme","pipo","pisy","pcgl","pcab","pied")
# speciesnames=c("Douglas Fir","Ponderosa Pine","Scotch Pine",
#                "White Spruce","Norway Spruce","Colorado Pinyon")

df <- read.csv(paste0(wdir,'species_data.csv'))
df <- df %>%
  mutate(tree_id = paste0(site_id, tree_id))
df <- df %>%
  mutate(pet.an = cwd.an + aet.an)

### Calculate site historic average aet, cwd and pet
site_clim <- df %>%
  filter(year < 1980) %>%                                                                      ## Switch to raw climate data rather than using tree dataframe
  group_by(site_id, year) %>%
  summarize(cwd.an = max(cwd.an, na.rm = TRUE), 
            aet.an = max(aet.an, na.rm = TRUE), 
            pet.an = max(pet.an, na.rm = TRUE)) %>%
  filter((cwd.an > -Inf) & (aet.an > -Inf)) %>%
  group_by(site_id) %>%
  summarize(cwd.ave = mean(cwd.an, na.rm = TRUE), 
            aet.ave = mean(aet.an, na.rm = TRUE), 
            pet.ave = mean(pet.an, na.rm = TRUE))

df <- merge(x = df, y = site_clim, by = "site_id", all.x = TRUE)

#remove sites with negative ages
df_trim <- df %>%
  filter(age>0)

#remove trees with a very short record (<5 years)
treeobs=df_trim %>%
  group_by(tree_id) %>%
  summarize(nyears=n())
df_trim=left_join(df_trim,treeobs, by = 'tree_id')
df_trim = df_trim %>%
  filter(nyears>5)

#remove sites with few trees
ntrees=df_trim%>%
  group_by(site_id,species_id)%>%
  summarize(ntrees=length(unique(tree_id)))
df_trim=left_join(df_trim,ntrees, by = c("site_id", "species_id"))
df_trim=df_trim%>%
  filter(ntrees>5)

complete_df=df %>%
  drop_na(c("cwd.an","aet.an","ring_width"))

site_summary <- complete_df %>% 
  select("site_id", "species_id", "species_name", "cwd.ave", "pet.ave", "aet.ave") %>%
  distinct() %>%
  as_tibble()

#### Run first stage ####
site_lm <- complete_df %>% 
  group_by(site_id, species_id) %>%
  do(fit_site = felm(ring_width ~ cwd.an+age+I(age^2)+I(age^3)+I(age^4)|tree_id|0|0, data = . ))
  # do(fit_site = felm(ring_width ~ cwd.an +aet.an+age+I(age^2)+I(age^3)+I(age^4)|tree_id|0|0, data = . ))

siteCoef=tidy(site_lm[[3]][[1]])%>%filter(term=="cwd.an")
for(i in 2:length(site_lm[[3]])){
  siteCoef=rbind(siteCoef,tidy(site_lm[[3]][[i]])%>%filter(term=="cwd.an"))
}
siteCoef$site_id <- site_lm$site_id
siteCoef$species_id <- site_lm$species_id
siteCoef <- siteCoef %>%
  left_join(site_summary, by = c("site_id","species_id"))

#remove one super extreme outlier
siteCoef_trimmed=siteCoef %>%
  group_by(species_id) %>%
  mutate(qhigh=quantile(estimate,0.99,na.rm=T),
         qlow=quantile(estimate,0.01,na.rm=T))
siteCoef_trimmed=siteCoef_trimmed %>%
  filter(estimate>qlow & estimate<qhigh)

### TEMP - KEEP ALL TO KEEP JUNIPER
siteCoef_trimmed <- siteCoef

# Calculate limits for each species
xlim.sp=site_summary %>%
  group_by(species_name) %>%
  mutate(cwd.qhigh = quantile(cwd.ave,0.95),
         cwd.qlow = quantile(cwd.ave,0.05),
         cwd.qmed = quantile(cwd.ave,0.5),
         pet.qhigh = quantile(pet.ave, 0.95), 
         pet.qlow = quantile(pet.ave, 0.05),
         pet.qmed = quantile(pet.ave,0.5),
         aet.qhigh = quantile(aet.ave, 0.95), 
         aet.qlow = quantile(aet.ave, 0.05),
         aet.qmed = quantile(aet.ave,0.5)) %>%
  select(species_name,cwd.qhigh,cwd.qlow,cwd.qmed,pet.qhigh,pet.qlow,pet.qmed, aet.qhigh, aet.qlow, aet.qmed) %>%
  distinct()

xlim.all <- site_summary %>%
  summarize(cwd.qhigh = quantile(cwd.ave,0.95),
            cwd.qlow = quantile(cwd.ave,0.05),
            cwd.qmed = quantile(cwd.ave,0.5),
            pet.qhigh = quantile(pet.ave, 0.95), 
            pet.qlow = quantile(pet.ave, 0.05),
            pet.qmed = quantile(pet.ave,0.5),
            aet.qhigh = quantile(aet.ave, 0.95), 
            aet.qlow = quantile(aet.ave, 0.05),
            aet.qmed = quantile(aet.ave,0.5))

### Generate summary plots of sample distributions  ####
# a=ggplot(siteCoef_trimmed,aes(x=cwd.ave,y=estimate,col=species))+theme_bw()+geom_point()+facet_wrap(~species)
# a=a+labs(x="Average Climatic Water Deficit",y="Marginal Effect of CWD")+scale_color_discrete(guide=FALSE)

kd <- kde2d(x = site_summary$cwd.ave, y = site_summary$aet.ave, n = 25, lims = c(c(0,750),c(100, 600)))
filled.contour(kd, color.palette = pal)

x11()
par(mfrow=c(4,5),mar=c(5,5,4,2))
for(i in 1:dim(xlim.sp)[1]){
  spec_name = xlim.sp$species_name[i]
  dat = siteCoef_trimmed %>% filter(species_name==spec_name)
  lim = xlim.sp %>% filter(species_name==spec_name)
  # kd <- kde2d(x = dat$cwd.ave, y = dat$pet.ave, n = 25, lims = c(c(lim$cwd.qlow, lim$cwd.qhigh), c(lim$pet.qlow, lim$pet.qhigh)))
  kd <- kde2d(x = dat$cwd.ave, y = dat$aet.ave, n = 25, 
              lims = c(c(xlim.all$cwd.qlow,xlim.all$cwd.qhigh),c(xlim.all$aet.qlow, xlim.all$aet.qhigh)))
  contour(kd)
  # kde2d(x, y, h, n = 25, lims = c(range(x), range(y)))
}

#### Run second stage of model  ####
dat = siteCoef_trimmed %>% 
  ungroup() %>%
  drop_na(std.error) %>%
  mutate(errorweights=1/std.error/sum(1/std.error,na.rm=T)) 
grandmodel <- lm(estimate~cwd.ave+I(cwd.ave^2) + aet.ave + I(aet.ave^2),weights=errorweights,data= dat)
# grandmodel <- lm(estimate~cwd.ave+I(cwd.ave^2), weights=errorweights,data= dat)

cwd.x <- seq(xlim.all$cwd.qlow, xlim.all$cwd.qhigh, length.out=1000)
aet.med <- xlim.all$aet.qmed

pred <- predict(grandmodel,data.frame(cwd.ave=cwd.x, aet.ave = aet.med),interval="prediction",level=0.90)
pred <- cbind(pred,cwd.x)
pred <- pred %>% as_tibble()
plot(cwd.x, pred$fit, type="l",xlab="Historic CWD",ylab="Estimated effect of CWD",las=1,lwd=2,col="#0dcfca",
     ylim = c(pred$lwr %>% min(), pred$upr %>% max()))
polygon(c(cwd.x,rev(cwd.x)),c(pred$lwr,rev(pred$upr)),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)

# aet.x <- seq(xlim.all$aet.qlow, xlim.all$aet.qhigh, length.out=1000)
# cwd.med <- xlim.all$cwd.qmed
# pred <- predict(grandmodel,data.frame(cwd.ave=cwd.med, aet.ave = aet.x),interval="prediction",level=0.90)
# pred <- cbind(pred,aet.x)
# pred <- pred %>% as_tibble()
# plot(aet.x, pred$fit, type="l",xlab="Historic AET",ylab="Marginal effect of AET",las=1,lwd=2,col="#0dcfca")


#### Run species specific second stage ####
grandmodels=siteCoef_trimmed %>%
  group_by(species_name) %>%
  drop_na(std.error) %>%
  mutate(errorweights=1/std.error/sum(1/std.error,na.rm=T)) %>%
  do(speciesmod=lm(estimate~cwd.ave+I(cwd.ave^2),weights=errorweights,data= .))
  # do(speciesmod=lm(estimate~cwd.ave + I(cwd.ave^2) + aet.ave + I(aet.ave^2),weights=errorweights,data= .))

# Review results
grandmodels %>% filter(species_name=="blue oak") %>% pull(speciesmod)

y=list()
predictions=data.frame()
for(i in 1:dim(grandmodels)[1]){
  species=as.data.frame(grandmodels)[i,1]
  cwd.x=seq(xlim.sp$cwd.qlow[which(xlim.sp$species_name==species)],xlim.sp$cwd.qhigh[which(xlim.sp$species_name==species)],length.out=1000)
  aet.x=seq(xlim.sp$aet.qlow[which(xlim.sp$species_name==species)],xlim.sp$aet.qhigh[which(xlim.sp$species_name==species)],length.out=1000)
  prd=predict(grandmodels[[2]][[i]],data.frame(cwd.ave=cwd.x),interval="prediction",level=0.90)
  # prd=predict(grandmodels[[2]][[i]],data.frame(cwd.ave=cwd.x, aet.ave=xlim.sp$aet.qmed[which(xlim.sp$species_name==species)]),interval="prediction",level=0.90)
  predictions=rbind(predictions,cbind(prd,cwd.x,aet.x))
}
predictions$species=rep(as.data.frame(grandmodels)[,1],each=length(cwd.x))
# x11()
# par(mfrow=c(5,4),mar=c(5,5,4,2))
for(i in 1:dim(xlim.sp)[1]){
  dat=predictions[which(predictions$species==xlim.sp$species_name[i]),]
  # plot(dat$pet.x,dat[,1],type="l",xlab="Average Climatic Water Deficit",ylab="",las=1,lwd=2,col="#0dcfca",ylim=c(min(dat$lwr),max(dat$upr)),main=dat$species[1])
  # polygon(c(dat$pet.x,rev(dat$pet.x)),c(dat[,2],rev(dat[,3])),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)
  plot(dat$cwd.x,dat[,1],type="l",xlab="Average Climatic Water Deficit",ylab="",las=1,lwd=2,col="#0dcfca",xlim=c(xlim.all$cwd.qlow, xlim.all$cwd.qhigh), ylim = c(-2e-3, 1e-3),main=dat$species[1])
  # ylim=c(min(dat$lwr),max(dat$upr))
  polygon(c(dat$cwd.x,rev(dat$cwd.x)),c(dat[,2],rev(dat[,3])),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)
  title(ylab="Effect of Water Deficit on Tree Growth", line=4)
  abline(h=0,lty=2,lwd=2)
}

