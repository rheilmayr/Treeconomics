# TODO:
#   Bootstrap first stage
#   Look at temperature and latitude in second stage regression (pure cwd in first stage)
#   Add measure of species distance from cwd mean - standardized metric of cwd.ave
#   add lagged climate data maybe up to 2 years
#   Hardwood vs softwood
#   Do you get snow?
#   PET in first stage (keep this separate from cwd model)
#   Try out single-stage approach - panel model with non-linear terms  
  


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

complete_df <- df_trim %>%
  drop_na(c("cwd.an","aet.an","pet.an", "ring_width"))

site_summary <- complete_df %>% 
  select("site_id", "species_id", "species_name", "cwd.ave", "pet.ave", "aet.ave") %>%
  distinct() %>%
  as_tibble()

#### Calculate standardized climate variables for each species ####
site_summary <- site_summary %>%
  group_by(species_id) %>%
  mutate(cwd.spstd = (cwd.ave - mean(cwd.ave)) / sd(cwd.ave),
         aet.spstd = (aet.ave - mean(aet.ave)) / sd(aet.ave),
         pet.spstd = (pet.ave - mean(pet.ave)) / sd(pet.ave))

to_merge <- site_summary %>%
  select("site_id", "species_id", "cwd.spstd", "aet.spstd", "pet.spstd")
complete_df <- complete_df %>%
  merge(y = to_merge, by = c("site_id", "species_id"), all.x = TRUE)

#### Run first stage ####
fit_mod <- function(d) {
  felm(ring_width ~ cwd.an + pet.an + age+I(age^2)+I(age^3)+I(age^4)|tree_id|0|0, data = d )
}

site_lm <- complete_df %>% 
  group_by(site_id, species_id) %>%
  nest() %>%
  mutate(mod = map(data, fit_mod),
         mod = map(mod, tidy)) %>%
  unnest(mod) %>%
  filter(term %in% c('cwd.an', 'pet.an'))

siteCoef <- site_lm %>%
  pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value")) %>%
  left_join(site_summary, by = c("site_id","species_id")) %>%
  select(-data)

#remove one super extreme outlier
siteCoef_trimmed <- siteCoef %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T),
         pet.qhigh=quantile(estimate_pet.an,0.99,na.rm=T),
         pet.qlow=quantile(estimate_pet.an,0.01,na.rm=T))
siteCoef_trimmed=siteCoef_trimmed %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         estimate_pet.an>pet.qlow & estimate_pet.an<pet.qhigh)


### Explore first stage ####
mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, data= siteCoef)
mod %>% summary()
mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd, data= siteCoef)
mod %>% summary()

siteCoef <- siteCoef %>%
  group_by("species_id") %>%
  mutate(cwd.q = ntile(cwd.ave, 5))
p <- siteCoef_trimmed %>% ggplot(aes(group=cwd.q, y=estimate_cwd.an)) + 
  geom_boxplot()
p



siteCoef_trimmed <- siteCoef_trimmed %>%
  group_by("species_id") %>%
  mutate(cwd.q = ntile(cwd.ave, 5),
         pet.q = ntile(pet.ave, 5),
         cwd.high = ntile(cwd.ave, 2)-1,
         pet.high = ntile(pet.ave, 2)-1,
         climate = paste0('cwd', cwd.high, '_pet', pet.high)) %>%
  ungroup()
siteCoef_trimmed <- siteCoef_trimmed %>%
  mutate(cwd.q = as.factor(ntile(cwd.spstd, 7)),
         pet.q = as.factor(ntile(pet.spstd, 7)),
         cwd.high = as.factor(ntile(cwd.spstd, 2)-1),
         pet.high = as.factor(ntile(pet.spstd, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high))
plotCoef <- siteCoef_trimmed %>%
  filter(species_id == "pcgl") %>%
    mutate(cwd.q = as.factor(ntile(cwd.spstd, 5)),
         pet.q = as.factor(ntile(pet.spstd, 5)),
         aet.q = as.factor(ntile(aet.spstd, 5)),
         cwd.high = as.factor(ntile(cwd.spstd, 2)-1),
         pet.high = as.factor(ntile(pet.spstd, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high))

p1 <- plotCoef %>% ggplot(aes(x=climate, y=estimate_cwd.an)) +
  geom_boxplot()
p2 <- plotCoef %>% ggplot(aes(x=climate, y=estimate_pet.an)) +
  geom_boxplot()
p1 / p2


p1 <- plotCoef %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot()
p2 <- plotCoef %>% ggplot(aes(x=cwd.q, y=estimate_pet.an)) +
  geom_boxplot()
p3 <- plotCoef %>% ggplot(aes(x=pet.q, y=estimate_cwd.an)) +
  geom_boxplot()
p4 <- plotCoef %>% ggplot(aes(x=pet.q, y=estimate_pet.an)) +
  geom_boxplot()

(p1 | p2) / (p3 | p4)


p <- plotCoef %>% ggplot(aes(x=cwd.spstd, y=pet.spstd, colour = estimate_cwd.an)) +
  geom_point() +
  xlim(c(-2,2)) +
  ylim(c(-2,2))
p

p <- plotCoef %>% ggplot(aes(x=cwd.spstd, y=pet.spstd, colour = estimate_pet.an)) +
  geom_point() +
  xlim(c(-2,2)) +
  ylim(c(-2,2))
p


p <- plotCoef %>% ggplot(aes(x=aet.ave, y=cwd.ave, colour = estimate_cwd.an)) +
  geom_point() 
p

+ scale_colour_brewer(palette = "Greens")

lm(estimate_cwd.an ~ cwd.spstd + I(cwd.spstd),weights=errorweights,data= .)




# ### TEMP - KEEP ALL TO KEEP JUNIPER
# siteCoef_trimmed <- siteCoef

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
  distinct() %>%
  ungroup()

xlim.all <- site_summary %>%
  ungroup() %>%
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

kd <- kde2d(x = site_summary$cwd.spstd, y = site_summary$pet.spstd, 
            n = 25, lims = c(c(-2.5,2.5),c(-2.5, 2.5)))
filled.contour(kd, color.palette = pal)
# points(site_summary$cwd.spstd, site_summary$pet.spstd)

x11()
par(mfrow=c(4,5),mar=c(5,5,4,2))
for(i in 1:dim(xlim.sp)[1]){
  spec_name = xlim.sp$species_name[i]
  dat = siteCoef_trimmed %>% filter(species_name==spec_name)
  lim = xlim.sp %>% filter(species_name==spec_name)
  # kd <- kde2d(x = dat$cwd.ave, y = dat$pet.ave, n = 25, lims = c(c(lim$cwd.qlow, lim$cwd.qhigh), c(lim$pet.qlow, lim$pet.qhigh)))
  kd <- kde2d(x = dat$cwd.ave, y = dat$aet.ave, n = 25, 
              lims = c(c(xlim.all$cwd.qlow,xlim.all$cwd.qhigh),c(xlim.all$aet.qlow, xlim.all$aet.qhigh)))
  filled.contour(kd, color.palette = pal)
  title(paste0(spec_name, dim(dat)[1]))
  # kde2d(x, y, h, n = 25, lims = c(range(x), range(y)))
}

#### Run second stage of model  ####
dat = siteCoef_trimmed %>% 
  ungroup() %>%
  drop_na(std.error_cwd.an) %>%
  mutate(errorweights=1/std.error_cwd.an/sum(1/std.error_cwd.an,na.rm=T)) 
grandmodel <- lm(estimate_cwd.an~cwd.ave + I(cwd.ave^2),weights=errorweights,data= dat)
# grandmodel <- lm(estimate~cwd.ave+I(cwd.ave^2), weights=errorweights,data= dat)

cwd.x <- seq(xlim.all$cwd.qlow, xlim.all$cwd.qhigh, length.out=1000)
aet.med <- xlim.all$aet.qmed
pet.med <- xlim.all$pet.qmed

pred <- predict(grandmodel,data.frame(cwd.ave=cwd.x, pet.ave = pet.med),interval="prediction",level=0.90)
pred <- cbind(pred,cwd.x)
pred <- pred %>% as_tibble()
plot(cwd.x, pred$fit, type="l",xlab="Historic CWD",ylab="Estimated effect of CWD",las=1,lwd=2,col="#0dcfca",
     ylim = c(pred$lwr %>% min(), pred$upr %>% max()))
points(dat$cwd.ave,dat$estimate,col=dat$species_id,pch=18)
polygon(c(cwd.x,rev(cwd.x)),c(pred$lwr,rev(pred$upr)),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)


pet.x <- seq(xlim.all$pet.qlow, xlim.all$pet.qhigh, length.out=1000)
cwd.med <- xlim.all$cwd.qmed

pred <- predict(grandmodel,data.frame(cwd.ave=cwd.med, pet.ave = pet.x),interval="prediction",level=0.90)
pred <- cbind(pred,pet.x)
pred <- pred %>% as_tibble()
plot(pet.x, pred$fit, type="l",xlab="Historic PET",ylab="Marginal effect of CWD",las=1,lwd=2,col="#0dcfca")
points(dat$pet.ave,dat$estimate,add=T)
polygon(c(pet.x,rev(pet.x)),c(pred$lwr,rev(pred$upr)),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)


#### Run species specific second stage ####
grandmodels=siteCoef_trimmed %>%
  group_by(species_name) %>%
  drop_na(std.error) %>%
  mutate(errorweights=1/std.error/sum(1/std.error,na.rm=T)) %>% ## LETS RECONSIDER THIS 
  # do(speciesmod=lm(estimate~cwd.ave+I(cwd.ave^2)+I(cwd.ave^3),weights=errorweights,data= .))
  do(speciesmod=lm(estimate~cwd.ave + I(cwd.ave^2) + pet.ave + I(pet.ave^2),weights=errorweights,data= .))

# Review results
grandmodels %>% filter(species_name=="blue oak") %>% pull(speciesmod)

y=list()
predictions=data.frame()
for(i in 1:dim(grandmodels)[1]){
  species=as.data.frame(grandmodels)[i,1]
  cwd.x=seq(xlim.sp$cwd.qlow[which(xlim.sp$species_name==species)],xlim.sp$cwd.qhigh[which(xlim.sp$species_name==species)],length.out=1000)
  pet.x=seq(xlim.sp$pet.qlow[which(xlim.sp$species_name==species)],xlim.sp$pet.qhigh[which(xlim.sp$species_name==species)],length.out=1000)
  prd=predict(grandmodels[[2]][[i]],data.frame(cwd.ave=cwd.x, pet.ave = xlim.sp$aet.qmed[which(xlim.sp$species_name==species)]),interval="prediction",level=0.90)
  # prd=predict(grandmodels[[2]][[i]],data.frame(cwd.ave=cwd.x, aet.ave=xlim.sp$aet.qmed[which(xlim.sp$species_name==species)]),interval="prediction",level=0.90)
  # prd=predict(grandmodels[[2]][[i]],data.frame(cwd.ave=xlim.sp$cwd.qmed[which(xlim.sp$species_name==species)], aet.ave=aet.x),interval="prediction",level=0.90)
  predictions=rbind(predictions,cbind(prd,cwd.x,pet.x))
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





##### 
#### Single-stage panel model #####
low <- complete_df %>% filter(cwd.dev<=-0.1)
mod_low <- felm(ring_width ~ cwd.an + I(cwd.an^2) + age + 
                  I(age^2) + I(age^3) + I(age^4) | site_id + year, 
                data = low)

med <- complete_df %>% filter((cwd.dev<0.1) & (cwd.dev>-0.1))
mod_med <- felm(ring_width ~ cwd.an + I(cwd.an^2) + age + 
                  I(age^2) + I(age^3) + I(age^4) | tree_id + year, 
                data = med)


high <- complete_df %>% filter(cwd.dev>=0.1)
mod_high <- felm(ring_width ~ cwd.an + I(cwd.an^2) + age + 
                   I(age^2) + I(age^3) + I(age^4) | site_id + year, 
                 data = high)

mod <- felm(ring_width ~ cwd.an + aet.an + 
              age + I(age^2) + I(age^3) + I(age^4) | tree_id + year, 
            data = high)