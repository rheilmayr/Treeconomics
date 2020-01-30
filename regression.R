#install.packages("plm")
#install.packages("lfe")

library(tidyverse)
library(lfe)
library(broom)
library(purrr)
library(ggiraphExtra)

### Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
#for Fran
wdir="C:/Users/fmoor/Google Drive/Treeconomics/Data/"
setwd(wdir)

### Read data
spp=c("psme","pipo","pisy","pcgl","pcab","pied")
speciesnames=c("Douglas Fir","Ponderosa Pine","Scotch Pine","White Spruce","Norway Spruce","Colorado Pinyon")

df <- read.csv(paste0(wdir,'topsixspecies_data.csv'))
df <- df %>%
  mutate(tree_id = paste0(site_id, tree_id))
df <- df %>%
  mutate(pet.an = cwd.an + aet.an)

### Calculate site average aet and cwd
site_clim <- df %>%
  filter(year < 1980) %>%                                                                      ## Switch to raw climate data rather than using tree dataframe
  group_by(site_id, year) %>%
  summarize(cwd.an = max(cwd.an, na.rm = TRUE), aet.an = max(aet.an, na.rm = TRUE), pet.an = max(pet.an, na.rm = TRUE)) %>%
  filter((cwd.an > -Inf) & (aet.an > -Inf)) %>%
  group_by(site_id) %>%
  summarize(cwd.ave = mean(cwd.an, na.rm = TRUE), aet.ave = mean(aet.an, na.rm = TRUE), pet.ave = mean(pet.an, na.rm = TRUE))

df <- merge(x = df, y = site_clim, by = "site_id", all.x = TRUE)

#remove one site with negative ages
df=df %>%
  filter(age>0)

#remove trees with a very short record (<5 years)
treeobs=df %>%
  group_by(tree_id) %>%
  summarize(nyears=n())
df=left_join(df,treeobs)
df=df %>%
  filter(nyears>5)

#remove sites with few trees
ntrees=df%>%
  group_by(site_id,species)%>%
  summarize(ntrees=length(unique(tree_id)))
df=left_join(df,ntrees)
df=df%>%
  filter(ntrees>5)

complete_df=df %>%
  drop_na(c("cwd.an","pet.an","ring_width"))

site_lm <- complete_df %>% 
  group_by(site_id,species) %>%
  do(fit_site = felm(ring_width ~ cwd.an+pet.an+age+I(age^2)+I(age^3)+I(age^4)|tree_id|0|0, data = . ))

siteCoef=tidy(site_lm[[3]][[1]])%>%filter(term==c("cwd.an", "pet.an"))
for(i in 2:length(site_lm[[3]])){
  siteCoef=rbind(siteCoef,tidy(site_lm[[3]][[i]])%>%filter(term==c("cwd.an", "pet.an")))
}
siteCoef$site_id=site_lm$site_id
siteCoef$species=site_lm$species
siteCoef=left_join(siteCoef,unique(complete_df[,which(colnames(complete_df)%in%c("site_id","species","cwd.ave", "pet.ave"))]))
#remove one super extreme outlier
siteCoef_trimmed=siteCoef %>%
  group_by(species) %>%
  mutate(qhigh=quantile(estimate,0.99,na.rm=T)) %>%
  mutate(qlow=quantile(estimate,0.01,na.rm=T))
siteCoef_trimmed=siteCoef_trimmed %>%
  filter(estimate>qlow&estimate<qhigh)
siteCoef_trimmed$species=recode_factor(siteCoef_trimmed$species,psme="Douglas Fir",pipo="Ponderosa Pine",pisy="Scotch Pine",pcgl="White Spruce",pcab="Norway Spruce",pied="Colorado Pinyon")
a=ggplot(siteCoef_trimmed,aes(x=cwd.ave,y=estimate,col=species))+theme_bw()+geom_point()+facet_wrap(~species)
a=a+labs(x="Average Climatic Water Deficit",y="Marginal Effect of CWD")+scale_color_discrete(guide=FALSE)

x11()
par(mfrow=c(2,3),mar=c(5,5,4,2))
for(i in 1:dim(xlim)[1]){
  spec_id = xlim$species[i]
  dat=siteCoef_trimmed %>% filter(species==spec_id)
  lim = xlim %>% filter(species==spec_id)
  # kd <- kde2d(x = dat$cwd.ave, y = dat$pet.ave, n = 25, lims = c(c(lim$cwd.qlow, lim$cwd.qhigh), c(lim$pet.qlow, lim$pet.qhigh)))
  kd <- kde2d(x = dat$cwd.ave, y = dat$pet.ave, n = 25, lims = c(c(0,750),c(100, 1000)))
  contour(kd)
  # kde2d(x, y, h, n = 25, lims = c(range(x), range(y)))
}

grandmodels=siteCoef_trimmed %>%
  group_by(species) %>%
  drop_na(std.error) %>%
  mutate(errorweights=1/std.error/sum(1/std.error,na.rm=T)) %>%
  # do(speciesmod=lm(estimate~cwd.ave+I(cwd.ave^2),weights=errorweights,data= .)) 
  do(speciesmod=lm(estimate~cwd.ave+I(cwd.ave^2) + pet.ave + I(pet.ave^2),weights=errorweights,data= .))

# Review results
grandmodels %>% filter(species=="Colorado Pinyon") %>% pull(speciesmod)

y=list()
xlim=siteCoef_trimmed %>%
  group_by(species) %>%
  mutate(cwd.qhigh=quantile(cwd.ave,0.95),
         cwd.qlow=quantile(cwd.ave,0.05),
         cwd.qmed=quantile(cwd.ave,0.5),
         pet.qhigh = quantile(pet.ave, 0.95), 
         pet.qlow = quantile(pet.ave, 0.05),
         pet.qmed = quantile(pet.ave,0.5)) %>%
  select(species,cwd.qhigh,cwd.qlow,cwd.qmed,pet.qhigh,pet.qlow,pet.qmed) %>%
  distinct()
predictions=data.frame()
for(i in 1:dim(grandmodels)[1]){
  species=as.data.frame(grandmodels)[i,1]
  cwd.x=seq(xlim$cwd.qlow[which(xlim$species==species)],xlim$cwd.qhigh[which(xlim$species==species)],length.out=1000)
  pet.x=seq(xlim$pet.qlow[which(xlim$species==species)],xlim$pet.qhigh[which(xlim$species==species)],length.out=1000)
  # predictions=rbind(predictions,cbind(predict(grandmodels[[2]][[i]],data.frame(cwd.ave=xlim$cwd.qmed[which(xlim$species==species)], pet.ave=pet.x),interval="prediction",level=0.90),cwd.x,pet.x))
  predictions=rbind(predictions,cbind(predict(grandmodels[[2]][[i]],data.frame(cwd.ave=cwd.x, pet.ave=xlim$pet.qmed[which(xlim$species==species)]),interval="prediction",level=0.90),cwd.x,pet.x))
  # predictions=rbind(predictions,cbind(predict(grandmodels[[2]][[i]],data.frame(cwd.ave=cwd.x, pet.ave=pet.x),interval="prediction",level=0.90),cwd.x,pet.x))
}
predictions$species=rep(as.data.frame(grandmodels)[,1],each=length(x))
x11()
par(mfrow=c(2,3),mar=c(5,5,4,2))
for(i in 1:dim(xlim)[1]){
  dat=predictions[which(predictions$species==xlim$species[i]),]
  # plot(dat$pet.x,dat[,1],type="l",xlab="Average Climatic Water Deficit",ylab="",las=1,lwd=2,col="#0dcfca",ylim=c(min(dat$lwr),max(dat$upr)),main=dat$species[1])
  # polygon(c(dat$pet.x,rev(dat$pet.x)),c(dat[,2],rev(dat[,3])),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)
  plot(dat$cwd.x,dat[,1],type="l",xlab="Average Climatic Water Deficit",ylab="",las=1,lwd=2,col="#0dcfca",ylim=c(min(dat$lwr),max(dat$upr)),main=dat$species[1])
  polygon(c(dat$cwd.x,rev(dat$cwd.x)),c(dat[,2],rev(dat[,3])),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)
  title(ylab="Effect of Water Deficit on Tree Growth", line=4)
  abline(h=0,lty=2,lwd=2)
}

