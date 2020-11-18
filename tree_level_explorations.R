#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 9/21/20
# Purpose: Assess tree-level evolution of drought sensitivity
#
# Input files:
# - 
#
# ToDo:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyr)
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(lfe)
library(broom)
library(purrr)
library(fixest)
library(dlnm)
library(tidyverse)
library(data.table)
# library(tidylog)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out\\dendro\\")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long.csv")) %>% 
  as_tibble()
dendro_df <- dendro_df %>% 
  select(-core_id) %>% 
  filter(year>=1900) %>% 
  filter(year<2020)
  # TODO: should be debugged in dendro prep - 1b


# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'out\\climate\\essentialcwd_data.csv')
cwd_df <- read_csv(cwd_csv)
cwd_df <- cwd_df %>% 
  rename("collection_id" = site)

# 3. Site summary data
site_df <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id) %>% 
  mutate(sp_id = str_to_lower(sp_id))

# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family) %>% 
  rename(sp_id = species_id)
site_df <- site_df %>% 
  left_join(sp_info, by = "sp_id")

# Merge site / species data back to dendro
dendro_df <- dendro_df %>%  
  left_join(site_df,by = "collection_id") %>% 
  mutate(ln_rwi = log(width))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site historic climate ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  group_by(collection_id, year) %>%
  summarise(aet.an = sum(aet),
            pet.an = sum((aet+cwd)),
            cwd.an = sum(cwd))


# Calculate site-level historic climate
hist_clim_df <- clim_df %>%
  group_by(collection_id) %>%
  filter(year<1980) %>%
  summarise(aet.ave = mean(aet.an),
            cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an),
            cwd.min = min(cwd.an))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create lagged cwd  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlags=30
for(i in 1:nlags){
  lagdf=clim_df %>%
    select(c(collection_id,year,cwd.an)) %>%
    mutate(year = year + i)
  colnames(lagdf)[3]=paste0("cwd_L",i)
  if(i==1) clim_lagged=left_join(clim_df,lagdf, by = c("collection_id", "year"))
  if(i>1) clim_lagged=left_join(clim_lagged,lagdf, by = c("collection_id", "year"))
  print(i)
}

dendro_lagged <- dendro_df %>% 
  left_join(clim_lagged, by = c("collection_id", "year"))

# fwrite(dendro_lagged,file="C:\\Users\\fmoore\\Desktop\\treedendrolagged.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Expore site-level DLNM  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keep_sites = dendro_df %>% 
  group_by(collection_id) %>% 
  tally() %>% 
  arrange(desc(n))

site = "NETH031"
site = "AUST115"
site = "RUSS223"
site = "CA667"
site = "MONG039"
site = "SWIT285"
site = "WV006"

i = 203
site = keep_sites[i,1] %>% pull()

gendat=dendro_lagged%>%
  filter(collection_id==site)
  
shock = gendat$cwd.an %>% quantile(.9, na.rm = TRUE)
at_vals = (shock - 50):(shock + 50)

lagged=data.frame(gendat$cwd.an,gendat[,grep("cwd_L",colnames(gendat))])
cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=2),arglag=list(knots=logknots(30, 3)))
gendat$id=interaction(gendat$collection_id,gendat$tree)
lagmod=lm(gendat$width~cblagged+gendat$pet.an,data=gendat)
crosspredict = crosspred(cblagged, lagmod, cen=0, at = at_vals, cumul = TRUE)
plot(crosspredict,
     var=shock,
     xlab=paste0("Lagged Effect of CWD=", as.integer(shock)),
     ylab="Ring Width Growth",main=paste("Site=",site),cumul=FALSE)






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Genus-specific distributed lag models  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock <- 300
at_vals = 0:1000
relgenus=c("Juniperus", "Pinus", "Pseudotsuga", "Abies", "Fagus", "Picea", "Larix", "Quercus", "Tsuga")
genlist=list()

i = 1
for(i in 1:length(relgenus)){
  gendat=dendro_lagged%>%
    filter(genus==relgenus[i])%>%
    filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
  # shock = gendat$cwd.an %>% quantile(.9, na.rm = TRUE)
  # at_vals = (shock - 50):(shock + 50)
  lagged=data.frame(gendat$cwd.an,gendat[,grep("cwd_L",colnames(gendat))])
  cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(30,4)))
  gendat$id=interaction(gendat$collection_id,gendat$tree)
  lagmod=lm(gendat$width~cblagged+gendat$pet.an,data=gendat)
  # crosspredict = crosspred(cblagged, lagmod, cen=0, at=0:cblim[i]*1, cumul = FALSE)
  genlist[[i]]=list(cblagged,lagmod)
  genlist[[i]][[3]]=crosspred(cblagged,lagmod,cen=0,at=at_vals,cumul=FALSE)
  print(i)
}

+x11()
par(mfrow=c(3,3))
for(i in 1:length(relgenus)){
  plot(genlist[[i]][[3]],
       var=300,
       xlab=paste0("Lagged Effect of CWD=", as.integer(shock)),
       ylab="Log Ring Width Growth",main=paste("Genus=",relgenus[i]),cumul=FALSE)
}

intlist=list()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Site-level distributed lag model  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sites <- dendro_lagged %>% 
  pull(collection_id) %>% 
  unique()
n_site = indices[i]
indices = c(1703, 302, 1252, 1367, 2028)
genlist=list()

for(i in 1:length(indices)){
  n_site <- indices[i]
  gendat=dendro_lagged%>%
    filter(collection_id==sites[n_site])%>%
    filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
  lagged=data.frame(gendat$cwd.an,gendat[,grep("cwd_L",colnames(gendat))])
  cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(30,4)))
  gendat$id=interaction(gendat$collection_id,gendat$tree)
  lagmod=lm(gendat$width~cblagged+gendat$pet.an,data=gendat)
  genlist[[i]]=list(cblagged,lagmod)
  genlist[[i]][[3]]=crosspred(cblagged,lagmod,cen=0,at=0:cblim[i]*1,cumul=TRUE)
  print(i)
}


x11()
par(mfrow=c(3,3))
for(i in 1:length(indices)){
  plot(genlist[[i]][[3]],var=300,xlab="Lagged Effect of CWD=800",ylab="Log Ring Width Growth",main=paste("Genus=",relgenus[i]),cumul=FALSE)
}
# subset_dat <- dendro_lagged%>%
#   filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T),
#          genus %in% c("Juniperus", "Pinus", "Pseudotsuga", "Abies", "Fagus", "Picea", "Larix"))
# lagged=data.frame(subset_dat$cwd.an,subset_dat[,grep("cwd_L",colnames(subset_dat))])
# cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(nlags,4)))
# subset_dat$id=interaction(subset_dat$collection_id,subset_dat$tree)
# lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
# prediction <- crosspred(cblagged,lagmod,cen=0,at=0:3000*1,cumul=TRUE)
# plot(prediction, var = 100, xlab = "Lagged Effect of CWD=800", ylab = "log Ring Width Growth")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lag interaction models  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i in 1:length(relgenus)){
  gendat=dendro_lagged%>%
    filter(genus==relgenus[i])%>%
    filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))%>%
    mutate(lag_5=rowMeans(select(.,cwd_L1,cwd_L2,cwd_L3,cwd_L4,cwd_L5)))%>%
    mutate(lag_6_10=rowMeans(select(.,cwd_L6,cwd_L7,cwd_L8,cwd_L9,cwd_L10)))
  gendat$id=interaction(gendat$collection_id,gendat$tree)
  intmod=felm(ln_rwi~cwd.an*lag_5+cwd.an*lag_6_10+pet.an|id|0|collection_id,data=gendat)
  intlist[[i]]=intmod
  print(i)
}

#try interactions model for whole dataset
dendro_lagged <- dendro_lagged %>%
  mutate(lag_5 = rowMeans(select(.,cwd_L1,cwd_L2,cwd_L3,cwd_L4,cwd_L5)),
         lag_2_5 = rowMeans(select(.,cwd_L2,cwd_L3,cwd_L4,cwd_L5)),
         lag_6_10 = rowMeans(select(.,cwd_L6,cwd_L7,cwd_L8,cwd_L9,cwd_L10)))
dendro_lagged$id=interaction(dendro_lagged$collection_id,dendro_lagged$tree)
intmod=felm(ln_rwi~cwd.an*lag_2_5+cwd.an*lag_6_10+pet.an|id|0|collection_id,data=dendro_lagged)
intmod=felm(ln_rwi~cwd.an + lag_5 + lag_6_10 + pet.an|id|0|collection_id,data=dendro_lagged)
summary(intmod)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Early life drought models  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
young_trees <- dendro_lagged %>% 
  group_by(collection_id, tree) %>% 
  summarise(min_year = min(year)) %>% 
  mutate(young = min_year>1901)

# flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv')) %>%
#   select(collection_id,tree,young)
young_trees <- young_trees %>%
  select(-min_year) %>% 
  right_join(dendro_lagged, by = c("collection_id", "tree")) %>%
  filter(young==TRUE)%>%
  select(collection_id,tree,year,cwd.an)%>%
  group_by(collection_id,tree)%>%
  filter(year<(min(year)+10))%>%
  summarize(earlylifecwd=mean(cwd.an))

youngtrees_df <- dendro_lagged %>%
  select(collection_id, tree, year, sp_id, ln_rwi, cwd.an, pet.an) %>% 
  inner_join(young_trees, by = c("collection_id", "tree"))

earlylifelist=list()
#interactions effect
for(i in 1:length(relgenus)){
  gendat=youngtrees_df%>%
    filter(genus == relgenus[i])%>%
    filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
  gendat$id=interaction(gendat$collection_id,gendat$tree)
  intmod=felm(ln_rwi~cwd.an*I(earlylifecwd>300)+pet.an-I(earlylifecwd>300)|id|0|collection_id,data=gendat)
  earlylifelist[[i]]=intmod
  print(i)
}

youngtrees_df$id=interaction(youngtrees_df$collection_id,youngtrees_df$tree)
youngtrees_df=youngtrees_df%>%
  filter(earlylifecwd<quantile(earlylifecwd,p=0.99,na.rm=T))%>%
  filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
earlylifemod=felm(ln_rwi~cwd.an*I(earlylifecwd>300)+pet.an-I(earlylifecwd>300)|id|0|collection_id,data=youngtrees_df)
earlylifemod=felm(ln_rwi~cwd.an*earlylifecwd + pet.an - earlylifecwd|id|0|collection_id,data=youngtrees_df)


youngtrees_df <- youngtrees_df %>% 
  left_join(hist_clim_df, by = "collection_id")

subset <- youngtrees_df %>% 
  filter(sp_id == "psme")
earlylifemod=felm(ln_rwi~cwd.an*earlylifecwd + cwd.an*cwd.ave + pet.an |id|0|collection_id,data=subset)
summary(earlylifemod)


earlylifemod=felm(ln_rwi~cwd.an*cwd.ave + pet.an |id|0|collection_id,data=subset)
summary(earlylifemod)


hist_clim_temp <- hist_clim_df %>% 
  left_join(site_df, by = "collection_id") 
hist_clim_temp <- hist_clim_temp %>% 
  group_by(sp_id) %>% 
  summarise(sp.cwd.median = median(cwd.ave, na.rm = TRUE)) %>% 
  right_join(hist_clim_temp, by = "sp_id") %>% 
  mutate(high_cwd_site = cwd.ave > sp.cwd.median)

cblim <- 1000
subset_dat <- dendro_lagged %>%
  left_join(hist_clim_temp %>% select(collection_id, high_cwd_site), by = "collection_id") %>% 
  filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T),
         sp_id =="pipo")
lagged=data.frame(subset_dat$cwd.an,subset_dat[,grep("cwd_L",colnames(subset_dat))])
cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(nlags,4)))
subset_dat$id=interaction(subset_dat$collection_id,subset_dat$tree)
lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
prediction <- crosspred(cblagged,lagmod,cen=0,at=0:3000*1,cumul=TRUE)
plot(prediction, var = 100, xlab = "Lagged Effect of CWD=800", ylab = "log Ring Width Growth")


subset_dat <- dendro_lagged %>%
  left_join(hist_clim_temp %>% select(collection_id, high_cwd_site), by = "collection_id") %>% 
  filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T),
         sp_id =="psme",
         high_cwd_site==0)
lagged=data.frame(subset_dat$cwd.an,subset_dat[,grep("cwd_L",colnames(subset_dat))])
cblagged=crossbasis(lagged,lag=c(0,nlags),argvar=list("bs",degree=3),arglag=list(knots=logknots(nlags,4)))
subset_dat$id=interaction(subset_dat$collection_id,subset_dat$tree)
lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
prediction <- crosspred(cblagged,lagmod,cen=0,at=0:3000*1,cumul=TRUE)
plot(prediction, var = 100, xlab = "Lagged Effect of CWD=800", ylab = "log Ring Width Growth")





calc_germ_clim <- function(site, germ_year){
  pre_germ_clim <- clim_df %>% 
    filter(collection_id == site,
           year < germ_year,
           year < 1970) %>% 
    summarise(pre_cwd = mean(cwd.an),
              pre_pet = mean(cwd.an))
  post_germ_clim <- clim_df %>% 
    filter(collection_id == site,
           year>germ_year,
           year<1970) %>% 
    summarise(post_cwd = mean(cwd.an),
              post_pet = mean(pet.an))
  germ_clim <- pre_germ_clim %>% 
    left_join(post_germ_clim, by = "collection_id")
  return(germ_clim)
}

young_trees <- young_trees %>% 
  filter(young==TRUE) %>% 
  mutate(germ_clim = map2(collection_id, min_year, calc_germ_clim))

youngtrees_df <- dendro_lagged %>%
  select(collection_id, tree, year, sp_id, ln_rwi, cwd.an, pet.an) %>% 
  inner_join(young_trees, by = c("collection_id", "tree"))

lagmod=felm(subset_dat$ln_rwi~cblagged+subset_dat$pet.an|id|0|collection_id,data=subset_dat)
         
         