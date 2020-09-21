#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Creat predictions of growth impacts from climate change
#
# Input files:
# - ss_mod: R model object saved from Second stage
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Dendrochronologies
dendro_dir <- paste0(wdir, "out\\dendro\\")
dendro_df <- read_csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  filter(year>1900) %>% 
  select(-core_id)

# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'out\\climate\\essentialcwd_data.csv')
cwd_df <- read_csv(cwd_csv)
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site))

# 3. Site summary data
site_df <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize and merge site historic climate ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  rename(collection_id = site_id) %>% 
  group_by(collection_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum((aet+cwd)))

dendro_df <- dendro_df %>% 
  left_join(clim_df, by = c("collection_id", "year"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clean up data  ------------------------------
# TODO: should be integrated into dendro prep - 1b
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#drop a couple spurious years
dendro_df=dendro_df%>%
  filter(year<2020)

dendro_df <- dendro_df %>% 
  rename(site = collection_id)

nlags=30
for(i in 1:nlags){
  lagdf=dendro_df %>%
    select(-c(width,pet.an,aet.an,site)) %>%
    mutate(year = year + i)
  colnames(lagdf)[5]=paste0("cwd_L",i)
  if(i==1) dendro_lagged=left_join(dendro_df,lagdf)
  if(i>1) dendro_lagged=left_join(dendro_lagged,lagdf)
  print(i)
}

# fwrite(dendro_lagged,file="C:\\Users\\fmoore\\Desktop\\treedendrolagged.csv")

#merge in species data
dendro_lagged=left_join(dendro_lagged,site_df,by = "collection_id")
dendro_lagged <- dendro_lagged %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id),
         genus = substr(species_id, 1, 2)) # TODO: need to correct this

#log ring width for dependent variable
dendro_lagged$ln_rwi=log(dendro_lagged$width)

relgenus=c("ju","pi","ps")
genlist=list()
cblim=c(1000,1400,900)

for(i in 1:length(relgenus)){
  gendat=dendro_lagged%>%
    filter(genus==relgenus[i])%>%
    filter(cwd.an<quantile(cwd.an,p=0.99,na.rm=T))
  lagged=data.frame(gendat$cwd.an,gendat[,grep("cwd_L",colnames(gendat))])
  cblagged=crossbasis(lagged,lag=c(0,30),argvar=list("bs",degree=3),arglag=list(knots=logknots(30,4)))
  gendat$id=interaction(gendat$collection_id,gendat$tree)
  lagmod=felm(gendat$ln_rwi~cblagged+gendat$pet.an|id|0|collection_id,data=gendat)
  genlist[[i]]=list(cblagged,lagmod)
  genlist[[i]][[3]]=crosspred(cblagged,lagmod,cen=0,at=0:cblim[i]*1,cumul=TRUE)
  print(i)
}

x11()
par(mfrow=c(2,2))
for(i in 1:length(relgenus)){
  plot(genlist[[i]][[3]],var=100,xlab="Lagged Effect of CWD=800",ylab="Log Ring Width Growth",main=paste("Genus=",relgenus[i]),cumul=FALSE)
}

intlist=list()
#interactions effect
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
dendro_lagged=dendro_lagged%>%
  mutate(lag_5=rowMeans(select(.,cwd_L1,cwd_L2,cwd_L3,cwd_L4,cwd_L5)))%>%
  mutate(lag_6_10=rowMeans(select(.,cwd_L6,cwd_L7,cwd_L8,cwd_L9,cwd_L10)))
dendro_lagged$id=interaction(dendro_lagged$collection_id,dendro_lagged$tree)
intmod=felm(ln_rwi~cwd.an*lag_5+cwd.an*lag_6_10+pet.an|id|0|collection_id,data=dendro_lagged)

#early life drought

flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv')) %>%
  select(collection_id,tree,young)
young_trees=left_join(dendro_df,flm_df)%>%
  filter(young==TRUE)%>%
  select(collection_id,tree,year,cwd.an)%>%
  group_by(collection_id,tree)%>%
  filter(year<(min(year)+10))%>%
  summarize(earlylifecwd=mean(cwd.an))

dendro_df=left_join(dendro_df,site_df,by = "collection_id")
dendro_df <- dendro_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id),
         genus = substr(species_id, 1, 2))
youngtrees_df=inner_join(dendro_df,young_trees)
youngtrees_df$ln_rwi=log(youngtrees_df$width)

earlylifelist=list()
#interactions effect
for(i in 1:length(relgenus)){
  gendat=youngtrees_df%>%
    filter(genus%in%relgenus)%>%
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

