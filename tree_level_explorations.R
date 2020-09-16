#get dendro_df data frame from first 79 lines of code from 4a

library(dlnm)
library(tidyverse)
library(data.table)
library(lfe)

#drop a couple spurious years
dendro_df=dendro_df%>%
  filter(year<2020)

#filter non-unique ids
index <- dendro_df %>% select(collection_id, year, tree, core)
dendro_df <- dendro_df[!duplicated(index),]

dendro_df <- dendro_df %>% 
  rename(site = collection_id)

nlags=30
for(i in 1:nlags){
  lagdf=dendro_df%>%
    select(-c(width,pet.an,aet.an,site))%>%
    mutate(year=year+i)
  colnames(lagdf)[5]=paste0("cwd_L",i)
  if(i==1) dendro_lagged=left_join(dendro_df,lagdf)
  if(i>1) dendro_lagged=left_join(dendro_lagged,lagdf)
  print(i)
}

# fwrite(dendro_lagged,file="C:\\Users\\fmoore\\Desktop\\treedendrolagged.csv")

#merge in species data
site_df <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id)
dendro_lagged=left_join(dendro_lagged,site_df,by = "collection_id")
dendro_lagged <- dendro_lagged %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id),
         genus = substr(species_id, 1, 2))

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

