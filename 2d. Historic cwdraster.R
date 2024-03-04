#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 7/10/20
# Purpose: Create gridded 0.5 degree cwd, pet and aet based on 1901-1980 climatological baseline (CRU Data)
#
# Input files:
#   cru_ts4.04.1901.2019.tmp.dat.nc: Monthly mean temp data from CRU.
#   cru_ts4.04.1901.2019.pre.dat.nc: Monthly precip data from CRU. Accessed from https://crudata.uea.ac.uk/cru/data/hrg/.
#   sr_cru_max.asc: Soil water capacity raster. Accessed from Wang-Erlandsson et al., 2016.
#   slopeaspectraster.Rdat: FRAN - WHERE DOES THIS GET CALCULATED?
# 
# Output files:
#   monthlycrubaseline_tas:
#   monthlycrubaseline_pr:
#   griddedbaselinecwddata.csv:
#   HistoricCWD_AETGrids.Rdat: 
#   HistoricCWD_AETGrids_Annual.Rdat: 
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load packages --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(plyr)
library(raster)
library(sp)
library(rgdal)
library(reshape2)
library(ncdf4)
library(foreach)
library(parallel)
library(doParallel)
library(data.table)
library(geosphere)
library(dplyr)
library(tidyr)
source("f_cwd_function_new.R")
library(zoo)
library(tidyverse)

max_clusters = 20

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wdir <- "remote/"
# franspath="C:/Users/fmoore/Box/Davis Stuff/Treeconomics"

## Load CRU temp and precip data
cruyears=rep(1901:2019,each=12);crumonths=rep(1:12,length(1901:2019))
tas=stack(paste0(wdir, "0_raw/CRUData/cru_ts4.04.1901.2019.tmp.dat.nc"))
pr=stack(paste0(wdir,"0_raw/CRUData/cru_ts4.04.1901.2019.pre.dat.nc"))

## Load soil water capacity raster
swc=raster(paste0(wdir,"0_raw/other data for cwd/sr_cru_max.asc"))
#convert swc from mm to cm
swc=swc/10

## Load topography raster
load(paste0(wdir,"0_raw/other data for cwd/slopeaspectraster.Rdat"))
terrainraster=crop(terrainraster,extent(swc))
#convert slope and aspect to degrees
terrainraster=terrainraster*180/pi


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate precip / temp baselines --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseyears=1901:1980
basemonths=rep(1:12,length(baseyears))
years=rep(baseyears,each=12)

tas=tas[[which(cruyears%in%baseyears)]];pr=pr[[which(cruyears%in%baseyears)]]

## Don't need to run if getting annual cwd values (not normals for baseline period)
# tas_months=calc(tas[[which(basemonths==1)]],function(x) mean(x, na.rm=T));pr_months=calc(pr[[which(basemonths==1)]],function(x) mean(x, na.rm=T))
# for(j in 2:12){
#   tas_months=stack(tas_months,calc(tas[[which(basemonths==j)]],function(x) mean(x, na.rm=T)))
#   pr_months=stack(pr_months,calc(pr[[which(basemonths==j)]],function(x) mean(x, na.rm=T)))
#   print(j)
# }
# 

# save(tas_months,pr_months,file=paste0(wdir,"in/CRUData/monthlycrubaseline.Rdat"))
#writeRaster(tas_months %>% mean(), file=paste0(wdir,"1_input_processed/climate/monthlycrubaseline_tas"), overwrite = TRUE)
#writeRaster(pr_months %>% sum(), file=paste0(wdir,"1_input_processed/climate/monthlycrubaseline_pr"), overwrite = TRUE)


# melt permanent site data down into a long data frame for cwd function
sitedata=cbind(as.data.frame(terrainraster),as.data.frame(swc))
sitedata=cbind(sitedata,coordinates(swc))
colnames(sitedata)=c("slope","aspect","swc","lon","lat")
sitedata$site=1:dim(sitedata)[1]

# fix a few spurious slope values
sitedata$slope[which(sitedata$slope<0)]=0 

# crop temp and precip data to swc extent
#tas_months=crop(tas_months,extent(swc));pr_months=crop(pr_months,extent(swc))
tas=crop(tas,extent(swc));pr=crop(pr,extent(swc))

# # set up parallel processing
# cl=makeCluster(4)
# clusterExport(cl,c("data","setorder"))
# registerDoParallel(cl)

tasdata=as.data.frame(tas);prdata=as.data.frame(pr)
tasdata$site=1:dim(tasdata)[1];prdata$site=1:dim(tasdata)[1]
# tasdata=melt(tasdata,id.vars="site");prdata=melt(prdata,id.vars="site")
# climatedata=cbind(tasdata,prdata[,3])
# colnames(climatedata)=c("site","month","temp","precip")
# 
# dat=merge(sitedata,climatedata)
# dat=dat[complete.cases(dat),]
# 
# test=cwd_function(site=as.factor(dat$site),slope=dat$slope,latitude=dat$lat,foldedaspect=dat$aspect,ppt=dat$precip,tmean=dat$temp,month=dat$month,soilawc=dat$swc,year=NULL,type="normal")
# fwrite(test,file=paste0(wdir,"1_input_processed/climate/griddedbaselinecwddata.csv"))

# #merge in lat longs from site data
# sitecrosswalk=data.frame(site_grid=unique(dat$site),site=1:length(unique(dat$site)))
# sitecrosswalk <- sitecrosswalk %>% mutate(site = as.character(site))
# test=merge(test,sitecrosswalk, by = "site")
# 
# test$site_grid=as.numeric(test$site_grid)
# test=merge(test[,c(3,28,29,30)],sitedata[,c(4:6)],by.x="site_grid",by.y="site")
# test$cells=cellFromXY(cwd_historic,test[,c(5,6)])
# 
# testmonths=split(test,test$month)
# 
# cwd_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
# aet_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
# 
# cwd_historic[testmonths[[1]]$cells]=testmonths[[1]]$cwd
# aet_historic[testmonths[[1]]$cells]=testmonths[[1]]$aet
# 
# for(j in 2:12){
#   cwdtemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc));aettemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
#   cwdtemp[testmonths[[j]]$cells]=testmonths[[j]]$cwd
#   aettemp[testmonths[[j]]$cells]=testmonths[[j]]$aet
#   cwd_historic=stack(cwd_historic,cwdtemp);aet_historic=stack(aet_historic,aettemp)
# }
# 
# save(cwd_historic,aet_historic,file=paste0(wdir,"1_input_processed/climate/HistoricCWD_AETGrids.Rdat"))

tasdata=reshape2::melt(tasdata,id.vars=c("site"));prdata=reshape2::melt(prdata,id.vars=c("site"))
tasdata$year=as.numeric(substr(tasdata$variable,2,5))
tasdata$month=as.numeric(substr(tasdata$variable,7,8))

climatedata=cbind(tasdata,prdata[,3])
climatedata=climatedata[,-which(colnames(climatedata)=="variable")]
colnames(climatedata)=c("site","temp","year","month","precip")

climatedata <- as_tibble(climatedata)
sitedata <- as_tibble(sitedata)

dat=inner_join(climatedata,sitedata,by="site")
dat=dat[complete.cases(dat),]

cl=makeCluster(max_clusters)
clusterExport(cl,c("data","setorder"))
registerDoParallel(cl)

sites=unique(dat$site)
#loop through 1000 sites at a time -> 61 groups of sites
sitegroup=data.frame(site=sites,sitegroup=c(rep(1:60,each=1000),rep(61,97)))
dat=left_join(dat,sitegroup)
sitegroups=1:61

for(y in 1:length(sitegroups)){
  groupdat=dat%>%filter(sitegroup==sitegroups[y])
  group_pet <- pet_function(site=groupdat$site, year=groupdat$year, month=groupdat$month,
                             slope=groupdat$slope, latitude=groupdat$lat, aspect=groupdat$aspect,
                             tmean=groupdat$temp)
  group_pet <- group_pet[,c("site", "month", "year", "petm")]
  groupdat <- merge(groupdat, group_pet, by = c("site", "month", "year"))
  
  groupdat <- groupdat %>% 
    as_tibble() %>% 
    rename(tmean = temp, latitude = lat) %>% 
    arrange(site, year, month) %>% 
    group_by(site) %>% 
    nest() %>% 
    mutate(data = map(.x = data, .f = pet_spei_function)) %>% 
    unnest(data)  %>% 
    rename(temp = tmean, lat = latitude) %>% 
    mutate(site = as.character(site)) %>% 
    as.data.table()
  
  group_cwd <- cwd_function(site=groupdat$site, year=groupdat$year, month=groupdat$month,
                           petm = groupdat$petm, tmean=groupdat$temp,  
                           ppt = groupdat$precip, soilawc = groupdat$swc)
  
  group_cwd <- group_cwd[,c("site", "month", "year", "cwd")]
  groupdat <- merge(groupdat, group_cwd, by = c("site", "month", "year"))
  
  group_cwd_spei <- cwd_function(site=groupdat$site, year=groupdat$year, month=groupdat$month,
                            petm = groupdat$pet_spei, tmean=groupdat$temp,  
                            ppt = groupdat$precip, soilawc = groupdat$swc)
  
  group_cwd_spei <- group_cwd_spei[,c("site", "month", "year", "cwd")]
  names(group_cwd_spei) <- c("site", "month", "year", "cwd_spei")
  groupdat <- merge(groupdat, group_cwd_spei, by = c("site", "month", "year"))
  
  fwrite(groupdat,file=paste0(wdir,"1_input_processed/climate/cwd_group",sitegroups[y],".csv"))
  print(y)
}

for(i in 1:length(sitegroups)){
  temp=fread(paste0(wdir,"1_input_processed/climate/cwd_group",sitegroups[i],".csv"))
  #annual totals
  temp[,water_year:=year]
  temp[(lat>=0) & (month>=10),water_year:=year+1] # Northern hemisphere water year is october through september
  temp[(lat<0) & (month>=7),water_year:=year+1] # Southern hemisphere water year is July through June
  
  temp=temp%>%group_by(site,water_year)%>%summarize(pet=sum(petm),cwd=sum(cwd),ppt = sum(precip), tmean = mean(temp), pet_spei = sum(pet_spei), cwd_spei = sum(cwd_spei))
  temp=left_join(temp,sitedata%>%select(lon,lat,site))
  
  if(i==1) cwdhist=temp
  if(i>1) cwdhist=rbind(cwdhist,temp)
}

cwdhist <- cwdhist %>% rename(year = water_year)

cwdgrid=pivot_wider(cwdhist[,c("lat", "lon", "year", "cwd")],id_cols=c(lon,lat),names_from=year,values_from=cwd)
petgrid=pivot_wider(cwdhist[,c("lat", "lon", "year", "pet")],id_cols=c(lon,lat),names_from=year,values_from=pet)
pptgrid=pivot_wider(cwdhist[,c("lat", "lon", "year", "ppt")],id_cols=c(lon,lat),names_from=year,values_from=ppt)
tmpgrid=pivot_wider(cwdhist[,c("lat", "lon", "year", "tmean")],id_cols=c(lon,lat),names_from=year,values_from=tmean)
cwdspgrid=pivot_wider(cwdhist[,c("lat", "lon", "year", "cwd_spei")],id_cols=c(lon,lat),names_from=year,values_from=cwd_spei)
petspgrid=pivot_wider(cwdhist[,c("lat", "lon", "year", "pet_spei")],id_cols=c(lon,lat),names_from=year,values_from=pet_spei)


cwd_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
pet_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
ppt_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
tmp_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
cwdsp_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
petsp_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))


cwdgrid$cells=cellFromXY(cwd_historic,as.matrix(cwdgrid[,c(1,2)]))
petgrid$cells=cellFromXY(pet_historic,as.matrix(petgrid[,c(1,2)]))
pptgrid$cells=cellFromXY(ppt_historic,as.matrix(pptgrid[,c(1,2)]))
tmpgrid$cells=cellFromXY(tmp_historic,as.matrix(tmpgrid[,c(1,2)]))
cwdspgrid$cells=cellFromXY(cwdsp_historic,as.matrix(cwdspgrid[,c(1,2)]))
petspgrid$cells=cellFromXY(petsp_historic,as.matrix(petspgrid[,c(1,2)]))

cwdgrid=as.data.frame(cwdgrid)
petgrid=as.data.frame(petgrid)
pptgrid=as.data.frame(pptgrid)
tmpgrid=as.data.frame(tmpgrid)
petspgrid=as.data.frame(petspgrid)
cwdspgrid=as.data.frame(cwdspgrid)

cwd_historic[cwdgrid$cells]=cwdgrid[,grep(baseyears[1],colnames(cwdgrid))]
pet_historic[petgrid$cells]=petgrid[,grep(baseyears[1],colnames(petgrid))]
ppt_historic[pptgrid$cells]=pptgrid[,grep(baseyears[1],colnames(pptgrid))]
tmp_historic[tmpgrid$cells]=tmpgrid[,grep(baseyears[1],colnames(tmpgrid))]
petsp_historic[petspgrid$cells]=petspgrid[,grep(baseyears[1],colnames(petspgrid))]
cwdsp_historic[cwdspgrid$cells]=cwdspgrid[,grep(baseyears[1],colnames(cwdspgrid))]

for(i in 2:length(baseyears)){
  cwdtemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  pettemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  ppttemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  tmptemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  petsptemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  cwdsptemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  
  cwdtemp[cwdgrid$cells]=cwdgrid[,grep(baseyears[i],colnames(cwdgrid))]
  pettemp[petgrid$cells]=petgrid[,grep(baseyears[i],colnames(petgrid))]
  ppttemp[pptgrid$cells]=pptgrid[,grep(baseyears[i],colnames(pptgrid))]
  tmptemp[tmpgrid$cells]=tmpgrid[,grep(baseyears[i],colnames(tmpgrid))]
  petsptemp[petspgrid$cells]=petspgrid[,grep(baseyears[i],colnames(petspgrid))]
  cwdsptemp[cwdspgrid$cells]=cwdspgrid[,grep(baseyears[i],colnames(cwdspgrid))]
  
  cwd_historic=addLayer(cwd_historic,cwdtemp)
  pet_historic=addLayer(pet_historic,pettemp)
  ppt_historic=addLayer(ppt_historic,ppttemp)
  tmp_historic=addLayer(tmp_historic,tmptemp)
  petsp_historic=addLayer(petsp_historic,petsptemp)
  cwdsp_historic=addLayer(cwdsp_historic,cwdsptemp)
  
  print(i)
}
save(cwd_historic, pet_historic, ppt_historic, tmp_historic, petsp_historic, cwdsp_historic, file=paste0(wdir,"1_input_processed/climate/HistoricCWD_AETGrids_Annual.Rdat"))


