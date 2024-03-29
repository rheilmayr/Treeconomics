#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 7/10/20
# Purpose: Create gridded 0.5 degree cwd, pet and aet based on 1901-1980 climatological baseline (CRU Data)
#
# Input files:
#   cru_ts4.04.1901.2019.tmp.dat.nc: Monthly mean temp data from CRU. FRAN - WHERE DOES THIS GET CALCULATED?
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
# TODO:
#   FRAN: Should we be calculating CWD / AET for each year first, then taking the average of that for the baseline?
#   FRAN: Would be great if you could annotate a bit more. Why are some datasets labeled test? What are "sites" mentioned below?
#   FRAN: Should we clean up outputs? At least move into a single directory?
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
source("f_cwd_function.R")

max_clusters = 20

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wdir <- "remote/"
# franspath="C:/Users/fmoore/Box/Davis Stuff/Treeconomics"

## Load CRU temp and precip data
cruyears=rep(1901:2019,each=12);crumonths=rep(1:12,length(1901:2019))
tas=stack(paste0(wdir, "in/CRUData/cru_ts4.04.1901.2019.tmp.dat.nc"))
pr=stack(paste0(wdir,"in/CRUData/cru_ts4.04.1901.2019.pre.dat.nc"))

## Load soil water capacity raster
swc=raster(paste0(wdir,"in/CMIP5 Data/other data for cwd/sr_cru_max.asc"))
#convert swc from mm to cm
swc=swc/10

## Load topography raster
load(paste0(wdir,"in/CMIP5 Data/other data for cwd/slopeaspectraster.Rdat"))
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

tas_months=calc(tas[[which(basemonths==1)]],function(x) mean(x, na.rm=T));pr_months=calc(pr[[which(basemonths==1)]],function(x) mean(x, na.rm=T))
for(j in 2:12){
  tas_months=stack(tas_months,calc(tas[[which(basemonths==j)]],function(x) mean(x, na.rm=T)))
  pr_months=stack(pr_months,calc(pr[[which(basemonths==j)]],function(x) mean(x, na.rm=T)))
  print(j)
}


# save(tas_months,pr_months,file=paste0(wdir,"in/CRUData/monthlycrubaseline.Rdat"))
writeRaster(tas_months %>% mean(), file=paste0(wdir,"in/CRUData/monthlycrubaseline_tas"), overwrite = TRUE)
writeRaster(pr_months %>% sum(), file=paste0(wdir,"in/CRUData/monthlycrubaseline_pr"), overwrite = TRUE)


# melt permanent site data down into a long data frame for cwd function
sitedata=cbind(as.data.frame(terrainraster),as.data.frame(swc))
sitedata=cbind(sitedata,coordinates(swc))
colnames(sitedata)=c("slope","aspect","swc","lon","lat")
sitedata$site=1:dim(sitedata)[1]

# fix a few spurious slope values
sitedata$slope[which(sitedata$slope<0)]=0 

# crop temp and precip data to swc extent
tas_months=crop(tas_months,extent(swc));pr_months=crop(pr_months,extent(swc))

# set up parallel processing
cl=makeCluster(4)
clusterExport(cl,c("data","setorder"))
registerDoParallel(cl)

tasdata=as.data.frame(tas_months);prdata=as.data.frame(pr_months)
colnames(tasdata)=1:12;colnames(prdata)=1:12
tasdata$site=1:dim(tasdata)[1];prdata$site=1:dim(tasdata)[1]
tasdata=melt(tasdata,id.vars="site");prdata=melt(prdata,id.vars="site")
climatedata=cbind(tasdata,prdata[,3])
colnames(climatedata)=c("site","month","temp","precip")

dat=merge(sitedata,climatedata)
dat=dat[complete.cases(dat),]

test=cwd_function(site=as.factor(dat$site),slope=dat$slope,latitude=dat$lat,foldedaspect=dat$aspect,ppt=dat$precip,tmean=dat$temp,month=dat$month,soilawc=dat$swc,year=NULL,type="normal")
fwrite(test,file=paste0(wdir,"/Data/griddedbaselinecwddata.csv"))

#merge in lat longs from site data
sitecrosswalk=data.frame(site_grid=unique(dat$site),site=1:length(unique(dat$site)))
test=merge(test,sitecrosswalk)

test$site_grid=as.numeric(test$site_grid)
test=merge(test[,c(3,28,29,30)],sitedata[,c(4:6)],by.x="site_grid",by.y="site")
test$cells=cellFromXY(cwd_historic,test[,c(5,6)])

testmonths=split(test,test$month)

cwd_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
aet_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))

cwd_historic[testmonths[[1]]$cells]=testmonths[[1]]$cwd
aet_historic[testmonths[[1]]$cells]=testmonths[[1]]$aet

for(j in 2:12){
  cwdtemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc));aettemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  cwdtemp[testmonths[[j]]$cells]=testmonths[[j]]$cwd
  aettemp[testmonths[[j]]$cells]=testmonths[[j]]$aet
  cwd_historic=stack(cwd_historic,cwdtemp);aet_historic=stack(aet_historic,aettemp)
}

save(cwd_historic,aet_historic,file=paste0(wdir,"in/CRUData/historic_raster/HistoricCWD_AETGrids.Rdat"))


tasdata=melt(tasdata,id.vars=c("site","year"));prdata=melt(prdata,id.vars=c("site","year"))
climatedata=cbind(tasdata,prdata[,4])
colnames(climatedata)=c("site","year","month","temp","precip")

climatedata=as.tbl(climatedata);sitedata=as.tbl(sitedata)

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
  test=cwd_function(site=as.factor(groupdat$site),slope=groupdat$slope,latitude=groupdat$lat,foldedaspect=groupdat$aspect,ppt=groupdat$precip,tmean=groupdat$temp,month=groupdat$month,soilawc=groupdat$swc,year=groupdat$year,type="annual")
  fwrite(test,file=paste0(wdir,"in/CRUData/historic_raster/cwd_group",sitegroups[y],".csv"))
  print(y)
}

for(i in 1:length(sitegroups)){
  temp=fread(paste0(wdir,"in/CRUData/historic_raster/cwd_group",sitegroups[i],".csv"))
  #annual totals
  temp=temp%>%group_by(site,year)%>%summarize(aet=sum(aet),cwd=sum(cwd))
  temp=left_join(temp,sitedata%>%select(lon,lat,site))
  
  if(i==1) cwdhist=temp
  if(i>1) cwdhist=rbind(cwdhist,temp)
}

cwdgrid=pivot_wider(cwdhist[,c(2,4,5,6)],id_cols=c(lon,lat),names_from=year,values_from=cwd)
aetgrid=pivot_wider(cwdhist[,c(2,3,5,6)],id_cols=c(lon,lat),names_from=year,values_from=aet)

cwd_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
aet_historic=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))

cwdgrid$cells=cellFromXY(cwd_historic,as.matrix(cwdgrid[,c(1,2)]))
aetgrid$cells=cellFromXY(aet_historic,as.matrix(aetgrid[,c(1,2)]))

cwdgrid=as.data.frame(cwdgrid);aetgrid=as.data.frame(aetgrid)

cwd_historic[cwdgrid$cells]=cwdgrid[,grep(baseyears[1],colnames(cwdgrid))]
aet_historic[aetgrid$cells]=aetgrid[,grep(baseyears[1],colnames(aetgrid))]

for(i in 2:length(baseyears)){
  cwdtemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  aettemp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
  cwdtemp[cwdgrid$cells]=cwdgrid[,grep(baseyears[i],colnames(cwdgrid))]
  aettemp[aetgrid$cells]=aetgrid[,grep(baseyears[i],colnames(aetgrid))]
  cwd_historic=addLayer(cwd_historic,cwdtemp)
  aet_historic=addLayer(aet_historic,aettemp)
  print(i)
}
save(cwd_historic,aet_historic,file=paste0(wdir,"in/CRUData/historic_raster/HistoricCWD_AETGrids_Annual.Rdat"))


