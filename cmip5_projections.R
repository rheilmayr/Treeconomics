#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 
# Purpose: Calculate CWD / PET based on CMIP5 data
#
# Input files:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
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
library(tidyverse)
source("cwd_function.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create CWD projections --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

tasfiles=list.files(paste0(wdir,"in/CMIP5 Data/tas"),pattern=".nc")
prfiles=list.files(paste0(wdir,"in/CMIP5 Data/pr"),pattern=".nc")

years=list(1970:2000,2045:2055,2090:2100)
yearmonths=rep(1861:2100,each=12)

for(i in 1:length(tasfiles)){
  tasfile=brick(paste0(wdir,"in/CMIP5 Data/tas/",tasfiles[i]))
  prfile=brick(paste0(wdir,"in/CMIP5 Data/pr/",prfiles[i]))
  for(k in 1:length(years)){
    tas_period=tasfile[[which(yearmonths%in%years[[k]])]]
    pr_period=prfile[[which(yearmonths%in%years[[k]])]]
    tas_months=calc(tas_period[[which(1:dim(tas_period)[3]%%12==1)]],function(x) mean(x, na.rm=T));pr_months=calc(pr_period[[which(1:dim(pr_period)[3]%%12==1)]],function(x) mean(x, na.rm=T))
    for(j in 2:12){
      tas_months=stack(tas_months,calc(tas_period[[which(1:dim(tas_period)[3]%%12==ifelse(j==12,0,j))]],function(x) mean(x, na.rm=T)))
      pr_months=stack(pr_months,calc(pr_period[[which(1:dim(pr_period)[3]%%12==ifelse(j==12,0,j))]],function(x) mean(x, na.rm=T)))
    }
    if(k==1){
      save(tas_months,file=paste0(wdir,"in/CMIP5 Data/monthly tas rasters/start/tas_",i,"_1970_2000.Rdat"))
      save(pr_months,file=paste0(wdir,"in/CMIP5 Data/monthly pr rasters/start/pr_",i,"_1970_2000.Rdat"))
    }
    if(k==2){
      save(tas_months,file=paste0(wdir,"in/CMIP5 Data/monthly tas rasters/mid/tas_",i,"_2045_2055.Rdat"))
      save(pr_months,file=paste0(wdir,"in/CMIP5 Data/monthly pr rasters/mid/pr_",i,"_2045_2055.Rdat"))
    }
    if(k==3){
      save(tas_months,file=paste0(wdir,"in/CMIP5 Data/monthly tas rasters/end/tas_",i,"_2090_2100.Rdat"))
      save(pr_months,file=paste0(wdir,"in/CMIP5 Data/monthly pr rasters/end/pr_",i,"_2090_2100.Rdat"))
    }
  }
  print(i)
}


#get soil, slope, latitude, elevation for cwd calculation

swc=raster(paste0(wdir,"in/CMIP5 Data/other data for cwd/sr_cru_max.asc"))
#convert swc from mm to cm
swc=swc/10

load(paste0(wdir,"in/CMIP5 Data/other data for cwd/slopeaspectraster.Rdat"))

terrainraster=crop(terrainraster,extent(swc))
#convert slope and aspect to degrees
terrainraster=terrainraster*180/pi

#melt permanent site data down into a long data frame for cwd function
sitedata=cbind(as.data.frame(terrainraster),as.data.frame(swc))
sitedata=cbind(sitedata,coordinates(swc))
colnames(sitedata)=c("slope","aspect","swc","lon","lat")
sitedata$site=1:dim(sitedata)[1]


#need to get correction of CMIP data using WorldClim 1970-2000 climatology
precipfiles=list.files(paste0(wdir,"in/WorldClim/Lower Resolution/precip"),pattern=".tif")
tmeanfiles=list.files(paste0(wdir,"in/WorldClim/Lower Resolution/tmean"),pattern=".tif")

for(i in 1:12){
  p=raster(paste0(wdir,"in/WorldClim/Lower Resolution/precip/",precipfiles[i]))
  tmean=raster(paste0(wdir,"in/WorldClim/Lower Resolution/tmean/",tmeanfiles[i]))

  if(i==1){precipclim=p;tempclim=tmean}
  if(i>1) {precipclim=stack(precipclim,p);tempclim=stack(tempclim,tmean)}
}

precipclim=crop(precipclim,extent(swc));tempclim=crop(tempclim,extent(swc))
precipclim=aggregate(precipclim,fact=3);tempclim=aggregate(tempclim,fact=3)


#get model-specific correction factors based on 1970-2000 climatology
tempmodelfiles=list.files(paste0(wdir,"in/CMIP5 Data/monthly tas rasters/start"),pattern=".Rdat")
precipmodelfiles=list.files(paste0(wdir,"in/CMIP5 Data/monthly pr rasters/start"),pattern=".Rdat")

pr_correction=list();tas_correction=list()
for(i in 1:length(tempmodelfiles)){
  load(paste0(wdir,"in/CMIP5 Data/monthly tas rasters/start/",tempmodelfiles[i]))
  load(paste0(wdir,"in/CMIP5 Data/monthly pr rasters/start/",precipmodelfiles[i]))
  tas_months=rotate(tas_months);pr_months=rotate(pr_months)
  tas_months=resample(tas_months,swc);pr_months=resample(pr_months,swc)
  
  #convert units - temperatre from kelvin to degrees, precip from kg per m2 per second to mm per month
  tas_months=tas_months-273.15; pr_months=pr_months*60*60*24*30
  
  #get grid-cell level correction based on 1970-2000 climatology
  pr_correction[[i]]=precipclim-pr_months
  tas_correction[[i]]=tempclim-tas_months
  print(i)
}

save(sitedata,pr_correction,tas_correction,file=paste0(wdir,"in/CMIP5 Data/other data for cwd/sitedata_climatologycorrection.Rdat"))

datfolder=paste0(wdir, "in/CMIP5 Data/")

period=c("start","mid","end")

swc=raster(paste0(wdir,"in/CMIP5 Data/other data for cwd/sr_cru_max.asc"))
load(paste0(wdir,"in/CMIP5 Data/other data for cwd/sitedata_climatologycorrection.Rdat"))
sitedata$slope[which(sitedata$slope<0)]=0 # fix a few suprious slope values

cl=makeCluster(20)
clusterExport(cl,c("data","setorder"))
registerDoParallel(cl)

for(i in 1:length(period)){
  print(period[i])
  tempmodelfiles=list.files(paste0(datfolder,"monthly tas rasters/",period[i]),pattern=".Rdat")
  precipmodelfiles=list.files(paste0(datfolder,"monthly pr rasters/",period[i]),pattern=".Rdat")
  for(j in 1:length(tempmodelfiles)){
    load(paste0(datfolder,"monthly tas rasters/",period[i],"/",tempmodelfiles[j]))
    load(paste0(datfolder,"monthly pr rasters/",period[i],"/",precipmodelfiles[j]))
    tas_months=rotate(tas_months);pr_months=rotate(pr_months)
    tas_months=resample(tas_months,swc);pr_months=resample(pr_months,swc)
    
    #convert units - temperatre from kelvin to degrees, precip from kg per m2 per second to mm per month
    tas_months=tas_months-273.15; pr_months=pr_months*60*60*24*30
    
    #apply grid-cell specific bias correction
    tas_months=tas_months+tas_correction[[j]];pr_months=pr_months+pr_correction[[j]]
    
    tasdata=as.data.frame(tas_months);prdata=as.data.frame(pr_months)
    colnames(tasdata)=1:12;colnames(prdata)=1:12
    tasdata$site=1:dim(tasdata)[1];prdata$site=1:dim(tasdata)[1]
    tasdata=reshape2::melt(tasdata,id.vars="site");prdata=reshape2::melt(prdata,id.vars="site")
    climatedata=cbind(tasdata,prdata[,3])
    colnames(climatedata)=c("site","month","temp","precip")
    
    dat=merge(sitedata,climatedata)
    dat=dat[complete.cases(dat),]

    test=cwd_function(site=as.factor(dat$site),slope=dat$slope,latitude=dat$lat,foldedaspect=dat$aspect,ppt=dat$precip,tmean=dat$temp,month=dat$month,soilawc=dat$swc,year=NULL,type="normal")
    fwrite(test,file=paste0(datfolder,"cwd calcs/",period[i],"_",j,".csv"))
    print(j)
  }
}



### --------Put CWD / AET Data into grids ----------------

cwddir=datfolder
cwdfiles=list.files(paste0(cwddir,"/cwd calcs/"),full.names=TRUE)

load(file=paste0(datfolder, "other data for cwd/sitedata_climatologycorrection.Rdat"))

climdat=as.data.frame(tas_correction[[1]][[1]])
climdat$site=1:dim(climdat)[1]
sitedata=merge(sitedata,climdat) 
sitedata=sitedata[complete.cases(sitedata),]

sitecrosswalk=data.frame(site_grid=unique(sitedata$site),site=1:length(unique(sitedata$site)))

period=c("start","end")

for(i in 1:length(period)){
  periodfiles=cwdfiles[grep(period[i],cwdfiles)]
  for(j in 1:length(periodfiles)){
    cwd_raster_temp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
    aet_raster_temp=raster(nrow=nrow(swc),ncol=ncol(swc),ext=extent(swc))
    
    cwddat=fread(periodfiles[j])
    cwddat=cwddat%>%
      select(site,month,aet,cwd)%>%
      group_by(site)%>%
      summarize(aet=sum(aet),cwd=sum(cwd))
    
    cwddat=merge(cwddat,sitecrosswalk)
    cwddat=merge(cwddat,sitedata[,c(1,5,6)],by.x="site_grid",by.y="site")
    cwddat$cells=cellFromXY(cwd_raster_temp,cwddat[,c(5,6)])
    
    cwd_raster_temp[cwddat$cells]=cwddat$cwd;aet_raster_temp[cwddat$cells]=cwddat$aet
    if(j==1){cwd_raster=cwd_raster_temp;aet_raster=aet_raster_temp}
    if(j>1){
      cwd_raster=stack(cwd_raster,cwd_raster_temp)
      aet_raster=stack(aet_raster,aet_raster_temp)
    }
  }
  print(period[i])
  save(aet_raster,cwd_raster,file=paste0(wdir,"/in/CMIP5 CWD/cmip5_cwdaet_",period[i],".Rdat"))
}



