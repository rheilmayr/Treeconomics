#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/6/20
# Purpose: Combine site-level weather, topography and swc for cwd calcs
#
# Input files:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(plyr)
library(raster)
library(sp)
library(rgdal)
library(reshape2)
library(seegSDM)
library(data.table)
library(ncdf4)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set working directory
wdir <- 'remote\\'

# 1. Load site topography 
sites=fread(paste0(wdir, 'out//dendro//site_summary_slopeaspect.csv'))
sites <- sites %>% 
  dplyr::rename(elevation_itrdb = elevation,
                elevation = demelevation)
plots=unique(data.frame(latitude=sites$latitude,longitude=sites$longitude,site_id=sites$collection_id))
plots=SpatialPointsDataFrame(coords=plots[,c(2,1)],data=as.data.frame(plots[,3]))

# 2. Define directories for CRU and WorldClim data 
cru_dir <- paste0(wdir,"in/CRUData/")
wclim_dir <- paste0(wdir,"in/WorldClim/")

# 3. Load soil water capacity data
swc <- raster(paste0(wdir,"in\\wang_swc\\sr_cru_max.asc"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get CRU tmin, tmax and precip --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
monthyears=data.frame(year=rep(1901:2019,each=12),month=rep(1:12,length(1901:2019)))

vars=c("pre","tmn", 'tmx')
for(i in 1:length(vars)){
  crudat=stack(paste0(cru_dir, "cru_ts4.04.1901.2019.",vars[i],".dat.nc"))
  NAvalue(crudat)=-999
  tempdat=extract(crudat,plots)
  nas=which(is.na(tempdat[,1]));napoints=nearestLand(plots@coords[nas,],crudat[[1]],max_distance = 100000)
  tempdat[nas,]=extract(crudat,napoints)
  
  temp=cbind(monthyears,t(tempdat))
  colnames(temp)=c("year","month",as.character(plots@data[,1]))
  temp=melt(temp,id.vars=c("year","month"),variable.name="site_id",value.name=vars[i])
  if(i==1) climdat=temp; if(i>1) climdat=cbind(climdat,temp[,4])
  print(i)
}

colnames(climdat)[5:6]=vars[2:3]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add downscaling correction ---------------------------------------------
# Uses 1970-2000 means and compares to higher resolution World Clim data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baselines <- climdat %>%
  filter(year>1969 & year<2001) %>%
  group_by(site_id, month) %>%
  dplyr::summarize(pre_baseline = mean(pre),
                   tmn_baseline = mean(tmn),
                   tmx_baseline = mean(tmx))

downscaled=list()
vars=c("prec","tmax","tmin")
for(i in 1:length(vars)){
  downscaled[[i]]=matrix(nrow=length(plots),ncol=12)
  for(j in 1:12){
    ind=ifelse(j<10,paste0("0",j),j)
    wcdat=raster(paste0(wclim_dir,vars[i],"/wc2.0_30s_",vars[i],"_",ind,".tif"))
    downscaled[[i]][,j]=extract(wcdat,plots)
    if(i==1&j==1){
      nas=which(is.na(downscaled[[i]][,j]))
      napoints=nearestLand(plots@coords[nas,],wcdat,max_distance = 100000)
    }
    downscaled[[i]][nas,j]=extract(wcdat,napoints)
    print(paste(i,j))
  }
  downscaled[[i]]=data.frame(downscaled[[i]])
  colnames(downscaled[[i]])=c(1:12)
  downscaled[[i]]=as.data.frame(downscaled[[i]]);downscaled[[i]]$site_id=plots@data[,1]
} 

dsdata=melt(downscaled[[1]],id.vars="site_id",variable.name="month",value.name=vars[1])
for(i in 2:3) dsdata=cbind(dsdata,melt(downscaled[[i]],id.vars="site_id",variable.name="month",value.name=vars[i])[,3])
colnames(dsdata)=c("site_id","month",vars)

baselines=merge(baselines,dsdata)
baselines$pre_correction=baselines$prec-baselines$pre_baseline
baselines$tmax_correction=baselines$tmax-baselines$tmx_baseline
baselines$tmin_correction=baselines$tmin-baselines$tmn_baseline

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add soil water capacity ---------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotswc=extract(swc,plots)
nas=which(is.na(plotswc))
napoints=nearestLand(plots@coords[nas,],swc,max_distance = 100000)
plotswc[nas]=extract(swc,napoints)
plotswc=data.frame(swc=plotswc,site_id=plots@data[,1])

sites <- sites %>% 
  dplyr::rename(site_id = collection_id) %>% 
  left_join(plotswc, by = c("site_id"))
sites <- sites %>% 
  left_join(climdat,by=c("site_id"))
sites <- sites %>% 
  left_join(baselines[,c(1:2,9:11)],by=c("month","site_id"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Save combined topography, weather, and swc data ------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(sites,file=paste0(wdir,"out/climate/sitedataforcwd.Rdat"))
fwrite(sites,file=paste0(wdir,"out/climate/sitedataforcwd.csv"))

       