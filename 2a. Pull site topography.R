#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 7/9/20
# Purpose: Pull ASTER topography
#
# Input files:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(raster)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set working directory
wdir <- 'remote//'

# 1. Load site latitude and longitude
site <- read_csv(paste0(wdir, 'out//dendro//site_summary.csv'))
site <- site %>% 
  rename(elevation_itrdb = elevation)

#2. Set Directory for DEM tiles
demfiles=list.files("E:\\DEM Data",full.names=TRUE)
demfiles=demfiles[grep("ASTGT",demfiles)]

#3. Read in DEM for each site and extract slope and aspect (in radians) for each point
for(i in 1:dim(site)[1]){
  #get dem granule for each site location
  lat=site$latitude[i];lon=site$longitude[i]
  demlat=ifelse(lat>0,floor(lat),-1*floor(lat))
  demlon=ifelse(lon>0,floor(lon),-1*floor(lon))
  demlat=formatC(demlat,width=2,flag="0");demlon=formatC(demlon,width=3,flag="0")
  demlatdir=ifelse(lat>0,"N","S");demlondir=ifelse(lon>0,"E","W")
  demname=paste0(demlatdir,demlat,demlondir,demlon)
  demfile=grep(demname,demfiles)
  if(length(demfile)==0){ #some sites don't have dems as on islands
    siteslopeaspect=c(NA,NA,NA)
  }
  if(length(demfile)>0){
    dem=raster(demfiles[demfile])
    ter=terrain(dem,opt=c("slope","aspect"))
    if(lon%%1==0) lon=lon+0.001;if(lat%%1==0) lat=lat+0.001 #small adjustment for points right on the edge of the dem granule
    siteslopeaspect=c(extract(ter,data.frame(lon=lon,lat=lat)),extract(dem,data.frame(lon=lon,lat=lat)))
  }
  if(i==1) slopeaspect=siteslopeaspect
  if(i>1) slopeaspect=rbind(slopeaspect,siteslopeaspect)
  print(i)
}

site=cbind(site,slopeaspect)
colnames(site)[13:15]=c("slope","aspect","demelevation")

write.csv(site, file=paste0(wdir, 'out//dendro//site_summary_slopeaspect.csv'))
