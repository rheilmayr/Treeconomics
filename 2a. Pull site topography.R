#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 7/9/20
# Purpose: Pull ASTER topography
#
# Input files:
#   site_summary.csv: Lat/Lon data for ITRDB sites. Created by "1b. Parse ITRDB.R"
#   site_summary_fia.csv: Lat/Lon data for FIA sites. Created by "1c. Parse Klesse.R"
#   Aster tiles: Elevation data tiles. Accessed from https://lpdaac.usgs.gov/products/astgtmv003/
# 
# Output files:
#   site_summary_slopeaspect.csv: Updated site summary with slope and aspect data.
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(raster)
select = dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set working directory
wdir <- 'remote//'

# 1. Load ITRDB site latitude and longitude
site_itrdb <- read_csv(paste0(wdir, 'out//dendro//site_summary.csv')) %>% 
  drop_na() %>% 
  mutate(plot_cn = collection_id) %>% 
  select(plot_cn, latitude, longitude, datasource)

# site <- site %>%
#   rename(elevation_itrdb = elevation)

# 2. Load FIA site latitude and longitude
site_fia <- read_csv(paste0(wdir, 'out//dendro//site_summary_fia.csv')) %>% 
  drop_na() %>% 
  mutate(datasource = "klesse_2018") %>% 
  select(plot_cn, latitude, longitude, datasource) %>% 
  distinct() # NOTE: Using plot_cn rather than collection_id since ~1/3 of collections are repeated locations but unique species.

# Combine sites
site <- rbind(site_itrdb, site_fia)

# 3. Set Directory for DEM tiles
demfiles=list.files("E:\\DEM Data",full.names=TRUE)
demfiles=demfiles[grep("ASTGT",demfiles)]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extract topography --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read in DEM for each site and extract slope and aspect (in radians) for each point
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Write data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site=cbind(site,slopeaspect)
colnames(site)[6:8]=c("slope","aspect","demelevation")

write.csv(site, file=paste0(wdir, 'out//dendro//site_summary_slopeaspect.csv'))
