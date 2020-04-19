#load in site lat longs
library(plyr)
library(raster)
library(sp)
library(rgdal)
library(reshape2)
library(ncdf4)
library(downloader)

franspath="C:/Users/fmoore/Box/Davis Stuff/Treeconomics"
plots=read.csv(paste0(franspath,"/Data/plot level records - final.csv"))
plotspoints=SpatialPoints(unique(plots[,c(2,1)]))
crs(plotspoints)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

tasfiles=list.files(paste0(franspath,"/Data/CMIP5 Data/tas"))
prfiles=list.files(paste0(franspath,"/Data/CMIP5 Data/pr"))

years=list(1970:2000,2045:2055,2090:2100)
yearmonths=rep(1861:2100,each=12)

for(i in 1:length(tasfiles)){
  tasfile=brick(paste0(franspath,"/Data/CMIP5 Data/tas/",tasfiles[i]))
  prfile=brick(paste0(franspath,"/Data/CMIP5 Data/pr/",prfiles[i]))
  for(k in 1:length(years)){
    tas_period=tasfile[[which(yearmonths%in%years[[k]])]]
    pr_period=prfile[[which(yearmonths%in%years[[k]])]]
    tas_months=calc(tas_period[[which(1:dim(tas_period)[3]%%12==1)]],function(x) mean(x, na.rm=T));pr_months=calc(pr_period[[which(1:dim(pr_period)[3]%%12==1)]],function(x) mean(x, na.rm=T))
    for(j in 2:12){
      tas_months=stack(tas_months,calc(tas_period[[which(1:dim(tas_period)[3]%%12==ifelse(j==12,0,j))]],function(x) mean(x, na.rm=T)))
      pr_months=stack(pr_months,calc(pr_period[[which(1:dim(pr_period)[3]%%12==ifelse(j==12,0,j))]],function(x) mean(x, na.rm=T)))
    }
    if(k==1){
      save(tas_months,file=paste0(franspath,"/Data/CMIP5 Data/monthly tas rasters/start/tas_",i,"_1970_2000.Rdat"))
      save(pr_months,file=paste0(franspath,"/Data/CMIP5 Data/monthly pr rasters/start/pr_",i,"_1970_2000.Rdat"))
    }
    if(k==2){
      save(tas_months,file=paste0(franspath,"/Data/CMIP5 Data/monthly tas rasters/mid/tas_",i,"_2045_2055.Rdat"))
      save(pr_months,file=paste0(franspath,"/Data/CMIP5 Data/monthly pr rasters/mid/pr_",i,"_2045_2055.Rdat"))
    }
    if(k==3){
      save(tas_months,file=paste0(franspath,"/Data/CMIP5 Data/monthly tas rasters/end/tas_",i,"_2090_2100.Rdat"))
      save(pr_months,file=paste0(franspath,"/Data/CMIP5 Data/monthly pr rasters/end/pr_",i,"_2090_2100.Rdat"))
    }
  }
  print(i)
}

#get soil, slope, latitude, elevation for cwd calculation

swc=raster(paste0(franspath,"/Data/CMIP5 Data/other data for cwd/sr_cru_max.asc"))
labs=read.delim(paste0(franspath,"/Data/CMIP5 Data/other data for cwd/demdatalinks.txt"))
for(i in 1:dim(labs)[1]){
  download.file(as.character(labs[i,]), destfile=paste0(franspath,"/Data/CMIP5 Data/other data for cwd/DEM Data/",substr(as.character(labs[i,]),79,85),".nc"), method="auto",mode="wb",extra="--load-cookies C:\\Users\\fmoore\\Desktop\\cookies.txt --save-cookies C:\\.urs_cookies --auth-no-challenge=on --keep-session-cookies --user=fmoore125 --ask-password --content-disposition")
}


testdem=raster(paste0(franspath,"/Data/CMIP5 Data/other data for cwd/DEM Data/ASTGTMV003_N22E010_dem.nc"))
