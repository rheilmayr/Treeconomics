#downscale CRU temperature and precip data based on WorldClim data
library(plyr)
library(raster)
library(sp)
library(rgdal)
library(reshape2)


franspath="C:/Users/fmoor/Box/Davis Stuff/Treeconomics"
climdat=read.csv(paste0(franspath,"/Data/cru_variables_1901-2016.csv"))
climdat=climdat[which(climdat$year%in%1970:2000),]
climdat$point=interaction(as.factor(climdat$latitude),as.factor(climdat$longitude),drop=TRUE)
baselines=ddply(climdat,.(month,point),summarize,pre_baseline=mean(pre),tmn_baseline=mean(tmn),tmx_baseline=mean(tmx),latitude=latitude[1],longitude=longitude[1])
baselines=baselines[,-2]
plots=SpatialPoints(unique(baselines[,6:5]))
crs(plots)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

downscaled=list()
vars=c("prec","tmax","tmin")
for(i in 1:length(vars)){
  downscaled[[i]]=matrix(nrow=length(plots),ncol=12)
  for(j in 1:12){
    ind=ifelse(j<10,paste0("0",j),j)
    wcdat=raster(paste0(franspath,"/Data/WorldClim Data for Downscaling/",vars[i],"/wc2.0_30s_",vars[i],"_",ind,".tif"))
    downscaled[[i]][,j]=extract(wcdat,plots)
    print(paste(i,j))
  }
  downscaled[[i]]=data.frame(downscaled[[i]])
  colnames(downscaled[[i]])=c(0:11)
  downscaled[[i]]$latitude=plots$latitude;downscaled[[i]]$longitude=plots$longitude
} 

dsdata=melt(downscaled[[1]],id.vars=c("latitude","longitude"))
for(i in 2:3) dsdata=cbind(dsdata,melt(downscaled[[i]],id.vars=c("latitude","longitude"))[,4])
colnames(dsdata)=c("latitude","longitude","month",vars)

baselines=merge(baselines,dsdata)
baselines$pre_correction=baselines$prec*10-baselines$pre_baseline #note CRU data is in tenths of a degree / mm while world clim is in degrees / mm 
baselines$tmax_correction=baselines$tmax*10-baselines$tmx_baseline
baselines$tmin_correction=baselines$tmin*10-baselines$tmn_baseline

write.csv(baselines,paste0(franspath,"/Data/downscalingcorrection.csv"))
