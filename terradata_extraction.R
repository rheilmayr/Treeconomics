library(data.table)
library(raster)
library(reshape2)
library(terra)
library(ncdf4)
library(tidyr)
wdir <- 'remote/'

vars=c("def","pet","ppt")

years=1958:2020

locs=read.csv("out/dendro/site_summary.csv")
lonlat=data.frame(lon=locs[,which(colnames(locs)=="longitude")],lat=locs[,which(colnames(locs)=="latitude")])

for(var in vars){
  for(i in years){
    print(paste(var,i))
    file=paste0("TerraClimate_Feb2024/",var,"/TerraClimate_",var,"_",i,".nc")
    dat=terra::rast(file)
    temp=extract(dat,lonlat,ID=FALSE)
    colnames(temp)=1:12;temp$collection_id=locs$collection_id
    if(i==1958){
      longdat=reshape2::melt(temp)
      longdat$year=i
      next
    }
    temp=reshape2::melt(temp);temp$year=i
    longdat=rbind(longdat,temp)
  }
  colnames(longdat)[2:3]=c("Month",var)
  fwrite(longdat,file=paste0("TerraClimate_Feb2024/itrdbsites_",var,".csv"))
}

climatologies=c("19611990","19812010")

for(i in climatologies){
  for (var in vars){
    dat=terra::rast(paste0("TerraClimate_Feb2024/TerraClimate",i,"_",var,".nc"))
    temp=extract(dat,lonlat,ID=FALSE)
    colnames(temp)=1:12;temp$collection_id=locs$collection_id
    if(var=="def"){
      longdat=reshape2::melt(temp)
      longdat$var=var
      next
    } 
    temp=reshape2::melt(temp)
    temp$var=var
    longdat=rbind(longdat,temp)
  }
  fwrite(longdat,file=paste0("TerraClimate_Feb2024/itrdbsites_defpetppt_",i,".csv"))
}

