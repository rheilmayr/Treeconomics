#get gridded 0.5 degree cwd and aet based on 1901-1980 climatolgoical baseline (CRU Data)
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

franspath="C:/Users/fmoore/Box/Davis Stuff/Treeconomics"

cruyears=rep(1901:2019,each=12);crumonths=rep(1:12,length(1901:2019))

tas=stack(paste0(franspath,"/Data/CRUData/cru_ts4.04.1901.2019.tmp.dat.nc"))
pr=stack(paste0(franspath,"/Data/CRUData/cru_ts4.04.1901.2019.pre.dat.nc"))

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

save(tas_months,pr_months,file=paste0(franspath,"/Data/CRUData/monthlycrubaseline.Rdat"))

#get soil, slope, latitude, elevation for cwd calculation

swc=raster(paste0(franspath,"/Data/CMIP5 Data/other data for cwd/sr_cru_max.asc"))
#convert swc from mm to cm
swc=swc/10

load(paste0(franspath,"/Data/CMIP5 Data/other data for cwd/slopeaspectraster.Rdat"))

terrainraster=crop(terrainraster,extent(swc))
#convert slope and aspect to degrees
terrainraster=terrainraster*180/pi

#melt permanent site data down into a long data frame for cwd function
sitedata=cbind(as.data.frame(terrainraster),as.data.frame(swc))
sitedata=cbind(sitedata,coordinates(swc))
colnames(sitedata)=c("slope","aspect","swc","lon","lat")
sitedata$site=1:dim(sitedata)[1]

sitedata$slope[which(sitedata$slope<0)]=0 # fix a few suprious slope values

#crop temp and precip data to swc extent
tas_months=crop(tas_months,extent(swc));pr_months=crop(pr_months,extent(swc))

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
fwrite(test,file=paste0(franspath,"/Data/griddedbaselinecwddata.csv"))

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

save(cwd_historic,aet_historic,file=paste0(franspath,"/Data/HistoricCWD_AETGrids.Rdat"))


#version with annual cwwd / pet instead of climatology
tas=crop(tas,extent(swc));pr=crop(pr,extent(swc))

for(i in 1:length(baseyears)){
  tempt=as.data.frame(dropLayer(tas,which(years!=baseyears[i])))
  tempp=as.data.frame(dropLayer(pr,which(years!=baseyears[i])))
  colnames(tempt)=1:12;colnames(tempp)=1:12
  tempt$year=rep(baseyears[i],dim(tempt)[1]);tempp$year=rep(baseyears[i],dim(tempp)[1])
  tempt$site=1:dim(tempt)[1];tempp$site=1:dim(tempp)[1]
  
  if(i==1){tasdata=tempt;prdata=tempp}
  if(i>1){
    tasdata=rbind(tasdata,tempt)
    prdata=rbind(prdata,tempp)
  }
  print(i)
}


tasdata=melt(tasdata,id.vars=c("site","year"));prdata=melt(prdata,id.vars=c("site","year"))
climatedata=cbind(tasdata,prdata[,4])
colnames(climatedata)=c("site","year","month","temp","precip")

climatedata=as.tbl(climatedata);sitedata=as.tbl(sitedata)

dat=inner_join(climatedata,sitedata,by="site")
dat=dat[complete.cases(dat),]

cl=makeCluster(20)
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
  fwrite(test,file=paste0(franspath,"/Data/Baseline CWD/cwd_group",sitegroups[y],".csv"))
  print(y)
}

for(i in 1:length(sitegroups)){
  temp=fread(paste0(franspath,"/Data/Baseline CWD/cwd_group",sitegroups[i],".csv"))
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
save(cwd_historic,aet_historic,file=paste0(franspath,"/Data/HistoricCWD_AETGrids_Annual.Rdat"))

cwd_function <- function(site,slope,latitude,foldedaspect,ppt,tmean,month,soilawc=200,year=NULL,type=c('normal','annual')) {
  
  # Packages
  require(data.table)
  require(geosphere)
  
  #R SCRIPT
  if (type =="annual"){
    data<-as.data.table(cbind(as.character(site),slope,latitude,foldedaspect,ppt,tmean,month,year))[,soilawc:=soilawc]
    colnames(data)<-c("site","slope","latitude","foldedaspect","ppt","tmean","month","year","soilawc")}
  
  if (type =="normal"){
    data<-as.data.table(cbind(site,slope,latitude,foldedaspect,ppt,tmean,month))[,soilawc:=soilawc]
    colnames(data)<-c("site","slope","latitude","foldedaspect","ppt","tmean","month","soilawc")
  }
  data$site<-as.character(data$site)
  data$slope<-as.numeric(as.character(data$slope))
  data$latitude<-as.numeric(as.character(data$latitude))
  data$foldedaspect<-as.numeric(as.character(data$foldedaspect))
  data$ppt<-as.numeric(as.character(data$ppt))
  data$tmean<-as.numeric(as.character(data$tmean))
  data$month<-as.numeric(as.character(data$month))
  
  ### for 30 year normal data only, we'll create iterations so that way the water storage can carry over from one year to the next
  if (type == "normal"){year<-c(rep(1,length(data$site)),rep(2,length(data$site)),rep(3,length(data$site)),rep(4,length(data$site)),rep(5,length(data$site)),rep(6,length(data$site)),rep(7,length(data$site)),rep(8,length(data$site)),rep(9,length(data$site)),rep(10,length(data$site)))
  data<-rbind(data,data,data,data,data,data,data,data,data,data)}
  data$year<-as.numeric(as.character(year))
  # calculate daylength
  
  daylength<-vector()
  datasite<-data[,.(latitude=mean(latitude)),by=.(site)]
  for (i in 1:length(datasite$latitude)){
    dl<-daylength(datasite$latitude[i],1:365)
    day<- tapply(dl, rep(1:12, c(31,28,31,30,31,30,31,31,30,31,30,31)), mean)
    site <- as.vector(rep(datasite$site[i],length(day)))
    month<-as.vector(c(1,2,3,4,5,6,7,8,9,10,11,12))
    join<-cbind(site,month,day)
    daylength<-rbind(daylength,join)
  }
  daylength<-as.data.frame(daylength)
  daylength$site<-as.character(daylength$site)
  daylength$month<-as.numeric(as.character(daylength$month))
  daylength$day<-as.numeric(as.character(daylength$day))
  data<-merge(data,daylength,by=c("site","month"))
  
  ##########data$year<-year
  # calculate all other variables
  data[,fm:=ifelse(tmean<0,0,ifelse(tmean>=6,1,tmean*0.166666666))]
  data[,rainm:=fm*ppt]
  data[,snowm:=(1-fm)*ppt]
  data<-setorder(data,site,year,month)
  sites<-unique(data$site)
  mergedata<-vector()
  
  mergedata=foreach(i=1:length(sites),.combine="rbind")%dopar%{
    packm<-vector()
    sitedat<-setorder(data,site,year,month)[site==sites[i],]
    for (j in 1:length(sitedat$site)){
      packmsite<-ifelse(j>1,(((1-sitedat$fm[j])^2)*sitedat$ppt[j])+((1-sitedat$fm[j])*packm[j-1]),(((1-sitedat$fm[j])^2)*sitedat$ppt[j]))
      packm<-c(packm,packmsite)
    }
    site<-as.character(sitedat$site)
    year<-as.numeric(as.character(sitedat$year))
    month<-as.numeric(as.character(sitedat$month))
    mergedat<-cbind(site,year,month,packm)
    #mergedata<-rbind(mergedata,mergedat)
    #print(paste0(i,"_1"))
  }
  mergedt<-as.data.table(mergedata)
  mergedt$year<-as.numeric(as.character(mergedt$year))
  mergedt$month<-as.numeric(as.character(mergedt$month))
  mergedt$site<-as.character(mergedt$site)
  data<-merge(data,mergedt,by=c("site","year","month"))
  data$packm<-as.numeric(as.character(data$packm))
  data[,meltm:=fm*snowm*data.table::shift(packm,1L,type="lag",fill=0)]
  data[,wm:=rainm+meltm]
  data[month==1 | month==3 |month==5|month==7|month==8|month==10|month==12,days:=31]
  data[month==4 | month==6 |month==9|month==11,days:=30]
  data[month==2,days:=28]
  data[,ea:=exp(((17.3*tmean)/(tmean+273.2)))*0.611]
  # convert slope, folded aspect, and latitude to radians
  data[,sloprad:=slope*0.0174532925]
  data[,afrad:=foldedaspect*0.0174532925]
  data[,latrad:=latitude*0.0174532925]
  #calculate heat load
  data[,heatload:=0.339+0.808*(cos(latrad)*cos(sloprad))-0.196*(sin(latrad)*sin(sloprad))-0.482*(cos(afrad)*sin(sloprad))]
  
  data[,petm:=ifelse(tmean<0,0,((((ea*tmean)/(tmean+273.3))*day*days*29.8)*heatload/10))]
  
  mergedata=foreach(i=1:length(sites),.combine="rbind")%dopar%{
    soilm<-vector()
    sitedat<-setorder(data,site,year,month)[site==sites[i],]
    for (j in 1:length(sitedat$site)){
      soilmsite<-ifelse(j>1,ifelse((sitedat$wm[j]-sitedat$petm[j]+soilm[j-1])<=0,0,ifelse((sitedat$wm[j]-sitedat$petm[j]+soilm[j-1])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]+soilm[j-1]),sitedat$soilawc[j])),ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,0,ifelse((sitedat$wm[j]-sitedat$petm[j])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]),sitedat$soilawc[j])))
      soilm<-c(soilm,soilmsite)
    }
    site<-as.character(sitedat$site)
    year<-as.numeric(as.character(sitedat$year))
    month<-as.numeric(as.character(sitedat$month))
    mergedat<-cbind(site,year,month,soilm)
  }
  mergedt<-as.data.table(mergedata)
  mergedt$year<-as.numeric(as.character(mergedt$year))
  mergedt$month<-as.numeric(as.character(mergedt$month))
  mergedt$site<-as.character(mergedt$site)
  
  data<-merge(data,mergedt,by=c("site","year","month"))
  data$soilm<-as.numeric(as.character(data$soilm))
  data[,soilm1:=data.table::shift(soilm,1L,type="lag",fill=0)]
  
  data[,deltsoil:=(soilm1*(1-(exp(-1*(petm-wm)/soilawc))))]
  data[,deltsoilwm:=ifelse(deltsoil>0,wm+deltsoil,wm)]
  data[,aet:=ifelse(deltsoilwm<petm,deltsoilwm,petm)]
  data[,cwd:=petm-aet]
  
  ##### for 800 m normal data then we subset to get just the last simulation
  if (type == "normal"){data<-data[year==10,]}
  return(data)}
