#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 1/26/24
# Purpose: Calculate CWD from known PET, precip and temperature variables
#
# Code adapated from Redmond, MD. 2022. CWD and AET function (Version V1.0.3). Zenodo. https://doi.org/10.5281/zenodo.6416352
# Relied upon the appendix from https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2699.2009.02268.x#b45
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# DATA INPUT: site (i.e. a unique id for a given location), latitude (decimal degrees), 
# ppt (mmh2o/m), tmean (C), petd(mmh2o/day), soilawc (in cm.; default is 200 cm), month, and year (if using 
# annual data only; default is null))

#### Note: you input each variable as a vector, but they need to all line up correctly 
# and be the same length, except for soil awc (that way if you don't have soil awc data 
# for each site, you can specify a specific value, such as soilawc = 300)

# PACKAGES THAT MUST BE INSTALLED BEFORE RUNNING THE SCRIPT: data.table and geosphere


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load required packages -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(data.table)
library(geosphere)
library(assertr)
library(parallel)
library(doParallel)
library(lubridate)
library(zoo)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define function -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_function <- function(site, month, latitude, ppt, tmean, petd, soilawc=200, year=NULL, type=c('normal','annual')) {

  # Packages
  require(data.table)
  require(geosphere)

  #R SCRIPT
  if (type =="annual"){
    data<-as.data.table(cbind(as.character(site),year,month,latitude,ppt,tmean,petd))[,soilawc:=soilawc]
    colnames(data)<-c("site","year","month","latitude","ppt","tmean","petd","soilawc")}

  if (type =="normal"){
    data<-as.data.table(cbind(as.character(site),month,latitude,ppt,tmean,petd))[,soilawc:=soilawc]
    colnames(data)<-c("site","month","latitude","ppt","tmean","petd","soilawc")
  }
  data$site<-as.character(data$site)
  data$latitude<-as.numeric(as.character(data$latitude))
  data$ppt<-as.numeric(as.character(data$ppt))
  data$petd<-as.numeric(as.character(data$petd))
  data$tmean<-as.numeric(as.character(data$tmean))
  data$month<-as.numeric(as.character(data$month))
  data$date <- as.yearmon(paste(data$year, data$month), "%Y %m")
  
  ### for 30 year normal data only, we'll create iterations so that way the water storage can carry over from one year to the next
  if (type == "normal"){year<-c(rep(1,length(data$site)),rep(2,length(data$site)),rep(3,length(data$site)),rep(4,length(data$site)),rep(5,length(data$site)),rep(6,length(data$site)),rep(7,length(data$site)),rep(8,length(data$site)),rep(9,length(data$site)),rep(10,length(data$site)))
  data<-rbind(data,data,data,data,data,data,data,data,data,data)}
  data$year<-as.numeric(as.character(year))
  
  ## Calculate monthly average daylength
  daylength<-vector()
  datasite<-data[,.(latitude=mean(latitude)),by=.(site)]
  for (i in 1:length(datasite$latitude)){
    dl<-daylength(datasite$latitude[i],1:365)
    day_l<- tapply(dl, rep(1:12, c(31,28,31,30,31,30,31,31,30,31,30,31)), mean)
    site <- as.vector(rep(datasite$site[i],length(day_l)))
    month<-as.vector(c(1,2,3,4,5,6,7,8,9,10,11,12))
    join<-cbind(site,month,day_l)
    daylength<-rbind(daylength,join)
  }
  daylength<-as.data.frame(daylength)
  daylength$site<-as.character(daylength$site)
  daylength$month<-as.numeric(as.character(daylength$month))
  daylength$day_l<-as.numeric(as.character(daylength$day_l))
  data<-merge(data,daylength,by=c("site","month"))

  ## Days in month
  data$days = days_in_month(data$date)
  
  ## Convert daily to monthly PET
  data$petm = data$days * data$petd
  
  ## Melt factor
  data[,fm:=ifelse(tmean<0,0,ifelse(tmean>=6,1,tmean*0.166666666))]
  
  ## Share of precipitation falling as rain
  data[,rainm:=fm*ppt]
  
  ## Share of precipitation falling as snow
  data[,snowm:=(1-fm)*ppt]
  
  ## Track snowpack
  data<-setorder(data,site,year,month)
  sites<-unique(data$site)
  # mergedata<-vector()

  mergedata=foreach(i=1:length(sites),.combine="rbind")%dopar%{
    sitedat<-setorder(data,site,year,month)[site==sites[i],]
    
    # Snowpack in first month
    packm <- ((1-sitedat$fm[1])^2) * sitedat$ppt[1]
    
    # Snowpack in later months
    for (j in 2:length(sitedat$site)){
      packmsite<-(((1-sitedat$fm[j])^2)*sitedat$ppt[j])+((1-sitedat$fm[j])*packm[j-1])
      packm<-c(packm,packmsite)
    }
    site<-as.character(sitedat$site)
    year<-as.numeric(as.character(sitedat$year))
    month<-as.numeric(as.character(sitedat$month))
    mergedat<-cbind(site,year,month,packm)
  }
  mergedt<-as.data.table(mergedata)
  mergedt$year<-as.numeric(as.character(mergedt$year))
  mergedt$month<-as.numeric(as.character(mergedt$month))
  mergedt$site<-as.character(mergedt$site)
  data<-merge(data,mergedt,by=c("site","year","month"))
  data$packm<-as.numeric(as.character(data$packm))
  data[,meltm:=fm*(snowm+data.table::shift(packm,1L,type="lag",fill=0))]
  
  ## Monthly water input to system
  data[,wm:=rainm+meltm]

  ## Soil moisture
  mergedata=foreach(i=1:length(sites),.combine="rbind")%dopar%{
    soilm<-vector()
    sitedat<-setorder(data,site,year,month)[site==sites[i],]
    
    # Soil moisture in year 1
    
    # Soil moisture in later years
    for (j in 1:length(sitedat$site)){
      soilmsite<-ifelse(j>1,ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,soilm[j-1]*((exp(-1*(sitedat$petm[j]-sitedat$wm[j])/sitedat$soilawc))),ifelse((sitedat$wm[j]-sitedat$petm[j]+soilm[j-1])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]+soilm[j-1]),sitedat$soilawc[j])),ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,0,ifelse((sitedat$wm[j]-sitedat$petm[j])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]),sitedat$soilawc[j])))
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
  
  data[,deltsoil:=soilm1-soilm]#change to soilm1 - soilm
  data[,deltsoilwm:=wm+deltsoil]
  data[,aet:=ifelse(wm>petm,petm,wm+deltsoil)]
  data[,cwd:=petm-aet]
  data[,cwb:=ppt-petm]
  
  ##### for 800 m normal data then we subset to get just the last simulation
  if (type == "normal"){data<-data[year==10,]}
  
  return(data)}

