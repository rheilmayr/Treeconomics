#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 1/26/24
# Purpose: Calculate CWD from known PET, precip and temperature variables
#
# 
# CWD original code provided by the Great Basin Landscape Ecology Lab: https://naes.unr.edu/weisberg/old_site/downloads/
# Original script says it adapts the method described by Dobrowski et al., 2012: https://doi.org/10.1111/gcb.12026
# But probably relied upon the appendix from https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2699.2009.02268.x#b45
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# DATA INPUT: site (i.e. a unique id for a given location), latitude (decimal degrees), 
# ppt, tmean,  soilawc (in cm.; default is 200 cm), month, and year (if using 
# annual data only; default is null))

#### Note: you input each variable as a vector, but they need to all line up correctly 
# and be the same length, except for soil awc (that way if you don't have soil awc data 
# for each site, you can specify a specific value, such as soilawc = 300)

# PACKAGES THAT MUST BE INSTALLED BEFORE RUNNING THE SCRIPT: data.table and geosphere

# Code adapated from Redmond, MD. 2022. CWD and AET function (Version V1.0.3). Zenodo. https://doi.org/10.5281/zenodo.6416352
# Relied upon the appendix from https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2699.2009.02268.x#b45


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load required packages -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(data.table)
library(geosphere)
library(assertr)
library(parallel)
library(doParallel)
library(lubridate)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define function -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_function <- function(site, year=NULL, month, latitude, ppt, tmean, pet, soilawc=200, type=c('normal','annual')) {

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
  data$foldedaspect <- abs(180 - abs(foldedaspect-225))  ## Were we calculating this incorrectly?
  data$ppt<-as.numeric(as.character(data$ppt))
  data$tmean<-as.numeric(as.character(data$tmean))
  data$month<-as.numeric(as.character(data$month))
  data$date <- as.yearmon(paste(data$year, data$month), "%Y %m")
  
  
  as.Date(paste(data$year, data$month, sep = "-"), "%Y-%M")

  ### for 30 year normal data only, we'll create iterations so that way the water storage can carry over from one year to the next
  if (type == "normal"){year<-c(rep(1,length(data$site)),rep(2,length(data$site)),rep(3,length(data$site)),rep(4,length(data$site)),rep(5,length(data$site)),rep(6,length(data$site)),rep(7,length(data$site)),rep(8,length(data$site)),rep(9,length(data$site)),rep(10,length(data$site)))
  data<-rbind(data,data,data,data,data,data,data,data,data,data)}
  data$year<-as.numeric(as.character(year))
  # calculate daylength

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

  ## Melt factor
  data[,fm:=ifelse(tmean<0,0,ifelse(tmean>=6,1,tmean*0.166666666))]
  
  ## Precipitation falling as rain
  data[,rainm:=fm*ppt]
  
  ## Precipitation falling as snow
  data[,snowm:=(1-fm)*ppt]
  
  ## Track snowpack
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
  
  ## Monthly water input to system
  data[,wm:=rainm+meltm]
  
  ## Days in month
  data$days = days_in_month(data$date)

  ## Soil moisture
  mergedata<-vector()
  for (i in 1:length(sites)){
    soilm<-vector()
    sitedat<-setorder(data,site,year,month)[site==sites[i],]
    for (j in 1:length(sitedat$site)){
      soilmsite<-ifelse(j>1,ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,soilm[j-1]*((exp(-1*(sitedat$petm[j]-sitedat$wm[j])/sitedat$soilawc))),ifelse((sitedat$wm[j]-sitedat$petm[j]+soilm[j-1])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]+soilm[j-1]),sitedat$soilawc[j])),ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,0,ifelse((sitedat$wm[j]-sitedat$petm[j])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]),sitedat$soilawc[j])))
      soilm<-c(soilm,soilmsite)
      site<-as.character(sitedat$site)
      year<-as.numeric(as.character(sitedat$year))
      month<-as.numeric(as.character(sitedat$month))
    }
    mergedat<-cbind(site,year,month,soilm)
    mergedata<-rbind(mergedata,mergedat)
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
  
  ##### for 800 m normal data then we subset to get just the last simulation
  if (type == "normal"){data<-data[year==10,]}
  
  return(data)}

