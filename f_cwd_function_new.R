# Source code: Redmond, MD. 2022. CWD and AET function (Version V1.0.3). Zenodo. https://doi.org/10.5281/zenodo.6416352
# Original code provided by the Great Basin Landscape Ecology Lab: https://naes.unr.edu/weisberg/old_site/downloads/
# Probably relied upon the appendix from https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2699.2009.02268.x#b45


# THIS FUNCTION CALCULATES CLIMATIC WATER DEFICIT USING EITHER MONTHLY ANNUAL DATA oR MEAN MONTHLY DATA

# DATA INPUT: site (i.e. a unique id for a given location - this vector should be  a character), slope (degrees), latitude (decimal degrees), folded aspect (degrees), ppt (monthly total precipitation, in mm.), tmean (mean monthly temperature, in degrees C),  soilawc (soil available water capacity in the top 150-200 cm of the soil, in mm.; default is 200 mm), month, and year (if using annual data only; default is null))

####note: you input each variable as a vector, but they need to all line up correctly and be the same length, except for soil awc (that way if you don't have soil awc data for each site, you can specify a specific value, such as soilawc = 300)


####note: if you are calculating CWD annually, the script runs such that the amount of water in the soil is at capacity initially (i.e. equal to soil available water capacity). 

# PACKAGES THAT MUST BE INSTALLED BEFORE RUNNING THE SCRIPT: data.table and geosphere

# EXAMPLE SCRIPTS:
###cwd_data<-cwd_function(site=data$site,slope=data$slope,latitude=data$latitude,foldedaspect=data$foldedaspect,ppt=data$ppt,tmean=data$tmean,month=data$month,year=data$year,type="annual")

###cwd_normal_data<-cwd_function(site=data$site,slope=data$slope,latitude=data$latitude,foldedaspect=data$foldedaspect,ppt=data$ppt,tmean=data$tmean,month=data$month,type="normal")

#example script with soil awc specified as one number across all sites:
### cwd_normal_data<-cwd_function(site=data$site,slope=data$slope,latitude=data$latitude,foldedaspect=data$foldedaspect,ppt=data$ppt,tmean=data$tmean,month=data$month,soilawc=300,type="normal")
#example script where I have unique soil awc data for each site:
### cwd_normal_data<-cwd_function(site=data$site,slope=data$slope,latitude=data$latitude,foldedaspect=data$foldedaspect,ppt=data$ppt,tmean=data$tmean,month=data$month,soilawc=data$soilawc,type="normal")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load required packages -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(data.table)
library(geosphere)
library(assertr)
library(parallel)
library(doParallel)
library(lubridate)
library(SPEI)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define function -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pet_spei_function <- function(data){
  # Applies SPEI's implementation (non heatload adjusted) of the thornthwaite equation
  data <- data %>% 
    arrange(year, month) 
  data$pet_spei <- thornthwaite(Tave = data$tmean, lat = data$latitude[1])
  return(data)
}



cwd_function <- function(site, year, month, petm, ppt, tmean, soilawc,normals=FALSE) {
  data<-as.data.table(cbind(as.character(site),year,month,petm,ppt,tmean,soilawc))
  colnames(data)<-c("site","year","month","petm","ppt","tmean","soilawc")
  data$site<-as.character(data$site)
  data$month<-as.numeric(as.character(data$month))
  data$year<-as.numeric(as.character(data$year))
  data$petm<-as.numeric(as.character(data$petm))
  data$ppt<-as.numeric(as.character(data$ppt))
  data$tmean<-as.numeric(as.character(data$tmean))
  data$soilawc<-as.numeric(as.character(data$soilawc))

  ### for 30 year normal data only, we'll create iterations so that way the water storage can carry over from one year to the next
  if (normals==TRUE){
    year<-c(rep(1,length(data$site)),rep(2,length(data$site)),rep(3,length(data$site)),rep(4,length(data$site)),rep(5,length(data$site)),rep(6,length(data$site)),rep(7,length(data$site)),rep(8,length(data$site)),rep(9,length(data$site)),rep(10,length(data$site)))
    data<-rbind(data,data,data,data,data,data,data,data,data,data)
    data$year=year
    }

  data[,fm:=ifelse(tmean<0,0,ifelse(tmean>=6,1,tmean*0.166666666))]
  data[,rainm:=fm*ppt]
  data[,snowm:=(1-fm)*ppt]
  data<-setorder(data,site,year,month)
  sites<-unique(data$site)
  mergedata<-vector()
  mergedata=foreach(s=1:length(sites),.combine="rbind")%dopar%{
    packm<-vector()
    sitedat<-setorder(data,site,year,month)[site==sites[s],]
    for (m in 1:length(sitedat$site)){
      packmsite<-ifelse(m>1,(((1-sitedat$fm[m])^2)*sitedat$ppt[m])+((1-sitedat$fm[m])*packm[m-1]),(((1-sitedat$fm[m])^2)*sitedat$ppt[m]))
      packm<-c(packm,packmsite)
    }
    site<-as.character(sitedat$site)
    year<-as.numeric(as.character(sitedat$year))
    month<-as.numeric(as.character(sitedat$month))
    meltm <- sitedat$fm *(sitedat$snowm + data.table::shift(packm,1L,type="lag",fill=0))
    mergedat<-cbind(site,year,month,packm,meltm)
    # mergedata<-rbind(mergedata,mergedat)
  }
  mergedt<-as.data.table(mergedata)
  mergedt$year<-as.numeric(as.character(mergedt$year))
  mergedt$month<-as.numeric(as.character(mergedt$month))
  mergedt$site<-as.character(mergedt$site)
  data<-merge(data,mergedt,by=c("site","year","month"))
  data$packm<-as.numeric(as.character(data$packm))
  data$meltm<-as.numeric(as.character(data$meltm))
  data[,wm:=rainm+meltm]
  
  mergedata<-vector()
  mergedata=foreach(s=1:length(sites),.combine="rbind")%dopar%{
    soilm<-vector()
    sitedat<-setorder(data,site,year,month)[site==sites[s],]
    for (j in 1:length(sitedat$site)){
      soilmsite<-ifelse(j>1,
                        ## Later periods. Three possible outcomes:
                        ## a) wm - petm <= 0:      soilm is 0 if swc is 0, otherwise soilm is a decay of previous soilm as a function of scale of potential deficit
                        ## b) (wm + soilm[j-1]) - petm < soilawc: soilm adds surplus
                        ## c) (wm + soilm[j-1]) - petm >= soilawc: soilm = soilawc
                        ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,
                               ifelse(sitedat$soilawc==0, 0, soilm[j-1]*((exp(-1*(sitedat$petm[j]-sitedat$wm[j])/sitedat$soilawc)))),
                               ifelse((sitedat$wm[j]-sitedat$petm[j]+soilm[j-1])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]+soilm[j-1]),sitedat$soilawc[j])),

                        ## First period - looks at single period water balance. Three possible outcomes:
                        ## a) wm - petm <= 0:      soilm = 0
                        ## b) wm - petm < soilawc: soilm = wm - petm
                        ## c) wm - petm >= soilawc: soilm = soilawc
                        ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,0,
                               ifelse((sitedat$wm[j]-sitedat$petm[j])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]),sitedat$soilawc[j])))
      
      soilm<-c(soilm,soilmsite)
    }
    site<-as.character(sitedat$site)
    year<-as.numeric(as.character(sitedat$year))
    month<-as.numeric(as.character(sitedat$month))
    mergedat<-cbind(site,year,month,soilm)
    # mergedata<-rbind(mergedata,mergedat)
  }
  mergedt<-as.data.table(mergedata)
  mergedt$year<-as.numeric(as.character(mergedt$year))
  mergedt$month<-as.numeric(as.character(mergedt$month))
  mergedt$site<-as.character(mergedt$site)
  
  data<-merge(data,mergedt,by=c("site","year","month"))
  data$soilm<-as.numeric(as.character(data$soilm))
  data[,soilm1:=data.table::shift(soilm,1L,type="lag",fill=0)]
  data[1,soilm1:=soilm] # No change in first period
  
  data[,deltsoil:=soilm1-soilm] # change to soilm1 - soilm
  data[,aet:=ifelse(wm>petm, petm, wm+deltsoil)]
  data[,cwd:=petm-aet]
  data<-setorder(data,site,year,month)
  
  if(normals==TRUE)data=data[which(data$year==10),]
  
  return(data)
}

pet_function <- function(site,year,month,slope,latitude,aspect,tmean) {
  data<-as.data.table(cbind(as.character(site),year,month,slope,latitude,aspect,tmean))
  colnames(data)<-c("site","year", "month", "slope","latitude","aspect","tmean")
  
  data$site<-as.character(data$site)
  data$slope<-as.numeric(as.character(data$slope))
  data$latitude<-as.numeric(as.character(data$latitude))
  data$aspect<-as.numeric(as.character(data$aspect))
  data$tmean<-as.numeric(as.character(data$tmean))
  data$month<-as.numeric(as.character(data$month))
  data$year<-as.numeric(as.character(year))
  data[,yearmonth:= as.yearmon(paste(year, month), "%Y %m")]
  data[,days:=days_in_month(yearmonth)]
  data<-setorder(data,site,year,month)
  
  data$foldedaspect <- ifelse(data$latitude>0,180 - ( data$aspect - 225), 180 - (data$aspect - 315)) # convert aspect (from radians to degrees) into folded aspect, if latitude <0 (southern hemisphere) maximum heat exposure in NW direction (315 degrees) else in SW direction (225 degrees).
  
  
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
  data<-setorder(data,site,year,month)
  
  data[,ea:=exp(((17.3*tmean)/(tmean+237.3)))*0.611] 
  # convert slope, folded aspect, and latitude to radians
  data[,sloprad:=slope*pi/180]
  data[,afrad:=foldedaspect*pi/180]
  data[,latrad:=latitude*pi/180]
  # calculate heat load
  data[,heatload:=0.339+0.808*(cos(latrad)*cos(sloprad))-0.196*(sin(latrad)*sin(sloprad))-0.482*(cos(afrad)*sin(sloprad))]
  
  # calculate pet
  data[,petm:=ifelse(tmean<0,0,((((ea)/(tmean+273.2))*day*days*29.8)*heatload))]
  
  return(data)
}

