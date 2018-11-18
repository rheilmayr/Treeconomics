
# THIS FUNCTION CALCULATES CLIMATIC WATER DEFICIT USING EITHER MONTHLY ANNUAL DATA oR MEAN MONTHLY DATA

# DATA INPUT: site (i.e. a unique id for a given location), slope (degrees), latitude (decimal degrees), folded aspect (degrees), ppt, tmean,  soilawc (in cm.; default is 200 cm), month, and year (if using annual data only; default is null))

####note: you input each variable as a vector, but they need to all line up correctly and be the same length, except for soil awc (that way if you don't have soil awc data for each site, you can specify a specific value, such as soilawc = 300)

# PACKAGES THAT MUST BE INSTALLED BEFORE RUNNING THE SCRIPT: data.table and geosphere

# EXAMPLE SCRIPTS:
franspath="C:/Users/fmoor/Box/Davis Stuff/Treeconomics"

###cwd_normal_data<-cwd_function(site=data$site,slope=data$slope,latitude=data$latitude,foldedaspect=data$foldedaspect,ppt=data$ppt,tmean=data$tmean,month=data$month,type="normal")
library(dplyr)
#example script with soil awc specified as one number across all sites:
### cwd_normal_data<-cwd_function(site=data$site,slope=data$slope,latitude=data$latitude,foldedaspect=data$foldedaspect,ppt=data$ppt,tmean=data$tmean,month=data$month,soilawc=300,type="normal")
#example script where I have unique soil awc data for each site:
data=read.csv(paste0(franspath,"/Data/181116-climate_soil_data_with_corrections.csv"))
#add corrections from World Clim to CWD to get downscaled variables
data$pre_corrected=data$pre+data$pre_correction
data$tmx_corrected=data$tmx+data$tmax_correction
data$tmn_corrected=data$tmn+data$tmin_correction

cols=c("pre_corrected","tmn_corrected","tmx_corrected","swc") #tmax, tmin and precip are in units of 1/10th of a degree / mm. swc is in units of mm and needs to be in units of cm
for(i in 1:length(cols)){
  col=which(colnames(data)==cols[i])
  data[,col]=data[,col]/10
}
data$tmean=(data$tmn_corrected+data$tmx_corrected)/2
cwd_data<-cwd_function(site=data$site_id,slope=data$slope,latitude=data$latitude,foldedaspect=data$aspect,ppt=data$pre_corrected,tmean=data$tmean,month=data$month,year=data$year,soilawc=data$swc,type="annual")

#cwd_normal_data<-cwd_function(site=data$site_id,slope=data$slope,latitude=data$latitude,foldedaspect=data$aspect,ppt=data$pre,tmean=data$tmean,month=data$month,soilawc=data$swc,type="normal")

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
for (i in 1:length(sites)){
  packm<-vector()
  sitedat<-setorder(data,site,year,month)[site==sites[i],]
  for (j in 1:length(sitedat$site)){
    packmsite<-ifelse(j>1,(((1-sitedat$fm[j])^2)*sitedat$ppt[j])+((1-sitedat$fm[j])*packm[j-1]),(((1-sitedat$fm[j])^2)*sitedat$ppt[j]))
    packm<-c(packm,packmsite)
    site<-as.character(sitedat$site)
    year<-as.numeric(as.character(sitedat$year))
    month<-as.numeric(as.character(sitedat$month))
  }
  mergedat<-cbind(site,year,month,packm)
  mergedata<-rbind(mergedata,mergedat)
  print(paste0(i,"_1"))
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

mergedata<-vector()
for (i in 1:length(sites)){
  soilm<-vector()
  sitedat<-setorder(data,site,year,month)[site==sites[i],]
  for (j in 1:length(sitedat$site)){
    soilmsite<-ifelse(j>1,ifelse((sitedat$wm[j]-sitedat$petm[j]+soilm[j-1])<=0,0,ifelse((sitedat$wm[j]-sitedat$petm[j]+soilm[j-1])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]+soilm[j-1]),sitedat$soilawc[j])),ifelse((sitedat$wm[j]-sitedat$petm[j])<=0,0,ifelse((sitedat$wm[j]-sitedat$petm[j])<sitedat$soilawc[j],(sitedat$wm[j]-sitedat$petm[j]),sitedat$soilawc[j])))
    soilm<-c(soilm,soilmsite)
    site<-as.character(sitedat$site)
    year<-as.numeric(as.character(sitedat$year))
    month<-as.numeric(as.character(sitedat$month))
  }
  mergedat<-cbind(site,year,month,soilm)
  mergedata<-rbind(mergedata,mergedat)
  print(paste0(i,"_2"))
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

write.csv(data,file=paste0(franspath,"/Data/cwd_data.csv"))
write.csv(data[,c(1:9,23,27:29)],file=paste0(franspath,"/Data/essentialcwd_data.csv"))
