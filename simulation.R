library(dplyr)
library(ggplot2)
library(broom)
library(purrr)
library(patchwork)

nsites=75
nyears=100
tempnoise=1.5 #standard deviation of monthly temperature around climatology, in degrees
seasonalamplitude=10 #difference between coldest and hottest months of year in degrees (future - make this vary based on average temp?)

annualtemps=runif(nsites,min=5,max=25)

tempdata=matrix(nrow=nsites,ncol=12*nyears)

#northern hemisphere - cold in beginning of year
seasonalcycle=cos(rep(seq(from=pi,to=3*pi,length.out=12),nyears))*seasonalamplitude

for(i in 1:nsites){
  temprecord=seasonalcycle+annualtemps[i]+rnorm(length(seasonalcycle),mean=0,sd=tempnoise)
  tempdata[i,]=temprecord
}

pet=function(temp,day=days,dl=DL){
  ea=0.611*exp(17.3*temp/(temp+237.2))
  return(29.8*day*dl*ea/(temp+273.2))
}

pet_wrong=function(temp,day=days,dl=DL){
  ea=0.611*exp(17.3*temp/(temp+237.2))*temp/10
  return(29.8*day*dl*ea/(temp+273))
}

plot(0:30,pet(0:30,day=rep(30,length(0:30)),dl=rep(12,length(0:30))),type="l",col="black",ylim=c(0,300))
lines(0:30,pet_wrong(0:30,day=rep(30,length(0:30)),dl=rep(12,length(0:30))),type="l",col="red")
legend("topleft",col=c("black","red"),legend=c("Correct","Error"),lwd=1)

days=matrix(nrow=nsites,data=rep(rep(c(31,28,31,30,31,30,31,31,30,31,30,31),nyears),each=nsites))
DL=matrix(nrow=nsites,data=rep(rep(c(5,6,7,9,10,11,13,14,12,10,8,6),nyears),each=nsites)) #could correlate with annual temps as well

petrecord_right=as.data.frame(pet(tempdata))
petrecord_wrong=as.data.frame(pet_wrong(tempdata))
colnames(petrecord_right)=rep(1:nyears,each=12);colnames(petrecord_wrong)=colnames(petrecord_right)
petrecord_right$id=1:nsites;petrecord_wrong$id=1:nsites
petrecord_right=pivot_longer(petrecord_right,cols=-id,names_to="year",values_to = "pet")
petrecord_wrong=pivot_longer(petrecord_wrong,cols=-id,names_to="year",values_to = "pet")
petrecord_right$month=rep(1:12,nsites*nyears);petrecord_wrong$month=rep(1:12,nsites*nyears)
petrecord_right$type="right";petrecord_wrong$type="wrong"
petall=rbind(petrecord_right,petrecord_wrong)
annualtemps=data.frame(id=1:nsites,meantemp=annualtemps)
petall=merge(petall,annualtemps,by="id")

#error across all sites
a=ggplot(petall%>%pivot_wider(names_from = type,values_from = pet),aes(x=right,y=wrong,col=meantemp))+geom_point()
a=a+geom_abline(slope=1,intercept = 0,col="darkred")+theme_classic()+labs(x="Correct PET",y="Wrong PET")+labs(title="All Sites All Months")
a

#time series for a few years for one site
b=ggplot(petall%>%filter(id==9&year%in%c(1:3)),aes(x=month,y=pet,col=type,group=interaction(year,type)))+geom_line()
b=b+theme_classic()+labs(title="One Site for Multiple Years")
a+b

##---------add rainfall and estimate annual cwd -----------------

#create cwd based on random average precip across sites and log norm precip across years within sites
meanprecip=runif(nsites,min=5,max=100)
precipdata=matrix(nrow=nsites,ncol=12*nyears)

for(i in 1:nsites){
  precipdata[i,]=meanprecip[i]*rlnorm(nyears*12,mean=0,sd=0.15)
}
precipdata=as.data.frame(precipdata)
colnames(precipdata)=rep(1:nyears,each=12)
precipdata$id=1:nsites
precipdata=pivot_longer(precipdata,cols=-id,names_to="year",values_to="ppt")
precipdata$month=rep(1:12,nsites*nyears)

cwddata=merge(petall,precipdata,by=c("id","year","month"))

#aggregate up to year
cwdall_year=cwddata%>%group_by(type,id,year)%>%dplyr::summarise(totalpet=sum(pet),totalppt=sum(ppt),totalcwd=max(0,sum(pet-ppt)),annualtemps=mean(meantemp))
meancwd=cwdall_year%>%group_by(id,type)%>%dplyr::summarise(sitemeancwd=mean(totalcwd))
cwdall_year=merge(cwdall_year,meancwd)

#summary plots of cwd
#error across all sites
plotdat=cwdall_year%>%select(id,year,totalcwd,type)%>%pivot_wider(names_from = type,values_from = totalcwd)
plotdat=merge(plotdat,cwdall_year%>%filter(type=="right")%>%select(id,year,sitemeancwd))
a=ggplot(plotdat,aes(x=right,y=wrong,col=sitemeancwd))+geom_point()
a=a+geom_abline(slope=1,intercept = 0,col="darkred")+theme_classic()+labs(x="CWD with Correct PET",y="CWD with Wrong PET",col="True Site Mean\nAnnual CWD")
a


#----------------relationship with outcome variable------------------------

#---------true relationship - no effect of average climate -----null effect
siteresponse=-0.05 #common effect of cwd on response variable across all sites (no spoiled tree or dry range edge sensitivity)
responseerror=0.5 #standard deviation of response error
peteffect=0.02 #benefit of higher pet - common across all sites

truedriver=cwdall_year%>%filter(type=="right")
y=siteresponse*truedriver$totalcwd+peteffect*truedriver$totalpet+rnorm(nsites*nyears,mean=0,sd=responseerror)
truedriver$y=y

data=merge(cwdall_year,truedriver[,c(2,3,9)],by=c("id","year"))

#fit regression model for each site for wrong and right pet
mods=data%>%
  nest(.by=c("id",'type','sitemeancwd',"annualtemps"))%>%
  dplyr::mutate(
    models=lapply(data,function(df) lm(y~totalpet+totalcwd,data=df)),
    tidied=map(models,tidy)
    )%>%
  unnest(tidied)%>%
  filter(term!="(Intercept)")

#compare regression coefficients visually
d=ggplot(mods%>%filter(term=="totalcwd")%>%select(type,annualtemps,estimate,id)%>%pivot_wider(names_from = type,values_from = estimate),aes(x=right,y=wrong,col=annualtemps))+geom_point()
d=d+theme_classic()+labs(title="True Effect: No Spoiled Tree or Dry Range Edge Effect")
d=d+geom_vline(xintercept = siteresponse,lty=2)+guides(col="none")+geom_hline(yintercept = siteresponse,lty=2)

##-----introduce true spoiled tree effect - correlation of slope with temperature

spoiledeffect=0.0006 #change in response per degree average cwd - smaller at hotter sites

truedriver=cwdall_year%>%filter(type=="right")
siteeffect=ifelse(siteresponse+spoiledeffect*(truedriver$sitemeancwd-mean(truedriver$sitemeancwd))>0,0,siteresponse+spoiledeffect*(truedriver$sitemeancwd-mean(truedriver$sitemeancwd)))
y=siteeffect*truedriver$totalcwd+peteffect*truedriver$totalpet+rnorm(nsites*nyears,mean=0,sd=responseerror)
truedriver$y=y

data_spoiled=merge(cwdall_year,truedriver[,c(2,3,9)],by=c("id","year"))

#fit regression model for each site for wrong and right pet
mods_spoiled=data_spoiled%>%
  nest(.by=c("id",'type','sitemeancwd',"annualtemps"))%>%
  dplyr::mutate(
    models=lapply(data,function(df) lm(y~totalpet+totalcwd,data=df)),
    tidied=map(models,tidy)
  )%>%
  unnest(tidied)%>%
  filter(term!="(Intercept)")

#compare regression coefficients visually
e=ggplot(mods_spoiled%>%filter(term=="totalcwd")%>%select(type,annualtemps,estimate,id)%>%pivot_wider(names_from = type,values_from = estimate),aes(x=right,y=wrong,col=annualtemps))+geom_point()
e=e+theme_classic()+labs(title="True Effect: Spoiled Tree")+geom_abline(slope=1,intercept=0)
e=e+guides(col="none")

##-----dry range sensitive - correlation of slope with temperature in opposite direction

truedriver=cwdall_year%>%filter(type=="right")
siteeffect=ifelse(siteresponse-spoiledeffect*(truedriver$sitemeancwd-mean(truedriver$sitemeancwd))>0,0,siteresponse-spoiledeffect*(truedriver$sitemeancwd-mean(truedriver$sitemeancwd)))
y=siteeffect*truedriver$totalcwd+peteffect*truedriver$totalpet+rnorm(nsites*nyears,mean=0,sd=responseerror)
truedriver$y=y

data_dryrange=merge(cwdall_year,truedriver[,c(2,3,9)],by=c("id","year"))

#fit regression model for each site for wrong and right pet
mods_dryrange=data_dryrange%>%
  nest(.by=c("id",'type','sitemeancwd',"annualtemps"))%>%
  dplyr::mutate(
    models=lapply(data,function(df) lm(y~totalpet+totalcwd,data=df)),
    tidied=map(models,tidy)
  )%>%
  unnest(tidied)%>%
  filter(term!="(Intercept)")
#compare regression coefficients visually
f=ggplot(mods_dryrange%>%filter(term=="totalcwd")%>%select(type,annualtemps,estimate,id)%>%pivot_wider(names_from = type,values_from = estimate),aes(x=right,y=wrong,col=annualtemps))+geom_point()
f=f+theme_classic()+labs(title="True Effect: Dry Range Edge")+geom_abline(slope=1,intercept=0)

d+e+f

###----- second stage regression under three true models ------

secondstage_null=mods%>%
  filter(term=="totalcwd")%>%
  nest(.by=type)%>%
  dplyr::mutate(
    models=lapply(data,function(df) lm(estimate~sitemeancwd,weights = 1/std.error,data=df)),
    tidied=map(models,tidy)
  )%>%
  unnest(tidied)%>%
  filter(term=="sitemeancwd")

secondstage_spoiled=mods_spoiled%>%
  filter(term=="totalcwd")%>%
  nest(.by=type)%>%
  dplyr::mutate(
    models=lapply(data,function(df) lm(estimate~sitemeancwd,weights = 1/std.error,data=df)),
    tidied=map(models,tidy)
  )%>%
  unnest(tidied)%>%
  filter(term=="sitemeancwd")

secondstage_dryrange=mods_dryrange%>%
  filter(term=="totalcwd")%>%
  nest(.by=type)%>%
  dplyr::mutate(
    models=lapply(data,function(df) lm(estimate~sitemeancwd,weights = 1/std.error,data=df)),
    tidied=map(models,tidy)
  )%>%
  unnest(tidied)%>%
  filter(term=="sitemeancwd")


