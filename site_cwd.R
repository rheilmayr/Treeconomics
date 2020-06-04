

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load packages --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("cwd_function.R")
library(dplyr)
library(naniar)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Pre-processed climate and soil data
data=read.csv(paste0(wdir,"CRU\\181116-climate_soil_data_with_corrections.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set NAs for precip and temp
data <- data %>%
  mutate(pre = na_if(pre, pre==-999),
         tmn = na_if(tmn, tmn==-999),
         tmx = na_if(tmx, tmx==-999))

# Add corrections from World Clim to CWD to get downscaled variables
data$pre_corrected=data$pre+data$pre_correction
data$tmx_corrected=data$tmx+data$tmax_correction
data$tmn_corrected=data$tmn+data$tmin_correction


# Unit conversions
# tmax, tmin and precip are in units of 1/10th of a degree / mm. swc is in units of mm and needs to be in units of cm
cols=c("pre_corrected","tmn_corrected","tmx_corrected","swc") 
for(i in 1:length(cols)){
  col=which(colnames(data)==cols[i])
  data[,col]=data[,col]/10
}
data$tmean=(data$tmn_corrected+data$tmx_corrected)/2 # IS IT VALID TO TAKE MIDPOINT BETWEEN MAX AND MIN TO CALCULATE TMEAN? 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Characterize missing data ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data %>% gg_miss_upset()
any_missing <- data %>%
  filter_all(any_vars(is.na(.)))
share_missing <- dim(any_missing)[1] / dim(data)[1]
share_missing
miss_var_summary(data)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate CWD and save data --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_data<-cwd_function(site=data$site_id,slope=data$slope,latitude=data$latitude,
                       foldedaspect=data$aspect,ppt=data$pre_corrected,
                       tmean=data$tmean,month=data$month,year=data$year,
                       soilawc=data$swc,type="annual")
write.csv(data,file=paste0(wdir,"/cwd_data.csv"))
write.csv(data[,c(1:9,23,27:29)],file=paste0(wdir,"/essentialcwd_data.csv"))