
library(tidyverse)
source("f_cwd_function.R")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Pre-processed climate and soil data
data=fread(paste0(wdir,"1_input_processed/climate/sitedataforcwd.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set NAs for precip and temp
data$pre = na_if(data$pre,-999)
data$tmn = na_if(data$tmn,-999)
data$tmx = na_if(data$tmx,-999)

# Add corrections from World Clim to CWD to get downscaled variables
data$pre_corrected=data$pre+data$pre_correction
data$tmx_corrected=data$tmx+data$tmax_correction
data$tmn_corrected=data$tmn+data$tmin_correction


# Unit conversions
data$swc=data$swc/10 #convert swc from mm to cm
data$tmean=(data$tmn_corrected+data$tmx_corrected)/2 
data$slope <- data$slope * 57.2958 # convert slope from radians to degrees
data$aspect <- data$aspect * 57.2958 # convert aspect from radians to degrees

data <- data %>% 
  select(site_id, slope, latitude, longitude, aspect, pre_corrected, tmean, month, year, swc)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Focus on example site --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_data <- data %>% 
  filter(site_id == "CA524") # PIPO site in CA

cwd_data <- cwd_function(site=test_data$site_id,
                       slope=test_data$slope,
                       latitude=test_data$latitude,
                       foldedaspect=test_data$aspect,
                       ppt=test_data$pre_corrected,
                       tmean=test_data$tmean,
                       month=test_data$month,
                       year=test_data$year,
                       soilawc=test_data$swc,
                       type="annual")

norm_data <- test_data %>% 
  group_by(site_id, month) %>%
  summarise(slope = mean(slope),
          latitude = mean(latitude),
          aspect = mean(aspect),
          pre_corrected = mean(pre_corrected),
          tmean = mean(tmean),
          swc = mean(swc)
  )

cwd_data <- cwd_function(site=norm_data$site_id,
                         slope=norm_data$slope,
                         latitude=norm_data$latitude,
                         foldedaspect=norm_data$aspect,
                         ppt=norm_data$pre_corrected,
                         tmean=norm_data$tmean,
                         month=norm_data$month,
                         soilawc=norm_data$swc,
                         type="normal")
