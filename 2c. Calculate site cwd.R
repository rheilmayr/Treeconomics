#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 7/10/20
# Purpose: Convert site topography, soil and weather to CWD/AET/PET
#
# Input files:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load packages --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("cwd_function.R")
library(dplyr)
library(naniar)
library(data.table)
library(tidyverse)
library(geosphere)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

# 1. Pre-processed climate and soil data
data=fread(paste0(wdir,"out\\climate\\sitedataforcwd_210620.csv"))
miss_var_summary(data)


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
# Calculate CWD and save data --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl=makeCluster(4)
clusterExport(cl,c("data","setorder"))
registerDoParallel(cl)

cwd_data<-cwd_function(site=data$site_id,slope=data$slope,latitude=data$latitude,
                       foldedaspect=data$aspect,ppt=data$pre_corrected,
                       tmean=data$tmean,month=data$month,year=data$year,
                       soilawc=data$swc,type="annual")
# fwrite(cwd_data,file=paste0(wdir,"out/climate/cwd_data_200620.csv"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Characterize missing data ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_data <- data %>% 
  select(site_id, slope, latitude, longitude, aspect, pre_corrected, tmean, month, year, swc)

input_data %>% gg_miss_upset()
any_missing <- input_data %>%
  filter_all(any_vars(is.na(.)))
share_missing <- dim(any_missing)[1] / dim(data)[1]
share_missing
miss_var_summary(data)


out_data <- cwd_data %>% 
  select(site, year, month, wm, petm, soilm, deltsoil, aet, cwd)

out_data %>% gg_miss_upset()
any_missing <- input_data %>%
  filter_all(any_vars(is.na(.)))
share_missing <- dim(any_missing)[1] / dim(data)[1]
share_missing
miss_var_summary(data)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Write out file ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_data <- cwd_data %>%
  select(site, year, month, aet, cwd)
fwrite(cwd_data,file=paste0(wdir,"out/climate/essentialcwd_data.csv"))
