#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 7/10/20
# Purpose: Convert site topography, soil and weather to CWD/AET/PET
#
# Input files:
#   sitedataforcwd.csv: Compiled site-level climate and soil data. Created by "2b. Pull site weather.R"
# 
# Output files:
#   essentialcwd_data.csv: Compiled annual CWD/PET data for each site.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load packages --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("f_cwd_function.R")
library(dplyr)
library(naniar)
library(data.table)
library(tidyverse)
library(geosphere)
library(tictoc)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Pre-processed climate and soil data
data=fread(paste0(wdir,"1_input_processed/climate/sitedataforcwd.csv"))
miss_var_summary(data)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set NAs for precip and temp
data$pre = na_if(data$pre,-999)
data$tmn = na_if(data$tmn,-999)
data$tmx = na_if(data$tmx,-999)
data$tmp = na_if(data$tmp,-999)
data$pet = na_if(data$pet,-999)

# # Add corrections from World Clim to CWD to get downscaled variables
# data$pre_corrected=data$pre+data$pre_correction
# data$pre_corrected=ifelse(data$pre_corrected < 0, 0, data$pre_corrected)
# data$tmx_corrected=data$tmx+data$tmax_correction
# data$tmn_corrected=data$tmn+data$tmin_correction

# Unit conversions
data$swc=data$swc/10 #convert swc from mm to cm
data$tmean=data$tmp
# data$slope <- data$slope * 57.2958 # convert slope from radians to degrees
# data$aspect <- data$aspect * 57.2958 # convert aspect from radians to degrees
data$petd <- data$pet
data$pre_corrected <- data$pre

data <- data[,c("site_id", "year", "month", "latitude", "longitude", 
                "pre_corrected", "tmean", "petd", "swc")]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate CWD and save data --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl=makeCluster(4)
clusterExport(cl,c("data","setorder"))
registerDoParallel(cl)

tic()
cwd_data <- cwd_function(site=data$site_id,
                         year=data$year,
                         month=data$month,
                         latitude=data$latitude,
                         ppt=data$pre_corrected,
                         tmean=data$tmean,
                         petd=data$petd,
                         soilawc=data$swc,
                         type="annual")
toc()
# fwrite(cwd_data,file=paste0(wdir,"out/climate/cwd_data_200620.csv"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Characterize missing data ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_data <- data %>% 
  select(site_id, latitude, longitude, pre_corrected, tmean, petd, month, year, swc)

input_data %>% gg_miss_upset()
# any_missing <- input_data %>%
#   filter_all(any_vars(is.na(.)))
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
  select(site, year, month, tmean, ppt, aet, cwd, pet = petm, cwb)

fwrite(cwd_data,file=paste0(wdir,"1_input_processed/climate/essentialcwd_data.csv"))

