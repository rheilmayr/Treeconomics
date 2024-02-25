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
source("f_cwd_function_new.R")
library(dplyr)
library(naniar)
library(data.table)
library(tidyverse)
library(geosphere)
library(zoo)
library(tictoc)

library(tidyverse)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic()
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
data$pet_cru = na_if(data$pet,-999)

miss_var_summary(data)


# Add corrections from World Clim to CWD to get downscaled variables
data$pre_corrected=data$pre*data$pre_correction
# data$pre_corrected=ifelse(data$pre_corrected < 0, 0, data$pre_corrected)
data$tmx_corrected=data$tmx+data$tmax_correction
data$tmn_corrected=data$tmn+data$tmin_correction

# Unit conversions
data$swc=data$swc/10 #convert swc from mm to cm
data$tmean=(data$tmn_corrected+data$tmx_corrected)/2 
data$slope <- data$slope * 57.2958 # convert slope from radians to degrees
data$aspect <- data$aspect * 57.2958 # convert aspect from radians to degrees

data <- data %>% 
  select(site = site_id, slope, latitude, longitude, aspect, pre_corrected, tmean, month, year, swc, pet_cru)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate PET --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold_data <- data
data <- hold_data
# data <- data %>% rename(site = site_id)
# data <- data[1:100000,]
  
# data <- data %>% as.data.table()

## Calculate PET using heatload adjusted Thornthwaite equation (using modified version of Redmond script)
tic()
pet_data <- pet_function(site=data$site, year=data$year, month=data$month,
                         slope=data$slope, latitude=data$latitude, aspect=data$aspect,
                         tmean=data$tmean)

pet_data <- pet_data[,c("site", "month", "year", "petm")]
data <- merge(data, pet_data, by = c("site", "month", "year"))


## Add comparison PET using SPEI package implementation of (non heatload adjusted) thornthwaite equation
data <- data %>% 
  as_tibble() %>% 
  arrange(site, year, month) %>% 
  group_by(site) %>% 
  nest() %>% 
  mutate(data = map(.x = data, .f = pet_spei_function)) %>% 
  unnest(data)


lm(petm~pet_spei, data = data) %>% summary()
lm(petm~pet_cru, data = data) %>% summary()
lm(pet_spei~pet_cru, data = data) %>% summary()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate CWD and save data --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl=makeCluster(8)
clusterExport(cl,c("data","setorder"))
registerDoParallel(cl)

cwd_data <- cwd_function(site=data$site, year=data$year, month=data$month,
                         petm = data$petm, tmean=data$tmean,  
                         ppt = data$pre_corrected, soilawc = data$swc)

cwd_data <- cwd_data[,c("site", "month", "year", "cwd")]
data <- data %>% as.data.table()
data <- merge(data, cwd_data, by = c("site", "month", "year"))

cwd_cru <- cwd_function(site=data$site, year=data$year, month=data$month,
                        petm = data$pet_cru, tmean=data$tmean,  
                        ppt = data$pre_corrected, soilawc = data$swc)
cwd_cru <- cwd_cru[,c("site", "month", "year", "cwd")]
names(cwd_cru) <- c("site", "month", "year", "cwd_cru")
data <- merge(data, cwd_cru, by = c("site", "month", "year"))



cwd_spei <- cwd_function(site=data$site, year=data$year, month=data$month,
                        petm = data$pet_spei, tmean=data$tmean,  
                        ppt = data$pre_corrected, soilawc = data$swc)
cwd_spei <- cwd_spei[,c("site", "month", "year", "cwd")]
names(cwd_spei) <- c("site", "month", "year", "cwd_spei")
data <- merge(data, cwd_spei, by = c("site", "month", "year"))
toc()


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Characterize missing data ----------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# input_data <- data %>% 
#   select(site_id, slope, latitude, longitude, aspect, pre_corrected, tmean, month, year, swc)
# 
# input_data %>% gg_miss_upset()
# any_missing <- input_data %>%
#   filter_all(any_vars(is.na(.)))
# share_missing <- dim(any_missing)[1] / dim(data)[1]
# share_missing
# miss_var_summary(data)
# 
# 
# out_data <- cwd_data %>% 
#   select(site, year, month, wm, petm, soilm, deltsoil, aet, cwd)
# 
# out_data %>% gg_miss_upset()
# any_missing <- input_data %>%
#   filter_all(any_vars(is.na(.)))
# share_missing <- dim(any_missing)[1] / dim(data)[1]
# share_missing
# miss_var_summary(data)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Write out file ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd_data <- cwd_data %>%
  select(site, year, month, tmean, ppt, aet, cwd, pet = petm, pet_cru, pet_spei)

fwrite(cwd_data,file=paste0(wdir,"1_input_processed/climate/essentialcwd_data.csv"))

