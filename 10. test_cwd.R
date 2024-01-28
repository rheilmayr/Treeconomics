
library(tidyverse)
source("f_cwd_function.R")
library(SPEI)
library(ncdf4)
library(sf)
library(tictoc)
library(terra)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Pre-processed climate and soil data
data=fread(paste0(wdir,"1_input_processed/climate/sitedataforcwd.csv"))

# 2. Reference CWD for California (https://www.sciencebase.gov/catalog/item/5f29c62d82cef313ed9edb39)
bcm_dir <- paste0(wdir, "0_raw/BCM/cwd1990/")
bcm_files <- list.files(bcm_dir, full.names = TRUE)
bcm_rast <- rast(bcm_files)
crs(bcm_rast) = crs("epsg:3310")

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
# data$tmx_corrected=data$tmx+data$tmax_correction
# data$tmn_corrected=data$tmn+data$tmin_correction


# Unit conversions
data$swc=data$swc/10 #convert swc from mm to cm
data$tmean=data$tmp
data$slope <- data$slope * 57.2958 # convert slope from radians to degrees
data$aspect <- data$aspect * 57.2958 # convert aspect from radians to degrees
data$petd <- data$pet
data$pre_corrected <- data$pre

data <- data[,c("site_id", "year", "month", "latitude", "longitude", 
                "pre_corrected", "tmean", "petd", "swc")]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Focus on example site --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_data <- data %>% 
  filter(site_id == "CA524") # PIPO site in CA
site_loc <- st_as_sf(test_data[1,], coords = c("longitude","latitude"), crs = 4326)

## Test annual cwd calculation
cwd_data <- cwd_function(site=test_data$site_id,
                         year=test_data$year,
                         month=test_data$month,
                         latitude=test_data$latitude,
                         ppt=test_data$pre_corrected,
                         tmean=test_data$tmean,
                         petd=test_data$petd,
                         soilawc=test_data$swc,
                         type="annual")


year = 90
start_month = (year - 1) * 12 + 1
end_month = start_month + 11
start_month
end_month
cwd_data$cwd[start_month:end_month] %>% sum()


## Compare to USGS CWD (within California)
site_loc <- site_loc %>% st_transform(crs(bcm_rast))
cwd_reference <- extract(bcm_rast, site_loc)

cwd_reference <- as_tibble(cwd_reference) %>% 
  select(cwd1990jan, cwd1990feb, cwd1990mar, cwd1990apr, cwd1990may, cwd1990jun, cwd1990jul,
         cwd1990aug, cwd1990sep, cwd1990oct, cwd1990nov, cwd1990dec)

sum(cwd_reference)
