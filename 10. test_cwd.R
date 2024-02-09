
library(tidyverse)
source("f_cwd_function_new.R")
library(SPEI)
library(zoo)
library(lubridate)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Pre-processed climate and soil data
data=fread(paste0(wdir,"1_input_processed/climate/sitedataforcwd.csv"))

# 2. Treecon climate data
clim_df <- read_csv(file=paste0(wdir,"1_input_processed/climate/essentialcwd_data.csv"))


# 3. Terraclimate reference data
tc_pet <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_pet.csv"))
tc_cwd <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_def.csv"))
tc_df <- tc_pet %>% 
  left_join(tc_cwd, by = c("collection_id", "Month", "year"))

# tc_df <- tc_df %>% 
#   mutate(date = as.yearmon(year, Month), "%Y %m",
#          days = days_in_month(date))

tc_df <- tc_df %>% 
  # mutate(pet = pet*days,
  #        def = def*days) %>% 
  group_by(collection_id, year) %>% 
  summarise(pet_tc = sum(pet),
            cwd_tc = sum(def))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compare data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clim_df <- clim_df %>%
  rename(collection_id = site) %>% 
  inner_join(tc_df, by = c("collection_id", "year"))

mod <- lm(cwd ~ cwd_tc, clim_df)
summary(mod)

mod <- lm(pet ~ pet_tc, clim_df)
summary(mod)

clim_df %>% 
  ggplot(aes(x = cwd_tc, y = cwd)) +
  geom_point(alpha = .3)


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
  select(site_id, year, month, latitude, longitude, pre_corrected, tmean, swc)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Focus on example site --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_data <- data %>% 
  filter(site_id == "CA524") # PIPO site in CA

## Test annual cwd calculation
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

year = 10
start_month = (year - 1) * 12 + 1
end_month = start_month + 11
start_month
end_month
cwd_data$pet[start_month:end_month] %>% sum()

## Compare to SPEI's implementation of PET calculation
pet <- thornthwaite(test_data$tmean, test_data$latitude[1])
pet[start_month:end_month] %>% sum()



# ## Test monthly norm cwd calculation
# norm_data <- test_data %>% 
#   group_by(site_id, month) %>%
#   summarise(slope = mean(slope),
#           latitude = mean(latitude),
#           aspect = mean(aspect),
#           pre_corrected = mean(pre_corrected),
#           tmean = mean(tmean),
#           swc = mean(swc)
#   )
# 
# cwd_data <- cwd_function(site=norm_data$site_id,
#                          slope=norm_data$slope,
#                          latitude=norm_data$latitude,
#                          foldedaspect=norm_data$aspect,
#                          ppt=norm_data$pre_corrected,
#                          tmean=norm_data$tmean,
#                          month=norm_data$month,
#                          soilawc=norm_data$swc,
#                          type="normal")
# cwd_data
