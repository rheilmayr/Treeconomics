
library(tidyverse)
source("f_cwd_function_new.R")
library(SPEI)
library(zoo)
library(lubridate)
library(fixest)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Pre-processed climate and soil data
site_df <- fread(paste0(wdir,"1_input_processed/climate/sitedataforcwd.csv"))
lat_df <- site_df %>% 
  select(collection_id = site_id, latitude) %>% 
  unique() %>% 
  as_tibble()


# 2. Treecon climate data
clim_df <- read_csv(file=paste0(wdir,"1_input_processed/climate/essentialcwd_data.csv"))


# 3. Terraclimate reference data
tc_pet <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_pet.csv"))
tc_cwd <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_def.csv"))
tc_month_df <- tc_pet %>% 
  left_join(tc_cwd, by = c("collection_id", "Month", "year")) %>% 
  rename(month = Month,
         pet_tc = pet,
         cwd_tc = def)

# tc_df <- tc_df %>% 
#   mutate(date = as.yearmon(year, Month), "%Y %m",
#          days = days_in_month(date))

tc_df <- tc_month_df %>% 
  # mutate(pet = pet*days,
  #        def = def*days) %>% 
  group_by(collection_id, year) %>% 
  summarise(pet_tc = sum(pet_tc),
            cwd_tc = sum(cwd_tc))

ave_site_clim <- read_rds(paste0(wdir, "2_output/climate/site_ave_clim.gz"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compare data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clim_df <- clim_df %>%
  rename(collection_id = site) %>% 
  inner_join(tc_month_df, by = c("collection_id", "year", "month"))

clim_df <- clim_df %>% 
  left_join(lat_df, by = "collection_id") 

calc_pet_thornthwaite_spei <- function(data){
  data <- data %>% 
    arrange(year, month) 
  data$pet_spei <- thornthwaite(Tave = data$tmean, lat = data$latitude[1])
  return(data)
}

clim_df <- clim_df %>% 
  group_by(collection_id) %>% 
  nest() %>% 
  mutate(data = map(.x = data, .f = calc_pet_thornthwaite_spei)) %>% 
  unnest()

clim_df <- clim_df %>% 
  ungroup()

mod <- lm(cwd_tc ~ cwd, clim_df)
summary(mod)

mod <- feols(cwd ~ cwd_tc | collection_id, clim_df)
summary(mod)

mod <- lm(pet ~ pet_tc, clim_df)
summary(mod)

mod <- feols(pet ~ pet_tc | collection_id, clim_df)
summary(mod)

mod <- lm(pet ~ pet_spei, clim_df)
summary(mod)

mod <- feols(pet ~ pet_spei | collection_id, clim_df)
summary(mod)

mod <- lm(pet_tc ~ pet_spei, clim_df)
summary(mod)


clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = pet_tc, x = pet)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")

clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = pet_spei, x = pet)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")

clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = pet_tc, x = pet_spei)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")


clim_df %>%
  sample_n(30000) %>% 
  ggplot(aes(y = cwd_tc, x = cwd)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")


an_clim_df <- clim_df %>% 
  group_by(collection_id, year) %>% 
  summarize(pet_tc = sum(pet_tc),
            cwd_tc = sum(cwd_tc),
            cwd = sum(cwd),
            pet = sum(pet),
            pet_spei = sum(pet_spei)) %>% 
  ungroup()

mod <- lm(cwd_tc ~ cwd, an_clim_df)
summary(mod)

mod <- feols(cwd_tc ~ cwd | collection_id, an_clim_df)
summary(mod)

mod <- lm(pet_tc ~ pet, an_clim_df)
summary(mod)

mod <- feols(pet_tc ~pet | collection_id, an_clim_df)
summary(mod)

mod <- lm(pet_spei ~ pet, an_clim_df)
summary(mod)

mod <- feols(pet ~ pet_spei | collection_id, an_clim_df)
summary(mod)

mod <- lm(pet_tc ~ pet_spei, an_clim_df)
summary(mod)

an_clim_df %>%
  sample_n(10000) %>% 
  ggplot(aes(y = cwd_tc, x = cwd)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")


an_clim_df %>%
  sample_n(10000) %>% 
  ggplot(aes(y = pet_tc, x = pet)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")

an_clim_df %>%
  sample_n(10000) %>% 
  ggplot(aes(y = pet_spei, x = pet)) +
  geom_point(alpha = .3) +
  geom_smooth() +
  geom_abline(intercept = 0, slope = 1, size = 0.5, color = "red")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set NAs for precip and temp
site_df$pre = na_if(site_df$pre,-999)
site_df$tmn = na_if(site_df$tmn,-999)
site_df$tmx = na_if(site_df$tmx,-999)

# Add corrections from World Clim to CWD to get downscaled variables
site_df$pre_corrected=site_df$pre*site_df$pre_correction
site_df$tmx_corrected=site_df$tmx+site_df$tmax_correction
site_df$tmn_corrected=site_df$tmn+site_df$tmin_correction


# Unit conversions
site_df$swc=site_df$swc/10 #convert swc from mm to cm
site_df$tmean=(site_df$tmn_corrected+site_df$tmx_corrected)/2 
site_df$ppt=site_df$pre_corrected
site_df$slope <- site_df$slope * 57.2958 # convert slope from radians to degrees
site_df$aspect <- site_df$aspect * 57.2958 # convert aspect from radians to degrees

site_df <- site_df %>% 
  select(site_id, year, month, latitude, longitude, aspect, slope, pre_corrected, tmean, swc)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Focus on example site --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_site <- "CA524"
test_data <- site_df %>% 
  filter(site_id == test_site) # PIPO site in CA

site=test_data$site_id
slope=test_data$slope
latitude=test_data$latitude
aspect=test_data$aspect
ppt=test_data$pre_corrected
tmean=test_data$tmean
month=test_data$month
year=test_data$year
soilawc=test_data$swc
type="annual"

## Test annual cwd calculation
cwd_data <- cwd_function(site=test_data$site_id,
                         slope=test_data$slope,
                         latitude=test_data$latitude,
                         aspect=test_data$aspect,
                         ppt=test_data$pre_corrected,
                         tmean=test_data$tmean,
                         month=test_data$month,
                         year=test_data$year,
                         soilawc=test_data$swc,
                         type="annual")

year = 15
start_month = (year - 1) * 12 + 1
end_month = start_month + 11
start_month
end_month
cwd_data$pet[start_month:end_month] %>% sum()


## Compare to SPEI's implementation of PET calculation
pet <- thornthwaite(test_data$tmean, test_data$latitude[1])
cwd_data$pet_spei <- pet
pet[start_month:end_month] %>% sum()

spei_compare_mod <- lm(pet ~ pet_spei, cwd_data)
summary(spei_compare_mod)


## Compare to Terraclimate values
cwd_data <- cwd_data %>%
  rename(collection_id = site) %>% 
  left_join(tc_month_df, by = c("collection_id", "year", "month"))

tc_compare_mod <- lm(pet ~ pet_tc, cwd_data)
summary(spei_compare_mod)





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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compare site averages --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod <- lm(cwd.ave.tc ~ cwd.ave, data = ave_site_clim)
summary(mod)

ave_site_clim %>% 
  ggplot(aes(x = cwd.ave, y = cwd.ave.tc)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth()


mod <- lm(ppt.ave.tc ~ ppt.ave, data = ave_site_clim)
summary(mod)
ave_site_clim %>% 
  ggplot(aes(x = ppt.ave, y = ppt.ave.tc)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth()


mod <- lm(pet.ave.tc ~ pet.ave, data = ave_site_clim)
summary(mod)

ave_site_clim %>% 
  ggplot(aes(x = pet.ave, y = pet.ave.tc)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth()
