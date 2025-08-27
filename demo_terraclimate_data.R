# 1. Add terraclimate raster data of historic climates
cwd_tc <- rast(paste0(wdir,"0_raw/TerraClimate/TerraClimate19611990_def.nc")) %>%
  sum()
pet_tc <- rast(paste0(wdir,"0_raw/TerraClimate/TerraClimate19611990_pet.nc")) %>%
  sum()
clim_tc <- rast(list("cwd" = cwd_tc, "pet" = pet_tc, "ppt" = ppt_tc))


# 2. Add terraclimate site-month-year data
tc_pet <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_pet.csv"))
tc_cwd <- read_csv(paste0(wdir,"0_raw/TerraClimate/itrdbsites_def.csv"))
site_clim_df_tc <- tc_pet %>%
  left_join(tc_cwd, by = c("collection_id", "Month", "year")) %>%
  left_join(tc_ppt, by = c("collection_id", "Month", "year")) %>%
  rename(month = Month,
         tc_pet = pet,
         tc_cwd = def)


# 3. Add water year
site_clim_df[,water_year:=year]
site_clim_df[(latitude>=0) & (month>=10),water_year:=year+1] # Northern hemisphere water year is october through september
site_clim_df[(latitude<0) & (month>=7),water_year:=year+1] # Southern hemisphere water year is July through June
site_clim_df <- site_clim_df %>% 
  as_tibble() %>% 
  select(-year) %>% 
  rename(year = water_year)


# 4. Calculate site-level annual climate
site_clim_df = site_clim_df %>%
  group_by(location_id, year) %>%
  summarise(tc_cwd.an = sum(tc_cwd),
            tc_pet.an = sum(tc_pet),
            .groups = "drop")


# 5. New example pull_clim function
# Pull and organize climate distribution for species
pull_clim <- function(spp_code, clim_raster){
  print(spp_code)
  
  # Pull relevant range map
  sp_range <- range_sf %>%
    filter(sp_code == spp_code)
  
  # Pull clim values
  clim_vals <- clim_raster %>% 
    mask(mask = sp_range, touches = TRUE) %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  
  return(clim_vals)
}

pull_clim_tc <- partial(.f = pull_clim, clim_raster = clim_tc)

clim_df_tc <- species_list %>%
  mutate(clim_vals = map(sp_code,.f = pull_clim_tc))
