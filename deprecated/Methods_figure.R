#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 7/6/20
# Purpose: Make a figure showing detrending
#
# Input files:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




# Define path
wdir <- 'NewRemote/'

# Read in data
dendro_dir <- paste0(wdir, "out/dendro/")
dendro_df <- read.csv(paste0(dendro_dir, "rwi_long.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            rwl = mean(rwl),
            .groups = "drop") %>% 
  as_tibble()


## Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)

dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')

# Species information and calling pipo
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))

pipo_list <- site_smry %>% 
  filter(species_id == "pipo") %>% 
  pull(collection_id)

## chosing example sites
ex_sites <- c("AZ026", "AZ050")

pipo_dat <- dendro_df %>% 
  filter(collection_id %in% ex_sites) %>% 
  data.frame()

## modifying dataframe into correct format for dplR
widepipo <- pipo_dat %>% 
  select(collection_id, rwl, year, tree) %>% 
  unite("newnames", collection_id,tree, sep="_", remove = F) %>% 
  select(newnames, year, rwl) %>% 
  pivot_wider(names_from = newnames, values_from = rwl)

## detrending
detrendpipo <- dplR::detrend(rwl=widepipo,method = "Spline", make.plot = T, verbose = FALSE)


## chosing plots



## estimating the chronology
chronpipo <- dplR::chron(detrendpipo, prefix = "CAM")
plot(chronpipo, add.spline=TRUE, nyrs=20)
