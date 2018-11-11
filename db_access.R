library(tidyverse)
library(dbplyr)
library(RSQLite)

# Define path
wdir = 'C:/Users/rheil/Google Drive/Treeconomics/Data/'
tree_db = paste0(wdir, 'tree_ring_data.db')

# Connect to database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)

# Identify all trees of a given species
spp = 'psme' #Whitebark = pial; douglas fir = 'psme'
tree_db = tbl(conn, 'trees')
tree_ids = tree_db %>%
  filter(species_id == doug_spp) %>%
  select('tree_id') %>%
  collect()

# Pull observations of identified trees
obs_db = tbl(conn, 'observations')
obs_db = obs_db %>%
  filter(tree_id %in% tree_ids$tree_id) %>%
  arrange(tree_id, desc(year))
obs = obs_db %>%
  collect()

# Export data
obs = select(obs, c(observation_id, tree_id, site_id, year, ring_width, first_year))
obs$age = obs$year - obs$first_year
write.csv(obs, file = paste0(wdir,spp,'_data.csv'), row.names = FALSE)


# obs_db = tbl(conn, 'observations')
# obs_db = obs_db %>%
#   filter(tree_id %in% c('BLK47B  ', 'PFL241d ')) %>%
#   arrange(tree_id, desc(year))
# obs = obs_db %>%
#   collect()
# 
# 
obs_db = tbl(conn, 'observations')
subset = obs_db %>%
  filter(tree_id == 'WR16A') %>%
  filter(site_id == 'wr20') %>%
  arrange(desc(year))

obs = subset %>% collect()
# 
# species_db = tbl(conn, 'species')
# spp = species_db %>%
#   filter(species_name == "douglas-fir")
# spp = spp %>% collect
