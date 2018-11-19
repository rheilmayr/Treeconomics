library(tidyverse)
library(dbplyr)
library(RSQLite)

# Define path
wdir = 'C:/Users/rheil/Google Drive/Treeconomics/Data/'
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
cwd_csv = paste0(wdir, 'essentialcwd_data.csv')

# Connect to database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)

# Identify all trees of a given species
spp = 'psme' #Whitebark = pial; douglas fir = 'psme'
tree_db = tbl(conn, 'trees')
tree_ids = tree_db %>%
  filter(species_id == spp) %>%
  select('tree_id') %>%
  collect()

# Pull observations of identified trees
obs_db = tbl(conn, 'observations_new')
obs_db = obs_db %>%
  filter(tree_id %in% tree_ids$tree_id) %>%
  arrange(tree_id, desc(year))
obs = obs_db %>%
  collect()

# Load and join cwd data
cwd_df = read.csv(cwd_csv, sep=',')
cwd_df.groups = cwd_df %>%
  group_by(site, year)
aet.annual = cwd_df.groups %>%
  summarise(aet.an = sum(aet))
cwd.annual = cwd_df.groups %>%
  summarise(cwd.an = sum(cwd))
clim.join = aet.annual %>%
  left_join(cwd.annual, by = c("site" = "site", "year" = "year"))
obs.join = obs %>%
  left_join(clim.join, by = c("site_id" = "site", "year" = "year"))

# Export data
obs.join = select(obs.join, c(observation_id, tree_id, site_id, year, ring_width, first_year, cwd.an, aet.an))
obs.join$age = obs.join$year - obs.join$first_year
write.csv(obs.join, file = paste0(wdir,spp,'_data.csv'), row.names = FALSE)

