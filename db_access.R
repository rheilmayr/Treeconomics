library(tidyverse)
library(dbplyr)
library(RSQLite)

# Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
# #for Fran
# wdir="C:/Users/fmoor/Google Drive/Treeconomics/Data/"
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
cwd_csv = paste0(wdir, 'essentialcwd_data.csv')

# Connect to database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)

# Open tables
tree_db = as.data.frame(tbl(conn, 'trees'))
spp_db = as.data.frame(tbl(conn,'species'))

# Pull all tree records for top 20 species
#collect trees from identified species
count_spp <- tree_db %>% 
  group_by(species_id) %>%
  tally() %>%
  arrange(-n)

spp = count_spp[1:20,] %>% pull(species_id)
spp_names = spp_db %>%
  filter(species_id %in% spp) %>%
  as_tibble()

# ID trees from selected species
tree_ids = tree_db %>%
  filter(species_id %in% spp) %>%
  select('tree_id') %>%
  collect()

# Pull observations of identified trees
obs_db = tbl(conn, 'observations_new')
obs_db = obs_db %>%
  filter(tree_id %in% local(tree_ids$tree_id)) %>%
  arrange(tree_id, desc(year))
obs = obs_db %>%
  collect()

# Load and combine cwd data
cwd_df = read.csv(cwd_csv, sep=',')
cwd_df.groups = cwd_df %>%
  group_by(site, year)
aet.annual = cwd_df.groups %>%
  summarise(aet.an = sum(aet))
cwd.annual = cwd_df.groups %>%
  summarise(cwd.an = sum(cwd))
clim.join = aet.annual %>%
  left_join(cwd.annual, by = c("site" = "site", "year" = "year"))

# Pull tree obsevations
obs_db = as.data.frame(tbl(conn, 'observations_new'))
obs = obs_db %>%
  filter(tree_id %in% pull(tree_ids)) %>%
  arrange(tree_id, desc(year))%>%
  collect()

# Add climate data
obs = obs %>%
  inner_join(clim.join, by = c("site_id" = "site", "year" = "year")) #note that we loose a few sites because we are missing CWD data - probably becaue they are on the coast and more sites because we don't have cwd data before 1900

# Add species labels
obs <- obs %>%
  inner_join(tree_db, by = c("tree_id", "site_id"))

obs <- obs %>%
  inner_join(spp_names, by = 'species_id')

# Write to shared directory
obs$age=obs$year-obs$first_year
write.csv(obs, file = paste0(wdir,'species_data.csv'), row.names = FALSE)


