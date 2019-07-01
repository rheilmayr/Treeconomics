library(tidyverse)
library(dbplyr)
library(RSQLite)

# Define path
wdir = 'C:/Users/rheil/Google Drive/Treeconomics/Data/'
#for Fran
wdir="C:/Users/fmoor/Google Drive/Treeconomics/Data/"
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
cwd_csv = paste0(wdir, 'essentialcwd_data.csv')

# Connect to database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)

# Pull all tree records for top 6 species
spp=c("psme","pipo","pisy","pcgl","pcab","pied")
speciesnames=c("Douglas Fir","Ponderosa Pine","Scotch Pine","White Spruce","Norway Spruce","Colorado Pinyon")

#collect trees from identified species
tree_db = as.data.frame(tbl(conn, 'trees'))
tree_ids = tree_db %>%
  filter(species_id %in% spp) %>%

  # Identify all trees of a given species
# Hypotheses: 1) quercus, sugarpine, cedar, white fir, juniper, redwood would exhibit spoiled tree
# High elevation tree species will exhibit less of an effect because high cwd is correlated with release from growth constraints 
species_ids = tree_db %>%
  select(species_id) %>%
  collect() %>%
  group_by(species_id) %>%
  count() %>%
  arrange(desc(n))

spp = 'pipo' #Whitebark = pial; douglas fir = 'psme'; ponderosa = 'pipo'
tree_db = tbl(conn, 'trees')

tree_ids = tree_db %>%
  filter(species_id == spp) %>%
  select('tree_id') %>%
  collect()

# Pull observations of identified trees
obs_db = tbl(conn, 'observations_new')
obs_db = obs_db %>%
  filter(tree_id %in% local(tree_ids$tree_id)) %>%
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

# Pull tree obsevations
obs_db = as.data.frame(tbl(conn, 'observations_new'))
obs_species=data.frame()
for(i in 1:length(spp)){
  trees=tree_ids %>%
    filter(species_id==spp[i]) %>%
    select(tree_id)
  obs = obs_db %>%
    filter(tree_id %in% trees$tree_id) %>%
    arrange(tree_id, desc(year))
  obs = obs %>%
    collect()
  obs.join = obs %>%
    inner_join(clim.join, by = c("site_id" = "site", "year" = "year")) #note that we loose a few sites because we are missing CWD data - probably becaue they are on the coast and more sites because we don't have cwd data before 1900
  obs.join$species=rep(spp[i],dim(obs.join)[1])
  obs_species=rbind(obs_species,obs.join)
  print(i)
}

obs_species$age=obs_species$year-obs_species$first_year
write.csv(obs_species, file = paste0(wdir,'topsixspecies_data.csv'), row.names = FALSE)


