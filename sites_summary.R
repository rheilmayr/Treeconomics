library(tidyverse)
library(dbplyr)
library(RSQLite)
library(rnaturalearthdata)


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
sites = as.data.frame(tbl(conn, 'sites')) %>%
  collect()

# Plot locations of sites
map <- sites %>% 
  ggplot(aes(x = longitude, y = latitude)) +  
  geom_hex(bins = 50) 

coast_sf <- ne_coastline(scale = "large", returnclass = "sf")
+ geom_sf(data = coast_sf, alpha = 0, size=0.5,fill="white")

sr <- "+proj=longlat +datum=WGS84 +no_defs"

  