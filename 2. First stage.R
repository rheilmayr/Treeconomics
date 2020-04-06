#devtools::install_github("tidyverse/tidyr") #if necessary for pivot_wider function
library(tidyr)
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(lfe)
library(broom)
library(purrr)


# add nobs
# fix joins to prevent duplicate species_id


# Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
# #for Fran
wdir="C:/Users/fmoore/Google Drive/Treeconomics/Data/"
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
cwd_csv = paste0(wdir, 'essentialcwd_data.csv')

# Connect to database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tables = dbListTables(conn)

# Open tables
tree_db = as.data.frame(tbl(conn, 'trees'))
spp_db = as.data.frame(tbl(conn,'species'))
site_db = as.data.frame(tbl(conn, "sites"))
obs_db = tbl(conn, 'observations_new')

# Define regression model for first stage
fit_mod <- function(d) {
  felm(ring_width ~ cwd.an + pet.an + age+I(age^2)+I(age^3)+I(age^4)|tree_id|0|tree_id+year, data = d )
}
getcov <- function(m){
  return(m$vcv[2,1])
}

# Load site climate data
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) %>%
  select(-site)

clim_df = cwd_df %>%
  group_by(site_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum(petm),
            temp.an = mean(tmean),
            ppt.an = sum(ppt))
hist_clim_df <- clim_df %>%
  group_by(site_id) %>%
  filter(year<1980) %>%
  summarise(aet.ave = mean(aet.an),
            cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an),
            temp.ave = mean(temp.an),
            ppt.ave = mean(ppt.an))

# Create list of species
sp_list <- spp_db %>%
  pull(species_id)
sp_id <- "psme"

df_list = list()
for (i in 1:length(sp_list)) {
  sp_id = sp_list[i]

  # ID trees from selected species
  tree_ids = tree_db %>%
    filter(species_id == sp_id) %>%
    select('tree_id') %>%
    collect()
  
  # Pull observations of identified trees
  obs = obs_db %>%
    filter(tree_id %in% local(tree_ids$tree_id)) %>%
    arrange(tree_id, desc(year)) %>%
    collect()
  
  # Add climate data
  obs = obs %>%
    inner_join(clim_df, by = c("site_id", "year")) #note that we loose a few sites because we are missing CWD data - probably becaue they are on the coast and more sites because we don't have cwd data before 1900
  obs = obs %>%
    left_join(hist_clim_df, by = c("site_id"))
  obs <- obs %>%
    inner_join(tree_db, by = c("tree_id", "site_id"))
  
  # Add tree age
  obs$age=obs$year-obs$first_year
  
  ###### Set up for first stage #####
  obs <- obs %>%
    mutate(tree_id = paste0(site_id, tree_id))
  
  # Remove trees with negative ages
  df_trim <- obs %>%
    filter(age>0)
  print(paste0(sp_id, " - dropping observations due to invalid age: ",(df_trim$age <=0) %>% sum()))
  
  
  # Remove trees with a very short record (<5 years)
  treeobs=df_trim %>%
    group_by(tree_id) %>%
    summarize(nyears=n())
  df_trim=left_join(df_trim,treeobs, by = 'tree_id')
  print(paste0(sp_id, " - dropping trees due to short record: ",(treeobs$nyears <=5) %>% sum()))
  df_trim = df_trim %>%
    filter(nyears>5)
  
  # Remove sites with few trees
  ntrees=df_trim%>%
    group_by(site_id,species_id)%>%
    summarize(ntrees=length(unique(tree_id)))
  df_trim=left_join(df_trim,ntrees, by = c("site_id", "species_id"))
  print(paste0(sp_id, " - dropping sites due to few trees: ",(ntrees$ntrees <=5) %>% sum()))
  df_trim=df_trim%>%
    filter(ntrees>5)
  
  
  # Remove NA and move on to next interation in loop if no observations remain
  complete_df <- df_trim %>%
    drop_na(c("cwd.an","aet.an","pet.an", "ring_width"))
  if (dim(complete_df)[1]==0) {
    print(paste0("No valid sites for: ", sp_id))
    next
  }
  
  # Count number of valid tree-year observations per site
  nobs <- complete_df %>%
    group_by(site_id, species_id) %>%
    summarize(nobs = n())
  
  # Run first stage
  site_lm <- complete_df %>% 
    group_by(site_id, species_id) %>%
    nest() %>%
    mutate(mod = map(data, fit_mod),
           cwd.pet.cov = map(mod, getcov),
           mod = map(mod, tidy)) %>%
    unnest(mod,cwd.pet.cov) %>%
    filter(term %in% c('cwd.an', 'pet.an'))
  
  siteCoef <- site_lm %>%
    pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value")) %>%
    select(-data)
  
  # Attach climate, nobs and ntrees for each site
  site_df <- obs %>%
    select(c("site_id", "species_id")) %>%
    distinct() %>%
    left_join(hist_clim_df, by = c("site_id")) %>%
    left_join(siteCoef, by = c("site_id", "species_id")) %>%
    left_join(nobs, by = c("site_id", "species_id")) %>%
    left_join(ntrees, by = c("site_id", "species_id"))
  
  # Calculate standardized historic climate relative to species mean / std
  site_df <- site_df %>%
    mutate(cwd.spstd = scale(cwd.ave)[,1],
           aet.spstd = scale(aet.ave)[,1],
           pet.spstd = scale(pet.ave)[,1],
           ppt.spstd = scale(ppt.ave)[,1],
           temp.spstd = scale(temp.ave)[,1])

  df_list[[i]] <- site_df
}

full_df = bind_rows(df_list)
write.csv(full_df, paste0(wdir, "first_stage.csv"))