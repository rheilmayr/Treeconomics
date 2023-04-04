#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Create predictions of growth impacts from climate change
#
# Input files:
# - ss_mod: R model object saved from Second stage
# - 
#
# ToDo:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(fixest)
library(patchwork)
library(tidyverse)
library(dtplyr)
library(prediction)
library(tictoc)
library(tmap)
library(tidylog)
library(broom)
library(purrr)
library(sf)
library(gstat)
library(units)
library(sjPlot)

select <- dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'out/first_stage/site_pet_cwd_std.csv')) 
  #select(-X1)

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))


# 3. Site information
site_df <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id)) %>%
  mutate(geog = str_replace_all(collection_id, "[:digit:]", ""))

flm_df <- flm_df %>% 
  left_join(ave_site_clim, by = c("collection_id")) %>%
  left_join(site_df, by = c("collection_id"))

# 4. Spatially proximate blocks of sites
load(file=paste0(wdir,"out/spatial_blocks.Rdat"))


# Identify and trim extreme outliers
cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
pet_est_bounds = quantile(flm_df$estimate_pet.an, c(0.01, 0.99),na.rm=T)
cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)

flm_df <- flm_df %>% 
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]))

trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na() %>%
  mutate(cwd_errorweights = 1 / (std.error_cwd.an))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper function to trim to geography --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study_species <- c("psme", "pipo")
study_geog <- c("UT", "NM", "AZ", "CO", "WY", "MT", "ID")


study_df <- trim_df %>%
  filter(geog %in% study_geog,
         species_id %in% study_species)
full_df <- trim_df %>%
  filter(species_id %in% study_species)


study_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, 
                data = study_df, weights = cwd_errorweights)
study_mod %>% summary()

study_df %>% ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1)


full_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, 
               data = full_df, weights = cwd_errorweights)
full_mod %>% summary()

full_df %>% ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) + geom_point()


## Eventually should probably shift to our full second stage structure that includes squared terms
study_mod <- lm(estimate_cwd.an ~ cwd.spstd + I(cwd.spstd**2) +  pet.spstd + I(pet.spstd**2), 
                data = study_df, weights = cwd_errorweights)
study_mod %>% summary()

full_mod <- lm(estimate_cwd.an ~ cwd.spstd + I(cwd.spstd**2) +  pet.spstd + I(pet.spstd**2), 
               data = full_df, weights = cwd_errorweights)
full_mod %>% summary()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot sensitivity --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
)


##=================================================================================================================
##             PONDEROSA PINE         
##==================================================================================================================

plot_df <- flm_df %>%
  filter(species_id == "pipo")

allpipo=plot_df %>%
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd, color="x", fill="y")) +
  geom_point()+
  geom_smooth(method="loess", span=1)+
  labs(y="Marginal effect of CWD", x="Historic CWD (SD)")+
  #ggtitle("Pinus ponderosa (PIPO)")+
  scale_fill_manual(values='#21908CFF')+
  scale_color_manual(values = '#440154FF')+
  guides(fill=F, color=F)

allpipo

map+allpipo+plot_annotation(tag_levels="A") & theme(
                plot.tag = element_text(face = 'bold', size=15, family ="Helvetica"),
                text=element_text(family ="Helvetica"))

# unique(plot_df$collection_id)
# 
# ## Look at specific site and its nearby sites
# site_id <- "AZ036"
# site_id <- "NM512"
# site_id <- "WY058"
# 
# proximate_ids <- block_list[site_id] %>% unlist(use.names = FALSE)
# 
# plot_df_small <- plot_df %>%
#   filter(collection_id %in% proximate_ids)
# 
# azpipo <- plot_df_small %>%
#   ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
#   geom_point()+
#   geom_smooth(method="loess", span=1)+
#   labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
#   ggtitle(paste("PIPO: < 400km of", site_id))
# 
# azpipo
# 
# nmpipo <- plot_df_small %>%
#   ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
#   geom_point()+
#   geom_smooth(method="loess", span=1)+
#   labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
#   ggtitle(paste("PIPO: < 400km of", site_id))
# 
# nmpipo
# 
# wypipo <- plot_df_small %>%
#   ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
#   geom_point()+
#   geom_smooth(method="loess", span=1)+
#   labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
#   ggtitle(paste("PIPO: < 400km of", site_id))
# 
# wypipo
# 
# pipolots <- allpipo | (azpipo/nmpipo/wypipo)
# pipolots + plot_annotation(tag_levels="A") & theme(
#                 plot.tag = element_text(face = 'bold', size=15, family ="Helvetica"),
#                 text=element_text(family ="Helvetica"))


# ##=================================================================================================================
# ##                PSME      
# ##==================================================================================================================
# plot_df <- flm_df %>%
#   filter(species_id == "psme")
# 
# allpsme=plot_df %>%
#   filter(estimate_cwd.an<1) %>% 
#   ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
#   geom_point()+
#   geom_smooth(method="loess", span=1)+
#   labs(y="Marginal effect of CWD", x="Historic CWD (SD)")+
#   ggtitle("Pseudotsuga menziesii (psme)")
# 
# allpsme
# unique(plot_df$collection_id)
# 
# ## Look at specific site and its nearby sites
# site_id <- "AZ558"
# site_id <- "MT154"
# site_id <- "CO052"
# 
# proximate_ids <- block_list[site_id] %>% unlist(use.names = FALSE)
# 
# plot_df_small <- plot_df %>%
#   filter(collection_id %in% proximate_ids)
# 
# azpsme <- plot_df_small %>%
#   ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
#   geom_point()+
#   geom_smooth(method="loess", span=1)+
#   labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
#   ggtitle(paste("PSME: < 400km of", site_id))
# 
# azpsme
# 
# mtpsme <- plot_df_small %>%
#   ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
#   geom_point()+
#   geom_smooth(method="loess", span=1)+
#   labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
#   ggtitle(paste("PSME: < 400km of", site_id))
# 
# mtpsme
# 
# copsme <- plot_df_small %>%
#   ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
#   geom_point()+
#   geom_smooth(method="loess", span=1)+
#   labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
#   ggtitle(paste("PSME: < 400km of", site_id))
# 
# copsme
# 
# psmelots <- allpsme | (azpsme/mtpsme/copsme)
# psmelots + plot_annotation(tag_levels="A") & theme(
#   plot.tag = element_text(face = 'bold', size=15, family ="Helvetica"),
#   text=element_text(family ="Helvetica"))

##=================================================================================================================
##                PIED     
##==================================================================================================================
plot_df <- flm_df %>%
  filter(species_id == "pied")

allpied=plot_df %>%
  filter(estimate_cwd.an<1) %>% 
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
  geom_point()+
  geom_smooth(method="loess", span=1)+
  labs(y="Marginal effect of CWD", x="Historic CWD (SD)")+
  ggtitle("Pinus edulis (PIED)")

allpied
unique(plot_df$collection_id)

## Look at specific site and its nearby sites
site_id <- "AZ102"
site_id <- "NM038"
site_id <- "UT534"

proximate_ids <- block_list[site_id] %>% unlist(use.names = FALSE)

plot_df_small <- plot_df %>%
  filter(collection_id %in% proximate_ids)

azpied <- plot_df_small %>%
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
  geom_point()+
  geom_smooth(method="loess", span=1)+
  labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
  ggtitle(paste("PIED: < 400km of", site_id))

azpied

nmpied <- plot_df_small %>%
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
  geom_point()+
  geom_smooth(method="loess", span=1)+
  labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
  ggtitle(paste("PIED: < 400km of", site_id))

nmpied

utpied <- plot_df_small %>%
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd)) +
  geom_point()+
  geom_smooth(method="loess", span=1)+
  labs(y="Marg. effect CWD", x="Historic CWD (SD)")+
  ggtitle(paste("PIED: < 400km of", site_id))

utpied

piedlots <- allpied | (azpied/nmpied/utpied)
piedlots + plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=15, family ="Helvetica"),
  text=element_text(family ="Helvetica"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Identify spatially proximate blocks of sites ---------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_list <- trim_df %>%
  pull(collection_id) %>%
  unique()
n_sites <- length(site_list)

site_dist=st_distance(site_points)
rownames(site_dist)=site_points$collection_id
colnames(site_dist)=site_points$collection_id
# save(site_dist,file=paste0(wdir,"out/site_distances.Rdat"))
# load(paste0(wdir,"out/site_distances.Rdat"))

dist_df <- as_tibble(site_dist) %>% 
  drop_units() 


threshold <- 900000
dist_df <- dist_df %>%
  lazy_dt() %>% 
  mutate(collection_id = names(dist_df)) %>% 
  # select(collection_id, site_list) %>% 
  # filter(collection_id %in% site_list) %>% 
  mutate(across(.cols = !collection_id, ~(.x < threshold))) %>% 
  # mutate(across(.cols = !collection_id, ~ifelse((.x < range), collection_id, "DROP"))) %>% 
  as_tibble()

block_list <- c()
for (site in site_list){
  block_sites <- dist_df %>% 
    filter(get(site) == TRUE) %>% 
    pull(collection_id)
  block_list[site] <- list(block_sites)
}
#save(block_list,file=paste0(wdir,"out/spatial_blocks_", as.character(threshold/1000), ".Rdat"))
load(file=paste0(wdir,"out/spatial_blocks.Rdat"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Simulate regional studies --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand_order <- flm_df %>%
  filter(species_id =="pipo") %>%
  pull(collection_id) %>%
  sample()

study_sites <- list()
included_sites <- list()

'%ni%' <- function(x,y)!('%in%'(x,y))

for (id in rand_order) {
  if (id %ni% included_sites) {
    study_sites <- study_sites %>% append(id) %>% unlist()
    proximate_ids <- block_list[id] %>% unlist(use.names = FALSE)
    included_sites <- included_sites %>% append(proximate_ids) %>% unlist()
  }
}

block_df <- tibble("collection_id" = study_sites)

estimate_block_ss <- function(id) {
  proximate_ids <- block_list[id] %>% unlist(use.names = FALSE)
  
  sp_id <- flm_df %>%
    filter(collection_id == id) %>%
    pull(species_id)
  proximate_ids <- flm_df %>%
    filter(species_id == sp_id,
           collection_id %in% proximate_ids) %>%
    pull(collection_id)
  
  
  if (length(proximate_ids)<1) {
    return(out = tibble(NaN))
  } else {
    reg_df <- flm_df %>%
      filter(collection_id %in% proximate_ids)
    
    reg_mod <- lm(estimate_cwd.an ~ cwd.spstd, data = reg_df)
    coef <- reg_mod %>% 
      tidy() %>%
      filter(term == "cwd.spstd")
    return(out = tibble(coef))
  }
}

block_df <- block_df %>%
  mutate(cwd_spstd_coef = map(collection_id, .f = estimate_block_ss))

block_df <- block_df %>%
  unnest(cwd_spstd_coef) %>%
  select(-`NaN`, -term) %>%
  drop_na()

block_df <- block_df %>%
  left_join(flm_df, by = "collection_id")

# block_df %>%
#   mutate(sig = p.value < 0.05) %>%
#   filter(estimate >-100, estimate < 100) %>%
#   ggplot(aes(x = pet.spstd, y = estimate, color = sig)) +
#   geom_point() +
#   # geom_smooth() +
#   theme_bw()

block_df <- block_df %>% 
  # filter(p.value < 0.05) %>%
  mutate(sig = p.value < 0.05,
         effect = ifelse(sig & (estimate > 0), "Drought-naive", ifelse(sig & (estimate < 0), "range-edge", "neither")),
         spoiled = estimate > 0)

block_df %>%
  group_by(effect) %>%
  tally()


block_df %>%
  mutate(sig = p.value < 0.05) %>%
  filter(estimate >-10, estimate < 10) %>%
  ggplot(aes(x = cwd.spstd, y = estimate, color=effect, label=collection_id)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_text(hjust=0, vjust=0)

block_df %>%
  mutate(sig = p.value < 0.05) %>%
  filter(estimate >-10, estimate < 10) %>%
  ggplot(aes(x = cwd.spstd, y = estimate, color=effect)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0, linetype="dashed")


block_df %>%
  filter(pet.spstd > -0.25, estimate < -0.5)


mod_df <- flm_df %>% 
  filter(species_id == "pied") 

mod <- lm( ~ pet.spstd + cwd.spstd, data = mod_df)
summary(mod)
tab_model(mod, show.ci = F, show.se = T)
