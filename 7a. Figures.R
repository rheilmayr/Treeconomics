#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run regressions to explore impact of historical climate on weather sensitivity
#
# Input files:
# - clim_niche.csv: Data documenting historic climate across each species range
#     Generated using species_niche.R
# - tree_ring_data_V2.db: Compiled database of ITRDB observations 
# - essentialcwd_data.csv: File detailing plot-level weather history
#
# ToDo:
# - Think through places where Monte carlo uncertainties can be added to figures
# - track down lost observations - currently dropping a lot due to NAN or failed RWI generation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
library(stringr)
library(raster)
library(rgdal)
library(viridis)
library(patchwork)
library(Hmisc)
library(latticeExtra)
library(prediction)
library(colorspace)
library(ggnewscale)
library(reshape2)
library(ggExtra)
library(ggsn)
library(maptools)
library(broom)
library(ggExtra)
library(extrafont)
library(marginaleffects)
library(tmap)
library(fixest)
library(forcats)
library(car)
librarian::shelf(ggplotify)

loadfonts(device = "win")
theme(family="Serif")

select <- dplyr::select
summarize <- dplyr::summarize

options(scipen=999)

theme_set(
  theme_bw(base_size = 25)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv"))

# 3. Site information
site_smry <- read_csv(paste0(wdir, 'out/dendro/site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id, latitude, longitude) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)
site_loc <- site_smry %>%
  select(collection_id, latitude, longitude)
flm_df <- flm_df %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


# 4. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_smry <- site_smry %>% 
  left_join(sp_info, by = c("species_id"))

trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()

# DNLM results
dnlm_results <- read_rds(paste0(wdir, "out/first_stage/dnlm_lagged_effects"))
# dnlm_results <- read_rds(paste0(wdir, "out/first_stage/dnlm_orig")) # Something has changed - old version creates positive pet effect, but gone after july tweaks to data...

# 5. Prediction rasters
rwi_list <- list.files(paste0(wdir, "out/predictions/pred_10000/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))

## 6. Second stage model
mod_df <- read_rds(paste0(wdir, "out/second_stage/ss_bootstrap.rds"))

# 7. Genus second stage model
genus_models <- readRDS(paste0(wdir, "out/second_stage/ss_conley_genus.rds"))

sp_models <- readRDS(paste0(wdir, "out/second_stage/ss_conley_species.rds"))

# 2. Species range maps
range_file <- paste0(wdir, 'in/species_ranges/merged_ranges.shp')
range_sf <- st_read(range_file)

# 8. Change in CWD
cmip_end <- load(paste0(wdir, 'in/CMIP5 CWD/cmip5_cwdaet_end.Rdat'))
cwd_cmip_end <- cwd_raster %>% mean()
rm(cwd_raster)
rm(aet_raster)

cmip_start <- load(paste0(wdir, 'in/CMIP5 CWD/cmip5_cwdaet_start.Rdat'))
cwd_cmip_start <- cwd_raster %>% mean()
rm(cwd_raster)
rm(aet_raster)
cwd_cmip_change <- cwd_cmip_end - cwd_cmip_start
cwd_cmip_pchange <- cwd_cmip_change / cwd_cmip_start

# cwdlist=list();aetlist=list()
# j=1
# for(i in c("start","end")){
#   load(paste0(wdir,"in/CMIP5 CWD/cmip5_cwdaet_",i,".Rdat"))
#   cwdlist[[j]]=cwd_raster;aetlist[[j]]=aet_raster
#   j=j+1
# }

# # 2. Species range maps
# range_file <- paste0(wdir, 'in//species_ranges//merged_ranges.shp')
# range_sf <- st_read(range_file) %>% 
#   filter(sp_code %in% (trim_df %>% pull(species_id) %>% unique()))
# 
# # 1. Historic climate raster
# clim_file <- paste0(wdir, 'in//CRUData//historic_raster//HistoricCWD_AETGrids_Annual.Rdat')
# load(clim_file)
# cwd_historic <- cwd_historic %>% mean(na.rm = TRUE)


# # 4. Site information
# site_smry <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
# site_smry <- site_smry %>% 
#   select(collection_id, sp_id, latitude, longitude) %>% 
#   mutate(species_id = tolower(sp_id)) %>% 
#   select(-sp_id)
# 
# # 5. Species information
# sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
# sp_info <- sp_info %>% 
#   select(species_id, genus, gymno_angio, family)
# site_smry <- site_smry %>% 
#   left_join(sp_info, by = c("species_id"))
# 
# # 5. Prediction rasters
# sp_predictions <- readRDS(paste0(wdir, "out/predictions/sq_sp_predictions.rds"))
# 
# 

# 2. Historic climate raster
clim_file <- paste0(wdir, 'in/CRUData/historic_raster/HistoricCWD_AETGrids.Rdat')
load(clim_file)
cwd_historic <- sum(cwd_historic)
names(cwd_historic) = "cwd"


# Average site conditions
ave_site_clim_df <- read_rds(paste0(wdir, "out/climate/site_ave_clim.gz"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define palettes ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_palette <- scale_colour_brewer(
  type = "seq",
  palette = 5,
  direction = 1,
  aesthetics = "colour"
)
















rwi_change <- sp_plot_dat %>% 
  filter(pet_hist > 1 & pet_hist < 1.5) %>% 
  ggplot(aes(x = cwd_hist, y = rwi_pred_change_mean)) +
  geom_hex(bins = 30) +
  scale_fill_viridis() +
  xlab("Mean CWD in 1970-2000\n(Deviation from species mean)") +
  ylab("RWI") +
  xlim(c(-2,3)) +
  labs(fill = "Count of\ngrid cells") +
  theme_bw(base_size = 12)
rwi_change

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Species-level changes in RWI  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_plot_dat <- sp_plot_dat %>% 
  mutate(rwi_change_bin = case_when(
    (rwi_pred_change_mean < -.5) ~ "1. RWI decline > 0.5",
    (rwi_pred_change_mean < -0.25) ~ "2. RWI decline between 0.25 and .5",
    (rwi_pred_change_mean < 0) ~ "3. RWI decline between 0 and 0.25",
    (rwi_pred_change_mean >= 0) ~ "4. RWI increase",
  ))
sp_plot_dat %>% 
  ggplot(aes(x = sp_code, fill = rwi_change_bin)) +
  geom_bar(position = "fill")


# transect_dat <- plot_dat %>% 
#   select(cwd.q, pet.q, rwi_change, rwi_change_psens, rwi_change_pclim) %>% 
#   pivot_longer(c(rwi_change, rwi_change_psens, rwi_change_pclim), names_to = 'scenario', values_to = "rwi_change") %>% 
#   mutate(scenario = fct_relevel(scenario, "rwi_change", "rwi_change_psens", "rwi_change_pclim"))
# 
# transect_1 <- transect_dat %>% 
#   filter(pet.q == 1) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-0.6, 0.2)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted")) +
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Constant shift in climate,\nvariable sensitivity",
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std above mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# 
# 
# transect_2 <- transect_dat %>% 
#   filter(pet.q == -1) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-0.4, 0.4)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted"))+
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Constant shift in climate,\nvariable sensitivity",
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std below mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#       legend.text = element_text(size=13),
#       legend.title = element_text(size=18),
#       legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# 
# 
# locator <- rwi_bin + 
#   theme_bw(base_size = 20)+
#   theme(legend.position = c(.18,.83),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank())+
#   geom_hline(yintercept = 1, size = 1) + 
#   geom_hline(yintercept = -1, size = 1)
# 
# locator | transect_1 / transect_2




# 
# lgd_pos <- c(.2, .8)
# cwd_sens_bin <- cwd_sens_bin +
#   theme(
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"))
# ggsave(paste0(wdir, "figures\\", "pres_cwd_sens.svg"), cwd_sens_bin, width = 9, height = 9)
# 
# pet_sens_bin <- pet_sens_bin +
#   theme(
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank())
# ggsave(paste0(wdir, "figures\\", "pres_pet_sens.svg"), pet_sens_bin, width = 9, height = 9)
# 

# cwd_change_bin <- cwd_change_bin +
#   theme(
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"))
# cwd_change_bin
# ggsave(paste0(wdir, "figures\\", "pres_cwd_change.svg"), cwd_change_bin, width = 9, height = 9)

# pet_change_bin <- pet_change_bin + 
#   theme(
#     legend.position = c(lgd_pos),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     plot.margin = margin(t=0, r=0, b=0, l=0, "cm"),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank())
# pet_change_bin
# ggsave(paste0(wdir, "figures\\", "pres_pet_change.svg"), pet_change_bin, width = 9, height = 9)
# 
# change_plot <- cwd_change_bin / pet_change_bin
# ggsave(paste0(wdir, "figures\\", "pred_full_b.svg"), change_plot, width = 9, height = 15)
# 
# 


# test +
#   scale_fill_continuous(type="viridis", name="Change in CWD\n1970-2000 to 2091-2100\n(mm per year)")
# 
# test %>% tmap_grob()
# 
# tm_shape(land) +
#   tm_raster(data = cwd_cmip_change, palette = terrain.colors(10)) +
# 
# map <-ggplot() +
#   geom_map(data = world, map = world,
#            aes(x=long, y=lat, map_id = region),
#            color = "black", fill = "lightgray", size = 0.1,alpha=.7) +
#   theme_bw(base_size =20) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_raster(data = cwd_cmip_change)
# 
# 
# 
#   geom_point(data = site_loc, aes(y = latitude, x = longitude), color = "#443A83FF")+
#   #geom_sf(data = flm_df, color = "#443A83FF", alpha = .2) +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   scale_x_continuous("Longitude", breaks = xr, labels = xlabels) +
#   scale_y_continuous("Latitude", breaks = yr, labels = ylabels)
# 
# map
# 
# 
# 
# 
# 
# # cwdmeanchange_end=mean(cwdlist[[3]])-mean(cwdlist[[1]])
# 
# 
# # ranges=st_read(paste0(wdir,"in/species_ranges/merged_ranges.shp"))
# r <- raster(nrows=nrow(cwd_cmip_change),ncols=ncol(cwd_cmip_change));extent(r)=extent(cwd_cmip_change)
# ranges_raster=rasterize(ranges_dissolve,r,field=1)
# test=rasterToPolygons(ranges_raster,dissolve = TRUE)
# testsf=st_as_sf(test)
# 
# plot(ranges)
# 
# world <- ne_download(scale = "medium", type="land",returnclass = "sf",category="physical")
# world=st_crop(world,extent(cwdmeanchange_end))
# world=st_transform(world,projection(cwdmeanchange_end))
# 
# library(reshape2)
# temp=as.matrix(cwd_cmip_change)
# colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# temp=melt(temp)
# colnames(temp)=c("lat","long","cwd")
# 
# area_grid=raster::area(cwdmeanchange_end)
# 
# #find global spatial average of cwd and pet changes - calculate areas of grid cells
# # cwdchanges=cellStats((mask(cwdlist[[3]],test)-mask(cwdlist[[1]],test))*area_grid,stat="sum")/cellStats(area_grid,stat="sum")
# # min=which.min(cwdchanges);max=which.max(cwdchanges)
# # mincwd=cwdlist[[3]][[min]]-cwdlist[[1]][[min]];maxcwd=cwdlist[[3]][[max]]-cwdlist[[1]][[max]]
# 
# world_raster=rasterize(world,r,field=1)
# world_raster=data.frame(temp[,1:2],melt(as.matrix(world_raster))[,3])
# colnames(world_raster)=c("lat","long","land")
# world_raster$land[which(world_raster$land==1)]="land"
# 
# a=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# a=a+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
# a=a+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# a=a+scale_fill_continuous(type="viridis",na.value="transparent",name="Change in CWD\n1970-2000 to 2091-2100\n(mm per year)")
# a=a+geom_raster()
# a=a+annotate("text",x=-150,y=-45,size=10,label="a)")
# # 
# # 
# # mincwd=raster::mask(mincwd,test)
# # temp=as.matrix(mincwd)
# # colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# # temp=melt(temp)
# # colnames(temp)=c("lat","long","cwd")
# # 
# # b=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# # b=b+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
# # b=b+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# # b=b+scale_fill_continuous(type="viridis",na.value="transparent",limits=c(-200,2700),guide=NULL)
# # b=b+annotate("text",x=-142,y=-45,size=8,label="Min")+geom_raster()
# # 
# # maxcwd=raster::mask(maxcwd,test)
# # temp=as.matrix(maxcwd)
# # colnames(temp)=xFromCol(cwdmeanchange_end);rownames(temp)=yFromRow(cwdmeanchange_end)
# # temp=melt(temp)
# # colnames(temp)=c("lat","long","cwd")
# # 
# # d=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# # d=d+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="lightgrey")
# # d=d+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# # d=d+scale_fill_continuous(type="viridis",na.value="transparent",limits=c(-200,2700),guide=NULL)
# # d=d+annotate("text",x=-142,y=-45,size=8,label="Max")+geom_raster()
# 
# #a/(b+d)
# 
# #get spatial averages in cwd and pet across species range for all three time points for all model runs
# 
# petlist=list()
# for(i in 1:3) petlist[[i]]=cwdlist[[i]]+aetlist[[i]]
# 
# cwdchange_dist=matrix(nrow=dim(cwdlist[[1]])[3],ncol=length(cwdlist))
# petchange_dist=matrix(nrow=dim(petlist[[1]])[3],ncol=length(petlist))
# 
# area_mask=mask(area_grid,test)
# 
# for(i in 1:length(cwdlist)){
#   cwd_temp=cellStats(cwdlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
#   pet_temp=cellStats(petlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
#   cwdchange_dist[,i]=cwd_temp;petchange_dist[,i]=pet_temp
#   print(i)
# }
# 
# colnames(cwdchange_dist)=c("Historic","Mid-Century","End-Century")
# colnames(petchange_dist)=c("Historic","Mid-Century","End-Century")
# cwdchange_dist=as.data.frame(cwdchange_dist);cwdchange_dist$Model=1:nrow(cwdchange_dist)
# petchange_dist=as.data.frame(petchange_dist);petchange_dist$Model=1:nrow(petchange_dist)
# 
# cwdchange_dist=cwdchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="CWD")
# petchange_dist=petchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="PET")
# 
# changedist=cbind(cwdchange_dist,petchange_dist[,3])
# changedist=changedist%>%pivot_longer(c("CWD","PET"),names_to="Variable",values_to="CWD_PET")
# changedist$Time=ordered(changedist$Time,c("Historic","Mid-Century","End-Century"))
# changedist$Variable=ordered(changedist$Variable,c("PET","CWD"))
# historic=changedist%>%group_by(Variable,Time)%>%dplyr::summarise(hist_mean=quantile(CWD_PET,0.5))
# changedist$CWD_PET[which(changedist$Time=="Historic")]=NA
# 
# e=ggplot(changedist,aes(x=Time,y=CWD_PET,col=Variable))+geom_boxplot(position="identity",width=0.15,inherit.aes=FALSE,aes(x=Time,y=CWD_PET,group=interaction(Variable,Time)),col="black",outlier.shape = NA,lwd=0.75)
# e=e+scale_y_continuous(limits=c(400,1500),name="CWD or PET (mm per yer)")
# e=e+geom_jitter(width=0.1)+theme_bw()
# e=e+geom_line(data=historic,aes(x=Time,y=hist_mean,group=Variable),col="black",lty=2)
# e=e+geom_point(data=historic%>%filter(Time=="Historic"),aes(y=hist_mean,col=Variable),x="Historic",size=3)
# e=e+labs(x="",color="")
# e=e+scale_color_manual(values=c("#972c25","#8e9169"))
# e=e+scale_x_discrete(labels=c("Historic"="Historic\n1970-2000","Mid-Century"="Mid-Century\n2045-2055","End-Century"="End-Century\n2090-2100"))
# e=e+annotate("text",x=0.7,y=470,size=10,label="b)")
# x11()
# a/e
# 
# ###Supplementary Figure Showing PET map with AET and PET changes
# 
# petmeanchange_end=mean(petlist[[3]])-mean(petlist[[1]])
# 
# petmeanchange_end=raster::mask(petmeanchange_end,test)
# temp=as.matrix(petmeanchange_end)
# colnames(temp)=xFromCol(petmeanchange_end);rownames(temp)=yFromRow(petmeanchange_end)
# temp=melt(temp)
# colnames(temp)=c("lat","long","cwd")
# 
# a=ggplot(temp,aes(x=long,y=lat,fill=cwd))
# a=a+geom_tile(data=world_raster%>%filter(!is.na(land)),inherit.aes=FALSE,aes(x=long,y=lat),fill="grey")
# a=a+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank(),panel.grid=element_blank())
# a=a+scale_fill_continuous(type="viridis",na.value="transparent",name="Change in PET\n1970-2000 to 2091-2100\n(mm per year)")
# a=a+geom_raster()
# a=a+annotate("text",x=-150,y=-45,size=10,label="a)")
# 
# #get model averages of AET changes
# aetchange_dist=matrix(nrow=dim(aetlist[[1]])[3],ncol=3)
# for(i in 1:length(cwdlist)){
#   aet_temp=cellStats(aetlist[[i]]*area_mask,stat="sum")/cellStats(area_mask,stat="sum")
#   aetchange_dist[,i]=aet_temp
#   print(i)
# }
# 
# colnames(aetchange_dist)=c("Historic","Mid-Century","End-Century")
# aetchange_dist=as.data.frame(aetchange_dist);aetchange_dist$Model=1:nrow(aetchange_dist)
# aetchange_dist=aetchange_dist%>%pivot_longer(c("Historic","Mid-Century","End-Century"),names_to="Time",values_to="AET")
# 
# changedist=cbind(petchange_dist,aetchange_dist[,3])
# changedist=changedist%>%pivot_longer(c("PET","AET"),names_to="Variable",values_to="PET_AET")
# changedist$Time=ordered(changedist$Time,c("Historic","Mid-Century","End-Century"))
# historic=changedist%>%group_by(Variable,Time)%>%dplyr::summarise(hist_mean=quantile(PET_AET,0.5))
# changedist$PET_AET[which(changedist$Time=="Historic")]=NA
# 
# e=ggplot(changedist,aes(x=Time,y=PET_AET,col=Variable))+geom_boxplot(position="identity",width=0.15,inherit.aes=FALSE,aes(x=Time,y=PET_AET,group=interaction(Variable,Time)),col="black",outlier.shape = NA,lwd=0.75)
# e=e+scale_y_continuous(limits=c(400,1500),name="PET or AET (mm per yer)")
# e=e+geom_jitter(width=0.1)+theme_bw()
# e=e+geom_line(data=historic,aes(x=Time,y=hist_mean,group=Variable),col="black",lty=2)
# e=e+geom_point(data=historic%>%filter(Time=="Historic"),aes(y=hist_mean,col=Variable),x="Historic",size=3)
# e=e+labs(x="",color="")
# e=e+scale_color_manual(values=c("#579473","#972c25"))
# e=e+scale_x_discrete(labels=c("Historic"="Historic\n1970-2000","Mid-Century"="Mid-Century\n2045-2055","End-Century"="End-Century\n2090-2100"))
# e=e+annotate("text",x=0.7,y=1400,size=10,label="b)")
# 
# x11()
# a/e
