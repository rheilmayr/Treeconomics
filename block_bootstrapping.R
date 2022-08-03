# library(tidyverse)
# library(sf)
# library(parallel)
# library(foreach)
# library(doParallel)

# 
# wdir<- 'remote\\'
# 
# site_df <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
# site_df <- site_df %>%
#   select(collection_id, sp_id, latitude, longitude)
# site_df <- site_df %>%
#   rename(species_id = sp_id) %>%
#   mutate(species_id = str_to_lower(species_id))
# 
# site_points=st_as_sf(site_df,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84
# +datum=WGS84 +no_defs +towgs84=0,0,0")
# 
# site_dist=st_distance(site_points)
# rownames(site_dist)=site_points$collection_id;colnames(site_dist)=site_points$collection_id
# save(site_dist,file=paste0(wdir,"out/site_distances.Rdat"))
load(paste0(wdir,"out/site_distances.Rdat"))

blockbootstrap_func=function(data,n_samps=n_mc,blockdist=363000,weights="cwd_errorweights",n_obs=n_sites,distmat=site_dist){
  #blockdist is a threshold distance (in m) so that for each site sampled all other sites within that radius are also included in the sample
  #default value is 50km
  site_list=data%>%
    group_by(collection_id)%>%
    dplyr::summarise(weights=unique(!!as.name(weights)))
  
  #bootstrap samples
  bootstrapsamp=foreach(i=1:n_samps,.combine="rbind",.packages = c("tibble"))%dopar%{
    samp=sample(1:n_obs,size=n_obs,prob=site_list$weights,replace=TRUE)
    temp=data.frame()
    j=1
    while(nrow(temp)<n_obs){
      #find sample site and other sites within the block (i.e. within the blockdist threshold)
      target=site_list[samp[j],1]
      block=colnames(distmat)[(which(as.numeric(distmat[which(rownames(distmat)==as.character(target)),])<blockdist))]
      iter_temp=sample(1:n_samps,size=1) #randomly sample from the firt stage draws
      temp=rbind(temp,data[which(data$collection_id%in%block&data$iter_idx==iter_temp),]) 
      j=j+1
    }
    tibble(iter_idx=i,data=temp)
  }
  bootstrapsamp=bootstrapsamp%>%
    group_by(iter_idx)%>%
    nest()
  return(bootstrapsamp)
  
}
