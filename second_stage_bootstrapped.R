library(MASS)
library(tidyverse)
library(lfe)
library(broom)
library(purrr)
library(patchwork)
library(ggpubr)
library(ggiraphExtra)
library(hexbin)
select <- dplyr::select


### Define path
wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
#for Fran
wdir="C:/Users/fmoore/Google Drive/Treeconomics/Data/"
#for Joan
# wdir = "~/Google Drive/Treeconomics/Data/"
setwd(wdir)

# Read data
full_df <- read_csv(paste0(wdir,'first_stage.csv')) %>%
  select(-X1)

# Remove extreme outliers
trim_df <- full_df %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.975,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.025,na.rm=T),
         pet.qhigh=quantile(estimate_pet.an,0.975,na.rm=T),
         pet.qlow=quantile(estimate_pet.an,0.025,na.rm=T)) %>%
  ungroup()
trim_df <- trim_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         estimate_pet.an>pet.qlow & estimate_pet.an<pet.qhigh)


##### Run second stage model #####
# Define model
ss_mod <- function(d) {
  d <- d %>% mutate(errorweights = sqrt(nobs / sum(nobs)))
  mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=d)
  return(mod)
}

# Subset to species with sufficient plots to characterize climate niche
sp_count <- trim_df %>%
  group_by(species_id) %>%
  summarise(sites_per_sp = n_distinct(site_id))
keep_sp <- sp_count %>%
  filter(sites_per_sp > 40) %>%
  pull(species_id)

# Run model
mod_dat <- trim_df %>%
  filter(species_id %in% keep_sp)
cwd.mod <- mod_dat %>% ss_mod()
cwd.mod %>% summary(robust=TRUE)

##get correct standard errors through sampling distribution of first stage estimates 
nsamp=1000
sampdat=matrix(nrow=nsamp,ncol=2)

sampfunc=function(dat){
  meandat=as.numeric(dat[which(names(dat)=="estimate_cwd.an")]);stddat=as.numeric(dat[which(names(dat)=="std.error_cwd.an")])
  return(rnorm(1,mean=meandat,sd=stddat))
}

for(i in 1:nsamp){
  #draw from parameter distributions from first stage
  newdat=mod_dat
  newdat$estimate_cwd.an=apply(mod_dat,MARGIN=1,sampfunc)
  sampmod=ss_mod(newdat)
  sampdat[i,]=coef(sampmod)[2:3]
  if(i%%100==0) print(i)
}
colnames(sampdat)=c("cwd.spstd","pet.spstd")
sampdat=as.data.frame(sampdat)

p1=ggplot(sampdat)+geom_histogram(aes(x=cwd.spstd),fill="#32718e")+theme_bw()+labs(x="Coefficient Value",y="Frequency",title="Effect of Standardized CWD on CWD-Response")
p1=p1+geom_vline(xintercept = cwd.mod$coefficients[2],lwd=2)
p2=ggplot(sampdat)+geom_histogram(aes(x=pet.spstd),fill="#ed6716")+theme_bw()+labs(x="Coefficient Value",y="Frequency",title="Effect of Standardized PET on CWD-Response")
p2=p2+geom_vline(xintercept = cwd.mod$coefficients[3],lwd=2)
p1+p2

cwd.spstd.std.error=sd(sampdat[,1]);pet.spstd.std.error=sd(sampdat[,2])


