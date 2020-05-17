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
# - add nobs
# - fix joins to prevent duplicate species_id
# - think through how to deal with CWD outliers
# - track down lost observations - currently dropping a lot due to NAN or failed RWI generation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote\\'

# 1. Historic species niche data
niche_csv <- paste0(wdir, 'clim_niche.csv')
niche_df <- read_csv(niche_csv) %>% 
  select(-X1) %>% 
  drop_na()

# 2. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'first_stage\\lm_cwd_pet.csv')) %>%
  select(-X1)

# Remove extreme outliers
flm_df <- flm_df %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T),
         pet.qhigh=quantile(estimate_pet.an,0.99,na.rm=T),
         pet.qlow=quantile(estimate_pet.an,0.01,na.rm=T)) %>%
  ungroup()
flm_df <- trim_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         estimate_pet.an>pet.qlow & estimate_pet.an<pet.qhigh)


# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) %>%
  select(-site)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Characterize site-level climate in relation to species range -------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  group_by(site_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum(petm),
            temp.an = mean(tmean),
            ppt.an = sum(ppt))


# Calculate site-level historic climate
hist_clim_df <- clim_df %>%
  group_by(site_id) %>%
  filter(year<1980) %>%
  summarise(aet.ave = mean(aet.an),
            cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an),
            temp.ave = mean(temp.an),
            ppt.ave = mean(ppt.an))


# Merge historic climate to flm data
flm_df <- flm_df %>% 
  inner_join(hist_clim_df, by = c("site_id"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species historic climate -------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niche_smry <- niche_df %>%
  group_by(sp_code) %>%
  summarize(sp_cwd_mean = mean(cwd),
            sp_cwd_sd = sd(cwd),
            sp_aet_mean = mean(aet),
            sp_aet_sd = sd(aet),
            sp_pet_mean = mean(pet),
            sp_pet_sd = sd(pet)) %>%
  ungroup()

# Merge chronology and species climate niche data
flm_df <- flm_df %>%
  inner_join(niche_smry, by = c("species_id" = "sp_code")) #note: currently losing a bunch of observations because we don't yet have range maps

# Calculate species-niche standardized climate
flm_df <- flm_df %>%
  mutate(aet.spstd = (aet.ave - sp_aet_mean) / sp_aet_sd,
         pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run second stage model --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, data = flm_df)
summary(mod)

mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd, data = flm_df)
summary(mod)


##### Run second stage model #####
# Define model
ss_mod <- function(d) {
  d <- d %>% mutate(errorweights = nobs / sum(nobs))
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
cwd.mod %>% summary()


#### Generate plots  ####
site_summary <- trim_df
## Summary plot of sample distribution
hex <- site_summary %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs)) +
  geom_hex()
hex

## Plot estimates by historic cwd / pet - would be awesome to combine these and make prettier
plot_dat <- trim_df
coef_plot1 <- 
  plot_dat %>% 
  ggplot(aes(x = cwd.spstd, y = pet.spstd, 
             z = estimate_pet.an, weight = nobs)) +
  stat_summary_hex() +
  xlim(c(-3, 4)) +
  #ylim(c(-1.5, 1.5))+
  scale_fill_gradientn (colours = c("darkblue","lightblue", "lightgrey")) +
  theme_bw()+
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme(legend.title = element_blank(),
        legend.position = "bottom")

coef_plot1 

coef_plot2 <- plot_dat %>% ggplot(aes(x = cwd.spstd, y = estimate_cwd.an, weight = nobs)) +
  stat_summary_bin(bins = 20)+
  theme_bw()+
  xlim(-3, 4)+
  ylab("Estimate CWD")+
  xlab("Deviation from mean CWD")


coef_plot2

coef_plot3 <- plot_dat %>% ggplot(aes(x = pet.spstd, y = estimate_cwd.an, weight = nobs)) +
  stat_summary_bin(bins = 20)+
  xlim(-3,4)+
  xlab("Deviation from mean PET")+
  ylab("Estimate CWD")+
  theme_bw()+
  rotate()

coef_plot3
coef_plot1/coef_plot2


ggarrange(coef_plot2, NULL,coef_plot1, coef_plot3, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = T)




## Plot prediction along one axis
cwd.x <- seq(-2, 2, length.out=1000)
pred <- predict(cwd.mod, data.frame(cwd.spstd=cwd.x, pet.spstd = 0),
                interval="prediction",level=0.90)
pred <- cbind(pred, cwd.x) %>%
  as_tibble()
plot(cwd.x, pred$fit, type="l", xlab="Historic CWD",ylab="Estimated effect of CWD",las=1,lwd=2,col="#0dcfca",
     ylim = c(pred$lwr %>% min(), pred$upr %>% max()))
polygon(c(cwd.x,rev(cwd.x)),c(pred$lwr,rev(pred$upr)),col=rgb(t(col2rgb("#0dcfca")),max=255,alpha=70),border=NA)


## Plot prediction along two axes
cwd.x <- seq(-2, 2, length.out=50)
pet.x <- seq(-2, 2, length.out=50)
grid.x <- expand.grid(cwd.x, pet.x)

pred <- predict(cwd.mod,data.frame(cwd.spstd=grid.x$Var1, pet.spstd = grid.x$Var2),
                interval="prediction",level=0.90)
pred <- cbind(pred,grid.x) %>% 
  as_tibble()

p <- pred %>% ggplot(aes(x = Var1, y = Var2, fill = fit)) +
  geom_tile()
p


## Boxplots of effects by cwd / pet bins
plot_dat <- trim_df %>%
  group_by("species_id") %>%
  mutate(cwd.q = as.factor(ntile(cwd.ave, 5)),
         pet.q = as.factor(ntile(pet.ave, 5)),
         cwd.high = as.factor(ntile(cwd.ave, 2)-1),
         pet.high = as.factor(ntile(pet.ave, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high)) %>%
  ungroup()

p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p2 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p3 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.003, 0.003)
p4 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.005, 0.005)

(p1 | p2) / (p3 | p4)
