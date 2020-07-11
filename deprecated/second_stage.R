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
# wdir="C:/Users/fmoore/Google Drive/Treeconomics/Data/"
#for Joan
# wdir = "~/Google Drive/Treeconomics/Data/"
setwd(wdir)

# Read data
full_df <- read_csv(paste0(wdir,'first_stage.csv')) %>%
  select(-X1)

# Remove extreme outliers
trim_df <- full_df %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T),
         pet.qhigh=quantile(estimate_pet.an,0.99,na.rm=T),
         pet.qlow=quantile(estimate_pet.an,0.01,na.rm=T)) %>%
  ungroup()
trim_df <- trim_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         estimate_pet.an>pet.qlow & estimate_pet.an<pet.qhigh)


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
