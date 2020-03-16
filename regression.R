# TODO:
#   Bootstrap first stage
#   Look at temperature and latitude in second stage regression (pure cwd in first stage)
#   add lagged climate data maybe up to 2 years
#   Do you get snow?


library(MASS)
library(tidyverse)
library(lfe)
library(broom)
library(purrr)
library(patchwork)
library(ggiraphExtra)
library(ggplot2)
select <- dplyr::select


### Define path
#wdir = 'D:/cloud/Google Drive/Treeconomics/Data/'
#for Fran
# wdir="C:/Users/fmoor/Google Drive/Treeconomics/Data/"
#for Joan
wdir = "~/Google Drive/Treeconomics/Data/"
setwd(wdir)

### Plotting parameters
pal=colorRampPalette(c('white','blue','yellow','red','darkred'))

### Read and prep data  ####
# spp=c("psme","pipo","pisy","pcgl","pcab","pied")
# speciesnames=c("Douglas Fir","Ponderosa Pine","Scotch Pine",
#                "White Spruce","Norway Spruce","Colorado Pinyon")

df <- read_csv(paste0(wdir,'species_data.csv'))
df <- df %>%
  mutate(tree_id = paste0(site_id, tree_id))
df <- df %>%
  mutate(pet.an = cwd.an + aet.an)

sp_details <- read.csv(paste0(wdir, 'species_list_jd.csv'))

### Calculate site historic average aet, cwd and pet  - NOTE - IS THIS WEIGHTING BY NUMBER OF OBS?
site_clim <- df %>%
  filter(year < 1980) %>%                                                                      ## Switch to raw climate data rather than using tree dataframe
  group_by(site_id, year) %>%
  summarize(cwd.an = max(cwd.an, na.rm = TRUE), 
            aet.an = max(aet.an, na.rm = TRUE), 
            pet.an = max(pet.an, na.rm = TRUE)) %>%
  filter((cwd.an > -Inf) & (aet.an > -Inf)) %>%
  group_by(site_id) %>%
  summarize(cwd.ave = mean(cwd.an, na.rm = TRUE), 
            aet.ave = mean(aet.an, na.rm = TRUE), 
            pet.ave = mean(pet.an, na.rm = TRUE))

df <- merge(x = df, y = site_clim, by = "site_id", all.x = TRUE)


#remove sites with negative ages
df_trim <- df %>%
  filter(age>0)

#remove trees with a very short record (<5 years)
treeobs=df_trim %>%
  group_by(tree_id) %>%
  summarize(nyears=n())
df_trim=left_join(df_trim,treeobs, by = 'tree_id')
df_trim = df_trim %>%
  filter(nyears>5)

#remove sites with few trees
ntrees=df_trim%>%
  group_by(site_id,species_id)%>%
  summarize(ntrees=length(unique(tree_id)))
df_trim=left_join(df_trim,ntrees, by = c("site_id", "species_id"))
df_trim=df_trim%>%
  filter(ntrees>5)

complete_df <- df_trim %>%
  drop_na(c("cwd.an","aet.an","pet.an", "ring_width"))


site_summary <- complete_df %>% 
  select("site_id", "species_id", "species_name", "cwd.ave", "pet.ave", "aet.ave") %>%
  distinct() %>%
  as_tibble()

nobs <- complete_df %>%
  group_by(site_id, species_id) %>%
  summarise(nobs = n())
            # ntrees = n_distinct(tree_id))

site_summary <- site_summary %>%
  merge(nobs, by = c('site_id', 'species_id'), all.x = TRUE) %>%
  merge(sp_details, by = 'species_id', all.x = TRUE)

#### Calculate standardized climate variables for each species ####
site_summary <- site_summary %>%
  group_by(species_id) %>%
  mutate(cwd.spstd = (cwd.ave - mean(cwd.ave)) / sd(cwd.ave),
         aet.spstd = (aet.ave - mean(aet.ave)) / sd(aet.ave),
         pet.spstd = (pet.ave - mean(pet.ave)) / sd(pet.ave))

to_merge <- site_summary %>%
  select("site_id", "species_id", "cwd.spstd", "aet.spstd", "pet.spstd")
complete_df <- complete_df %>%
  merge(y = to_merge, by = c("site_id", "species_id"), all.x = TRUE)

#### Run first stage ####
fit_mod <- function(d) {
  felm(ring_width ~ cwd.an + pet.an + age+I(age^2)+I(age^3)+I(age^4)|tree_id|0|0, data = d )
}

site_lm <- complete_df %>% 
  group_by(site_id, species_id) %>%
  nest() %>%
  mutate(mod = map(data, fit_mod),
         mod = map(mod, tidy)) %>%
  unnest(mod) %>%
  filter(term %in% c('cwd.an', 'pet.an'))

siteCoef <- site_lm %>%
  pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value")) %>%
  left_join(site_summary, by = c("site_id","species_id")) %>%
  select(-data)

#remove one super extreme outlier
siteCoef_trimmed <- siteCoef %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T),
         pet.qhigh=quantile(estimate_pet.an,0.99,na.rm=T),
         pet.qlow=quantile(estimate_pet.an,0.01,na.rm=T)) %>%
  ungroup()
siteCoef_trimmed=siteCoef_trimmed %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         estimate_pet.an>pet.qlow & estimate_pet.an<pet.qhigh)

# Add weighting variable
siteCoef_trimmed <- siteCoef_trimmed %>%
  mutate(errorweights = nobs / sum(nobs)) 


# Calculate limits
xlim.all <- site_summary %>%
  ungroup() %>%
  summarize(cwd.qhigh = quantile(cwd.ave,0.95),
            cwd.qlow = quantile(cwd.ave,0.05),
            cwd.qmed = quantile(cwd.ave,0.5),
            pet.qhigh = quantile(pet.ave, 0.95), 
            pet.qlow = quantile(pet.ave, 0.05),
            pet.qmed = quantile(pet.ave,0.5),
            aet.qhigh = quantile(aet.ave, 0.95), 
            aet.qlow = quantile(aet.ave, 0.05),
            aet.qmed = quantile(aet.ave,0.5))


#### Run second stage of model  ####
cwd.mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd + factor(species_id), weights=errorweights, data= siteCoef_trimmed)
cwd.mod %>% summary()

pet.mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd + factor(species_id), weights=errorweights, data= siteCoef_trimmed)
pet.mod %>% summary()

## Run second stage by subgroup
ss_mod <- function(d) {
  d <- d %>% mutate(errorweights = nobs / sum(nobs)) 
  lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=d)
}

ss_lm <- siteCoef_trimmed %>% 
  group_by(biome) %>%    ### NOTE: can change this to run model by family, genus, biome, conifer_broad, gymno_angio
  nest() %>%
  mutate(mod = map(data, ss_mod),
         mod = map(mod, tidy)) %>%
  unnest(mod) %>%
  filter(term %in% c('(Intercept)', 'cwd.spstd'))


#### Generate plots  ####
## Summary plot of sample distribution
hex <- site_summary %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs)) +
  geom_hex()
hex

## Plot estimates by historic cwd / pet - would be awesome to combine these and make prettier
plot_dat <- siteCoef_trimmed
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

library(ggpubr)

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
plot_dat <- siteCoef_trimmed %>%
  group_by("species_id") %>%
  mutate(cwd.q = as.factor(ntile(cwd.ave, 5)),
         pet.q = as.factor(ntile(pet.ave, 5)),
         cwd.high = as.factor(ntile(cwd.ave, 2)-1),
         pet.high = as.factor(ntile(pet.ave, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high)) %>%
  ungroup()

p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot()
p2 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_pet.an)) +
  geom_boxplot()
p3 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_cwd.an)) +
  geom_boxplot()
p4 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_pet.an)) +
  geom_boxplot()

(p1 | p2) / (p3 | p4)
