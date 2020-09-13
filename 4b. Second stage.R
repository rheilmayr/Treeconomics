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
# - Mechanisms: Think through this model in more detail
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(MASS)
library(tidyverse)
library(broom)
library(purrr)
library(patchwork)
library(patchwork)
library(dbplyr)
library(RSQLite)
library(lmtest)
library(sandwich)
library(prediction)
library(Hmisc)
library(hexbin)

library(ggpubr)
library(ggiraphExtra)
library(modi)
library(margins)


select <- dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote\\'

# 1. Historic species niche data
niche_csv <- paste0(wdir, 'out/climate/clim_niche.csv')
niche_df <- read_csv(niche_csv) %>% 
  select(-X1)

# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'out/climate/essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) %>%
  select(-site)

# 3. Tree-level regressions
flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\site_log_pet_cwd.csv')) %>%
  select(-X1)

# 4. Site information
site_df <- read_csv(paste0(wdir, 'out\\dendro\\site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id)
flm_df <- flm_df %>% 
  left_join(site_df, by = "collection_id")
flm_df <- flm_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# 5. Species information
sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
flm_df <- flm_df %>% 
  left_join(sp_info, by = "species_id")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Characterize site-level climate in relation to species range -------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate site-level annual climate
clim_df = cwd_df %>%
  group_by(site_id, year) %>%
  summarise(aet.an = sum(aet),
            cwd.an = sum(cwd),
            pet.an = sum(aet+cwd)) %>% 
  drop_na()

# Calculate site-level historic climate
hist_clim_df <- clim_df %>%
  group_by(site_id) %>%
  filter(year<1980) %>%
  summarise(aet.ave = mean(aet.an),
            cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an))

# Merge historic climate to flm data
flm_df <- flm_df %>% 
  inner_join(hist_clim_df, by = c("collection_id" = "site_id"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species historic climate -------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niche_df <- niche_df %>% 
  select(sp_code, sp_pet_mean = pet_mean, sp_pet_sd = pet_sd, sp_cwd_mean = cwd_mean, sp_cwd_sd = cwd_sd)

# Merge chronology and species climate niche data
flm_df <- flm_df %>%
  inner_join(niche_df, by = c("species_id" = "sp_code")) #note: currently losing ~7k trees (10%) because we don't have range maps

# Calculate species-niche standardized climate
flm_df <- flm_df %>%
  mutate(pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
         cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)

# Calculate tree-species-niche standardized climate
flm_df <- flm_df %>%
  mutate(pet.trstd = (pet.trmean - sp_pet_mean) / sp_pet_sd,
         cwd.trstd = (cwd.trmean - sp_cwd_mean) / sp_cwd_sd)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Trim data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add weighting based on inverse of first stage variance
flm_df <- flm_df %>% 
  mutate(errorweights = 1 / (std.error_cwd.an^2),
         pet_errorweights = 1 / (std.error_pet.an^2))

# Identify and trim extreme outliers
flm_df <- flm_df %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T)) %>%
  ungroup()

trim_df <- flm_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         abs(cwd.spstd)<5,
         abs(pet.spstd)<5) %>% 
  drop_na()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin <- -2.5
xmax <- 2.5

### Summary plot of sample distribution
hex <- flm_df %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs)) +
  geom_hex() +
  xlim(xmin, xmax) +
  ylim(xmin, xmax) +
  labs(fill = "Number of tree-year\nobservations") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  theme_bw(base_size = 22) +
  coord_fixed()
hex
# ggsave(paste0(wdir, 'figures\\clim_density.svg'), plot = hex)

### Binned plot of cwd sensitivity
plot_dat <- trim_df %>%
  filter(((abs(cwd.spstd)<3) & (abs(pet.spstd<3)))) %>%
  drop_na()
nbins = 8
plot_dat <- plot_dat %>% 
  mutate(cwd.q = as.numeric(ntile(cwd.spstd, nbins)),
         pet.q = as.numeric(ntile(pet.spstd, nbins)))

cwd.quantiles = quantile(plot_dat$cwd.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)
pet.quantiles = quantile(plot_dat$pet.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)

group_dat <- plot_dat %>% 
  group_by(cwd.q, pet.q) %>% 
  dplyr::summarize(wvar = wtd.var(estimate_cwd.an, errorweights, na.rm = TRUE),
            wsd = sqrt(wvar),
            wmean = weighted.mean(estimate_cwd.an, errorweights, na.rm = TRUE),
            n = n(),
            error = qt(0.975, df = n-1)*wsd/sqrt(n),
            lower = wmean - error,
            upper = wmean + error) %>% 
  filter(n>4)


binned_margins <- group_dat %>% 
  ggplot(aes(x = cwd.q, y = pet.q, fill = wmean)) +
  geom_tile() +
  # xlim(c(-3, 4)) +
  #ylim(c(-1.5, 1.5))+
  scale_fill_gradientn (colours = c("darkblue","lightblue")) +
  # scale_fill_viridis_c() +
  theme_bw(base_size = 22)+
  ylab("Deviation from mean PET")+
  xlab("Deviation from mean CWD")+
  theme(legend.position = "right") +
  labs(fill = "Marginal effect\nof CWD") +
  scale_x_continuous(labels = cwd.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
  scale_y_continuous(labels = pet.quantiles, breaks = seq(0.5, nbins+0.5, 1)) +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  coord_fixed()

binned_margins 
# ggsave(paste0(wdir, 'figures\\binned_margins.svg'), plot = binned_margins)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Primary second stage model --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run regression and cluster s.e.

# trim_df <- trim_df %>% 
#   filter(gymno_angio=="gymno")

mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=trim_df)
cluster_vcov <- vcovCL(mod, cluster = trim_df$collection_id)
coeftest(mod, cluster = trim_df$collection_id)
saveRDS(mod, paste0(wdir, "out\\second_stage\\ss_mod.rds"))



pet_mod <- lm(estimate_pet.an ~ cwd.spstd + pet.spstd, weights = pet_errorweights, data=flm_df)
cluster_vcov <- vcovCL(pet_mod, cluster = flm_df$collection_id)
coeftest(pet_mod, cluster = flm_df$collection_id)
saveRDS(pet_mod, paste0(wdir, "out\\second_stage\\pet_ss_mod.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot marginal effects ---------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predictions <- prediction(mod, at = list(cwd.spstd = seq(xmin, xmax, .1)), vcov = cluster_vcov, calculate_se = T) %>% 
  summary() %>% 
  rename(cwd.spstd = "at(cwd.spstd)")

margins_plot <- ggplot(predictions, aes(x = cwd.spstd)) + 
  geom_line(aes(y = Prediction)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
  geom_line(aes(y = upper), linetype = 3) +
  geom_line(aes(y = lower), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") + 
  theme_bw(base_size = 22) +
  scale_y_continuous(labels = scales::scientific)

# margins_plot

histogram <- ggplot(plot_dat, aes(x = cwd.spstd)) + 
  geom_histogram(bins = 40, alpha=0.2, fill = "darkblue") +
  xlim(c(xmin, xmax)) +
  theme_bw(base_size = 22) + 
  ylab("# trees") +
  theme(aspect.ratio = 0.3,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

histogram / margins_plot
=ggsave(paste0(wdir, 'figures\\cwd_margins.svg'), plot = margins_plot)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Robustness tests --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Squared terms
sq_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd + I(cwd.spstd^2) + I(pet.spstd^2), weights = errorweights, data=trim_df)
coeftest(sq_mod, vcov = vcovCL, cluster = trim_df$collection_id)
sq_cluster_vcov <- vcovCL(sq_mod, cluster = trim_df$collection_id)

# Drop PET
nopet_mod <- lm(estimate_cwd.an ~ cwd.spstd, weights = errorweights, data=trim_df)
coeftest(nopet_mod, vcov = vcovCL, cluster = trim_df$collection_id)
nopet_cluster_vcov <- vcovCL(nopet_mod, cluster = trim_df$collection_id)

# Don't trim dataset
allobs_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=flm_df)
coeftest(allobs_mod, vcov = vcovCL, cluster = flm_df$collection_id)
allobs_cluster_vcov <- vcovCL(allobs_mod, cluster = flm_df$collection_id)

# Limit to gymnosperms
subset_df <- trim_df %>% filter(gymno_angio=="gymno")
conifer_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=subset_df)
coeftest(conifer_mod, vcov = vcovCL, cluster = subset_df$collection_id)
conifer_cluster_vcov <- vcovCL(conifer_mod, cluster = subset_df$collection_id)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate mechanism - tree or site --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Focus on trees from genera with effect, and for which we have full weather history
trim_df <- flm_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         abs(cwd.spstd)<5,
         abs(pet.spstd)<5,
         genus %in% c("ju", "pi", "ps"),
         young==TRUE) %>% 
  drop_na()

mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trstd + pet.spstd + pet.trstd, weights = errorweights, data=trim_df)
coeftest(mod, vcov = vcovCL, cluster = trim_df$collection_id)

# Could include all trees from sites with any variation in tree-level drought history (including trees born <1901) 
young_sites <- flm_df %>% 
  filter(young == T) %>% 
  pull(collection_id) %>% 
  unique()
trim_df <- flm_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
         abs(cwd.spstd)<5,
         abs(pet.spstd)<5,
         genus %in% c("ju", "pi", "ps"),
         collection_id %in% young_sites) %>% 
  drop_na()
mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trstd + pet.spstd + pet.trstd, weights = errorweights, data=trim_df)
coeftest(mod, vcov = vcovCL, cluster = trim_df$collection_id)

# Include collection fixed effects
mod <- feols(estimate_cwd.an ~ cwd.trstd + pet.trstd | collection_id, weights = trim_df$errorweights, data=trim_df)
summary(mod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate variation by genus  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predict_cwd_effect <- function(trim_df){
  trim_df <- trim_df %>%
    filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
           abs(cwd.spstd)<5,
           abs(pet.spstd)<5) %>%
    drop_na()
  gen_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=trim_df)
  cluster_vcov <- vcovCL(gen_mod, cluster = trim_df$collection_id)
  print(coeftest(gen_mod, vcov = vcovCL, cluster = trim_df$collection_id))
  predictions <- prediction(gen_mod, at = list(cwd.spstd = seq(-3, 3, .5)), vcov = cluster_vcov, calculate_se = T) %>% 
    summary() %>% 
    rename(cwd.spstd = "at(cwd.spstd)") 
  return(predictions)
}

## Plot marginal effects by gymno / angio
grp_freq <- flm_df %>% 
  group_by(gymno_angio) %>% 
  summarise(n_collections = n_distinct(collection_id)) %>% 
  arrange(desc(n_collections))

grp_keep <- grp_freq %>% 
  filter(n_collections>50) %>% 
  pull(gymno_angio)

grp_df <- flm_df %>%
  filter(gymno_angio %in% grp_keep) %>% 
  group_by(gymno_angio) %>% 
  nest() %>%
  mutate(prediction = map(data, predict_cwd_effect))

grp_df <- grp_df %>% 
  unnest(prediction) %>% 
  select(-data)

grp_df$gymno_angio <- grp_df$gymno_angio %>% recode("angio" = "Angiosperm", "gymno" = "Gymnosperm")

margins_plot <- ggplot(grp_df, aes(x = cwd.spstd)) + 
  geom_line(aes(y = Prediction)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
  facet_wrap(~gymno_angio, scales = "free", ncol = 1) +
  geom_line(aes(y = upper), linetype = 3) +
  geom_line(aes(y = lower), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") + 
  theme_bw(base_size = 22) +
  scale_y_continuous(labels = scales::comma)
margins_plot
ggsave(paste0(wdir, 'figures\\gymno_angio.svg'), plot = margins_plot)


## Plot marginal effects by genus
genus_freq <- flm_df %>% 
  group_by(genus) %>% 
  summarise(n_collections = n_distinct(collection_id)) %>% 
  arrange(desc(n_collections))

genus_keep <- genus_freq %>% 
  filter(n_collections>50) %>% 
  pull(genus)

genus_df <- flm_df %>% 
  filter(genus %in% genus_keep) %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(prediction = map(data, predict_cwd_effect))

genus_key <- sp_info %>% 
  select(genus, gymno_angio) %>% 
  unique()

genus_df <- genus_df %>% 
  unnest(prediction) %>% 
  select(-data) %>% 
  left_join(genus_key, by = "genus")

margins_plot <- ggplot(genus_df, aes(x = cwd.spstd, fill = gymno_angio)) + 
  geom_line(aes(y = Prediction)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
  facet_wrap(~genus, scales = "free") +
  geom_line(aes(y = upper), linetype = 3) +
  geom_line(aes(y = lower), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") + 
  theme_bw(base_size = 22) +
  scale_y_continuous(labels = scales::scientific)
margins_plot









#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Deprecated code below this  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
margins(mod, at  = list(pet.spstd = -2:2))
cplot(mod, 'pet.spstd', what = "prediction")
persp(mod, 'cwd.spstd', 'pet.spstd')
img <- image(mod, 'cwd.spstd', 'pet.spstd')




group_dat <- plot_dat %>% 
  group_by(cwd.q) %>% 
  summarize(wvar = weighted.var(estimate_cwd.an, errorweights),
            wsd = sqrt(wvar),
            wmean = weighted.mean(estimate_cwd.an, errorweights),
            n = n(),
            error = qt(0.975, df = n-1)*wsd/sqrt(n),
            lower = wmean - error,
            upper = wmean + error)


ggplot(group_dat, aes(x=cwd.q, y=wmean)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) 


pet_dat <- plot_dat %>% 
  group_by(pet.q) %>% 
  summarize(wvar.cwd = weighted.var(estimate_cwd.an, errorweights),
            wsd.cwd = sqrt(wvar),
            wmean = weighted.mean(estimate_cwd.an, errorweights),
            n = n(),
            error = qt(0.975, df = n-1)*wsd/sqrt(n),
            lower = wmean - error,
            upper = wmean + error)

ggplot(group_dat, aes(x=pet.q, y=wmean)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) 


p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p1



## Plot estimates by historic cwd / pet - would be awesome to combine these and make prettier
plot_dat <- flm_df
coef_plot1 <- 
  plot_dat %>% 
  ggplot(aes(x = cwd.spstd, y = pet.spstd, 
             z = estimate_pet.an, weight = errorweights)) +
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

coef_plot2 <- plot_dat %>% ggplot(aes(x = cwd.spstd, y = estimate_cwd.an, weight = errorweights)) +
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



## Bootstrap distributions from first stage
N <- 100

draw_est_smpl <- function(dat){
  var1 <- as.numeric(dat[which(names(dat)=="std.error_cwd.an")])^2
  var2 <- as.numeric(dat[which(names(dat)=="std.error_pet.an")])^2
  cov <- as.numeric(dat[which(names(dat)=="pet_cwd_cov")])^2
  vcov <- matrix(c(var1, cov, cov, var2), 2)
  mu <- c(data$estimate_cwd.an, data$estimate_pet.an)
  smpl <- mvrnorm(N, mu = mu, Sigma = vcov )
  colnames(smpl)=c("cwd.est.smpl","pet.est.smpl")
  smpl <- smpl %>% as_tibble()
  return(smpl)
}

draw_df <- flm_df %>% 
  rowwise %>% 
  do(draw = draw_est_smpl(.))
draw_df <- cbind(flm_df, draw_df) %>% 
  unnest(draw)
  

## Boxplots of effects by cwd / pet bins


sig_df <- flm_df %>% 
  filter(p.value_cwd.an<0.05)

plot_dat <- flm_df %>%
  # filter(genus %in% c("pi", "ps", "ju")) %>% 
  group_by("species_id") %>%
  mutate(cwd.q = as.factor(ntile(cwd.ave, 5)),
         pet.q = as.factor(ntile(pet.ave, 5)),
         cwd.high = as.factor(ntile(cwd.ave, 2)-1),
         pet.high = as.factor(ntile(pet.ave, 2)-1),
         climate = paste0('cwd', cwd.high, '_pet', pet.high)) %>%
  ungroup()



plot_dat %>% 
  mutate(inv_std = 1/ std.error_cwd.an) %>% 
  group_by(cwd.q) %>% 
  summarize(wmean = weighted.mean(estimate_cwd.an, inv_std, na.rm = TRUE))

plot_dat %>% 
  group_by(cwd.q) %>% 
  summarize(mean = mean(estimate_cwd.an, na.rm = TRUE))


p1 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p2 <- plot_dat %>% ggplot(aes(x=cwd.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p3 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_cwd.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)
p4 <- plot_dat %>% ggplot(aes(x=pet.q, y=estimate_pet.an)) +
  geom_boxplot() +
  ylim(-0.0025, 0.0025)

(p1 | p2) / (p3 | p4)

