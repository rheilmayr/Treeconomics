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
library(patchwork)
library(hexbin)
library(dbplyr)
library(RSQLite)
library(modi)
library(margins)
library(fixest)


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


# 2. Site-specific weather history
cwd_csv <- paste0(wdir, 'essentialcwd_data.csv')
cwd_df <- read.csv(cwd_csv, sep=',')
cwd_df <- cwd_df %>% 
  mutate("site_id" = as.character(site)) %>%
  select(-site)


# 3. Site-level regressions
# flm_df <- read_csv(paste0(wdir, 'first_stage\\tree_log_pet_cwd.csv')) %>%
#   select(-X1)

flm_df <- read_csv(paste0(wdir, 'first_stage\\log_cwd_pet.csv')) %>%
  select(-X1)


# Remove extreme outliers
flm_df <- flm_df %>%
  group_by(species_id) %>%
  mutate(cwd.qhigh=quantile(estimate_cwd.an,0.99,na.rm=T),
         cwd.qlow=quantile(estimate_cwd.an,0.01,na.rm=T)) %>%
  ungroup()
flm_df <- flm_df %>%
  filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh)


# Connect to database
tree_db = paste0(wdir, 'tree_ring_data_V2.db')
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, tree_db)
tree_db = as.data.frame(tbl(conn, 'trees'))
sp_site_index <- tree_db %>% 
  select(species_id, site_id) %>% 
  distinct()
dbDisconnect(conn)
n_sp_sites <- sp_site_index %>% 
  group_by(species_id) %>% 
  tally()
freq_species <- n_sp_sites %>% 
  filter(n > 10) %>% 
  pull(species_id)

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
            ppt.an = sum(ppt)) %>% 
  drop_na()


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


# Calculate species niche based on ITRDB sites
clim_df <- clim_df %>% 
  left_join(sp_site_index, by = c("site_id"))

sp_clim_df <- clim_df %>% 
  group_by(species_id) %>% 
  filter(year<1980) %>% 
  summarise(pet.sp.ave = mean(pet.an),
            cwd.sp.ave = mean(cwd.an),
            pet.sp.sd = sd(pet.an),
            cwd.sp.sd = sd(cwd.an))

flm_df <- flm_df %>% 
  left_join(sp_clim_df, by = c("species_id"))


flm_df <- flm_df %>% 
  mutate(pet.spstd = (pet.ave - pet.sp.ave) / pet.sp.sd,
         cwd.spstd = (cwd.ave - cwd.sp.ave) / cwd.sp.sd)

# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Summarize species historic climate -------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# niche_smry <- niche_df %>%
#   group_by(sp_code) %>%
#   summarize(sp_cwd_mean = mean(cwd),
#             sp_cwd_sd = sd(cwd),
#             sp_aet_mean = mean(aet),
#             sp_aet_sd = sd(aet),
#             sp_pet_mean = mean(pet),
#             sp_pet_sd = sd(pet)) %>%
#   ungroup()
# 
# # Merge chronology and species climate niche data
# flm_df <- flm_df %>%
#   inner_join(niche_smry, by = c("species_id" = "sp_code")) #note: currently losing a bunch of observations because we don't yet have range maps
# 
# # Calculate species-niche standardized climate
# flm_df <- flm_df %>%
#   mutate(aet.spstd = (aet.ave - sp_aet_mean) / sp_aet_sd,
#          pet.spstd = (pet.ave - sp_pet_mean) / sp_pet_sd,
#          cwd.spstd = (cwd.ave - sp_cwd_mean) / sp_cwd_sd)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run second stage model --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flm_df <- flm_df %>% 
  mutate(errorweights = 1 / (std.error_cwd.an^2))

trim_df <- flm_df %>% 
  filter(species_id %in% freq_species) %>% 
  drop_na()

mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd + cwd.spstd:pet.spstd, weights=errorweights, data=flm_df)
summary(mod)

saveRDS(mod, paste0(wdir, "second_stage\\ss_mod.rds"))


trim_df <- flm_df %>% 
  filter(species_id == "psme") %>% 
  drop_na()

mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=trim_df)
summary(mod)

saveRDS(mod, paste0(wdir, "second_stage\\psme_mod.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generate plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_summary <- flm_df
## Summary plot of sample distribution
hex <- site_summary %>% ggplot(aes(x = cwd.spstd, y = pet.spstd, weight = nobs)) +
  geom_hex() +
  xlim(-2.5, 2.5) +
  ylim(-2.5, 2.5) +
  labs(fill = "Number of tree-year\nobservations") +
  ylab("Historic PET\n(Deviation from species mean)") +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  theme_bw(base_size = 22) +
  coord_fixed()
hex
ggsave(paste0(wdir, 'figures\\clim_density.svg'), plot = hex)



## Binned plot of cwd sensitivity
plot_dat <- trim_df %>%
  filter(abs(cwd.spstd)<3,abs(pet.spstd<3)) %>% 
  drop_na()
nbins = 8
plot_dat <- plot_dat %>% 
  mutate(cwd.q = as.numeric(ntile(cwd.spstd, nbins)),
         pet.q = as.numeric(ntile(pet.spstd, nbins)))

cwd.quantiles = quantile(plot_dat$cwd.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)
pet.quantiles = quantile(plot_dat$pet.spstd, probs = seq(0, 1, 1/nbins), names = TRUE) %>% round(2)

group_dat <- plot_dat %>% 
  group_by(cwd.q, pet.q) %>% 
  summarize(wvar = weighted.var(estimate_cwd.an, errorweights, na.rm = TRUE),
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
ggsave(paste0(wdir, 'figures\\binned_margins.svg'), plot = binned_margins)


## Modeled margins plots
cdat <- cplot(mod, 'cwd.spstd', what = "prediction", draw = FALSE)
cdat <- cdat %>%
  filter(abs(xvals)<2.5)
margins_plot <- ggplot(cdat, aes(x = xvals)) + 
  geom_line(aes(y = yvals)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
  geom_line(aes(y = upper), linetype = 3) +
  geom_line(aes(y = lower), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") + 
  theme_bw(base_size = 22) +
  scale_y_continuous(labels = scales::scientific)
# +
#   geom_point(data = plot_dat, aes(x = cwd.spstd, y = estimate_cwd.an)) +
#   ylim(c(-1.5e-4, 0)) +
#   xlim(c(-2.5, 2.5))




margins_plot
ggsave(paste0(wdir, 'figures\\cwd_margins.svg'), plot = margins_plot)


histogram <- ggplot(plot_dat, aes(x = cwd.spstd)) + 
  geom_histogram(bins = 40, alpha=0.2, fill = "darkblue") +
  xlim(c(-2.5, 2.5)) +
  theme_bw(base_size = 22) + 
  ylab("") +
  scale_x_continuous(labels = c(""), breaks = c(0)) +
  theme(aspect.ratio = 0.3,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


histogram / margins_plot
  
  
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

