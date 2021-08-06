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
library(tidylog)
library(fixest)


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

# 3a. Site-level regressions
flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_pet_cwd_std.csv')) %>%
  select(-X1)
flm_df <- flm_df %>%
  mutate(collection_id = tree_id) %>% 
  separate(collection_id, c("collection_id", "tree_n"), "_")


# # 3b. Tree-level regressions
# tree_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd.csv')) %>%
#   select(-X1)

# ggplot(flm_df %>% filter(abs(estimate_cwd.an)<0.01, cwd.min>5), aes(x = estimate_cwd.an)) + 
#   geom_histogram(bins = 100, alpha=0.4)
# dim(flm_df %>% filter(abs(estimate_cwd.an)<100))[1] / dim(flm_df)[1] 
# dim(tree_df %>% filter(abs(estimate_cwd.an)<100))[1] / dim(tree_df)[1] 


# old_flm_df <- read_csv(paste0(wdir, 'out\\first_stage\\tree_log_pet_cwd_old.csv')) %>%
#   select(-X1) %>%
#   select(-c(year, aet.an, cwd.an, pet.an)) %>%
#   unique()
# old_sites <- old_flm_df %>% select(collection_id) %>% unique() %>% pull()

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
            pet.ave = mean(pet.an),
            cwd.min = min(cwd.an))

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Trim data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add weighting based on inverse of first stage variance
flm_df <- flm_df %>% 
  mutate(errorweights = 1 / (std.error_cwd.an),
         errorweights2 = sqrt(ntrees),
         pet_errorweights = 1 / (std.error_pet.an),
         int_errorweights = 1 / (std.error_intercept))

# Identify and trim extreme outliers
cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)


flm_df <- flm_df %>% 
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (cwd.spstd<cwd_spstd_bounds[1]) |
           (cwd.spstd>cwd_spstd_bounds[2]) |
           (pet.spstd<pet_spstd_bounds[1]) |
           (pet.spstd>pet_spstd_bounds[2]))

flm_df %>% write.csv(paste0(wdir, "out/first_stage/tree_pet_cwd_std_augmented.csv"))

# flm_df <- flm_df %>%
#   # group_by(species_id) %>%
#   mutate() %>%
#   ungroup()

trim_df <- flm_df %>% 
  filter(outlier==0) %>% 
  drop_na()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate mechanism - tree or site --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod_df = trim_df
mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = errorweights, data=mod_df)
cluster_vcov <- vcovCL(mod, cluster = mod_df$species_id)
coeftest(mod, vcov = vcovCL, cluster = mod_df$species_id)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate mechanism - tree or site --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Focus on trees for which we have full weather history
young_df <- trim_df %>%
  filter(young==TRUE) %>%
  drop_na()

mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trmean + pet.spstd + pet.trmean, data=young_df)
cluster_vcov <- vcovCL(mod, cluster = young_df$collection_id)
coeftest(mod, vcov = vcovCL, cluster = young_df$collection_id)


mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trmean + pet.spstd + pet.trmean, weights = errorweights, data=young_df)
cluster_vcov <- vcovCL(mod, cluster = young_df$collection_id)
coeftest(mod, vcov = vcovCL, cluster = young_df$collection_id)

# Could include all trees from sites with any variation in tree-level drought history (including trees born <1901)
young_sites <- trim_df %>%
  filter(young == T) %>%
  pull(collection_id) %>%
  unique()
youngsite_df <- trim_df %>%
  filter(collection_id %in% young_sites) %>%
  drop_na()
mod <- lm(estimate_cwd.an ~ cwd.spstd + cwd.trmean + pet.spstd + pet.trmean, weights = errorweights, data=youngsite_df)
coeftest(mod, vcov = vcovCL, cluster = youngsite_df$collection_id)

# Include collection fixed effects
mod <- feols(estimate_pet.an ~ cwd.trmean + pet.trmean | collection_id, data=young_df)
summary(mod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Investigate variation by genus  ----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# predict_cwd_effect <- function(trim_df){
#   trim_df <- trim_df %>%
#     filter(estimate_cwd.an>cwd.qlow & estimate_cwd.an<cwd.qhigh,
#            abs(cwd.spstd)<5,
#            abs(pet.spstd)<5) %>%
#     drop_na()
#   gen_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=trim_df)
#   cluster_vcov <- vcovCL(gen_mod, cluster = trim_df$collection_id)
#   print(coeftest(gen_mod, vcov = vcovCL, cluster = trim_df$collection_id))
#   predictions <- prediction(gen_mod, at = list(cwd.spstd = seq(-3, 3, .5)), vcov = cluster_vcov, calculate_se = T) %>% 
#     summary() %>% 
#     rename(cwd.spstd = "at(cwd.spstd)") 
#   return(predictions)
# }
# # ## Plot marginal effects by gymno / angio
# grp_freq <- flm_df %>% 
#   group_by(gymno_angio) %>% 
#   summarise(n_collections = n_distinct(collection_id)) %>% 
#   arrange(desc(n_collections))
# 
# grp_keep <- grp_freq %>% 
#   filter(n_collections>50) %>% 
#   pull(gymno_angio)
# 
# grp_df <- flm_df %>%
#   filter(gymno_angio %in% grp_keep) %>% 
#   group_by(gymno_angio) %>% 
#   nest() %>%
#   mutate(prediction = map(data, predict_cwd_effect))
# 
# grp_df <- grp_df %>% 
#   unnest(prediction) %>% 
#   select(-data)
# 
# grp_df$gymno_angio <- grp_df$gymno_angio %>% recode("angio" = "Angiosperm", "gymno" = "Gymnosperm")
# 
# margins_plot <- ggplot(grp_df, aes(x = cwd.spstd)) + 
#   geom_line(aes(y = Prediction)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   facet_wrap(~gymno_angio, scales = "free", ncol = 1) +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") + 
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::comma)
# margins_plot
# ggsave(paste0(wdir, 'figures\\gymno_angio.svg'), plot = margins_plot)
# 
# ## Plot marginal effects by gymno / angio
# grp_freq <- flm_df %>% 
#   group_by(gymno_angio) %>% 
#   summarise(n_collections = n_distinct(collection_id)) %>% 
#   arrange(desc(n_collections))
# 
# grp_keep <- grp_freq %>% 
#   filter(n_collections>50) %>% 
#   pull(gymno_angio)
# 
# grp_df <- flm_df %>%
#   filter(gymno_angio %in% grp_keep) %>% 
#   group_by(gymno_angio) %>% 
#   nest() %>%
#   mutate(prediction = map(data, predict_cwd_effect))
# 
# grp_df <- grp_df %>% 
#   unnest(prediction) %>% 
#   select(-data)
# 
# grp_df$gymno_angio <- grp_df$gymno_angio %>% recode("angio" = "Angiosperm", "gymno" = "Gymnosperm")
# 
# margins_plot <- ggplot(grp_df, aes(x = cwd.spstd)) + 
#   geom_line(aes(y = Prediction)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill = "darkblue") +
#   facet_wrap(~gymno_angio, scales = "free", ncol = 1) +
#   geom_line(aes(y = upper), linetype = 3) +
#   geom_line(aes(y = lower), linetype = 3) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   xlab("Historic CWD\n(Deviation from species mean)") + 
#   ylab("Predicted sensitivity to CWD") + 
#   theme_bw(base_size = 22) +
#   scale_y_continuous(labels = scales::comma)
# margins_plot
# ggsave(paste0(wdir, 'figures\\gymno_angio.svg'), plot = margins_plot)


predict_cwd_effect <- function(trim_df){
  trim_df <- trim_df %>%
    filter(abs(cwd.spstd)<3,
           abs(pet.spstd)<3) %>%
    drop_na()
  gen_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights=errorweights, data=trim_df)
  cluster_vcov <- vcovCL(gen_mod, cluster = trim_df$collection_id)
  mod_cl <- tidy(coeftest(gen_mod, vcov = vcovCL, cluster = trim_df$collection_id))
  pred_min <- trim_df$cwd.spstd %>% quantile(0.025)
  pred_max <- trim_df$cwd.spstd %>% quantile(0.975)
  predictions <- prediction(gen_mod, at = list(cwd.spstd = seq(pred_min, pred_max, 0.1)), vcov = cluster_vcov, calculate_se = T) %>% 
    summary() %>% 
    rename(cwd.spstd = "at(cwd.spstd)")
  out <- tibble(predictions = list(predictions), model = list(mod_cl))
  return(out)
}



## Plot marginal effects by genus
genus_freq <- mod_df %>% 
  group_by(genus) %>% 
  summarise(n_collections = n_distinct(collection_id),
            min_cwd = min(cwd.spstd),
            max_cwd = max(cwd.spstd),
            range_cwd = max_cwd - min_cwd) %>% 
  arrange(desc(n_collections))

genus_keep <- genus_freq %>% 
  filter(n_collections>2) %>%
  pull(genus)


genus_key <- sp_info %>% 
  select(genus, gymno_angio) %>% 
  unique()

genus_df <- mod_df %>% 
  filter(genus %in% genus_keep) %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(prediction = map(data, predict_cwd_effect)) %>% 
  unnest(prediction) %>% 
  left_join(genus_key, by = "genus")


genus_coefs <- genus_df %>%
  unnest(model) %>% 
  filter(term == "cwd.spstd") %>% 
  select(genus, estimate, p.value)

genus_coefs <- genus_coefs %>% 
  mutate(lab = paste0("Slope: ", as.character(format(estimate, digits = 3, scientific = F)), 
                      "\nP value: ", as.character(format(p.value, digits = 3, scientific = T))))

genus_predictions <- genus_df %>% 
  unnest(predictions) %>% 
  select(-data, -model) %>% 
  left_join(genus_freq, by = "genus") %>% 
  left_join(genus_coefs, by = "genus")

saveRDS(genus_predictions, paste0(wdir, "out/second_stage/genus_mods.rds"))

# histogram <- ggplot(flm_df %>% filter(genus %in% genus_keep), aes(x = cwd.spstd)) + 
#   facet_wrap(~genus) +
#   geom_histogram(aes(fill = gymno_angio), bins = 40, alpha=0.4) +
#   xlim(c(-3, 3)) +
#   theme_bw(base_size = 22) + 
#   ylab("# sites") +
#   xlab("Historic CWD\n(Deviation from species mean)")
# histogram


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





sp <- "pipo"
sp_data <- flm_df %>% filter(species_id == sp)
sp_data %>% filter(cwd.spstd<1) %>% 
  ggplot(aes(x = cwd.spstd, y = estimate_cwd.an)) + 
  geom_point() + 
  geom_point(data = sp_data %>% filter(collection_id=="CO559"), color = 'red', alpha = 1, size = 5) + 
  geom_smooth(method = lm) +
  theme_bw(base_size = 25) +
  ylab("Sensitivity to CWD (\U03B2)")+
  xlab("Historic CWD\n(Deviation from species mean)") +
  geom_hline(yintercept = 0)

sp_data %>% arrange(estimate_cwd.an) %>% select(collection_id, estimate_cwd.an, p.value_cwd.an, cwd.ave, cwd.spstd)
(sp_data %>% arrange(desc(estimate_cwd.an)) %>% select(collection_id, estimate_cwd.an, p.value_cwd.an, cwd.ave, cwd.spstd))[10:20,]
