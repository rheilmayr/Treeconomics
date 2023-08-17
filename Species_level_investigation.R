library(tidyverse)
library(broom)
library(feols)

### Define path
wdir <- 'remote/'

# 1. Site-level regressions
flm_df <- read_csv(paste0(wdir, "out/first_stage/site_pet_cwd_std_augmented.csv")) 
trim_df <- flm_df %>% filter(outlier == 0)

# 2. Species niche
niche_df <- read_csv(paste0(wdir, "out//climate//clim_niche.csv")) %>%
  select(-...1, species_id = sp_code)



sp_reg <- function(sp_data){
  sp_mod <- lm(estimate_cwd.an ~ cwd.spstd + pet.spstd, weights = cwd_errorweights, data = sp_data)
  sp_mod <- tidy(sp_mod)
  cwd_result <- sp_mod %>% filter(term == "cwd.spstd")
  return(cwd_result)
}

sp_nest <- trim_df %>% 
  group_by(species_id) %>% 
  mutate(species_n = n()) %>%
  filter(species_n > 10) %>% 
  nest()

sp_regs <- sp_nest %>% 
  mutate(mod = map(data, .f = sp_reg)) %>% 
  unnest(mod) %>% 
  select(-data, term)

sp_regs <- sp_regs %>% 
  left_join(niche_df, by = "species_id")

sp_regs %>% 
  filter(p.value <= 0.05)
sp_regs %>% 
  filter(p.value > 0.05)


sp_regs %>%
  filter(p.value <= 0.05) %>% 
  ggplot(aes(x = cwd_mean, y = estimate, color =species_id)) +
  geom_point() +
  geom_smooth()

trim_df %>% 
  filter(species_id == "abla") %>% 
  ggplot(aes(y = estimate_cwd.an, x = cwd.spstd))+
  geom_point() +
  geom_smooth(method = "lm")

test_data <- trim_df %>% filter(species_id == "abla")
test_mod <- lm("estimate_cwd.an ~ cwd.spstd + pet.spstd", data = test_data, weights = cwd_errorweights)
summary(test_mod)
mod <- lm(estimate ~ pet_mean + cwd_mean, data = sp_regs)
summary(mod)

mod <- lm(estimate ~ pet_mean + cwd_mean, data = sp_regs %>% filter(p.value <=0.05))
summary(mod)

mod <- lm(estimate ~ pet_mean + cwd_mean, data = sp_regs, weights = (1/std.error))
summary(mod)

mod <- lm(estimate ~ pet_sd + cwd_sd, data = sp_regs)
summary(mod)

mod <- lm(estimate ~ pet_sd + cwd_sd, data = sp_regs)
summary(mod)

### Genus
trim_df <- fs_spl %>% 
  filter(outlier==0) %>% 
  drop_na()

run_ss_conley <- function(data, outcome = "cwd_coef"){
  if (outcome == "cwd_coef") {
    error_weights = data$cwd_errorweights
  } else if (outcome == "pet_coef") {
    error_weights = data$pet_errorweights
  } else if (outcome == "int_coef") {
    error_weights = data$int_errorweights
  }
  formula <- as.formula(paste(outcome, " ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)"))
  mod <- feols(formula, data=data, weights = error_weights,
               vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
  
  # formula <- as.formula(paste(outcome, " ~ cwd.spstd + I(cwd.spstd**2) + pet.spstd + I(pet.spstd**2)"))
  # mod <- lm(formula, data=data, weights = error_weights)
  return(mod)
}

# run_ss_conley <- function(data, outcome = "estimate_cwd.an"){
#   formula <- as.formula(paste(outcome, " ~ cwd.spstd + pet.spstd")) ## ADD BACK WEIGHTING HERE?
#   mod <- feols(formula, data = data, weights = data$cwd_errorweights, 
#         vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
#   return(mod)
# }

genus_lookup <- trim_df %>% select(collection_id, genus)
genus_freq <- trim_df %>% 
  group_by(genus) %>% 
  summarise(n_collections = n_distinct(collection_id),
            min_cwd = min(cwd.spstd),
            max_cwd = max(cwd.spstd),
            range_cwd = max_cwd - min_cwd) %>% 
  arrange(desc(n_collections))

gymno_key <- sp_info %>% 
  select(genus, gymno_angio) %>% 
  unique()



# N = 10: 16 genera; 26: 12 genera; 50: 9 genera
genus_keep <- genus_freq %>%
  top_n(12,n_collections) %>% 
  pull(genus)


genus_df <- trim_df %>% 
  mutate(int_coef = estimate_intercept,
         pet_coef = estimate_pet.an, 
         cwd_coef = estimate_cwd.an) %>% 
  filter(genus %in% genus_keep) %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(model_estimates = map(data, run_ss_conley))

genus_df$model_estimates



### Binned plot of cwd sensitivity
gen_marg_fx_df <- function(gen_mod, gen_data){
  # cwd_min = 
  cwd_inc <- 0.1
  at_pet <- 0
  cwd_min <- gen_data %>% pull(cwd.spstd) %>% min()
  cwd_max <- gen_data %>% pull(cwd.spstd) %>% max()
  gen_pred <- predictions(gen_mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = seq(cwd_min,cwd_max,cwd_inc)))
  return(gen_pred)
}


pull_coefs <- function(gen_mod, gen_dat){
  median_cwd <- gen_dat %>% pull(cwd.spstd) %>% median()
  lht <- linearHypothesis(gen_mod, c(paste0('cwd.spstd + ', as.character(median_cwd), ' * I(cwd.spstd^2) = 0')))
  pvalue <- lht$`Pr(>Chisq)`[2] %>% round(digits = 3)
  coefs = gen_mod$coefficients
  me <- (coefs['cwd.spstd'] + median_cwd * 2 * coefs['I(cwd.spstd^2)']) %>% round(digits = 3)
  
  # coef_table <- gen_mod %>% coeftable()
  # coef <- coef_table[2,1]  %>% round(digits = 3)
  # p <- coef_table[2, 4] %>% round(digits = 3)
  n <- gen_mod$nobs
  
  label <-  paste0("   n sites: ", n, ";\n   slope: ", me, ";\n   p value: ", pvalue, "\n")
  return(label)
}


coef_labels <- genus_df %>%
  mutate(labels = map2(model_estimates, data, pull_coefs)) %>%
  select(genus, labels) %>%
  unnest(labels) %>%
  arrange(genus)


genus_predictions <- genus_df %>% 
  mutate(predictions = map2(model_estimates, data, gen_marg_fx_df))

genus_predictions <- genus_predictions %>% 
  unnest(predictions) %>% 
  select(-data, -model_estimates) %>% 
  arrange(genus)

gen_plot <- genus_predictions %>% 
  # filter(genus %in% genus_keep) %>%
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  # geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkblue") +
  facet_wrap(~genus, scales = "free", ncol = 3) +
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)) +
  geom_line(aes(y = conf.low), linetype = 3) +
  geom_line(aes(y = conf.high), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") +
  xlim(c(-3, 3))


gen_plot <- gen_plot +
  geom_text(data = coef_labels, aes(label = labels, x = -Inf, y = -Inf),
            hjust = 0, vjust = 0)

gen_plot











run_ss_conley <- function(data, outcome = "cwd_coef"){
  if (outcome == "cwd_coef") {
    error_weights = data$cwd_errorweights
  } else if (outcome == "pet_coef") {
    error_weights = data$pet_errorweights
  } else if (outcome == "int_coef") {
    error_weights = data$int_errorweights
  }
  formula <- as.formula(paste(outcome, " ~ cwd.spstd + pet.spstd"))
  mod <- feols(formula, data=data, weights = error_weights,
               vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
  
  # formula <- as.formula(paste(outcome, " ~ cwd.spstd + I(cwd.spstd**2) + pet.spstd + I(pet.spstd**2)"))
  # mod <- lm(formula, data=data, weights = error_weights)
  return(mod)
}





sp_freq <- flm_df %>% 
  group_by(species_id) %>% 
  tally()

# N = 10: 16 genera; 26: 12 genera; 50: 9 genera
sp_keep <- sp_freq %>%
  top_n(25,n) %>% 
  pull(species_id)


sp_df <- trim_df %>% 
  mutate(int_coef = estimate_intercept,
         pet_coef = estimate_pet.an, 
         cwd_coef = estimate_cwd.an) %>% 
  filter(species_id %in% sp_keep) %>% 
  group_by(species_id) %>% 
  nest() %>% 
  mutate(model_estimates = map(data, run_ss_conley))

sp_df$model_estimates



### Binned plot of cwd sensitivity
gen_marg_fx_df <- function(gen_mod, gen_data){
  # cwd_min = 
  cwd_inc <- 0.1
  at_pet <- 0
  cwd_min <- gen_data %>% pull(cwd.spstd) %>% min()
  cwd_max <- gen_data %>% pull(cwd.spstd) %>% max()
  gen_pred <- predictions(gen_mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = seq(cwd_min,cwd_max,cwd_inc)))
  return(gen_pred)
}


pull_coefs <- function(gen_mod, gen_dat){
  median_cwd <- gen_dat %>% pull(cwd.spstd) %>% median()
  lht <- linearHypothesis(gen_mod, c(paste0('cwd.spstd + ', as.character(median_cwd), ' * I(cwd.spstd^2) = 0')))
  pvalue <- lht$`Pr(>Chisq)`[2] %>% round(digits = 3)
  coefs = gen_mod$coefficients
  me <- (coefs['cwd.spstd'] + median_cwd * 2 * coefs['I(cwd.spstd^2)']) %>% round(digits = 3)
  
  # coef_table <- gen_mod %>% coeftable()
  # coef <- coef_table[2,1]  %>% round(digits = 3)
  # p <- coef_table[2, 4] %>% round(digits = 3)
  n <- gen_mod$nobs
  
  label <-  paste0("   n sites: ", n, ";\n   slope: ", me, ";\n   p value: ", pvalue, "\n")
  return(label)
}


sp_predictions <- sp_df %>% 
  mutate(predictions = map2(model_estimates, data, gen_marg_fx_df))

sp_predictions <- sp_predictions %>% 
  unnest(predictions) %>% 
  select(-data, -model_estimates) %>% 
  arrange(species_id)

sp_plot <- sp_predictions %>% 
  # filter(genus %in% genus_keep) %>%
  ggplot(aes(x = cwd.spstd)) + 
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
  # geom_ribbon(aes(ymin=cwd_ci_min, ymax=cwd_ci_max), alpha=0.2, fill = "darkblue") +
  facet_wrap(~species_id, scales = "free", ncol = 5) +
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)) +
  geom_line(aes(y = conf.low), linetype = 3) +
  geom_line(aes(y = conf.high), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Historic CWD\n(Deviation from species mean)") + 
  ylab("Predicted sensitivity to CWD") +
  xlim(c(-3, 3))


sp_plot <- sp_plot +
  geom_text(data = coef_labels, aes(label = labels, x = -Inf, y = -Inf),
            hjust = 0, vjust = 0)

sp_plot
