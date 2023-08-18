#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 8/7/23
# Purpose: Create figure S4 illustrating variability across genera
#
# Input files:
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
wdir <- 'remote/'

# 1. Genus second stage model
genus_models <- readRDS(paste0(wdir, "2_output/second_stage/ss_conley_genus.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define style  ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_palette <- scale_colour_brewer(
  type = "seq",
  palette = 5,
  direction = 1,
  aesthetics = "colour"
)

theme_set(
  theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Marginal effects by genera ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


coef_labels <- genus_models %>%
  mutate(labels = map2(model_estimates, data, pull_coefs)) %>%
  select(genus, labels) %>%
  unnest(labels) %>%
  arrange(genus)


genus_predictions <- genus_models %>% 
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
  facet_wrap(~genus, scales = "free", ncol = 2) +
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
  xlim(c(-4, 3))


gen_plot <- gen_plot +
  geom_text(data = coef_labels, aes(label = labels, x = -Inf, y = -Inf),
            hjust = 0, vjust = 0)

gen_plot
ggsave(paste0(wdir, '3_results/figures/FigS4_genus_margins_nonlinear.svg'), gen_plot, width = 10, height = 11, units = "in")

