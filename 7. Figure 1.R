#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 4/2/23
# Purpose: Generate intro hypotheses figure
#
# Input files:
#
# ToDo:
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(patchwork)
library(viridis)
library(forcats)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wdir <- 'remote/'

# 1. Prediction rasters
rwi_list <- list.files(paste0(wdir, "2_output/predictions/sp_rwi/"), pattern = ".gz", full.names = TRUE)
sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))

theme_set(
  theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text=element_text(family ="Helvetica"))
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prep dataframe --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins <- seq(-2.5,2.5,0.2)
labels <- seq(-2.4,2.4,0.2)
bin_length <- length(bins)
dn_sens <- seq(-0.1, 0, (0.1) / (bin_length - 2))
cons_sens <- rep(-0.05, times = bin_length - 1)
res_sens <- seq(0, -0.1, (-0.1) / (bin_length - 2))

exposure_df <- sp_predictions %>% 
  filter(abs(cwd_hist) < 2.5) %>% 
  mutate(cwd_bin = cut(cwd_hist, bins)) %>% 
  mutate(cwd_change = cwd_cmip_end_mean - cwd_cmip_start_mean) %>% 
  group_by(cwd_bin) %>% 
  summarise(cwd_change = mean(cwd_change),
            cwd_end = mean(cwd_cmip_end_mean),
            cwd_start = mean(cwd_cmip_start_mean)) %>% 
  mutate(labels = labels)

exposure_df <- exposure_df %>% 
  mutate(dn = dn_sens,
         cons = cons_sens,
         res = res_sens)

# Sensitivity
sens_df <- exposure_df  %>% 
  pivot_longer(cols = c("dn", "cons", "res")) %>% 
  select(cwd_bin, labels, name, sens = value) %>% 
  mutate(name = fct_relevel(name, "cons", "res", "dn"))

# Exposure
exp_df <- exposure_df %>% 
  select(cwd_bin, labels, cwd_change)

# Vulnerability
vuln_df <- exposure_df  %>% 
  mutate(dn = (dn * cwd_end) - (dn * cwd_start),
         cons = (cons * cwd_end) - (cons * cwd_start),
         res = (res * cwd_end) - (res * cwd_start)) %>% 
  pivot_longer(cols = c("dn", "cons", "res")) %>% 
  select(cwd_bin, labels, name, rwi_change = value) %>% 
  mutate(name = fct_relevel(name, "cons", "res", "dn"))




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

discrete_pal <- c("#1e9c89","darkgrey", "#472e7c")
labels = c("H0: Consistent","H1: Drought-naive", "H2: Dry-range")

panel1 <- sens_df %>% 
  ggplot(aes(x = labels, y = sens, fill = name, group = name, color = name, linetype = name)) +
  geom_smooth() +
  scale_fill_manual(values = discrete_pal,labels=c("H0: Consistent","H1: Dry-range ", "H2: Drought-naïve"))+
  scale_color_manual(values = discrete_pal, labels=c("H0: Consistent","H1: Dry-range ", "H2: Drought-naïve"))+
  scale_linetype_manual(values=c("solid", "dotted", "dashed"), labels=c("H0: Consistent","H1: Dry-range ", "H2: Drought-naïve"))+
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5))+
  ylab("Sensitivity (Δ RWI / Δ CWD)")+
  xlab("Standardized aridity\n(Deviation from species mean)")+
  guides(fill="none", color="none", linetype = "none")+
  ggtitle("Sensitivity")

panel1

panel2 <- exp_df %>% 
  ggplot(aes(x = labels, y = cwd_change)) +
  geom_col(position = "dodge", alpha=.4)+
  ylab("Change in CWD")+
  xlab("Standardized aridity\n(Deviation from species mean)")+
  ggtitle("Exposure")+
  theme(plot.title = element_text(hjust = 0.5))


panel3 <- vuln_df %>%
  ggplot(aes(x = labels, y = rwi_change, fill = name, group = name, linetype = name)) +
  geom_col(position = "dodge", width =.2, alpha=0.4)+
  geom_smooth(aes(color=name), se=F, method="gam")+
  scale_fill_manual(values = discrete_pal,labels=c("H0: Consistent","H1: Dry-range ", "H2: Drought-naïve"))+
  scale_color_manual(values = discrete_pal, labels=c("H0: Consistent","H1: Dry-range ", "H2: Drought-naïve"))+
  scale_linetype_manual(values=c("solid", "dotted", "dashed"), labels=c("H0: Consistent","H1: Dry-range ", "H2: Drought-naïve"))+
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5), legend.position = "bottom")+
  ylab("Change in RWI")+
  xlab("Standardized aridity\n(Deviation from species mean)")+
  ggtitle("Vulnerability")
panel3

  

full_fig <- panel1 + panel2 + panel3 + 
  plot_layout(guides = "collect") &  
  theme(legend.position = 'bottom') & 
  plot_annotation(tag_levels="A") & 
  theme(plot.tag = element_text(face = 'bold', size=12))

full_fig

ggsave(paste0(wdir, '3_results/figures/Fig1_conceptual.svg'), plot = full_fig, bg= 'transparent', width = 10, height = 4)

#dims 15x6
