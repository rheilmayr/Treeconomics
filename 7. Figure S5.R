#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Create specification chart illustrating robustness of results
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(patchwork)
library(viridis)
source("f_spec_chart_function.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wdir <- 'remote/'

# 1. Specification comparisons
specs <- read_rds(paste0(wdir, "2_output/second_stage/robustness_specs.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create spec chart --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Create figure
highlight_n <- 1

## Define label structure
labels <- list("Detrending" = c("Spline", "Negative binomial", "Autoregressive"),
               "Energy control" = c("PET", "Mean temperature"),
               "First stage" = c("Annual RWI", "Cumulative\ndynamic lag"),
               "Model structure" = c("Two stage", "Single stage", "Squared\nCWD and PET", "Linear\nCWD and PET", "Temporally\ncorrelated error"),
               "Trimming and weighting" = c("Weight by\ninverse of s.e.", "Trim\noutliers in y", "Trim\noutliers in X"))


svg(paste0(wdir, '3_results/figures/FigS5_robustness.svg'), width = 9, height = 14)
robustness_fig <- schart(specs, labels, highlight=highlight_n, order = "asis", 
                         heights = c(.4,.6),
                         n=c(1, 2, 1, 1, 3, 3), ci = c(.95), 
                         ylab = "Slope of line relating sites' historic\nCWD to CWD's impact on growth\n(evaluated at median historic CWD)",
                         col.est=c("grey80", "dodgerblue4"),
                         col.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         bg.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         lwd.symbol=1)

text(x = 1, y = -.008, label = "1",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 3, y = -.008, label = "2",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 4, y = -.008, label = "3",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 6, y = -.008, label = "4",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 8, y = -.008, label = "5",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 10, y = -.008, label = "6",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 11, y = -.008, label = "7",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 12, y = -.008, label = "8",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 14, y = -.008, label = "9",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 15, y = -.008, label = "10",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 16, y = -.008, label = "11",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

dev.off()

