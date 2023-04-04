#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 4/29/19
# Purpose: Creates a symlink from code directory to a directory storing project data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(R.utils)

# Define the path to your local code directory
code_dir <- "D:/dev/Treeconomics/"

# Define the path to your local google drive Treeconomics\\Data directory 
data_dir <- "J:/.shortcut-targets-by-id/10TtqG9P3BY70rcYp-WACmO38J5zBeflA/Treeconomics/Data/"

createLink(paste0(code_dir, 'remote'), data_dir, overwrite = FALSE)

