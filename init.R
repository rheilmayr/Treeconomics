#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 4/29/19
# Purpose: Creates a symlink from code directory to a directory storing project data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Define the path to your local code directory
code_dir <- 'D:\\dev\\Treeconomics\\'

# Define the path to your local google drive Treeconomics\\Data directory 
data_dir <- 'D:\\cloud\\Google Drive\\Treeconomics\\Data\\'

library(R.utils)
createLink(paste0(code_dir, 'remote\\'), data_dir, overwrite = FALSE)