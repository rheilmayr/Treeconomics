#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 4/29/19
# Purpose: Creates a symlink from code directory to a directory storing project data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(R.utils)

# Define the path to your local code directory
code_dir <- "/Users/treelife/Documents/Treeconomics/Treecon/"

# Define the path to your local google drive Treeconomics\\Data directory 
data_dir <- "/Users/treelife/Google Drive/Treeconomics/Data/"

  
createLink(paste0(code_dir, 'remote'), data_dir, overwrite = FALSE)

