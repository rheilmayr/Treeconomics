#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 4/29/19
# Purpose: Creates a symlink from code directory to a directory storing project data
# Note: Should be run with administrator priviledges
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Build compute environment -----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(renv)
renv::restore()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create symbolic link to data directory -----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(R.utils)

# Define the path to your local code directory
code_dir <- "D:/dev/Treeconomics/"

# Define the root directory where you've unzipped the replication dataset.
# The replication dataset can be downloaded from XX.
# After unzipping, the root directory should contain the "1_input_processed" directory.
data_dir <- "G:/.shortcut-targets-by-id/10TtqG9P3BY70rcYp-WACmO38J5zBeflA/Treeconomics/Data/replication - downscale/"

createLink(paste0(code_dir, 'remote'), data_dir, overwrite = FALSE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add output directory structure -----------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir.create(file.path("remote/2_output/"), showWarnings = TRUE)
dir.create(file.path("remote/2_output/climate/"), showWarnings = TRUE)
dir.create(file.path("remote/2_output/first_stage/"), showWarnings = TRUE)
dir.create(file.path("remote/2_output/second_stage/"), showWarnings = TRUE)
dir.create(file.path("remote/2_output/predictions/"), showWarnings = TRUE)

dir.create(file.path("remote/3_results/"), showWarnings = TRUE)
dir.create(file.path("remote/3_results/figures/"), showWarnings = TRUE)
dir.create(file.path("remote/3_results/figures/methods_panels/"), showWarnings = TRUE)
dir.create(file.path("remote/3_results/tables/"), showWarnings = TRUE)


