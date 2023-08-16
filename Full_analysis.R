

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Preparation -------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Note: New users should edit this file to create a symbolic link pointing
## from the code repository, to the directory where data is being stored.

source("0_init.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Steps 1-2: Process data and standardize climate and weather data -------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Note: Raw data for these scripts are large and often require a cluster to 
## process efficiently. As a result, we  do not provide the input data in our
## data repository, and comment out these initial processing scripts.
## However, annotations in each script direct interested readers towards the original
## source for each dataset, and code allows readers to see how data were processed.

# source("1a. Pull ITRDB.R")
# source("1b. Parse ITRDB.R")
# source("1c. Parse FIA.R")
# source("2a. Pull site topography.R")
# source("2b. Pull site weather.R")
# source("2c. Calculate site cwd.R")
# source("2d. Historic cwdraster.R")
# source("2e. cmip5 cwd projections.R")
# source("3a. Collate ranges.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 3: Estimate site-level sensitivity to CWD and PET -----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Run first stage regressions, estimating site-level sensitivity to CWD and PET
source("4a. First stage.R")

## Run alternate version of first stage, using dynamic lag, nonlinear model to
## measure cumulative growth impact over 15 years
source("4b. DLNM first stage.R")

## Re-run first stage comparing FIA tree ring data against comparable ITRDB data.
## Note: Access to FIA data is not public. As a result, we are unable to share the 
## necessary input data to run script 4c

# source("4c. First stage FIA.R") 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 4: Characterize heterogeneity in sensitivity ----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Run second stage regressions, estimating heterogeneity in site-level sensitivity
source("5a. Second stage.R")

## Run alternate specifications of second stage model
source("5b. Robustness.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Steps 5-6: Predict sensitivity and changes in RWI ----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Use second stage regression results to generate rasters of predicted
## sensitivities and changes in RWI through 2100
source("6. Prediction.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generate figures -------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Fig. 1: Illustration of hypotheses
source("7. Figure 1.R")

## Fig. 2: Illustration of methods using a single species
source("7. Figure 2.R")

## Fig. 3: Key sensitivity patterns
source("7. Figure 3.R")

## Fig. 4: Changes in exposure
source("7. Figure 4.R")

## Fig. 5: Integrated results illustrating impacts on RWI (vulnerability)
source("7. Figure 5.R")

## Fig. S1: Summary of data
source("7. Figure S1.R")

## Fig. S2: Comparison of climate variation in ITRDB and species range
source("7. Figure S2.R")

## Fig. S3: First stage results
source("7. Figure S3.R")

## Fig. S4: Second stage results by genera
source("7. Figure S4.R")

## Fig. S5: Robustness of second stage results
source("7. Figure S5.R")

## Fig. S6: Comparison of FIA and ITRDB results
## Note: Access to FIA data is not public. As a result, we are unable to share the 
## necessary input data to run script "7. Figure S6.R"
source("7. Figure S6.R")

## Fig. S7: Prediction comparison across hypothesized patterns of sensitivity
source("7. Figure S7.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generate tables --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Table 1: Summary of species analyzed
source("8a. Table S1.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate stats reported in paper --------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Calculate all statistics reported in the paper
source("8b. Paper stats.R")

## Summarize why sites drop out of analysis (described in methods)
source("8c. Missing data.R")

