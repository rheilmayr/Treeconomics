The code in this repository reproduces the analysis presented in Heilmayr, Dudney and Moore, 2023 (https://www.science.org/doi/10.1126/science.adi1071)

# Data download
Pre-processed data can be downloaded from the Harvard Dataverse: https://doi.org/10.7910/DVN/DZEXQN

# Analysis
Once data have been downloaded, the "Full_analysis.R" provides a full walkthrough of the integrated analysis.

# Outline of analysis
## Environment configuration
- "0. init.R" - creates symbolic link to point towards file directory


## 1. Parse and detrend ITRDB data
- "1a. Pull ITRDB.R" - Pulls ITRDB files from NOAA ftp
- "1b. Parse ITRDB.R" - Reads each rwl file, detrends series, and combines into master dendro file
- "1c. Parse Klesse.R" - Reads FIA rw data from Klesse et al., 2018. Note that input data is confidential and can't be provided publicly


## 2. Prep climate data
- "2a. Pull site topography.R" - Extracts aspect, slope and elevation from ASTER data
- "2b. Pull site weather.R" - Extracts and downscales precip and temp history for each site
- "2c. Calculate site CWD.R" - Combines topography, swc and weather to calculate annual CWD / AET history for each site
- "2d. Historic cwdraster.R" - Calculates historic CWD / AET rasters using CRU data
- "2e. CMIP5 projections.R" - Calculates future CWD / AET rasters using CMIP5 data


## 3. Calculate species niches
- "3a. Collate ranges.R" - Combine individual species range maps into master file
- "3b. Species niche.R" - Calculate historic weather niche for each species. Apply species-level standardization to past and future weather and climate data


## 4. First stage regressions
- "4a. First stage.R" - Runs weather sensitivity regressions for each site
- "4b. DNLM first stage.R" - Re-runs first stage using dnlm model which is used as a robustness check


## 5. Second stage regressions
- "5a. Second stage.R" - Analyzes heterogeneity in site-level weather sensitivity across historic climate
- "5b. Robustness.R" - Runs robustness checks on second stage regressions


## 6. Prediction
- "6. Prediction.R" - Predicts sensitivity and RWI changes


## 7. Generate figures
- "7. Figure XX.R" - Generate the individual figures in the paper and SI


## 8. Generate tables and stats
- "8a. Table S1.R" - Generate supplementary table 1
- "8b. Paper stats.R" - Calculate all statistics reported in paper
- "8c. Missing data.R" - Track missing data


## Helper functions
- "f_cwd_function.R" - Functions to process monthly weather and soil data to generate annual PET, CWD and AET.
- "f_spec_chart_function.R" - Function to run specification chart for robustness results

