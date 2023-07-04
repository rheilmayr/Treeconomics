The code in this repository reproduces the analysis conducted in Heilmayr, Dudney and Moore, 2023. 

# Processing workflow
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


## 7. Generate figures, tables and stats
- "7a. Figures.R" - Generates most figures for paper
- "7b. Conceptual figure.R" - Generates figure 1 illustrating conceptual link between sensitivity, exposure and vulnerability
- "7c. Species level figures.R" - Generates figure 2 which illustrates methods using PIPO
- "7d. Paper tables.R" - Generates supplemental table summarizing species included in analysis 
- "7e. Paper stats.R" - Reproduces all stats from the paper


## Helper functions
- "f_cwd_function.R" - Functions to process monthly weather and soil data to generate annual PET, CWD and AET.
- "f_spec_chart_function.R" - Function to run specification chart for robustness results




## Data to provide (to run starting from script 3b)
#   merged_ranges.shp:
#   HistoricCWD_AETGrids.Rdat: Rasters describing historic CWD and AET
#     Generated using historic_cwdraster.R
#   monthlycrubaseline_tas:
#   cmip5_cwdaet_start.Rdat:
#   cmip5_cwdaet_end.Rdat:
#   essentialcwd_data.csv:
#   site_summary.csv:
#   rwi_long.csv: Directory containing processed RWI data from "1b. Parse ITRDB.R"
#   species_gen_gr.csv: Annotated data about species.

