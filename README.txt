Processing workflow:

# Environment configuration
"0. init.R"

# 1. Parse and detrend ITRDB data
"1a. Pull ITRDB.R" - Pulls ITRDB files from NOAA ftp
"1b. Parse ITRDB.R" - Reads each rwl file, detrends series, and combines into master dendro file

# 2. Prep climate data
"cwd_function.R" - Functions to process monthly weather and soil data to generate annual PET, CWD and AET.
"2a. Pull site topography.R" - Extracts aspect, slope and elevation from ASTER data
"2b. Pull site weather.R" - Extracts and downscales precip and temp history for each site
"2c. Calculate site CWD.R" - Combines topography, swc and weather to calculate annual CWD / AET history for each site
"2d. Historic cwdraster" - 
"2e. CMIP5 projections" -
"downscaling"


# 3. Calculate species niches
"3a. Collate_ranges.R" - Combine individual species range maps into master file
"3b. Species niche.R" - Calculate historic weather niche for each species

# 4. First stage regressions
"4a. First stage.R" - Calculates weather sensitivity regression for each site
"4b. Second stage.R" - Analyzes tree-level weather sensitivity based on historic conditions
"4c. Missing data.R" - Tracks data erosion to determine which steps are leading to data drops

# 5. Second stage regressions

# 6. Prediction
"6. Prediction.R"



# 7. Generate figures, tables and stats
