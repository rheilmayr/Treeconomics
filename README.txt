Processing workflow:

1. Parse ITRDB data
"1-a. Pull ITRDB.R" - Pulls ITRDB files from NOAA ftp
"1-b. Parse ITRDB.R" - Reads each rwl file, detrends series, and combines into master dendro file

2. Prep climate data
"cwd_function.R" - Functions to process monthly weather and soil data to generate annual PET, CWD and AET.
"2-a. Pull site topography.R" - Extracts aspect, slope and elevation from ASTER data
"2-b. Pull site weather.R" - Extracts and downscales precip and temp history for each site
"2-c. Calculate site CWD.R" - Combines topography, swc and weather to calculate annual CWD / AET history for each site

historic_cwdraster
cmpi5_projections


3. Calculate species niches
"3-a. Collate_ranges.R" - Combine individual species range maps into master file
"3-b. Species niche.R" - Calculate historic weather niche for each species

4. Analysis