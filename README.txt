Processing workflow:

0. Parse ITRDB data
0a. XXXX: Script that converts rwl files into tree ring database

1. Prep climate data
1a. XXXX: Missing script that pulls worldclim and cru data
1b. downscaling.R: Downscales CRU data using worldclim
1c. XXXX: Missing script that combines downscaling correction, CRU data, and soil data into "181116-climate_soil_data_with_corrections.csv" file that serves as input to cwd_function.
1d. cwd_function.R: Processes monthly weather and soil data to generate annual PET, CWD and AET.

2. Run analysis
2a. first_stage.R: 
2b. second_stage.R: