#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 2026-02-13
# Purpose: Pull TerraClimate data for FIA site locations
#
# Input files:
#   site_summary_fia.csv: FIA site metadata with lat/lon coordinates
#
# Output:
#   site_climate_fia.csv: Monthly TerraClimate data (tmax, tmin, ppt, def, pet)
#                         for all FIA sites from 1958-present
#
# Data source: TerraClimate (http://www.climatologylab.org/terraclimate.html)
# THREDDS server: http://thredds.northwestknowledge.net:8080/thredds/catalog.html
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(ncdf4)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define paths
wdir <- 'remote/'
dendro_dir <- paste0(wdir, '1_input_processed/dendro/')
output_dir <- paste0(wdir, '0_raw/TerraClimate/')

# Load site locations
site_df <- read_csv(paste0(dendro_dir, "site_summary_fia.csv"))

# Get unique site locations
sites <- site_df %>%
  distinct(collection_id, LAT, LON) %>%
  filter(!is.na(LAT), !is.na(LON))

cat(sprintf("Loaded %d unique FIA sites\n", nrow(sites)))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define extraction function --------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extract_terraclimate_point <- function(lon, lat, variable) {
  # Purpose: Extract TerraClimate time series for a single point location
  # Inputs:
  #   lon: Longitude (decimal degrees, -180 to 180)
  #   lat: Latitude (decimal degrees, -90 to 90)
  #   variable: TerraClimate variable name (e.g., "tmax", "tmin", "ppt", "def", "pet")
  # Returns:
  #   Data frame with columns: year, month, value

  # Build TerraClimate THREDDS URL
  base_url <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/"
  nc_url <- paste0(base_url, "agg_terraclimate_", variable, "_1958_CurrentYear_GLOBE.nc")

  # Try to open NetCDF connection
  nc <- tryCatch({
    nc_open(nc_url)
  }, error = function(e) {
    warning(sprintf("Failed to open NetCDF for variable %s: %s", variable, e$message))
    return(NULL)
  })

  if (is.null(nc)) {
    return(NULL)
  }

  # Ensure cleanup on function exit
  on.exit(nc_close(nc))

  # Extract coordinate arrays
  nc_lon <- ncvar_get(nc, "lon")
  nc_lat <- ncvar_get(nc, "lat")

  # Find closest grid cell to input coordinates
  # TerraClimate resolution is 1/24 degree (~4km)
  lon_diff <- abs(nc_lon - lon)
  lat_diff <- abs(nc_lat - lat)

  lon_idx <- which.min(lon_diff)
  lat_idx <- which.min(lat_diff)

  # Check if match is within tolerance (1/48 degree = ~2km)
  if (lon_diff[lon_idx] > 1/48 || lat_diff[lat_idx] > 1/48) {
    warning(sprintf("No grid cell found within tolerance for lon=%.4f, lat=%.4f", lon, lat))
    return(NULL)
  }

  # Extract time series for this location
  # start: c(lon_index, lat_index, first_time)
  # count: c(1, 1, all_times)
  time_dim <- nc$dim$time$len

  values <- ncvar_get(nc, variable,
                      start = c(lon_idx, lat_idx, 1),
                      count = c(1, 1, time_dim))

  # Get time dimension and convert to year/month
  time_vals <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value

  # TerraClimate time is "days since 1900-01-01"
  # Convert to date
  origin_date <- as.Date("1900-01-01")
  dates <- origin_date + time_vals

  # Create output data frame
  result <- tibble(
    year = as.integer(format(dates, "%Y")),
    month = as.integer(format(dates, "%m")),
    value = as.numeric(values)
  )

  return(result)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extract TerraClimate data for all sites -------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define variables to extract
variables <- c("tmax", "tmin", "ppt", "def", "pet")

cat("Starting TerraClimate data extraction...\n")
cat(sprintf("Processing %d sites × %d variables = %d extractions\n",
            nrow(sites), length(variables), nrow(sites) * length(variables)))

# Initialize storage for results
all_climate_data <- list()
extraction_count <- 0
failed_extractions <- 0

# Loop through each site
for (i in 1:nrow(sites)) {
  site <- sites[i, ]

  # Progress indicator every 50 sites
  if (i %% 50 == 0) {
    cat(sprintf("Progress: %d/%d sites (%.1f%%) - %d failures\n",
                i, nrow(sites), 100 * i / nrow(sites), failed_extractions))
  }

  # Extract all variables for this site
  site_data <- NULL

  for (var in variables) {
    extraction_count <- extraction_count + 1

    # Extract time series for this variable
    var_data <- extract_terraclimate_point(
      lon = site$LON,
      lat = site$LAT,
      variable = var
    )

    if (is.null(var_data)) {
      failed_extractions <- failed_extractions + 1
      next
    }

    # Rename value column to variable name
    var_data <- var_data %>%
      rename(!!var := value)

    # Join with existing site data
    if (is.null(site_data)) {
      site_data <- var_data
    } else {
      site_data <- site_data %>%
        left_join(var_data, by = c("year", "month"))
    }
  }

  # Add collection_id and store
  if (!is.null(site_data)) {
    site_data <- site_data %>%
      mutate(collection_id = site$collection_id, .before = year)

    all_climate_data[[i]] <- site_data
  }
}

cat(sprintf("\nExtraction complete: %d/%d successful\n",
            extraction_count - failed_extractions, extraction_count))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine and save output -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine all site data into single data frame
climate_df <- bind_rows(all_climate_data)

cat(sprintf("Combined data: %d rows, %d unique sites\n",
            nrow(climate_df),
            length(unique(climate_df$collection_id))))

# Verify output structure
cat("\nOutput structure:\n")
cat(sprintf("  Columns: %s\n", paste(names(climate_df), collapse = ", ")))
cat(sprintf("  Date range: %d-%d to %d-%d\n",
            min(climate_df$year), min(climate_df$month),
            max(climate_df$year), max(climate_df$month)))

# Data quality checks
cat("\nData quality summary:\n")
cat(sprintf("  Missing values: tmax=%d, tmin=%d, ppt=%d, def=%d, pet=%d\n",
            sum(is.na(climate_df$tmax)),
            sum(is.na(climate_df$tmin)),
            sum(is.na(climate_df$ppt)),
            sum(is.na(climate_df$def)),
            sum(is.na(climate_df$pet))))

# Basic range checks
cat(sprintf("  Temp range: %.1f to %.1f °C\n",
            min(climate_df$tmin, na.rm = TRUE),
            max(climate_df$tmax, na.rm = TRUE)))
cat(sprintf("  Precip range: %.1f to %.1f mm\n",
            min(climate_df$ppt, na.rm = TRUE),
            max(climate_df$ppt, na.rm = TRUE)))

# Save to CSV
output_file <- paste0(output_dir, "site_climate_fia.csv")
write_csv(climate_df, output_file)

cat(sprintf("\nData saved to: %s\n", output_file))
cat("Done!\n")


