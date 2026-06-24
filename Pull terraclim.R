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
library(tidylog)
library(ncdf4)
library(progress)


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
  distinct(plot_cn, latitude, longitude) %>%
  filter(!is.na(latitude), !is.na(longitude))

cat(sprintf("Loaded %d unique FIA sites\n", nrow(sites)))


itrdb_df <- read_csv(paste0(dendro_dir, "site_summary.csv"))
itrdb_df <- itrdb_df %>%
  select(site_id = collection_id, latitude, longitude) %>%
  mutate(source = "ITRDB") %>%
  unique()

all_sites <- sites %>%
  select(site_id = plot_cn, latitude, longitude) %>%
  mutate(source = "FIA") %>%
  rbind(itrdb_df)

all_sites %>% pull(site_id) %>% unique()

all_sites %>% write_csv(file = paste0(dendro_dir, "all_site_locations.csv"))

all_sites <- read_csv(paste0(dendro_dir, "all_site_locations.csv"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Checkpoint configuration -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Configuration parameters
checkpoint_batch_size <- 1  # Save after every N sites (1 = save after each site)
max_retries <- 3             # Number of retries for failed extractions
retry_delay <- 5             # Seconds to wait between retries

# Define checkpoint directory
checkpoint_dir <- paste0(output_dir, ".checkpoints/")
progress_file <- paste0(checkpoint_dir, "progress.txt")

# Create checkpoint directory if it doesn't exist
if (!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir, recursive = TRUE)
  cat("Created checkpoint directory:", checkpoint_dir, "\n")
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper functions for checkpointing ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
get_completed_sites <- function(progress_file) {
  if (!file.exists(progress_file)) {
    return(character(0))
  }

  # Read progress file, skip comment lines
  lines <- readLines(progress_file)
  lines <- lines[!grepl("^#", lines)]  # Remove comments

  if (length(lines) == 0) {
    return(character(0))
  }

  # Extract plot_cn values (first column before comma)
  completed <- sapply(strsplit(lines, ","), `[`, 1)
  completed <- na.omit(completed)

  return(completed)
}

record_site_completion <- function(plot_cn, progress_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("%s,%s\n", plot_cn, timestamp)

  # Append to progress file (atomic operation)
  cat(line, file = progress_file, append = TRUE)
}

save_site_checkpoint <- function(site_data, plot_cn, checkpoint_dir) {
  checkpoint_file <- paste0(checkpoint_dir, "site_", plot_cn, ".csv")
  write_csv(site_data, checkpoint_file)
}

extract_with_retry <- function(lon, lat, variable, max_retries = 3, retry_delay = 5) {
  for (attempt in 1:max_retries) {
    result <- tryCatch({
      extract_terraclimate_point(lon, lat, variable)
    }, error = function(e) {
      if (attempt < max_retries) {
        cat(sprintf("RETRY %d/%d (waiting %ds)... ", attempt, max_retries - 1, retry_delay))
        flush.console()
        Sys.sleep(retry_delay)
      }
      return(NULL)
    })

    if (!is.null(result)) {
      return(result)
    }
  }

  # All retries failed
  return(NULL)
}


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


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Connection pooling function (extract all sites for one variable) ------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extract_all_sites_for_variable <- function(sites_df, variable, max_retries = 3) {
#   # Purpose: Extract TerraClimate time series for ALL sites for a single variable
#   #          Uses connection pooling - opens ONE NetCDF connection and extracts all sites
#   # Inputs:
#   #   sites_df: Data frame with columns: plot_cn, latitude, longitude
#   #   variable: TerraClimate variable name (e.g., "tmax", "tmin", "ppt", "def", "pet")
#   #   max_retries: Number of retry attempts for opening NetCDF connection
#   # Returns:
#   #   Data frame with columns: plot_cn, year, month, [variable]

#   # Build TerraClimate THREDDS URL
#   base_url <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/"
#   nc_url <- paste0(base_url, "agg_terraclimate_", variable, "_1958_CurrentYear_GLOBE.nc")

#   # Open NetCDF connection ONCE for this variable
#   nc <- NULL
#   for (attempt in 1:max_retries) {
#     nc <- tryCatch({
#       nc_open(nc_url)
#     }, error = function(e) {
#       if (attempt < max_retries) {
#         cat(sprintf("  Connection failed, retrying (%d/%d)...\n", attempt, max_retries))
#         Sys.sleep(5)
#       }
#       return(NULL)
#     })
#     if (!is.null(nc)) break
#   }

#   if (is.null(nc)) {
#     warning(sprintf("Failed to open NetCDF for variable %s after %d attempts", variable, max_retries))
#     return(NULL)
#   }

#   on.exit(nc_close(nc))

#   # Extract coordinate arrays (once)
#   nc_lon <- ncvar_get(nc, "lon")
#   nc_lat <- ncvar_get(nc, "lat")
#   time_vals <- ncvar_get(nc, "time")

#   # Convert time to year/month
#   origin_date <- as.Date("1900-01-01")
#   dates <- origin_date + time_vals
#   years <- as.integer(format(dates, "%Y"))
#   months <- as.integer(format(dates, "%m"))

#   # Loop through sites and extract from the SAME connection
#   results_list <- vector("list", nrow(sites_df))
#   total_sites <- nrow(sites_df)

#   # Progress indicators every 10% or every 50 sites, whichever is smaller
#   progress_interval <- min(50, max(1, floor(total_sites / 10)))

#   for (i in 1:nrow(sites_df)) {
#     site <- sites_df[i, ]

#     # Show progress periodically
#     if (i %% progress_interval == 0 || i == total_sites) {
#       cat(sprintf("  Processing site %d/%d (%.1f%%)...\n",
#                   i, total_sites, 100 * i / total_sites))
#       flush.console()
#     }

#     # Find closest grid cell
#     lon_diff <- abs(nc_lon - site$longitude)
#     lat_diff <- abs(nc_lat - site$latitude)
#     lon_idx <- which.min(lon_diff)
#     lat_idx <- which.min(lat_diff)

#     # Check tolerance
#     if (lon_diff[lon_idx] > 1/48 || lat_diff[lat_idx] > 1/48) {
#       warning(sprintf("No grid cell found for site %s", site$plot_cn))
#       next
#     }

#     # Extract time series
#     time_dim <- nc$dim$time$len
#     values <- ncvar_get(nc, variable,
#                         start = c(lon_idx, lat_idx, 1),
#                         count = c(1, 1, time_dim))

#     # Create result
#     results_list[[i]] <- tibble(
#       plot_cn = site$plot_cn,
#       year = years,
#       month = months,
#       !!variable := as.numeric(values)
#     )
#   }

#   # Combine all sites for this variable
#   bind_rows(results_list)
# }


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Extract TerraClimate data for all sites -------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read completed sites
completed_sites <- get_completed_sites(progress_file)
cat(sprintf("Found %d already completed sites\n", length(completed_sites)))

# Filter to only unprocessed sites
sites_to_process <- sites %>%
  filter(!plot_cn %in% completed_sites) 
# %>%
#   head(10)  # Test with just 50 sites (remove this line for full run)

cat(sprintf("Processing %d remaining sites (out of %d total)\n",
            nrow(sites_to_process), nrow(sites)))
flush.console()

if (nrow(sites_to_process) == 0) {
  cat("All sites already completed! Skipping to merge step.\n")
  flush.console()
} else {
  # Define variables to extract
  variables <- c("tmax", "tmin", "ppt", "def", "pet")

  cat("\n=== Starting TerraClimate extraction ===\n")
  cat(sprintf("Processing %d sites × %d variables = %d total extractions\n",
              nrow(sites_to_process), length(variables),
              nrow(sites_to_process) * length(variables)))
  cat("Strategy: Extract all variables for each site, then save checkpoint immediately\n")
  cat("This allows resuming from any point if interrupted\n\n")
  flush.console()

  # Progress tracking
  total_sites <- nrow(sites_to_process)
  sites_processed <- 0
  failed_sites <- 0

  # OUTER LOOP: Process each site (one at a time)
  for (i in 1:nrow(sites_to_process)) {
    site <- sites_to_process[i, ]

    # Progress indicator
    cat(sprintf("[%d/%d] Processing site %s (%.1f%% complete)\n",
                i, total_sites, site$plot_cn,
                100 * i / total_sites))
    flush.console()

    # Extract all variables for this site
    site_data <- NULL
    variables_extracted <- 0

    for (var in variables) {
      cat(sprintf("  - Extracting %s... ", var))
      flush.console()

      # Extract time series for this variable with error handling
      var_data <- tryCatch({
        extract_terraclimate_point(
          lon = site$longitude,
          lat = site$latitude,
          variable = var
        )
      }, error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        flush.console()
        return(NULL)
      })

      if (is.null(var_data)) {
        cat("FAILED\n")
        flush.console()
        next
      }

      cat("OK\n")
      flush.console()

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

      variables_extracted <- variables_extracted + 1
    }

    # Save checkpoint for this site immediately (if we got any data)
    if (!is.null(site_data) && variables_extracted > 0) {
      # Add plot_cn as the identifier
      site_data <- site_data %>%
        mutate(plot_cn = site$plot_cn, .before = year)

      # Save checkpoint
      save_site_checkpoint(site_data, site$plot_cn, checkpoint_dir)
      record_site_completion(site$plot_cn, progress_file)

      sites_processed <- sites_processed + 1
      cat(sprintf("  ✓ Checkpoint saved (%d/%d variables extracted)\n",
                  variables_extracted, length(variables)))
      flush.console()
    } else {
      failed_sites <- failed_sites + 1
      cat(sprintf("  ✗ Failed to extract any data for this site\n"))
      flush.console()
    }

    # Periodic summary every 10 sites
    if (i %% 10 == 0 || i == total_sites) {
      cat(sprintf("\n--- Progress Summary ---\n"))
      cat(sprintf("  Sites processed: %d/%d (%.1f%%)\n",
                  i, total_sites, 100 * i / total_sites))
      cat(sprintf("  Successful: %d, Failed: %d\n",
                  sites_processed, failed_sites))
      cat(sprintf("  Checkpoints saved: %d\n\n", sites_processed))
      flush.console()
    }
  }

  cat(sprintf("\n=== Extraction Complete ===\n"))
  cat(sprintf("Successfully processed: %d/%d sites\n",
              sites_processed, total_sites))
  cat(sprintf("Failed: %d sites\n", failed_sites))
  cat(sprintf("Checkpoints saved in: %s\n", checkpoint_dir))
  flush.console()
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Merge checkpoint files into final output ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Merging checkpoint files ===\n")
flush.console()

# Get all checkpoint CSV files
checkpoint_files <- list.files(checkpoint_dir,
                                pattern = "^site_.*\\.csv$",
                                full.names = TRUE)

cat(sprintf("Found %d checkpoint files to merge\n", length(checkpoint_files)))
flush.console()

if (length(checkpoint_files) == 0) {
  stop("ERROR: No checkpoint files found!\n",
       "This means either:\n",
       "  1. The extraction failed before creating checkpoints\n",
       "  2. All sites were skipped (already completed) but no previous checkpoints exist\n",
       "Check the output above for error messages.")
}

# Read and combine all checkpoint files
climate_df <- checkpoint_files %>%
  map_df(read_csv, show_col_types = FALSE)

cat(sprintf("Combined data: %d rows, %d unique sites\n",
            nrow(climate_df),
            length(unique(climate_df$plot_cn))))

# Verify output structure
cat("\nOutput structure:\n")
cat(sprintf("  Columns: %s\n", paste(names(climate_df), collapse = ", ")))
cat(sprintf("  Date range: %d-%02d to %d-%02d\n",
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

# Save final CSV
output_file <- paste0(output_dir, "site_climate_fia.csv")
write_csv(climate_df, output_file)

cat(sprintf("\nFinal data saved to: %s\n", output_file))

# Cleanup message
cat("\n=== Cleanup ===\n")
cat("Checkpoint files preserved in:", checkpoint_dir, "\n")
cat("To clean up checkpoints and start fresh, delete the .checkpoints directory\n")
cat("Done!\n")


