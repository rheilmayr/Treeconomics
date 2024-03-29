#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/7/2020
# Purpose: Iterate through downloading ITRDB data from NOAA ftp server
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import packages --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(RCurl)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull ITRDB data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'

## Note - Downloaded on July 5, 2020; ITRDB v.7.22
continents <- list('northamerica/mexico', 'northamerica/canada', 'northamerica/usa')
for (continent in continents) {
  url = "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/measurements/"
  url <- paste0(url, continent, '/')
  filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames <- strsplit(filenames, "\r\n")
  filenames = unlist(filenames)
  for (filename in filenames) {
    download.file(paste(url, filename, sep = ""), paste(wdir, "raw_in/itrdb/rwi/", filename,
                                                  sep = ""))
    Sys.sleep(2)
  }
}

