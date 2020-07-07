
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import packages --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(RCurl)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull ITRDB data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
wdir <- 'remote\\'



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


test <- as_tibble(filenames)
i <- 5
filenames[i]
filenames <- filenames[i:length(filenames)]

