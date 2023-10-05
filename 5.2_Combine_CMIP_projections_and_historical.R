# Combine CMIP6 and OISST data
		# Written by Dave S
			# June 2020

# Make available the helper functions -------------------------------------

source("/Users/davidschoeman/Dropbox/Documents/Research_Projects/Multivariate_MPAs/MV_MPA_Code/MWH_helpers.R")

# Folders -----------------------------------------------------------------

	tos <- "/Volumes/CMIP6_Regridded/tosRegrid"
	hist <- "/Volumes/CMIP6_Regridded/histRegrid"
	mhw <- "/Volumes/CMIP6_MHWs"


# Combine CMIP6 historical data and projections --------------------------

	hh <- dir(hist, pattern = "tos_", full.names = TRUE, recursive = TRUE) # Get all instances of files for 2015
	registerDoMC(cores = detectCores()-2)
		plyr::ldply(hh, .fun = spliceCMIP, .parallel = TRUE)

		
