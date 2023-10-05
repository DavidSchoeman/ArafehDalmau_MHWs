# Prep CMIP6 data for combining with older data
		# Written by Dave S
			# June 2020

# Make available the helper functions -------------------------------------

	source("MWH_helpers.R")

# Folders -----------------------------------------------------------------

	tos <- "/Volumes/CMIP6_Regridded/tosRegrid"

# Compute daily anomalies relative to monthly means for CMIP6 2015 --------		
	# Do this in parallel	
	xx <- dir(tos, pattern = "_2015.nc", full.names = TRUE, recursive = TRUE) # Get all instances of files for 2015
		registerDoMC(cores = detectCores()-2)
			plyr::ldply(xx, .fun = pllYr, .parallel = TRUE)
		
# Goto 5.1_Combine_CMIP_OISST or 5.2_Combine_CMIP_hist, depending on what needs doing