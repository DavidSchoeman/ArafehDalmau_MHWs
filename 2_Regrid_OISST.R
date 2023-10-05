# Regrid CMIP6 data
	# Written by Dave S
		# May 2020

# Make available the helper functions -------------------------------------

	source("MWH_helpers.R")


# Folders -----------------------------------------------------------------

	oisst_d <- "/Volumes/CMIP6_Raw/OISST_day"
	oisst_y <- "/Volumes/CMIP6_Raw/OISST_yr"
	dRegrid <- "/Volumes/CMIP6_Regridded/oisstRegrid"


# Process OISST into the same format as tos -------------------------------

	# First, combine daily OISST data into years - otherwise it falls over
		yrs <- 1983:2015
		registerDoMC(cores = detectCores()-2)
			plyr::ldply(yrs, .fun = mOISST, .parallel = TRUE)
	# Next, merge all annual files into a single netCDF
		dO <- paste0("cdo mergetime ", oisst_y,"/OISST* OISST_05.nc")
			system(dO) # Merge years
	# Next, reproject to same grid and clean up
		dO <- paste0("cdo remapdis,global_1 OISST_05.nc ", dRegrid, "/OISST.nc")
			system(dO) # Distance-weighted mean interpolation
		file.remove("OISST_05.nc") # Clean up
	# Fially, make OISST into a 365-day calendar
		dO <- paste0("cdo delete,month=2,day=29 ", dRegrid, "/OISST.nc ", dRegrid, "/OISST_365.nc")
			system(dO)

# Goto 3_Prep_OISST_CMIP6.R