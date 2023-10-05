# Combine CMIP6 and OISST data
		# Written by Dave S
			# June 2020

# Make available the helper functions -------------------------------------

	source("MWH_helpers.R")

# Folders -----------------------------------------------------------------

	tos <- "/Volumes/CMIP6_Regridded/tosRegrid"
	oisst <- "/Volumes/CMIP6_Regridded/oisstRegrid"

	# OISST_365.nc has 1° OISST data for 1983 - 2015
	# OISST_2015.nc has 1° data for 2015, only
	# OISSTmeanYear.nc has 1° daily data of monthly means
	# Add OISSTmeanYear.nc to CMIP6 2015 anomalies (*CMIP6_2015ANOM.nc), then add those to overall anomalies (*20150101-21001231_anom.nc)
	# Next, multiply OISST with a mask from CMIP6
	# Finally, mergetime OISST with CMIP6 to form single time series	

# Adjust CMIP6 projections to 2015 OISST --------------------------

	xx <- dir(tos, pattern = "CMIP6_2015ANOM.nc", full.names = TRUE, recursive = TRUE) # Get all instances of files for 2015
	registerDoMC(cores = detectCores()-2)
		plyr::ldply(xx, .fun = adjCMIP_OISST, .parallel = TRUE)
	

# Convert the 1983-2014 OISST data to a "classic" netCDF, with the correct variable names --------

		r <- stack("/Volumes/CMIP6_Regridded/oisstRegrid/OISST_365.nc") # Get the data
		raster2netCDF4(r, yrSel = 1983:2014, pth = "/Volumes/CMIP6_Regridded/oisstRegrid", ncName = "OISST_cmip_dummy.nc") # Write it to a netCDF format compatible with CMIP6
		d <- dir(oisst, pattern = "OISST_cmip_dummy", full.names = TRUE)
		if(ncvar_get(nc_open(d), var = "lat")[1] > 0) { # Ensure that the latitudes run in the same direction as CMIP6 data
			system(paste0("cdo invertlat ", d, " /Volumes/CMIP6_Regridded/oisstRegrid/OISST_cmip_format.nc"))
			file.remove(d)
			} else {
				file.rename(d, "/Volumes/CMIP6_Regridded/oisstRegrid/OISST_cmip_format.nc")
				}
		
		
# Merge CMIP6 and OISST data -----------------------------

	xx <- dir(tos, pattern = "_OISST_ADJUSTED.nc", full.names = TRUE, recursive = TRUE) # Do this for all OISST-ajusted CMIP6 files
	registerDoMC(cores = detectCores()-2)
		plyr::ldply(xx, .fun = mergeCMIP_OISST, .parallel = TRUE)	
		

#******* Still to DO...create a mask and resolve?
		
		
# Wrapper to merge OISST and OISST-adjusted CMIP6 data ------------------
		
		maskCMIP_OISST <- function(d) {
			pth <- unlist(strsplit(d, "/tos_Oday"))[1]
			r <- stack(d)[[1:365]]
			r[][!is.na(r[])] <- 1 # Make all non-zero values be 1 - this keeps the next step manageable
			msk <- stackApply(r, rep(1, nlayers(r)), prod, na.rm = FALSE) # Will identify where NAs are across the year
			names(msk) <- "X2015.01.01" # The next stage needs a date
			raster2netCDF4(msk, yrSel = 2015, pth = pth, ncName = "/mask.nc") # Make a masking netCDF
			f <- dir(pth, pattern = "FIN", full.names = TRUE)
			for(i in 1:length(f)) {
				system(paste0("cdo div ", f[i], " -gec,0.5 mask.nc ", unlist(strsplit(f[i], ".nc"))[1], "masked.nc")) # Mask each 118-year stack
				}
			}
# Goto 6_Break_into_Blocks.R