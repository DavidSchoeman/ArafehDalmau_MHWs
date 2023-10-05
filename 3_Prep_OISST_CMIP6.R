# Prep OISST for combining with CMIP6 data
		# Written by Dave S
			# June 2020

# Make available the helper functions -------------------------------------

	source("MWH_helpers.R")

# Folders -----------------------------------------------------------------

	oisst <- "/Volumes/CMIP6_Regridded/oisstRegrid"


# Make daily netCDF of monthly means for OISST 2015 --------
	
	iFile <- paste0(oisst, "/OISST_365.nc")
	oFile1 <- paste0(unlist(strsplit(iFile, "_365"))[1], "_2015.nc")
	dO <- paste("cdo selyear,2015/2015", iFile, oFile1, sep = " ")
		system(dO)
	# Get monthly means
		dO <- paste0("cdo monmean OISST2015_1.nc ", oisst, "/OISST2015_monmean.nc") 
			system(dO)
	# Break the netCDF into monthly pieces
		paste0("cdo splitmon OISST2015_monmean.nc OISST2015_mon.nc")
		dO <- paste0("cdo splitmon ", oisst, "/OISST2015_monmean.nc ", oisst, "/OISST2015_mon.nc")
			system(dO)
	# Make each of these into a monthly stack
		m <- dir(oisst, pattern = "OISST2015_mon.nc", full.names = TRUE)
		for(i in 1:length(m)) {
			r <- stack(m[i])
			if(substr(names(r), 10, 11) == 16) {
				for(j in 1:15) {
					system(paste0("cdo shifttime,-", j,"days ", m[i], " dayOut-", j,".nc"))
					system(paste0("cdo shifttime,+", j,"days ", m[i], " dayOut+", j,".nc"))
				}
			} else {
				if(substr(names(r), 10, 11) == 15) {
					for(j in 1:14) {
						system(paste0("cdo shifttime,-", j,"days ", m[i], " dayOut-", j,".nc"))
						system(paste0("cdo shifttime,+", j,"days ", m[i], " dayOut+", j,".nc"))
					}
					system(paste0("cdo shifttime,+", j+1,"days ", m[i], " dayOut+", j+1,".nc"))
				} else {
					for(j in 1:13) {
						system(paste0("cdo shifttime,-", j,"days ", m[i], " dayOut-", j,".nc"))
						system(paste0("cdo shifttime,+", j,"days ", m[i], " dayOut+", j,".nc"))
					}
					system(paste0("cdo shifttime,+", j+1,"days ", m[i], " dayOut+", j+1,".nc"))
				}
			}	
			file.copy(m[i], "dayOut0.nc") # Make sure that the original file is included
			dO <- paste0("cdo mergetime dayOut*.nc ", oisst, "/mOISST_", unlist(strsplit(m[i], ".nc"))[2], ".nc")
				system(dO) # Combine them all together again
			file.remove(dir(pattern = "dayOut"))
		}
		dO <- paste0("cdo mergetime ", oisst, "/mOISST_*.nc ", oisst, "/OISSTmeanYear.nc") # Combine into a year
			system(dO)
		file.remove(dir(oisst, pattern = "mOISST_", full.names = TRUE)) # Clean up


