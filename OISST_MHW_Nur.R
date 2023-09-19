# Breaking combined data into blocks for processing
	# Written by Dave S
		# November 2020

# Depends on OISST data

# Make available the helper functions -------------------------------------
	source("Nur_MWH_helpers.R")


# Folders -----------------------------------------------------------------
	oFold <- "/Volumes/CMIP6_Regridded/NurPts" # A temporary output path...this is going to generate a LOT of files!
	iFold <- "/Volumes/CMIP6_MHWs"
	
# Read coordinates, find extent, square it out, then extract points -------------

	dat <- read.csv("Sites_MHWs.csv") # Data provided by Nur on 17 November 2020
	r <- rasterFromXYZ(dat[,c(3, 2, 4)], crs = projection(raster()))
	r <- crop(raster(), extent(r) + 2)
		r[] <- 1:ncell(r)
	E <- round(extent(r) + 2) # Extent with a little buffer
	toDo <- dat[,c(3, 2)] %>% 
		rename(x = Long, y = Lat)
	
# Get time series for relevant pixels --------------------------------------------------------------

	# Crop then merge the OISST data
		ncs <- dir("/Volumes/Data/OISST", pattern = "avhrr-only-v2.", full.names = TRUE)
		registerDoMC(cores = detectCores()-2)	
			plyr::ldply(ncs, .fun = clpOISST, .parallel = TRUE)
		ncs <- dir("/Volumes/Data/OISST", pattern = "Clipped", full.names = TRUE)
		for(j in 198:201) {
			nn <- ncs[grep(paste0("avhrr-only-v2.", j), ncs)]
			dO <- paste0("cdo mergetime ", paste(nn, collapse = " "), " /Volumes/Data/OISST_Combo/OISST", j, "0s.nc")
				system(dO)
			}
		nnn <- dir("/Volumes/Data/OISST_Combo", pattern = "OISST", full.names = TRUE)
		dO <- paste0("cdo mergetime ", paste(nnn, collapse = " "), " /Volumes/Data/OISST_Combo/OISST_Cropped_Combo.nc")
			system(dO)
		file.remove(nnn)
		file.remove(ncs)

	# For each subset netCDF, determine the relevant pixels and extract time series
		getOPix("/Volumes/Data/OISST_Combo/OISST_Cropped_Combo.nc")
		
		
# Go back to pixels and compute MHW stats ---------------------------------

		SST_files <- dir("/Volumes/Data/OISST_Combo", pattern = ".csv", full.names = TRUE, recursive = TRUE) # An output path...this is going to generate a LOT of files!
			SST_files <- SST_files[-grep("XY.csv", SST_files)] # Just the pixels
		nc <- nc_open("/Volumes/Data/OISST_Combo/OISST_Cropped_Combo.nc")
		Dates <- data.frame(Date = as.Date(nc$dim$time$vals, origin = "1978-01-01")) # Need it as a data frame
		registerDoMC(cores = detectCores()-2)	
			out <- plyr::ldply(SST_files, .fun = purrr_O_wraper, .parallel = TRUE)
		saveRDS(out, file = "OISST_Events.Rda")

# Compute summary stats ---------------------------------------------------

		mhw_O_sum("OISST_Events.Rda")
			

# Recompile to original cells and average ----------------------------------------------------------
	
	d <- read.csv("/Volumes/Data/OISST_Combo/XY.csv") %>% 
		bind_cols(dat) %>% 
		select(x, y, Pixel)		
	readRDS("/Users/davidschoeman/Dropbox (GCERG)/Documents/Student Documents/Ongoing Students/Nur/MHW_Nur/OISST_Stats.Rda") %>% 
		filter(year == 2014 | year == 2015) %>% 
		select(x, y, year, Cum_Int) %>% 
		pivot_wider(names_from = year, values_from = Cum_Int) %>% 
		left_join(d) %>% 
		write.csv("Recomputed_OISST_Cum_Int.csv", row.names = FALSE)
		
		
		
