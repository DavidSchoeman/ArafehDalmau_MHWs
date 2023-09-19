# Breaking combined data into blocks for processing
	# Written by Dave S
		# November 2020

# Depends on initial processing from /Users/davidschoeman/Dropbox (GCERG)/Documents/Research_Projects/Multivariate_MPAs/MV_MPA_Code

# Make available the helper functions -------------------------------------
	source("Nur_MWH_helpers.R")


# Folders -----------------------------------------------------------------
	oFold <- "/Volumes/CMIP6_Regridded/Nur" # A temporary output path...this is going to generate a LOT of files!
	iFold <- "/Volumes/CMIP6_MHWs"
	
# Read shape, find extent, square it out, then extract points -------------

	br <- raster() # A base raster of 1° resolution
	pa <- st_read("/Users/davidschoeman/Dropbox (GCERG)/Documents/Student Documents/Ongoing Students/Nur/MHW Planning Paper/SCB_Plannig_12_11_20/SCB_Planning_02.shp") %>% 
		st_transform(pa, crs = raster::projection(br)) %>%  # Project lonlat
		st_zm() # Remove z-dimension
	e <- round(extent(pa) + 2) # Extent with a little buffer
	pa1 <- rasterize(pa, raster::disaggregate(crop(br, e), 10)) # A finer raster version of pa
	# Reaggregate the finer raster to 1°, then extract points that correspond with MPA overlap
	toDoPts <- raster::aggregate(pa1, 10, sum) %>% 
		rasterToPoints() %>% 
		data.frame() # Points including weighting by cell
	toDo <- toDoPts %>% 
		dplyr::select(-layer) # Points without weighting
	
	
# Get pixels --------------------------------------------------------------

	# # Subset for just the years you want, and the area you want
		f <- dir(iFold, pattern = "18500101-21001231.nc", recursive = TRUE, full.names = TRUE)	
		registerDoMC(cores = detectCores()-2)	
			plyr::ldply(f, .fun = ptPixel, .parallel = TRUE)


# Go back to pixels and compute MHW stats ---------------------------------

	o1 <- dir(oFold, full.names = TRUE) # All subfolders
	o2 <- unlist(lapply(o1, function(x) {dir(x, full.names = TRUE)})) # All subfolders within those
	for(o in o2) { # Work one folder at a time to avoid crashing
		SST_files <- dir(o, pattern = ".csv", full.names = TRUE, recursive = TRUE) # An output path...this is going to generate a LOT of files!
		d <- as.Date("1983/01/01"):as.Date("2100/12/31") # Make a dat string****change the dates here, if needed
		d <- as.Date(d, origin = "1970/01/01") # Write them as dates
		d <- d[-grep("-02-29", d)] # Remove leap days to match format of netCDFs
		nd <- length(d) # How many dates are there? Need this as a check to kill cells that don't have past AND future values
		Dates <- data.frame(Date = d) # Need it as a data frame
		registerDoMC(cores = detectCores()-2)	
			out <- plyr::ldply(SST_files, .fun = purrr_wraper, .parallel = TRUE)
		nm <-unlist(strsplit(o, "/ssp"))[2] # Make a name for the output
		nm <- unlist(strsplit(nm, "/"))
		nm <- paste0("ssp", nm[1], "_", nm[2])
		eval(parse(text = paste0("saveRDS(out, file = '", o, "/MHWout_", nm, ".Rda')")))
		rm(out, SST_files)
	}	
			

# Compute summary stats ---------------------------------------------------

	f <- dir(oFold, pattern = "MHWout_", recursive = TRUE, full.names = TRUE)
	registerDoMC(cores = detectCores()-10) # Reduce number of cores to prevent memory overload	
		plyr::ldply(f, .fun = mhw_sum, .parallel = TRUE)			
			

# Compute trends ----------------------------------------------------------

		# Loop through files
		f <- dir(oFold, pattern = "MHWstatsTS_", recursive = TRUE, full.names = TRUE)	
		# f <- f[-grep("EC-Earth3", f)] # EC-Earth3 crashes the routine at the moment *** CHECK AND CORRECT DOWN THE LINE
		for(i in f) {
			cat(paste("\nStarting...", i, "\n"))
			eval(parse(text = paste0("event_sum <- readRDS('", i, "')")))		
			# Run the functions using plyr
			registerDoMC(cores = detectCores()-2)
			maxYear <- plyr::ddply(event_sum, c("x", "y"), mxYr, .parallel = TRUE) %>% 
				dplyr::rename(yrPerm = V1)  
			nTrend <- plyr::ddply(event_sum, c("x", "y"), lin_n, .parallel = TRUE) %>%
				dplyr::rename(nslope = slope, np = p)
			iTrend <- plyr::ddply(event_sum, c("x", "y"), lin_int, .parallel = TRUE) %>%
				dplyr::rename(islope = slope, ip = p)
			ciTrend <- plyr::ddply(event_sum, c("x", "y"), lin_cumi, .parallel = TRUE) %>%
				dplyr::rename(cislope = slope, cip = p)
			dTrend <- plyr::ddply(event_sum, c("x", "y"), lin_dur, .parallel = TRUE) %>%
				dplyr::rename(dslope = slope, dp = p)
			cdTrend <- plyr::ddply(event_sum, c("x", "y"), lin_cdur, .parallel = TRUE) %>%
				dplyr::rename(cdslope = slope, cdp = p)
			# Join, then rasterise
			out <- left_join(nTrend, iTrend) %>%
				left_join(., ciTrend) %>%
				left_join(., dTrend) %>%
				left_join(., cdTrend) %>%
				left_join(., maxYear)
			nm <- unlist(strsplit(i, "/MHWstatsTS")) # Start to create a name for the output
			nm <- paste0(nm[1], "/HMWtrendDF", nm[2]) # Make a name for the output
			eval(parse(text = paste0("saveRDS(out, file = '", nm, "')")))
			cat(paste("\n", i, "Completed\n"))
		}	
		
		# Write data frame to raster and save
		f <- dir(oFold, pattern = "HMWtrendDF_", recursive = TRUE, full.names = TRUE)	
		for(i in f) {
			cat(paste("\nStarting...", i, "\n"))
			eval(parse(text = paste0("df <- readRDS('", i, "')")))		
			nm <- gsub("HMWtrendDF", "HMWtrendRasterStack", i) # Make a name for the output
			eval(parse(text = paste0("r <- rasterFromXYZ(df, crs = projection(raster()))")))
			eval(parse(text = paste0("saveRDS(r, file = '", nm, "')")))
			cat(paste("\n", i, "Completed\n"))
		}

		
		
		
		
		
		
		