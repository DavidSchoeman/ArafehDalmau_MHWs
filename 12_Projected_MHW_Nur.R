# Similar to OISST, but computing projected cumulative MHW stats 
	# Written by Dave S
		# November 2020

# Depends on initial processing from /Users/davidschoeman/Dropbox (GCERG)/Documents/Research_Projects/Multivariate_MPAs/MV_MPA_Code

# Make available the helper functions -------------------------------------
	source("Nur_MWH_helpers.R")


# Folders -----------------------------------------------------------------
	oFold <- "/Volumes/CMIP6_Regridded/NurPts" # A temporary output path...this is going to generate a LOT of files!
	iFold <- "/Volumes/CMIP6_MHWs"
	
# Read coordinates, find extent, square it out, then extract points -------------

	dat <- read.csv("Sites_MHWs.csv") # Data provided by Nur on 17 November 2020
	# r <- rasterFromXYZ(dat[,c(3, 2, 4)], crs = projection(raster()))
	# r <- crop(raster(), extent(r) + 2)
	# 	r[] <- 1:ncell(r)
	# E <- round(extent(r) + 2) # Extent with a little buffer
	# toDo <- xyFromCell(r, raster::extract(r, dat[,c(3,2)])) %>% 
	# 	as.data.frame()
	toDo <- dat[,c(3, 2)] %>% 
		rename(x = Long, y = Lat)
	
# Get time series for relevant pixels --------------------------------------------------------------

	# Subset for just the years you want, and the area you want
		f <- dir(iFold, pattern = "18500101-21001231.nc", recursive = TRUE, full.names = TRUE)	
		registerDoMC(cores = detectCores()-2)	
			plyr::ldply(f, .fun = getExt, .parallel = TRUE)
		file.remove(dir(oFold, pattern = "temp.nc", full.names = TRUE, recursive = TRUE))

	# For each subset netCDF, determine the relevant pixels and extract time series
		ff <- dir(oFold, pattern = "_19830101-21001231.nc", recursive = TRUE, full.names = TRUE)	
		registerDoMC(cores = detectCores()-2)	
			plyr::ldply(ff, .fun = getPix, .parallel = TRUE)
		
		
# Go back to pixels and compute MHW stats ---------------------------------

	o1 <- dir(paste0(oFold, "/pts"), full.names = TRUE) # All subfolders
	o2 <- unlist(lapply(o1, function(x) {dir(x, full.names = TRUE)})) # All subfolders within those
	for(o in o2) { # Work one folder at a time to avoid crashing
		SST_files <- dir(o, pattern = ".csv", full.names = TRUE, recursive = TRUE) # An output path...this is going to generate a LOT of files!
			SST_files <- SST_files[-grep("XY.csv", SST_files)] # Just the pixels
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
			

# Recompile to original cells and average ----------------------------------------------------------
	
	o1 <- dir(paste0(oFold, "/pts"), full.names = TRUE) # All subfolders
	o2 <- unlist(lapply(o1, function(x) {dir(x, full.names = TRUE)})) # All subfolders within those
	for(o in o2) {
		sXY <- dir(o, pattern = "XY.csv", full.names = TRUE)
		d <- read.csv(sXY) %>% 
			bind_cols(dat) %>% 
			select(x, y, Pixel)
		nm <- dir(o, pattern = "MHWstatsTS_", full.names = TRUE)
		sDat <- readRDS(nm) %>%
			dplyr::filter(year > 2020) %>% 
			left_join(d) %>% 
			arrange(Pixel)
		eval(parse(text = paste0("saveRDS(sDat, file = '", nm, "')")))
		}
	for(i in o1) {
		f <- dir(i, pattern = "MHWstatsTS_", recursive = TRUE, full.names = TRUE)
		d <- purrr::map_dfr(f, readRDS) %>%  # Read and append the RDAs
			group_by(Pixel, year) %>% 
			summarise(mdCumIntensity = median(Cum_Int)) %>% 
			as.data.frame()
		ssp <- unlist(strsplit(i, "/pts/"))[2]
		eval(parse(text = paste0("write.csv(d, '", ssp, "_CumMHWInt.csv')")))
		}


# Messing with plots ------------------------------------------------------

	ggplot(read.csv("/Users/davidschoeman/Dropbox (GCERG)/Documents/Student Documents/Ongoing Students/Nur/MHW_Nur/ssp126_CumMHWInt.csv"), aes(x = year, y = mdCumIntensity, colour = Pixel)) +
		geom_line() +
		geom_hline(aes(yintercept = 600))
