# Regrid CMIP6 tos 
	# Written by Dave Schoeman
		# April 2020


# Make available the helper functions -------------------------------------

	source("MWH_helpers.R")


# Set up variables - change stuff here, if needed -------------------------

	# For projections
		dRaw <- "/Volumes/CMIP6_Raw/tosRaw"
			dRawP <- "/Volumes/CMIP6_Raw/tosRaw_processed"
		dRegrid <- "/Volumes/CMIP6_Regridded/tosRegrid"
		cmipVar <- "tos"
		cmipFreq <- "Oday"
	
	# For historical
		# dRaw <- "/Volumes/CMIP6_Raw/hist_dump_processed"
		# dRawP <- "/Volumes/CMIP6_Raw/tosRaw_hist"
		# dRegrid <- "/Volumes/CMIP6_Regridded/histRegrid"
		# cmipVar <- "tos"
		# cmipFreq <- "Oday"	
	
# Get details of models and SSPs --------
	
	# ...only for 365/366-day calendars, and deleting 29 Febs, where they exist
	dFiles <- dir(dRaw, recursive = TRUE, full.names = TRUE, pattern = ".nc") # What files am I working with?
		fDat <- unlist(lapply(dFiles, function(x) unlist(strsplit(x, paste0(cmipVar, "_", cmipFreq, "_")))[2]))
		mdls <- unlist(lapply(fDat, function(x) unlist(strsplit(x, "_"))[1]))
		ssps <- unlist(lapply(fDat, function(x) unlist(strsplit(x, "_"))[2]))
			ms <- unique(paste0(mdls, "_", ssps, "_"))
		fList <- lapply(ms, function(x) {grep(x, dFiles)})
		fList <- data.frame(mod = mdls, ssp = ssps) %>% 
			distinct()
		fList <- unique(paste0(mdls, "_", ssps, "_"))
	# Outputs a list of model/ssps to be processed

	
# Remap the raw files in parallel -----------------------------------------

	# For files starting 20150101 and ending 21001231, and ALSO having a 365-day or Gregorian calendar, remap to standard 1° global grid and 365-day calendar
	system.time({
		registerDoMC(cores = detectCores()-2)
		plyr::ldply(fList, .fun = go_remap, .parallel = TRUE)
		})
	# Output for each model/SSP, regridded combined netCDFs, ending with *_365.nc

	
# Convert remapped netCDFs to anomailes in parallel - NOT needed for "historical" -------------------------------------

	system.time({
		dFiles <- dir(dRegrid, recursive = TRUE, full.names = TRUE, pattern = "_365.nc") # What SSP folders do I have?
		registerDoMC(cores = detectCores()-2)	
			plyr::ldply(dFiles, .fun = mk_anom, .parallel = TRUE) # Regrids netCDFs in parallel
		})		
	# Output for each model/SSP, *_2015.nc (netCDF of SST for just 2015) and *_anom.nc (netCDF of daily anomalies relative to days of 2015)

# Goto 2_Regrid_OISST.R
