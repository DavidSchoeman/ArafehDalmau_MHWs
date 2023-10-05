# Breaking combined data into blocks for processing
	# Written by Dave S
		# June 2020

# Make available the helper functions
	source("MWH_helpers.R")

# Folders
	oFold <- "/Volumes/CMIP6_MHW_Data/temp_csv" # A temporary output path...this is going to generate a LOT of files!
	iFold <- "/Volumes/CMIP6_MHWs"
	mFold <- "/Volumes/CMIP6_Regridded/MHWncs"

# Subset for just the years you want, then break world into 6 big blocks
	f <- dir(iFold, pattern = "18500101-21001231.nc", recursive = TRUE, full.names = TRUE)	
	registerDoMC(cores = detectCores()-2)	
		plyr::ldply(f, .fun = bWorld1, .parallel = TRUE)
	
# Break each big block into smaller blocks - doing this in two steps uses multiple cores more efficiently - I think
	ff <- dir(oFold, pattern = "BBLOCK_", recursive = TRUE, full.names = TRUE)	
		registerDoMC(cores = detectCores()-2)	
		plyr::ldply(ff, .fun = bWorld2, .parallel = TRUE)
	
# For each masked netCDF, extract the boxes and write them to the appropriate folder
	fff <- dir(oFold, pattern = "BLOCK_", full.names = TRUE, recursive = TRUE)
	mods <- unique(unlist(lapply(fff, getBlockMod)))
	for(l in mods) {
		ffff <- fff[grepl(l, fff)]
		nFile <- dir(mFold, pattern = l, recursive = TRUE, full.names = TRUE)[1]
		dO <- paste0("cdo timmean ", nFile, " Mask.nc")
			system(dO)
		r <- stack("Mask.nc")
		pts <- as.data.frame(rasterToPoints(r)) # # First get just the cells we want
		names(pts)[3] <- "sst"
		pts <- na.omit(pts)[,1:2] # This gives us coordinates where there are ssts
		registerDoMC(cores = detectCores()-2)	
			plyr::ldply(ffff, .fun = pPixel, .parallel = TRUE)
		file.remove("Mask.nc") # Remove the mask
		file.remove(ffff) # Remove the blocks
		}

# Goto 7_Extract_SST_and_Compute_Events.R
	