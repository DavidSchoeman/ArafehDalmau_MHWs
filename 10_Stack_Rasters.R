# MHW trends for global analysis
	# Written by Dave S
		# July 2020

# Note that some models failed to yield useful estimates, so we will ignore those

# Make available the helper functions
	source("MWH_helpers.R")

# Loop through files
	oFold <- "/Volumes/CMIP6_MHW_Data/mhwRasters/"
	iFold <- "/Volumes/CMIP6_MHW_Data/temp_csv" # A temporary output path...this is going to generate a LOT of files!
	f <- dir(iFold, pattern = "HMWtrendRasterStack", recursive = TRUE, full.names = TRUE)	
		# f <- f[-grep("EC-Earth3|CNRM", f)] # EC-Earth3 crashes the routine at the moment *** CHECK AND CORRECT DOWN THE LINE
	ssps <- unique(unlist(lapply(f, getSSP)))
	for(i in ssps) {
		ff <- f[grep(i, f)] # Just the files for ssp i
		mods <- unique(unlist(lapply(ff, getBlockMod)))
		w <- raster() # A 1° world raster to act as a template
		for(k in c("cislope", "cdslope", "yrPerm"))	{ # Cumulative impact, cumulative 
			l <- which(names(raster::stack(readRDS(ff[grep(j, ff)]))) == k)
			for(j in mods) {
				if(j == mods[1]) {
					r <- extend(raster::stack(readRDS(ff[grep(j, ff)]))[[l]], w)
					} else {
						s <- extend(raster::stack(readRDS(ff[grep(j, ff)]))[[l]], w)
						r <- stack(r, s)
						}
				}
			names(r) <- mods # Name the layers
			# medianVal <- calc(r, median) # I removed na.rm = TRUE in order to ensure that all cells were computed with the same number of reps
			# minVal <- calc(r, min)
			# maxVal <- calc(r, max)
			# nMods <- calc(r, function(x) {sum(!is.na(x))}) # How many non-NAs do I have?
			medianVal <- calc(r, median, na.rm = TRUE) # I removed na.rm = TRUE in order to ensure that all cells were computed with the same number of reps
			minVal <- calc(r, min, na.rm = TRUE)
			maxVal <- calc(r, max, na.rm = TRUE)
			nMods <- calc(r, function(x) {sum(!is.na(x))}) # How many non-NAs do I have?
			rr <- raster::stack(medianVal, minVal, maxVal, nMods)
			names(rr) <- c("medianVal", "minVal", "maxVal", "nMods")
			out <- stack(rr, r)
			nm <- paste0(oFold, i, "_", k, ".Rda")
			eval(parse(text = paste0("saveRDS(out, file = '", nm, "')")))
			cat(paste("\n", nm, "Completed\n"))
			}
		}

# Data pre-processing complete	