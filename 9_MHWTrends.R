# MHW trends for global analysis
	# Written by Dave S
		# March 2020
			# Modified July 2020

# Note that instead of simple linear regression, we instead use a Siegel repeated medians, which is far more robust to autocorrelation
	# "Siegel repeated medians is more complicated. For each point, the slopes between it and the others are calcuated (resulting n-1 slopes) and the median is taken. 
	# This results in n medians and median from this medians is slope estimator. Intercept is calculated in similar way, for more information please take a look in function source.
	# The breakdown point of Theil-Sen method is about 29%, Siegel extended it to 50%, so these regression methods are very robust. 
	# Additionally, if the errors are normally distributed and no outliers are present, the estimators are very similar to classic least squares.
	
# Make available the helper functions
	source("MWH_helpers.R")

# Loop through files
	iFold <- "/Volumes/CMIP6_MHW_Data/temp_csv" # A temporary output path...this is going to generate a LOT of files!
		f <- dir(iFold, pattern = "MHWstatsTS", recursive = TRUE, full.names = TRUE)	
	f <- f[-grep("EC-Earth3", f)] # EC-Earth3 crashes the routine at the moment *** CHECK AND CORRECT DOWN THE LINE
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
	f <- dir(iFold, pattern = "HMWtrendDF", recursive = TRUE, full.names = TRUE)	
	for(i in f) {
		cat(paste("\nStarting...", i, "\n"))
		eval(parse(text = paste0("df <- readRDS('", i, "')")))		
		nm <- gsub("HMWtrendDF", "HMWtrendRasterStack", i) # Make a name for the output
		eval(parse(text = paste0("r <- rasterFromXYZ(df, crs = projection(raster()))")))
		eval(parse(text = paste0("saveRDS(r, file = '", nm, "')")))
		cat(paste("\n", i, "Completed\n"))
		}
		
# Goto MergeRasters.R

ii <- 1

i <- f[ii]
cat(paste("\nStarting...", i, "\n"))
eval(parse(text = paste0("df <- readRDS('", i, "')")))		
nm <- gsub("HMWtrendDF", "HMWtrendRasterStack", i) # Make a name for the output
eval(parse(text = paste0("r <- rasterFromXYZ(df, crs = projection(raster()))")))
ii <- ii + 1
		
		# Loop through files
		iFold <- "/Volumes/CMIP6_MHW_Data/temp_csv" # A temporary output path...this is going to generate a LOT of files!
		f <- dir(iFold, pattern = "MHWstatsTS", recursive = TRUE, full.names = TRUE)	
			f <- f[-grep("EC-Earth3", f)] # EC-Earth3 crashes the routine at the moment *** CHECK AND CORRECT DOWN THE LINE
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

# Goto 10_Stack_Rasters.R