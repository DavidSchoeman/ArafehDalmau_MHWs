# MHW summaries by year for global analysis
	# Written by Dave S
		# March 2020

# Make available the helper functions
	source("MWH_helpers.R")

# Summaries:
	iFold <- "/Volumes/CMIP6_MHW_Data/temp_csv" # A temporary output path...this is going to generate a LOT of files!
	f <- dir(iFold, pattern = "MHWout", recursive = TRUE, full.names = TRUE)
	registerDoMC(cores = detectCores()-10) # Reduce number of cores to prevent memory overload	
		plyr::ldply(f, .fun = mhw_sum, .parallel = TRUE)

# Goto 9_MHWTrends.R
