# Extracting data for computing MHW stats
	# Written by Dave S
		# March 2020

# Make available the helper functions
	source("MWH_helpers.R")

# Paths and names
	oFold <- "/Volumes/CMIP6_MHW_Data/temp_csv" # A temporary output path...this is going to generate a LOT of files!
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
	
	# Goto 8_SummariseMHW.R
