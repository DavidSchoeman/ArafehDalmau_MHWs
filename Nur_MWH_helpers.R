# Helper functions for processing CMIP6 data and computing MHW stats
	# Written by Dave S
		# April 2020

	# Packages
		library(plyr)
		library(ncdf4)
		library(lubridate)
		library(raster)
		library(SOAR)
		library(tidyverse)
		library(heatwaveR)
		library(doMC)
		library(mblm)
		library(stringr)
		library(PCICt)
		library(sf)
		library(tmap)
		library(tmaptools)


# Function to extract single-coordinate time series from netCDF -----------------------
		
		ptPixel <- function(i, sy = "19830101", ey = "21001231", ee = E) {
			ssp <- getDeets(i)$ssp # What ssp
			mod <- getDeets(i)$model
			ofile <- gsub(getDeets(i)$syr, sy, i)
			regDir <- paste0(oFold, "/pts/", ssp, "/", mod)
			if(!isTRUE(file.info(regDir)$isdir)) dir.create(regDir, recursive = TRUE) # If the folder doesn't exist, create it	
			ofile <- gsub(unlist(strsplit(ofile, "/tos"))[1], regDir, ofile)
			dO <- paste0("cdo sellonlatbox,", paste0(as.vector(ee), collapse = ","), " ", i, " ", unlist(strsplit(ofile, "tos_"))[1], "temp.nc")
			system(dO)
			dO <- paste0("cdo selyear,", substr(sy, 1, 4), "/", substr(ey, 1, 4), " ", unlist(strsplit(ofile, "tos_"))[1], "temp.nc", " ", ofile)
			system(dO)
			nc <- nc_open(ofile)
				for(j in 1:nrow(toDo)) {
					out <- ncvar_get(nc, varid = "tos",
													 start= c(which(nc$dim$lon$vals == toDo[j,1]), # look for long
													 				 which(nc$dim$lat$vals == toDo[j,2]),  # look for lat
													 				 1),
													 count = c(1,1,-1)) #count '-1' means 'all values along that dimension'
					if(sum(out, na.rm = TRUE) != 0) {
						pName <- paste0(toDo[j,1], "_", toDo[j,2])
						pName1 <- pName # Take a copy so that file names can include a "-"
						if(length(grep("-" , pName)) > 0) {pName <- gsub("-", "neg", pName)}
						if(grepl("\\d", substr(pName, 1, 1))) {pName <- paste0("pos", pName)} # Avoid starting an object name with a number
						eval(parse(text = paste0(pName, " <- out")))
						eval(parse(text = paste0("write.csv(", pName, ", '", unlist(strsplit(ofile, "tos_"))[1], pName1, ".csv', row.names = FALSE)")))
					}
				}
			file.remove(paste0(unlist(strsplit(ofile, "tos_"))[1], "temp.nc"))
		}
		
		

# Function to just get the subsetted, cropped raster for an extent --------

		getExt <- function(i, sy = "19830101", ey = "21001231", ee = E) {
			ssp <- getDeets(i)$ssp # What ssp
			mod <- getDeets(i)$model
			ofile <- gsub(getDeets(i)$syr, sy, i)
			regDir <- paste0(oFold, "/pts/", ssp, "/", mod)
			if(!isTRUE(file.info(regDir)$isdir)) dir.create(regDir, recursive = TRUE) # If the folder doesn't exist, create it	
			ofile <- gsub(unlist(strsplit(ofile, "/tos"))[1], regDir, ofile)
			dO <- paste0("cdo sellonlatbox,", paste0(as.vector(ee), collapse = ","), " ", i, " ", unlist(strsplit(ofile, "tos_"))[1], "temp.nc")
				system(dO)
			dO <- paste0("cdo selyear,", substr(sy, 1, 4), "/", substr(ey, 1, 4), " ", unlist(strsplit(ofile, "tos_"))[1], "temp.nc", " ", ofile)
				system(dO)
		}


# Function to extract time series for all non-NA pixels -------------------

	getPix <- function(i, sy = "19830101", ey = "21001231", ee = E) {
		regDir <- unlist(strsplit(i, "tos_Oday_"))[1]
		r <- stack(i)[[1]] # The first layer of the netCDF we're working with
		XY <- getNearestPixel(r, toDo) # Get coordinates of the nearest cell with data in it
		eval(parse(text = paste0("write.csv(XY, ", "'", regDir, "XY.csv', row.names = FALSE)")))
		XY <- unique(XY)
		nc <- nc_open(i)
		for(j in 1:nrow(XY)) {
			out <- ncvar_get(nc, varid = "tos",
											 start= c(which(nc$dim$lon$vals == XY[j,1]), # look for long
											 				 which(nc$dim$lat$vals == XY[j,2]),  # look for lat
											 				 1),
											 count = c(1,1,-1)) #count '-1' means 'all values along that dimension'
			if(sum(out, na.rm = TRUE) != 0) {
				pName <- paste0(XY[j,1], "_", XY[j,2])
				pName1 <- pName # Take a copy so that file names can include a "-"
				if(length(grep("-" , pName)) > 0) {pName <- gsub("-", "neg", pName)}
				if(grepl("\\d", substr(pName, 1, 1))) {pName <- paste0("pos", pName)} # Avoid starting an object name with a number
				eval(parse(text = paste0(pName, " <- out")))
				eval(parse(text = paste0("write.csv(", pName, ", '", regDir, pName1, ".csv', row.names = FALSE)")))
			}
		}
	}
	
	
# Function to get details from file name ----------------------------------
		
	getDeets <- function(x) {
		bits <- unlist(strsplit(x, "/"))
		fName <- bits[length(bits)]
		bits1 <- unlist(strsplit(fName, "_"))
		yrs <- bits1[length(bits1)]
		yrs <- unlist(strsplit(yrs, ".nc"))
		yrs <- unlist(strsplit(yrs, "-"))
		ssp <- bits[grep("ssp", bits)[1]]
		return(list(model = bits1[!grepl("tos|Oday|.nc", bits1)], ssp = ssp, syr = yrs[1], eyr = yrs[2]))
	}
		
	# The wrapper function to compute baseline climatology and then de --------
	
	# Based on https://cran.r-project.org/web/packages/heatwaveR/vignettes/gridded_event_detection.html
	event_only <- function(df){
		# First calculate the climatologies
		clim <- ts2clm(data = df, x = Date, y = SST, climatologyPeriod = c("1983-01-01", "2012-12-31")) # As per Oliver et al. and Smale et al.
		# Then the events
		event <- detect_event(data = clim, x = Date, y = SST)
		# Last, we return only the event dataframe of results
		return(event$event)
	}
	
	
# Compute MHW stats -------------------------
	
	# The 'purrr' wrapper function to pass to 'dplyr'
	purrr_wraper <- function(file_name){
		sst <- read.csv(file_name)
		if(length(na.omit(sst$x)) == nd) { # If ALL values are available for the cell
			bits <- unlist(strsplit(unlist(strsplit(file_name, paste0(o, "/")))[2], ".csv"))
			bits <- unlist(strsplit(bits, "_"))
			sstDat <- data.frame(Date = Dates$Date,
													 x = as.numeric(bits[1]),
													 y = as.numeric(bits[2]),
													 SST = round(sst$x, 2))
			event <- event_only(sstDat) %>% 
				mutate(x = sstDat$x[1], y = sstDat$y[1]) %>% 
				dplyr::select(x, y, everything())
		}
	}

	purrr_O_wraper <- function(file_name){
		sst <- read.csv(file_name)
			bits <- unlist(strsplit(unlist(strsplit(file_name, "/Volumes/Data/OISST_Combo/"))[2], ".csv"))
			bits <- unlist(strsplit(bits, "_"))
			sstDat <- data.frame(Date = Dates$Date,
													 x = as.numeric(bits[1]),
													 y = as.numeric(bits[2]),
													 SST = round(sst$x, 2))
			event <- event_only(sstDat) %>% 
				mutate(x = sstDat$x[1], y = sstDat$y[1]) %>% 
				dplyr::select(x, y, everything())
	}	
	
# Wrapper to create times series by pixel --------------------------------------------------------------------
	
	mhw_sum <- function(i) {
		eval(parse(text = paste0("MHW_result <- readRDS('", i, "')")))		
		event_sum <- MHW_result %>%
			mutate(year = year(date_start)) %>%
			group_by(x, y, year) %>%
			dplyr::summarise(n = n_distinct(event_no),
											 Peak_Int = max(intensity_max),
											 Cum_Int = sum(intensity_cumulative),
											 Duration = mean(as.numeric(date_end-date_start)),
											 Cum_Dur = sum(as.numeric(date_end-date_start)))
		nm <- gsub("MHWout", "MHWstatsTS", i)
		eval(parse(text = paste0("saveRDS(event_sum, file = '", nm, "')")))
	}

	mhw_O_sum <- function(i) {
		eval(parse(text = paste0("MHW_result <- readRDS('", i, "')")))		
		event_sum <- MHW_result %>%
			mutate(year = year(date_start)) %>%
			group_by(x, y, year) %>%
			dplyr::summarise(n = n_distinct(event_no),
											 Peak_Int = max(intensity_max),
											 Cum_Int = sum(intensity_cumulative),
											 Duration = mean(as.numeric(date_end-date_start)),
											 Cum_Dur = sum(as.numeric(date_end-date_start)))
		eval(parse(text = paste0("saveRDS(event_sum, file = 'OISST_Stats.Rda')")))
	}		
	
# Functions for the statistical analysis of trends in MHWs ----------------
	
	mxYr <- function(ev) {return(max(ev$year, na.rm = TRUE))}
	lin_n <- function(ev) {
		mod1 <- mblm(n ~ year, data = filter(ev, Cum_Dur <= 365))
		# extract slope coefficient and its p-value
		tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
										 p = summary(mod1)$coefficients[2,4])
		return(tr)
	}
	lin_int <- function(ev) {
		mod1 <- mblm(Peak_Int ~ year, data = filter(ev, Cum_Dur <= 365))
		# extract slope coefficient and its p-value
		tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
										 p = summary(mod1)$coefficients[2,4])
		return(tr)
	}
	lin_cumi <- function(ev) {
		mod1 <- mblm(Cum_Int ~ year, data = filter(ev, Cum_Dur <= 365))
		# extract slope coefficient and its p-value
		tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
										 p = summary(mod1)$coefficients[2,4])
		return(tr)
	}			
	lin_dur <- function(ev) {
		mod1 <- mblm(Duration ~ year, data = filter(ev, Cum_Dur <= 365))
		# extract slope coefficient and its p-value
		tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
										 p = summary(mod1)$coefficients[2,4])
		return(tr)
	}
	lin_cdur <- function(ev) {
		mod1 <- mblm(Cum_Dur ~ year, data = filter(ev, Cum_Dur <= 365))
		# extract slope coefficient and its p-value
		tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
										 p = summary(mod1)$coefficients[2,4])
		return(tr)
	}
	

# Extract value from nearest coastal cell ---------------------------------

	getNearestPixel <- function(r, xy, type = "land") { # Where r is a raster layer, xy are SpatialPoints (or just a matrix), all projected lonlat, and type is either "land" (default) or "sea"
		if(type == "land") {
			coast <- raster::boundaries(r, type = "inner")
		} else {
			coast <- raster::boundaries(r, type = "outer")			
		}
		coast[coast[] == 0] <- NA # Replace all 0s in coast with NAs
		cmat <- data.frame(xyFromCell(r, 1:ncell(r)), Temp = r[])
		cmat <- na.omit(cmat)
		tsp <- SpatialPoints(cmat[,1:2])
		if(class(xy) == "SpatialPoints") {
			XY <- xy@coords
		} else {
			XY <- xy
		}
		cXY <- tsp[apply(XY, 1, FUN = function(x) {which.min(pointDistance(tsp, x, lonlat = TRUE))})] %>% 
			as.data.frame()
		return(cXY)
	}
	

# Clipping OISST ----------------------------------------------------------

	clpOISST <- function(i, e = E) {
		nm <- paste0(unlist(strsplit(i, ".nc")), "_Clipped.nc")
		dO <- paste0("cdo sellonlatbox,", paste0(as.vector(e), collapse = ","), " ", i, " ", nm)
		system(dO)
	}

	
# OISST Function to extract time series for nearest non-NA pixels -------------------
	
	getOPix <- function(i, sy = "19830101", ey = "21001231", ee = E) {
		regDir <- "/Volumes/Data/OISST_Combo/"
		r <- stack(i)[[1]] # The first layer of the netCDF we're working with
		XY <- getNearestPixel(r, toDo) # Get coordinates of the nearest cell with data in it
		eval(parse(text = paste0("write.csv(XY, ", "'", regDir, "XY.csv', row.names = FALSE)")))
		XY <- unique(XY)
		nc <- nc_open(i)
		for(j in 1:nrow(XY)) {
			out <- ncvar_get(nc, varid = "sst",
											 start= c(which(nc$dim$lon$vals == XY[j,1]), # look for long
											 				 which(nc$dim$lat$vals == XY[j,2]),  # look for lat
											 				 1, 1),
											 count = c(1,1,1,-1)) #count '-1' means 'all values along that dimension'
			if(sum(out, na.rm = TRUE) != 0) {
				pName <- paste0(XY[j,1], "_", XY[j,2])
				pName1 <- pName # Take a copy so that file names can include a "-"
				if(length(grep("-" , pName)) > 0) {pName <- gsub("-", "neg", pName)}
				if(grepl("\\d", substr(pName, 1, 1))) {pName <- paste0("pos", pName)} # Avoid starting an object name with a number
				eval(parse(text = paste0(pName, " <- out")))
				eval(parse(text = paste0("write.csv(", pName, ", '", regDir, pName1, ".csv', row.names = FALSE)")))
			}
		}
	}	
	
# Old ptPixel -------------------------------------------------------------

	ptPixel <- function(i, sy = "19830101", ey = "21001231") {
		ssp <- getDeets(i)$ssp # What ssp
		mod <- getDeets(i)$model
		ofile <- gsub(getDeets(i)$syr, sy, i)
		regDir <- paste0(oFold, "/pts/", ssp, "/", mod)
		if(!isTRUE(file.info(regDir)$isdir)) dir.create(regDir, recursive = TRUE) # If the folder doesn't exist, create it	
		ofile <- gsub(unlist(strsplit(ofile, "/tos"))[1], regDir, ofile)
		dO <- paste0("cdo sellonlatbox,", paste0(as.vector(e), collapse = ","), " ", i, " ", unlist(strsplit(ofile, "tos_"))[1], "temp.nc")
		system(dO)
		dO <- paste0("cdo selyear,", substr(sy, 1, 4), "/", substr(ey, 1, 4), " ", unlist(strsplit(ofile, "tos_"))[1], "temp.nc", " ", ofile)
		system(dO)
		nc <- nc_open(ofile)
		for(j in 1:nrow(toDo)) {
			out <- ncvar_get(nc, varid = "tos",
											 start= c(which(nc$dim$lon$vals == toDo[j,1]), # look for long
											 				 which(nc$dim$lat$vals == toDo[j,2]),  # look for lat
											 				 1),
											 count = c(1,1,-1)) #count '-1' means 'all values along that dimension'
			if(sum(out, na.rm = TRUE) != 0) {
				pName <- paste0(toDo[j,1], "_", toDo[j,2])
				pName1 <- pName # Take a copy so that file names can include a "-"
				if(length(grep("-" , pName)) > 0) {pName <- gsub("-", "neg", pName)}
				if(grepl("\\d", substr(pName, 1, 1))) {pName <- paste0("pos", pName)} # Avoid starting an object name with a number
				eval(parse(text = paste0(pName, " <- out")))
				eval(parse(text = paste0("write.csv(", pName, ", '", unlist(strsplit(ofile, "tos_"))[1], pName1, ".csv', row.names = FALSE)")))
			}
		}
		file.remove(ofile, paste0(unlist(strsplit(ofile, "tos_"))[1], "temp.nc"))
	}
	