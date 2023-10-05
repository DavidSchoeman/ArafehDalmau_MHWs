# ArafehDalmau_MHWs
Analysis of marine heatwaves for Arafeh-Dalmau et al. (2023) One Earth.

# Caveats and Dependencies
This code is provided "as is" it runs start-to-end on my machine, but I cannot guarantee the same for others, as I have not tested it in an external environment.
Much of the code requires system calls to CDO (Climate Data Operators), described here: https://code.mpimet.mpg.de/projects/cdo/wiki. These system calls work well on MacOS, but I have not tested them in the Windows environment.

# Original data sources
OISST data were downloaded from NASA: https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
CMIP6 sea-surface temperature projections were downloaded in 2020/2021 from the ESGF website: https://aims2.llnl.gov/search 

# Notes on code development
Code was developed over the period 2020–2023. Some packages might now be deprecated or are nearing deprecation (e.g., raster). As a result, I have subsequently improved some workflows for a follow-on paper and I will update this repo with that improved code once that paper is published.

# How to use the code
After downloading the requisite data, the user should run the code chunks sequentially from "1_tos_Regrid_Anom" to "XXX". Each code chunk should produce output (files written to disk) that is used as input in the subsequent chunk(s).
Note that the folders/directories named in this code pertain to my own setup, and these paths will need to be edited to suit the user's own machine.

# Fair warning
Raw and processed data files (mainly in netCDF format) require several terrabytes of storage space. Do not attempt the workflow before provisioning your machine accordingly.
