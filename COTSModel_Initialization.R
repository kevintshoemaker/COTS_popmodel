##########################
#  Script for preparing the R workspace for modeling COTS population...  
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
#
#  NOTE: this is the only script that refers to specific users and file structures. 
# 
#  07 April 2015 -- started scripting

##########################


####################
#  GLOBAL PARAMETERS  (USER SPECIFIED PARAMS)
####################

NREPS <- 1      
NYEARS <- 10
NSEASONS <- 2
SEASONS <- c("summer","winter")

VERBOSE <- TRUE        # flag whether functions should return detailed information
DEBUG <- TRUE          # flag whether to output debug files etc. 

# save global params
setwd(RDATA_DIRECTORY)
save(NREPS,NYEARS,NSEASONS,SEASONS,file="GlobalParams.RData")


#############################
#  LOAD PACKAGES
#############################
# note: 'loadPackage' should install the package from CRAN automatically if it is not already installed
loadPackages()   # load all packages into the global environment

###############################
#        LOAD PROJECTION FOR READING IN SPATIAL DATA
###############################

projection <- "+proj=longlat +datum=WGS84"   #"+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

#############
#  END SCRIPT
#############