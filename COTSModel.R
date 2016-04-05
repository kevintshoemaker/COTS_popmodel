##########################
#  Script for modeling COTS outbreaks in the Great Barrier Reef 
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting

##########################


#######################
#   CLEAR THE WORKSPACE
#######################

rm(list=ls())

#########################
# SET PROJECT DIRECTORIES
#########################

USER = "KEVIN"

if(USER=="KEVIN") BASE_DIRECTORY <- "C:\\Users\\Kevin\\Dropbox\\CoTS_Model"             # NOTE: this should link to the Dropbox folder with shared project resources	                                                                        
if(USER=="KEVIN") CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\COTS_popmodel"              # NOTE: code directory should be your local copy of the GitHub repository

if(USER=="SAM") BASE_DIRECTORY <- ## ****FILL IN*****
if(USER=="SAM") CODE_DIRECTORY <- ## ****FILL IN*****

SPATIALDATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Spatial Layers",sep="")                          # directory for storing relevant spatial data (ASC, SHP files)
if(is.na(file.info(SPATIALDATA_DIRECTORY)[1,"isdir"])) dir.create(SPATIALDATA_DIRECTORY)

DATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Data",sep="")                                 # directory for storing data (CSV files)
if(is.na(file.info(DATA_DIRECTORY)[1,"isdir"])) dir.create(DATA_DIRECTORY)

FIGURES_DIRECTORY <- paste(BASE_DIRECTORY,"\\Figures\\RawFigures",sep="")               # directory for storing raw figures generated from R
if(is.na(file.info(FIGURES_DIRECTORY)[1,"isdir"])) dir.create(FIGURES_DIRECTORY)

RESULTS_DIRECTORY <- paste(BASE_DIRECTORY,"\\results",sep="")                           # directory for storing relevant results
if(is.na(file.info(RESULTS_DIRECTORY)[1,"isdir"])) dir.create(RESULTS_DIRECTORY)

RDATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\R_Workspaces",sep="")                        # directory for storing .RData files (R workspaces and data objects)
if(is.na(file.info(RDATA_DIRECTORY)[1,"isdir"])) dir.create(RDATA_DIRECTORY)

####################
#  GLOBAL PARAMETERS  (USER SPECIFIED PARAMS)
####################

NREPS <- 1      
NYEARS <- 10
NSEASONS <- 2
SEASONS <- c("summer","winter")


VERBOSE<-TRUE        # flag whether functions should return detailed information

      # save global params
setwd(RDATA_DIRECTORY)
save(NREPS,TIMESTEPS,file="GlobalParams.RData")

####################
#  LOAD FUNCTIONS
#################### 

setwd(CODE_DIRECTORY)
source("COTSModel_Utilityfunctions.R")     # load utility functions, e.g., for loading packages etc. 

#############################
#  LOAD PACKAGES
#############################
                           # note: loadPackage should install the package from CRAN automatically if it is not already installed
   
loadPackage("lhs")            # for latin hypercube sampling
loadPackage("RCurl")          # for loading source code from github
loadPackage("raster")         # for managing raster data
loadPackage("rgdal")

###############################
#        LOAD SPATIAL DATA
###############################

projection <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

	##########
	# Reef ID Layer

setwd(SPATIALDATA_DIRECTORY)
reef_ID <- readGDAL("reefs.asc",p4s=projection)
reef_ID <- raster(reef_ID)

plot(reef_ID)

###############################


###############################
# RUN COTS MODEL
###############################

initializeCOTSabund(,...)      # initialize the COTS abundance object (for year 0) 
initializeCoralCover(,...)    # initialize the coral cover object (for year 0)

for(year in 1:NYEARS){                  # loop through years
  for(season in SEASONS){               # loop through seasons
    doCOTSDispersal(season,COTSabund,...)
    doCOTSDemography(season,COTSabund,CoralCover...)
    doCoralDispersal(season,...)
    doCoralDisturbance(season,COTSabund,...)           # coral disturbance processes, including from COTS
    
    collectResults(year,season,COTSabund,CoralCover)   # collect results for analysis and visualization
  }
}







