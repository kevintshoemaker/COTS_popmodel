##########################
#  Script for modeling COTS outbreaks in the Great Barrier Reef 
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
# 
#  05 April 2015 -- started scripting

##########################

########################
#  TODO LIST

#  - how to initialize the COTS population??? [KTS developed a crude method using reef-level densities...]

########################

#######################
#   CLEAR THE WORKSPACE 
#######################

rm(list=ls())

#########################
# SET PROJECT DIRECTORIES (this should be the only place where local directories should be referenced)
#########################

USER = "KEVIN"

if(USER=="KEVIN") BASE_DIRECTORY <- "C:\\Users\\Kevin\\Dropbox\\CoTS_Model"             # NOTE: this should link to the Dropbox folder with shared project resources	                                                                        
if(USER=="KEVIN") CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\COTS_popmodel"              # NOTE: code directory should be your local copy of the GitHub repository

if(USER=="SAM") BASE_DIRECTORY ##<- ****FILL IN*****
if(USER=="SAM") CODE_DIRECTORY ##<- ****FILL IN*****

SPATIALDATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Spatial Layers",sep="")                          # directory for storing relevant spatial data (ASC, SHP files)
if(is.na(file.info(SPATIALDATA_DIRECTORY)[1,"isdir"])) dir.create(SPATIALDATA_DIRECTORY)

DATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\Data",sep="")                                 # directory for storing data (CSV files)
if(is.na(file.info(DATA_DIRECTORY)[1,"isdir"])) dir.create(DATA_DIRECTORY)

ENVDATA_DIRECTORY <- paste(DATA_DIRECTORY,"\\Environmental",sep="")                                 # directory for storing data (CSV files)
if(is.na(file.info(ENVDATA_DIRECTORY)[1,"isdir"])) dir.create(ENVDATA_DIRECTORY)

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

VERBOSE <- TRUE        # flag whether functions should return detailed information
DEBUG <- TRUE          # flag whether to output debug files etc. 

      # save global params
setwd(RDATA_DIRECTORY)
save(NREPS,NYEARS,NSEASONS,SEASONS,file="GlobalParams.RData")

####################
#  LOAD FUNCTIONS FROM SOURCE CODE (use your local GitHub repository)
#################### 

setwd(CODE_DIRECTORY)
source("COTSModel_Utilityfunctions.R")   # load utility functions, e.g., for loading packages etc. 
source("COTSModel_COTSfunctions.R")      # load functions for implementing COTS demography and dispersal
source("COTSModel_Coralfunctions.R")     # load functions for implemeting coral growth/recovery and dispersal
source("COTSModel_GISfunctions.R")     # load functions for implemeting coral growth/recovery and dispersal

#############################
#  LOAD PACKAGES
#############################
                           # note: 'loadPackage' should install the package from CRAN automatically if it is not already installed
loadPackages()   # load all packages into the global environment

###############################
#        LOAD SPATIAL DATA
###############################

projection <- "+proj=longlat +datum=WGS84"   #"+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"


###############################
# READ IN ENVIRONMENTAL COVARIATES
###############################

setwd(ENVDATA_DIRECTORY)
envData <- read.table("ENV_dataGBR.txt",header=T,sep="\t")
nrow(envData)
    # generate a template grid
head(envData)

minx <- min(envData$x)-0.005  # farthest west
miny <- min(envData$y)-0.005  # farthest south
maxx <- max(envData$x)+0.005  # farthest east
maxy <- max(envData$y)+0.005  # farthest north

studyRegion <- extent(minx,maxx,miny,maxy)
plot(studyRegion)

template <- raster(ext=studyRegion,resolution=0.01,vals=NA)   # empty raster
plot(template)

length(template@data@values)    # 1,335,372 cells

xcells <- (maxx-minx)/0.01
ycells <- (maxy-miny)/0.01
totcells <- xcells*ycells

   ## rasterize an environmental variable.
reefID <- as.numeric(envData$REEF_ID)
ndx <- !is.na(reefID)
reefIDraster <- rasterize(envData[,c('x','y')][ndx,], template, field=reefID[ndx])   # memory limitations

plot(reefIDraster)

##########
# Reef ID Layer: this layer is used for generating many key parameters and many other maps

   #reef_ID <- ReadRaster("reefs.asc",projection,plot=T)

GetReefData(reefIDraster)   # harvest key data from the reef_ID layer

slotNames(reefIDraster)
?raster

# read in reefs shapefile

?readOGR

reefdatadir <- "C:\\Users\\Kevin\\Dropbox\\CoTS_Model\\Spatial Layers\\Masked reefs"
setwd(reefdatadir)
GBR_polygons <- readOGR(reefdatadir,layer="GBR_reefs")

plot(GBR_polygons)

###############################
# MAKE UP DATA INPUTS
###############################

initDensityA <- round(rnorm(NREEFS,10000,1000))    # density per km2 of reef habitat
initDensityS <- round(rnorm(NREEFS,1000,100))

###############################
# RUN COTS MODEL
###############################

COTSabund <- initializeCOTSabund(,...)      # initialize the COTS abundance object (for year 0) 
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







