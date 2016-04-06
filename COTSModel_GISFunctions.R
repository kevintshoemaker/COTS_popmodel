##############################################
###  GIS functions for the COTS population model
#################################### 


#####################
##  ReadRaster: generic function for reading in ascii raster layers for COTS
#####################

ReadRaster <- function(rastername,projection=projection,plot=F){
  setwd(SPATIALDATA_DIRECTORY)
  newraster <- readGDAL(rastername,p4s=projection)
  newraster <- raster(newraster)
  if(plot) plot(newraster)
  return(newraster)
}


#####################
##  GetReefData: use the reef ID layer to get global data about the study area- e.g., number of reefs etc... 
#####################

GetReefData <- function(reef_ID){
  UNIQUEREEFIDS <<- unique(reef_ID@data@values)[-which(is.na(unique(reef_ID@data@values)))]
  NREEFS <<- length(UNIQUEREEFIDS)
}


