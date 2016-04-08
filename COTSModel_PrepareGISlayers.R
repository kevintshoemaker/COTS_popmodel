##########################
#  R script: generate key GIS layers for modeling COTS population in the GBR 
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
#  
#  08 April 2015 -- started scripting

##########################


###############################
#  PREPARE KEY GIS LAYERS
###############################

################
### create "reefraster" layer: represents percent of each grid cell with reef coverage. 
################

reefraster <- ReadRaster("reefraster2.asc",projection=projection,plot=F)   # uses "ReadRaster" utility function 
reefraster <- reclassify(reefraster,rcl=c(NA,NA,0, -Inf,0.5,0, 0.6,Inf,1))  # temporary layer: reclassify to binary
compute_percentReef <- function(t,na.rm=TRUE) sum(t,na.rm)  ## aggregation function
reefpercent <- aggregate(reefraster, fact=10, fun=compute_percentReef, expand=TRUE, na.rm=TRUE) 
reefpercent <- reefpercent-1    # for some reason, the result ends up between 1 and 101- reformulate for proper percent
plot(reefpercent)
setwd(SPATIALDATA_DIRECTORY)
writeRaster(reefpercent,filename="reefPercentRaster.asc",format="ascii",overwrite=T)   # write to file

rm(reefpercent) # remove from memory (save space!)
rm(reefraster)


################
### Create template raster for the GBR (proper extent and resolution)
################

setwd(ENVDATA_DIRECTORY)    # first load the xyz data for each population of interest 
PopData <- read.table("ENV_dataGBR.txt",header=T,sep="\t")
NPOPS <- nrow(PopData)
# generate a template grid
head(PopData)


minx <- min(PopData$x)-0.005  # farthest west
miny <- min(PopData$y)-0.005  # farthest south
maxx <- max(PopData$x)+0.005  # farthest east
maxy <- max(PopData$y)+0.005  # farthest north
studyRegion <- extent(minx,maxx,miny,maxy)
cat(sprintf("Western margin is %s degrees, southern margin is %s degrees, eastern margin is %s degrees, and northern margin is %s degrees",minx,miny,maxx,maxy))

plot(studyRegion, main="Study Region (blank)")

    # save study region extent (as "extent" object)
setwd(SPATIALDATA_DIRECTORY)
save(studyRegion,file="studyRegion.RData")


template <- raster(ext=studyRegion,resolution=0.01,vals=0)   # create empty raster
#plot(template)
setwd(SPATIALDATA_DIRECTORY)
writeRaster(template,filename="templateRaster.asc",format="ascii",overwrite=T)   # write to file

totalCells <- length(template@data@values)    # 1,335,372 cells
cat(sprintf("total cells in study region raster is %s",totalCells))


################
### Create Reef ID layer for the GBR
################

## rasterize an XYS variable.
reefID <- as.numeric(PopData$REEF_ID)
ndx <- !is.na(reefID)
reefIDraster <- rasterize(PopData[,c('x','y')][ndx,], template, field=reefID[ndx])   # memory limitations

setwd(SPATIALDATA_DIRECTORY)
writeRaster(template,filename="reefIDRaster.asc",format="ascii",overwrite=T)   # write to file
plot(reefIDraster)

rm(reefIDraster)
