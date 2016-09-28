##########################
#  R script: generate key GIS layers for modeling COTS population in the GBR 
#
#  Authors: Kevin Shoemaker, Sam Matthews, Camille Mellin, Damien Fordham
#  
#  08 April 2015 -- started scripting
#
#  26 September 2016 -- redefined GIS layers to comply with GBRMPA MArineBioregions_WGS84 Shapefile
#
#         URGENT: NEED TO OBTAIN ENTIRE NERP DATASET AND EXTRACT BASED ON NEW DEFINITIONS

##########################


###############################
#  PREPARE KEY GIS LAYERS
###############################

################
### create "reefraster" layer: represents percent of each grid cell with reef coverage. 

### UPDATE 28 September 2009 -- Convert to using the MarineBioregions_WGS84 Shapefile as Basis
################

### Create template raster for the GBR (proper extent and resolution)
################

setwd(ENVDATA_DIRECTORY)    # first load the xyz data for each population of interest 
EnvData <- read.table("ENV_dataGBR.txt",header=T,sep="\t")
NPOPS <- nrow(EnvData)
# generate a template grid
head(EnvData)


minx <- min(EnvData$x)-0.005  # farthest west
miny <- min(EnvData$y)-0.005  # farthest south
maxx <- max(EnvData$x)+0.005  # farthest east
maxy <- max(EnvData$y)+0.005  # farthest north
studyRegion <- extent(minx,maxx,miny,maxy)
cat(sprintf("Western margin is %s degrees, southern margin is %s degrees, eastern margin is %s degrees, and northern margin is %s degrees",minx,miny,maxx,maxy))

plot(studyRegion, main="Study Region (blank)")

# save study region extent (as "extent" object)
setwd(SPATIALDATA_DIRECTORY)
save(studyRegion,file="studyRegion.RData")


template <- raster(ext=studyRegion,resolution=0.01,vals=0)   # create empty raster
template001 <- raster(ext=studyRegion,resolution=0.001,vals=0)   # create empty raster at 0.001 res
#plot(template)
setwd(SPATIALDATA_DIRECTORY)
writeRaster(template,filename="templateRaster.asc",format="ascii",overwrite=T)   # write to file

##################
## create reefraster from GBRMPA Classification at 0.001 resolution
##################
reefshape <- readShapefile("MarineBioregions_WGS84", projection = projection, plot=T)
crs(reefshape) <- crs(projection)
reefraster001 <- rasterize(reefshape, template001, field = NA)
setwd(SPATIALDATA_DIRECTORY)
writeRaster(reefraster,filename="reefraster001.asc",format="ascii",overwrite=T)   # write to file

summary(reefshape)
##################
## Create ReefPercent Raster (aggregate reefraster to 0.01 resolution to calculate percent reef in each 0.01 cell)
##################

reefraster01 <- reclassify(reefraster,rcl=c(NA,NA,0, -Inf,0.5,0, 0.501,Inf,1))  # temporary layer: reclassify to binary
compute_percentReef <- function(t,na.rm=TRUE) sum(t,na.rm)  ## aggregation function
reefpercent <- aggregate(reefraster, fact=10, fun=compute_percentReef, expand=TRUE, na.rm=TRUE) 
reefpercent <- reefpercent-1    # for some reason, the result ends up between 1 and 101- reformulate for proper percent
plot(reefpercent)
setwd(SPATIALDATA_DIRECTORY)
writeRaster(reefpercent,filename="reefPercentRaster.asc",format="ascii",overwrite=T)   # write to file
rm(reefpercent)


##################
## Create ReefBioregion ReefID raster files  
##################

reefbioregion <- rasterize(reefshape, template, field = "BIOREGION")
writeRaster(reefbioregion,filename="reefbioregion.asc",format="ascii",overwrite=T)

reefID <- rasterize(reefshape, template, field = "OBJECTID")
writeRaster(reefID,filename="reefID.asc",format="ascii",overwrite=T)

##################
## create a XYZ grid from reefshape keeping all fields
##################
reefpercent.df <- data.frame(data.frame(x = coordinates(reefpercent)[,1], y = coordinates(reefpercent)[,2], reefpercent = reefpercent@data@values))
reefpercent.df <- subset(reefpercent.df, reefpercent>0)
coords <- reefpercent.df[,1:2]
sp <- SpatialPoints(coords = coords, proj4string = CRS(projection))

#extract vals from shapefile for every pixel that contains reef
vals <- over(sp, reefshape)
reefs <- cbind(coordinates(sp), vals, reefpercent = reefpercent.df[,3])

#create subset of pixels that contain reef but aren't assigned attributes
reef.NA <- subset(reefs, is.na(REEF_ID))
#creat subset that do have attributes
reef.YES <- subset(reefs, !is.na(REEF_ID))
#convert both the spatial points dataframes to determine the geographic distance between them
coordinates(reef.NA) <- ~x+y; coordinates(reef.YES) <- ~x+y
Gdist <- gDistance(reef.NA,reef.YES, byid = T)
#min.d <- apply(Gdist, 2, which.min) #creates vector of the minimum distances which we use to index the cords 
#min.gd <- apply(Gdist, 2, min)



##################
#Check how all shapefiles, rasters and environmental data line up
##################
envcoords <- EnvData[,2:3]
envpoints <- SpatialPoints(coords = envcoords, proj4string = CRS(projection))

plot(reefpercent, xlim=c(150.2,150.5), ylim=c(-22.4,-22.2))
#plot(reefraster001, xlim=c(150.2,150.5), ylim=c(-22.4,-22.2), add = TRUE)
plot(reefshape, xlim=c(150.2,150.5), ylim=c(-22.4,-22.2), add = TRUE)
points(sp, xlim=c(150.2,150.5), ylim=c(-22.4,-22.2), pch=20, col="blue", cex=0.5)
points(envpoints, xlim=c(150.2,150.5), ylim=c(-22.4,-22.2), pch=17, col="red", cex=0.5)

# Points lie directly in the centre of the grid cell
# WE ARE MISSING HUGE CHUNKS OF DATA, AND SOME DATA IS ON THE LAND

#merge env data with reef data

data <- dplyr::inner_join(reefs, EnvData, by = c("x", "y"))


# there is a severe mismatch between environmental data and reef percent we should 25026 cases 


#our reef ID is not based off the environmental variables and needs to be the XYZ grid from camille



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
writeRaster(reefIDraster,filename="reefIDRaster.asc",format="ascii",overwrite=T)   # write to file
plot(reefIDraster)

rm(reefIDraster)
