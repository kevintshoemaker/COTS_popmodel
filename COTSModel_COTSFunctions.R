########################
# Functions for modeling COTS demography and dispersal
########################
#   Authors: Kevin Shoemaker, Sam Matthews
#
#   28 September 2016 - Updated initializeCoTSabund to matrix form
#                     - started scripting mortality, fecundity, dispersal functions
#
#   THINGS TO DO:
#         1. Set initial adult densities (use 1996 to line up with coral model) - check dates on report year
#         2. Write out params of m,f and d functions
#         3. Do we want to have CoTS pops at the reef level


###################
#  BASIC VITAL RATE AND GROWTH FUNCTIONS 
###################

   # converts diameter (in mm) to mass (in grams)

COTS_MassFromDiam <- function(Diam){
  Mass <- 6.29*10^-5*Diam^2.929
  return(Mass)
}

   # converts female mass to total larval fecundity 
COTS_FecFromMass <- function(Mass){
  Fec <- 558*Mass^1.439
  return(Fec)
}

COTS_StableStage <- c(0.9803, 0.0171, 0.0026)   # very approximate stable stage distribution (J1, J2, Adult: see below for back-of-the-envelope calculation)

###################
# initializeCOTSabund
##########
# OBJECTIVE:
#    generate an object for storing the COTS abundance in each pixel. 
# PARAMS:
#    - reefmap: raster template for the study region: NA outside of reefs, reef ID value within reefs 
#    - PopData: data frame containiing PIXEL ID's, percent reef cover and environ vriable
#    - COTSInterp: txt file containing interpolated values of COTS Manta tow, giving a value for CoTS density
#    - Year: Which year we are using as our starting values
#    - Detectability: detectability of adult CoTS from MAnata tow surveys
#    - stagenames: vector of stagenames eg J1, J2, A1
#    - nstages: number of stages
#    - nreefs: number of reefs in simulation
#    - npops: number of separate populations - initially using every reef pixel as a reef
# RETURNS:
#    - COTSabund: spatially-structured and stage-structured COTS abundance
#           COTSabund$J_1: vector representing spatially structured abundance of Juvenile stage 1 individuals
#           COTSabund$J_2: vector representing spatially structured abundance of Juvenile stage 2 individuals
#           COTSabund$A: vector representing spatially structured abundance of reproductive adult individuals
#           COTSabund$S: vector representing spatially structured abundance of senile adult individuals
#           NOTE: larvae are not considered explicitly here. 
###################

setwd(DATA_DIRECTORY)
PopData <- read.table("PopData.txt", header = TRUE, sep = "\t")
COTSInterp <- read.table("COTS_Interpolated.txt", header = TRUE, sep = "\t")

stagenames <- c('J_1', 'J_2', 'A')

initializeCOTSabund <- function(PopData, COTSInterp, Year, stagenames, COTS_StableStage){
  
  ### set NA Values in interpolation to 0
  COTSInterp[is.na(COTSInterp)] <- 0
  
  npops <- length(PopData[,1])
  nstages <- length(stagenames)
  
  ### Set up the COTS abundance object
  COTSabund <- matrix(0,nrow=npops, ncol=nstages)
  colnames(COTSabund) <- stagenames
  
  ### Set up reference for year
  colname <- paste('X', Year, sep="")
  
  ### Update abundances based from interpolated manta tow data
  COTSabund[,'A'] <- COTSInterp[,colname] * 1500 * (PopData$reefpercent/100)   #need to multiply by function from observations to density
  COTSabund[,'J_2'] <- COTSabund[,'A'] * (COTS_StableStage[2]/COTS_StableStage[3])
  COTSabund[,'J_1'] <- COTSabund[,'A'] * (COTS_StableStage[1]/COTS_StableStage[3])
  return(COTSabund)
}

initCOTS <- initializeCOTSabund(PopData = PopData, COTSInterp = COTSInterp, Year=1996, stagenames, COTS_StableStage = COTS_StableStage)

    
###################
# CoTS_StageTransition
##########
# OBJECTIVE: Transition all individuals through life stages
#    
# PARAMS: 
#     - COTSabund: Matrix of CoTS abundances for which to transition
#     - COTSmort: named Vector of natural mortality rates for each stage
#     - COTSremain: named Vector of proportions of individuals that remain in current life stage
#    
# RETURNS: 
#     - newCOTSabund: COTS abund updated after 6 month time step
#     
###################

COTSmort <- c(0.7,0.4,0.1)
names(COTSmort) <- stagenames

COTSremain <- c(0.02,0.2,1) # proportion remianing in each life stage --> we can make this a function of resource
names(COTSremain) <- stagenames

COTS_StageTransition <- function(COTSabund, COTSmort, COTSremain) {
  
  #Set up matrices
  newCOTSabund <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(newCOTSabund) <- colnames(COTSabund)
  COTS_Mort <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(COTS_Mort) <- colnames(COTSabund)
  COTS_Remain <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(COTS_Remain) <- colnames(COTSabund)
  COTS_Trans <- matrix(0,nrow=nrow(COTSabund), ncol=ncol(COTSabund))
  colnames(COTS_Trans) <- colnames(COTSabund)
  
  # apply mortality
  COTS_Mort <- sweep(COTSabund,MARGIN=2,COTSmort,`*`)
  # update abundance
  newCOTSabund <- COTSabund - COTS_Mort
  
  # number of COTS remaining and transitioning for each stage based on post-mortality abundaces
  COTS_Remain <- sweep(newCOTSabund,MARGIN=2,COTSremain,`*`)
  COTS_Trans <- sweep(newCOTSabund,MARGIN=2,1-COTSremain,`*`)
  
  # update newCOTSabund
  newCOTSabund[, 'J_1'] <- COTS_Remain[,'J_1']
  newCOTSabund[, 'J_2'] <- COTS_Remain[,'J_2'] + COTS_Trans[,'J_1']
  newCOTSabund[, 'A'] <- COTS_Remain[,'A'] + COTS_Trans[,'J_2']
  return(newCOTSabund)
}

t1 <- COTS_StageTransition(COTSabund = initCOTS, COTSmort = COTSmort, COTSremain = COTSremain)

head(initCOTS)
head(t1)

###################
# CoTS_Fecundity
##########
# OBJECTIVE: 1. Assume a size distribution amongst adults and then sample diameters
#            2. Convert diameter to mass 
#            3. Convert mass to eggs
#            4. Multiply total eggs by connectivity matrix
#    
# PARAMS: 
#     - COTSabund: standard abundance matrix for the time step
#     - mean: mean of the size distribution of COTS adults
#     - sd: standard deviation of size distribution
#     - npops: number of pixels
#    
# RETURNS:
#     - TotalLarvae: Total larvae produced per pixel
#     
###################

COTS_Fecundity <- function(COTSabund, mean, sd, npops) {
    
  ### Intitialize matrix to store total eggs
  COTS_Eggs <- vector(mode = "numeric", length = npops)
    for (r in 1:npops) {
        Sizes <- rnorm(COTSabund[r,'A'], mean, sd)
        COTS_Eggs[r] <- sum(COTS_FecFromMass(COTS_MassFromDiam(Sizes)))
    }
    return(COTS_Eggs)
}


EGGS <- COTS_Fecundity(t1, 35, 10, npops)


###################
# CoTS_Dispersal
##########
# OBJECTIVE: To disperse larvae throughout the system based upon a connectivity matrix
#    
# PARAMS:
#     - nLarvae: vector of number of larvae produced for each pixel
#     - Conn: 13577x13577 Connectivity matrix containing proportion of original larvae reching every other cell
#    
# RETURNS:
#     - 
#     
###################


###################
# CoTS_Settlement
##########
# OBJECTIVE: To disperse larvae throughout the system based upon a connectivity matrix
#    
# PARAMS:
#     - 
#    
# RETURNS:
#     - 
#     
###################


####################
### COTS SANDBOX: for testing, etc.
####################


#find geographic distances between all sites

Pop1 <- PopData
coordinates(Pop1) <- ~x+y
Gdist <- gDistance(Pop1, Pop1, byid = T)

Gdist[1:10,1:10]
# scale geographic distances between 0-1

# 1/((1-GDistNorm) - 0.5)


m <- matrix(1, 3,3)
v <- 1:3
m*v
   ## build a very rough population transition matrix... 

#  typical reproductive female is 300 mm in diameter

Mass <- COTS_MassFromDiam(300)    # 1132 grams
Fec <- COTS_FecFromMass(Mass)     # each female produces approx 14 million larvae!

  # assume that maybe 0.0001 of these larvae establish on a reef
Fec <- Fec*0.0001

TransMat <- matrix(c(0,0.03,0,0,0,0,0.2,0,Fec/2,0,0.3,0.1,Fec/(2*6),0,0,0.6),nrow=4)


lambda(TransMat)         # strong positive growth rate: 3.57
stable.stage(TransMat)   #stable age distribution

#### take away: stable stage distribution is heavily biased towards juveniles 

