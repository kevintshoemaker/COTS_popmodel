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

### for testing:
### set small extent using crop function for testing

## NB Substantial coral mortality at 1,500 CoTS/km2 (=0.22/manta tow)
1500*(100/22)
initializeCOTSabund <- function(reefmap, PopData, COTSInterp, Year, stagenames, nstages, nreefs,npops,initDensityA,initDensityS){
  
     ### Set up the COTS abundance object
  COTSabund <- matrix(0,nrow=npops, ncol=nstages)
  colnames(COTSabund) <- stagenames
  nreefs <- length(unique(PopData$REEF_ID))
  npops <- length(PopData$PIXELID)
  
    ### set the initial abundances: 
  r=1
  for(r in 1:nReefs){
    thisReefID <- reefIDs[r]
    #mask <- reclassify(reefmap,rcl=c(NA,NA,0, -Inf,thisReefID-0.5,0, thisReefID-0.4,thisReefID+0.4,1, thisReefID+0.5,Inf,0))@data@values
    COTSabund[,'A'] <- COTSInterp$Year * 1500 * (ReefPercent/100)   #need to multiply by function from observations to density
    #COTSabund[,'A'] <- COTSabund[,'A'] + (mask*initDensityA[r])#use interpolated layers as starting point
    #COTSabund[,'S'] <- COTSabund[,'S'] + (mask*initDensityS[r])
    #totAdult <- initDensityA[r] + initDensityS[r]  
    #densJ2 <- totAdult * (COTS_StableStage[2]/COTS_StableStage[3])#correcting based on stable stage ratios
    #densJ1 <- totAdult * (COTS_StableStage[1]/COTS_StableStage[3])
    COTSabund[,'J_2'] <- COTSabund[,'A'] * (COTS_StableStage[2]/COTS_StableStage[3])
    COTSabund[,'J_1'] <- COTSabund[,'A'] * (COTS_StableStage[1]/COTS_StableStage[3])
  }
}

### I NEED TO REDO INTERPOLATIONS TO PUT INTO THE REGULAR GRID
    
###################
# CoTS_Mortality
##########
# OBJECTIVE:
#    
# PARAMS:
#    
# RETURNS: 
#     - newCOTSabund: COTS abund updated for stage specific mortality
#     
###################

###################
# CoTS_Fecundity
##########
# OBJECTIVE:
#    
# PARAMS:
#    
# RETURNS:
#     - TotalLarvae: Total larvae produced per grid cell
#     
###################

###################
# CoTS_Dispersal
##########
# OBJECTIVE:
#    
# PARAMS:
#    
# RETURNS:
#     - 
#     
###################

####################
### COTS SANDBOX: for testing, etc.
####################


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

