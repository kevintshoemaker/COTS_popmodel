########################
# Functions for modeling COTS demography and dispersal
########################



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
#    - initDensityA: for every reef in the study area, a vector of initial adult densities
#    - initDensityS: for every reef in the study area, a vector of initial senile adult densities
# RETURNS:
#    - COTSabund: spatially-structured and stage-structured COTS abundance
#           COTSabund$J_1: vector representing spatially structured abundance of Juvenile stage 1 individuals
#           COTSabund$J_2: vector representing spatially structured abundance of Juvenile stage 2 individuals
#           COTSabund$A: vector representing spatially structured abundance of reproductive adult individuals
#           COTSabund$S: vector representing spatially structured abundance of senile adult individuals
#           NOTE: larvae are not considered explicitly here. 
###################

### for testing:


initializeCOTSabund <- function(reefIDs, initDensityA, initDensityS){
  veclength <- length(reefIDs)
  
     ### Set up the COTS abundance object
  COTSabund <- list()
  COTSabund[['J_1']] <- numeric(veclength)
  COTSabund[['J_2']] <- numeric(veclength)
  COTSabund[['A']] <- numeric(veclength)
  COTSabund[['S']] <- numeric(veclength)
  
    ### set the initial abundances: 
  nReefs <- length(uniqueReefIDs)
  r=1
  for(r in 1:nReefs){
    thisReefID <- uniqueReefIDs[r]
    mask <- reclassify(reefmap,rcl=c(NA,NA,0, -Inf,thisReefID-0.5,0, thisReefID-0.4,thisReefID+0.4,1, thisReefID+0.5,Inf,0))@data@values
    COTSabund[['A']] <- COTSabund[['A']] + (mask*initDensityA[r])
    COTSabund[['S']] <- COTSabund[['S']] + (mask*initDensityS[r])
    totAdult <- initDensityA[r] + initDensityS[r]  
    densJ2 <- totAdult * (COTS_StableStage[2]/COTS_StableStage[3])
    densJ1 <- totAdult * (COTS_StableStage[1]/COTS_StableStage[3])
    COTSabund[['J_1']] <- COTSabund[['J_1']] + (mask*densJ1)
    COTSabund[['J_2']] <- COTSabund[['J_2']] + (mask*densJ2)
  }
}





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

