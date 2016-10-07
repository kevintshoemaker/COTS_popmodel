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
# CoTS_Fertilisation ----
###################
# OBJECTIVE: Determine the fertilisation success(%) at different densities of COTS (between 0-150,000)
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
#Define Fertilisation by distance fucntion
#####

Fert.data <- data.frame(Dist = c(0,2,4,8,16,32,64,100), PercFert = c(90,86.5, 71.8,71.9,41.5,26.8,20.5,5.8))

m1 <- as.formula(PercFert ~ p * exp(k * Dist)) #standard 
m2 <- as.formula(PercFert ~ p * exp(k * Dist) + q) #standard 
m3 <- as.formula(PercFert ~ p * exp(k * Dist) + (92-p)) #fixed intercept of 92%

em <-function(x,p,k,q) {(p*exp(k*x)) + q}
em.fixed <-function(x,p,k,f) {(p*exp(k*x)) + (f-p)}
em2<-function(x,p,k) {(p*exp(k*x))}


nls1 <- nls(m1,start=list(p=80,k=-0.05), data = Fert.data)
nls2 <- nls(m2,start=list(p=90,k=-0.05, q=10), data = Fert.data, control = list(maxiter=500))
nls3 <- nls(m3,start=list(p=80,k=-0.05), data = Fert.data)

###################
#Retrieve best fit model parameters
#############
BestFitPars <- function(nls.object){
  confints <- confint(nls.object)
  bestpars <- nls.object$m$getPars()
  upperpars <- confints[,2] 
  lowerpars <- confints[,1]
  pars <- data.frame(bestpars,upperpars,lowerpars)
  return(pars)
}

BestFitPars(nls1)
BestFitPars(nls2)
BestFitPars(nls3)

###################
# Plot Fertilisation Function 
############

Data <- Fert.data
nls.object <- nls1

nlsCIplot <- function(nls.object, Data) {
  
  (confints <- confint(nls.object))
  (bestpars <- nls.object$m$getPars())
  upperpars <- confints[,2] 
  lowerpars <- confints[,1]
  (pars <- data.frame(bestpars,upperpars,lowerpars))
  
  fit.data <- data.frame(x=seq(0,200,len=100), 
                         best=NA, upper=NA, lower=NA)
  fit.data$best <- em(fit.data$x, bestpars[1], bestpars[2], bestpars[3])
  
  fit.data$upper <- em(fit.data$x, upperpars[1], upperpars[2], upperpars[3])
  fit.data$lower <- em(fit.data$x, lowerpars[1], lowerpars[2], lowerpars[3])
  ggplot(fit.data, aes(y=best, x=x)) +
    geom_line() + theme_classic() +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) +
    geom_point(data=Data, aes(x=Dist,y=PercFert)) +
    xlab("Distance") +
    ylab("Percentage Eggs Fertilised") +
    ggtitle("Fertilisation Function")
}

nlsCIplot2 <- function(nls.object, Data) {
  
  (confints <- confint(nls.object))
  (bestpars <- nls.object$m$getPars())
  upperpars <- confints[,2] 
  lowerpars <- confints[,1]
  (pars <- data.frame(bestpars,upperpars,lowerpars))
  
  fit.data <- data.frame(x=seq(0,200,len=100), 
                         best=NA, upper=NA, lower=NA)
  fit.data$best <- em2(fit.data$x, bestpars[1], bestpars[2])
  
  fit.data$upper <- em2(fit.data$x, upperpars[1], upperpars[2])
  fit.data$lower <- em2(fit.data$x, lowerpars[1], lowerpars[2])
  ggplot(fit.data, aes(y=best, x=x)) +
    geom_line() + theme_classic() +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) +
    geom_point(data=Data, aes(x=Dist,y=PercFert)) +
    xlab("Distance") +
    ylab("Percentage Eggs Fertilised") +
    ggtitle("Fertilisation Function")
}

nlsCIplot(nls2, Fert.data)
nlsCIplot2(nls1, Fert.data)

(bestpars <- BestFitPars(nls1))

###################
# Fertilisation by Density ------
###############################
# This section will determine the Von Bertalanffy Growth Parameters for 10 different sex ratios and plot the results

geomSeries <- function(base, max) {
  base^(0:floor(log(max, base)))
}


Densities <- geomSeries(2, 10000)
Densities <- c(Densities, 1536, 2560, 3072, 3584, 5120, 6144,7168,9216)
Densities <- sort(Densities)
SexRatio <- c('M' = 0.2, 'F' = 0.8)
SexRatios <- list(c('M' = 0.1, 'F' = 0.9), c('M' = 0.2, 'F' = 0.8), c('M' = 0.3, 'F' = 0.7),
                  c('M' = 0.4, 'F' = 0.6), c('M' = 0.5, 'F' = 0.5), c('M' = 0.6, 'F' = 0.4),
                  c('M' = 0.7, 'F' = 0.3), c('M' = 0.8, 'F' = 0.2), c('M' = 0.9, 'F' = 0.1))

expParams <- bestpars[,1]

sexes <- seq(0.1,0.9, 0.1)
sexes <- paste("M", sexes, sep="")  
# I will have to convert fertilisation based on reefpercent

FertVsDensity <- function (SexRatio, Densities, expParams) {
  # Place all adults amongst landscape
  props <- data.frame(Dens=Densities, Female = NA, Male = NA)
  fertbydens <- list()
  
  for (i in 1:length(Densities)){
    coordsx <- runif(Densities[i], min = 0, max = 1000)
    coordsy <- runif(Densities[i], min = 0, max = 1000)
    coords <- data.frame(x = coordsx, y = coordsy)
    # randomly assign each coordinate pair as male or female
    coords$sex <- sample(c('M', 'F'), length(coordsx), replace = T, prob = SexRatio)
    #determine number of males and females
    props[i,2:3] <- prop.table(table(coords$sex))[1:2]*Densities[i]
  
    if(sum(is.na(props[i,2:3]))==1) {
      fertbydens[[i]] <- rep(NA, Densities[i])
      next
    } else {
    
      # now for each female, we calculate the distance and the probability that her eggs were NOT
      # fertilised by that male
    
      Females <- as.matrix(coords[coords[,'sex']=='F',][,1:2])
      Males <- as.matrix(coords[coords[,'sex']=='M',][,1:2])
      # Each row has the distance between each female and male
      DistMat <- SpatialTools::dist2(Females, Males)
      # convert distance into probability of not fertalising
      em2 <- function(x,p,k) {(1- (p*exp(k*x))/100)}
      DistMat.NotFert <- matrix(sapply(DistMat, FUN = em2, expParams[1], expParams[2]), 
                       nrow = length(Females[,1]), ncol = length(Males[,1]))
      #Multiply across each row to find out the probability that egss from that female were not fertilised
      # by any male
      NotFert <- apply(DistMat.NotFert, MARGIN = 1, FUN = prod)
      fertbydens[[i]] <- 1-NotFert
    }
  }
  # name the list by the desnity of the population
  names(fertbydens) <- Densities
  #convert into a dataframe for plotting
  dens <- rep(Densities, each = 1, times = props[,'Female'])
  fert <- unlist(fertbydens)
  fert.df <- data.frame(dens,fert)
  return(list(Props = props, FertbyDens.df = fert.df, FertbyDens = fertbydens))
}

f1 <- FertVsDensity(Densities, SexRatio, expParams) 

l1 <- lapply(SexRatios, FUN = FertVsDensity, expParams=expParams, Densities=Densities)

names(l1) <- sexes 
ggplot(data=f1[[2]], aes(x=dens, y=fert)) +
  geom_point() +
  geom_smooth() +
  theme_classic()
head(l1[[1]][[2]])
ggplot(data=l1[[1]][[2]], aes(x=dens, y=fert)) +
  geom_point() +
  geom_smooth() +
  theme_classic()
ggplot(data=l1[[5]][[2]], aes(x=dens, y=fert)) +
  geom_point() +
  geom_smooth() +
  theme_classic()
ggplot(data=l1[[7]][[2]], aes(x=dens, y=fert)) +
  geom_point() +
  geom_smooth() +
  theme_classic()
ggplot(data=l1[[9]][[2]], aes(x=dens, y=fert)) +
  geom_point() +
  geom_smooth() +
  theme_classic()

###################
# Von Bertalanffy Growth
#################
#
# This function takes the data created by the FertVsDensity function to define the parameters of the 
# Von Bertalanffy Growth for each of our potential densities
#
#
FvD <- l1

SSQ <- function(theta, x) {
  Linf <- theta[1]
  K <- theta[2]
  t0 <- theta[3]
  epsilon <- rep(0, length(dens))
  lpred <- rep(0, length(dens))
  for (i in 1:length(dens)) {
    lpred[i] <- Linf * (1 - exp(-K * (dens[i] - t0)))
    epsilon[i] <- (fert[i] - lpred[i])^2
  }
  ssq <- sum(na.omit(epsilon))
  return(ssq)
}

VonBertalannfyGrowth <- function(FvD){
    Model <- list()
    Parameters <- data.frame(SexRatio=names(FvD), Linf = NA, K = NA, t0 = NA)
    for (i in 1:length(Fvd.list)) {
      fert <- FvD[[i]][[2]][,2]
      dens <- FvD[[i]][[2]][,1]
      
      theta <- c(1, 0.001, 0.1)
      
      out <- optim(theta, fn = SSQ, method = "BFGS", x = na.omit(dens), hessian = TRUE)
      out$V <- solve(out$hessian)  #solve the hessian
      out$S <- sqrt(diag(out$V))  #Standard Error
      out$R <- out$V/(out$S %o% out$S)  #Correlation
      Model[[i]] <- out
      Parameters[i,2:4] <- out$par
    }
    return(list(Models=Model, Parameters=Parameters))
}


VBG.Models <- VonBertalannfyGrowth(FvD)
VBG.Models[[2]]

fert <- FvD[[5]][[2]][,2]
dens <- FvD[[5]][[2]][,1]

theta <- c(1, 0.1, 0.1)

out <- optim(theta, fn = SSQ, method = "BFGS", x = na.omit(dens), hessian = TRUE)
out$par

##################
## VBGPlot
############

VBGPlot <- function(SR, FvD, Parameters) {
  #Sex Ratio is integer from 1 to 9 relating to proportion of Males --> i.e 1 = 0.1M, 0.9F Sex Ratio
  data <- FvD[[SR]][[2]]

  fit <- Parameters[SR,2] * (1 - exp(-Parameters[SR, 3] * (data$dens - Parameters[SR,4])))
  fit.data <- data.frame(dens=data$dens, fit=fit)

  ggplot(data=data, aes(x=dens, y=fert)) +
    geom_point() +
    geom_smooth() +
    geom_line(data=fit.data, aes(x=dens, y=fit), col="green", size=1) +
    labs(x=expression(CoTS~Density~(km^{-2})),
         y="% of Eggs Fertilised") +
    ggtitle(paste(SR,"M:",10-SR, "F", sep="")) +
    theme_classic() 
  
}

Parameters <- data.frame(SexRatio=names(FvD), Linf = NA, K = NA, t0 = NA)
Parameters[2,2:4] <- out$par
VBGPlot(5, FvD, Parameters = VBG.Models[[2]])
VBGPlot(2, FvD, Parameters)

# SOmething is wrong with my optimising function

EvenSexRatio <- out$par


###################
# CoTS_Dispersal ----
###################
# OBJECTIVE: Take in the amount of eggs produced as well as the fecundity by density function to determin
#     how many larvae to disperse across our landscape via our connectivity matrix
#    
# PARAMS:
#     - nLarvae: vector of number of larvae produced for each pixel
#     - Conn: 13577x13577 Connectivity matrix containing proportion of original larvae reching every other cell
#    
# RETURNS:
#     - 
#     
###################

  #############
  # Build Dispersal Matrix (Pdist)----
  #############

  # Import coords of all of our sites ..NB for now this is the Env_Data
  setwd(DATA_DIRECTORY)
  PopData <- read.table("PopData.txt", header = TRUE, sep = "\t") 
  Coords <- PopData[,c('x','y')]
  coordinates(Coords) <- ~x+y
  Gdist <- gDistance(Coords, Coords, byid = T)
  Gdist[1:10,1:10]
  
  # Convert to km
  Gdist <- Gdist*100
  
  # Limit the distance matrix to ~500km. i.e set >500km to NA
  Gdist[Gdist>500] <- 0
  Gdist.Sp <- Matrix(Gdist, sparse=T)
  
  
  #function to ignore NA's when summing
  plus <- function(x) {
    if(all(is.na(x))){
      c(x[0],NA)} else {
        sum(x,na.rm = TRUE)}
  }
  # Assume probability to be the inverse of distance --> need cumulative distrubtion function
  Pdist <- apply(Gdist.Sp, 2, function(x) ((1/x)/plus(1/x)))

#############
#help
CoTS_Dispersal <- 
  # Assume 30% self recruitment {REFERENCE}
  
  # Assume some random proportion of larvae get washed out to sea {REFERENCE}
  
  # Convert between distance and 0-1 proportion for remaining eggs



###############################
# CoTS_Settlement
###############################
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


x = matrix(rnorm(20), ncol=4)
rownames(x) = paste("X", 1:nrow(x), sep=".")
y = matrix(rnorm(12), ncol=4)
rownames(y) = paste("Y", 1:nrow(y), sep=".")


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

