
####################################
###  Utility functions for the COTS population model
#################################### 


########################################
## GENERIC FUNCTION FOR INSTALLING/LOADING PACKAGES FROM CRAN
########################################

loadPackage <- function(pkg){

  if(pkg %in% rownames(installed.packages()) == FALSE) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)

}

loadPackages <- function(){
  loadPackage("lhs")            # for latin hypercube sampling
  loadPackage("RCurl")          # for loading source code from github
  loadPackage("raster")         # for managing raster data
  loadPackage("rgdal")          # for reading and writing all sorts of spatial data   
  loadPackage("popbio")         # for doing basic matrix population modeling
}

source_github <- function(baseurl,scriptname) {
  # load package
  loadPackage(RCurl)
 
  # read script lines from website
  url <- sprintf("%s%s",baseurl,scriptname)
  script <- getURL(url, ssl.verifypeer = FALSE)
  
  script <- gsub("\r\n", "\n", script)     # get rid of carriage returns (not sure why this is necessary...)
 
  # parse lines and evaluate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}


 




