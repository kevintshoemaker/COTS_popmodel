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