#PMACS_kickoff.R is reserved for setting/parsing basic parameters such as:
# - bushmanJobID
# - blatStartPort (set for deprecation)
# - cleanup
# - codeDir
# All other 'initiation' logic is done in programFlow.R/startIntSiteCaller()

args <- commandArgs(trailingOnly = FALSE)

### SET RUN PARAMETERS HERE ###
bushmanJobID <- "intSiteValidation" #allows simultaneous processing of datasets - make sure to use unique BLAT ports!
blatStartPort <- 5560 #this can get a bit weird since spawning a bunch of blat threads could result in conflicts with other processes
cleanup <- TRUE
### END RUN PARAMETERS ###

#eventually parse everything as arguments
codeDir <- file.path(sub("--codeDir=", "", grep("--codeDir=", args, value=T)))
if(length(codeDir) == 0 ) codeDir <- "/home/yinghua/projects/intSiteCaller"
stopifnot(length(list.files(path=codeDir, pattern="[intSiteLogic|programFlow].R"))>=2)

source(paste0(codeDir, "/programFlow.R")) #codeDir is only calculated above

startIntSiteCaller()
