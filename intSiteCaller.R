#check if environment is suitable for running intSiteCaller
#command line stuff
commandLinePrograms <- c("blat", "python")
programsPresent <- !sapply(sprintf("which %s > /dev/null 2>&1", commandLinePrograms), system)
if(any(!programsPresent)){
  stop(paste(commandLinePrograms[!programsPresent]), " is not available")
}

#R packages
rPackages <- c("ShortRead", "hiReadsProcessor", "GenomicRanges",
               "rtracklayer", "BSgenome", "argparse")
#presence of individual BSgenome packages (ex. hg18, hg19) is checked by
#get_reference_genome called from postTrimReads
rPackagesPresent <- is.element(rPackages, installed.packages()[,1])
if(any(!rPackagesPresent)){
  stop(paste(rPackages[!rPackagesPresent]), " is not available")
}

library("argparse", quietly=T)

#define args
parser <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')

parser$add_argument("-j", "--jobID", type="character", nargs=1, default="intSiteCallerJob",
                    help="Unique name by which to identify this intance of intSiteCaller [default: %(default)s]")
parser$add_argument("-c", "--codeDir", type="character", nargs=1, default=".",
                    help="Directory where intSiteCaller code is stored, can be relative or absolute [default: %(default)s]")
parser$add_argument("-p", "--primaryAnalysisDir", type="character", default=".",
                    help="Location of primary analysis directory, can be relative or absolute [default: %(default)s]")
parser$add_argument("-C", "--cleanup", type="logical", nargs="?", const=TRUE, default=FALSE,
                    help="Remove temporary files upon successful execution of intSiteCaller")

parsedArgs <- parser$parse_args(commandArgs(trailingOnly = TRUE))

#source is necessary so that processMetadata() is available
#parsedArgs$codeDir can be given in absolute path OR relative path from intSiteCaller.R
#so it's ok to just do paste0(codeDir, "/programFlow.R") here, in intSiteCaller.R
source(paste0(parsedArgs$codeDir, "/programFlow.R")) 

processMetadata()
