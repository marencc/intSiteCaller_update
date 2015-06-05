library(stringr)
library(plyr)
options(stringsAsFactors = FALSE)

testrunlog <- "testrun.log"
testrunmd5 <- "testrun.md5"
if( any(file.exists(testrunlog, testrunmd5)) ) stop(
    "Previous run trace left, run git clean -df first" )

cmd <- sprintf('Rscript ../../intSiteCaller.R -c=../.. > %s 2>&1', testrunlog)
message(cmd)
if( system(cmd)!=0 ) stop(cmd, " not executed")

jobid_lines <- grep("Job.*<\\d+>", readLines(testrunlog), value=TRUE)
stopifnot( length(jobid_lines)==1 )

start_jobid <- as.integer(str_match(jobid_lines, "<(\\d+.)>")[2])
stopifnot( length(start_jobid)==1 )

still_running <- TRUE
minutes <- 0
while( still_running ) {
    mybjobsid <- as.integer(system("bjobs | grep -v JOBID | awk '{print $1}'", intern=TRUE))
    still_running <- any( mybjobsid >= start_jobid )
    if ( still_running ) {
        Sys.sleep(60)
        minutes <- minutes + 1
    }
    message("Running: ", minutes)
}
message("Run stopped after: ", minutes)


cmd <- sprintf("ls */*fa */*RData | sort | xargs md5sum > %s", testrunmd5)
message(cmd)
system(cmd)

oldmd5 <- read.table("intSiteValidation.md5", )
colnames(oldmd5) <- c("oldmd5", "file")
oldmd5 <- arrange(oldmd5, file)
newmd5 <- read.table(testrunmd5)
colnames(newmd5) <- c("newmd5", "file")
newmd5 <- arrange(newmd5, file)

if( !all(newmd5$file==oldmd5$file) ) stop(
             message("File names produced are not same, check the md5 files"))

allmd5 <- merge(oldmd5, newmd5, by="file")
if( any(allmd5$oldmd5!=allmd5$newmd5) ) {
    message("The following files have different md5")
    print(allmd5[allmd5$oldmd5!=allmd5$newmd5,]) }


if( all(allmd5$oldmd5==allmd5$newmd5) ) {
    message("New run produced same files as old run")
    files.remove(testrunlog, testrunmd5) }

q(save="no")

