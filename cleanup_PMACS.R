cleanup = get(load("cleanup.RData"))

if(cleanup){
  system("mv indexPath.RData indexPath", ignore.stderr=T)
  system("rm *.RData", ignore.stderr=T)
  system("rm Data/*.fasta", ignore.stderr=T)
  system("rm */hits.R*.RData", ignore.stderr=T)
  system("rm */p*.fa*", ignore.stderr=T)
  system("rm */keys.RData", ignore.stderr=T)
  system("rm */*Status.RData", ignore.stderr=T)
  system("rm -r logs", ignore.stderr=T)
  system("rm -r Data/demultiplexedReps", ignore.stderr=T)
  system("mv indexPath indexPath.RData", ignore.stderr=T)
}
