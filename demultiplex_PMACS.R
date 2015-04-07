library("ShortRead")

I1 = readFasta(list.files("Data", pattern="correctedI1-.", full.names=T)) #readFasta("Data/correctedI1.fasta")

metadata = read.csv("allRepsBCs.csv")

I1 = I1[as.vector(sread(I1)) %in% metadata$I7.Index.Sequence]
#indexSeqs = as.vector(sread(I1))
#validIndexSeqs=subset(indexSeqs, indexSeqs %in% metadata$I7.Index.Sequence)
#reps = unlist(sapply(as.character(sread(I1)), function(x){metadata[metadata$I7.Index.Sequence %in% x, "Replicate"]}))
reps = metadata[match(as.character(sread(I1)), metadata$I7.Index.Sequence), "Replicate"]

#only necessary if using native data - can parse out description w/ python
I1Names =  sapply(strsplit(as.character(ShortRead::id(I1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!

rm(I1)

suppressWarnings(dir.create("Data/demultiplexedReps"))

#R1 

R1 = readFastq("Data/Undetermined_S0_L001_R1_001.fastq.gz")
R1Names = sapply(strsplit(as.character(ShortRead::id(R1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
names(R1Names) = NULL

R1 = R1[match(I1Names, R1Names)]
R1 = split(R1, reps)
for (i in 1:length(R1)){writeFastq(R1[[i]], paste0("Data/demultiplexedReps/&", names(R1[i]), "_S0_L001_R1_001.fastq.gz"), mode="w")}

rm(R1, R1Names)


#R2

R2 = readFastq("Data/Undetermined_S0_L001_R2_001.fastq.gz")
R2Names = sapply(strsplit(as.character(ShortRead::id(R2)), " "), "[[", 1) #for some reason we can't dynamically set name/id on ShortRead!
names(R2Names) = NULL

R2 = R2[match(I1Names, R2Names)]
R2 = split(R2, reps)

for (i in 1:length(R2)){writeFastq(R2[[i]], paste0("Data/demultiplexedReps/&", names(R2[i]), "_S0_L001_R2_001.fastq.gz"), mode="w")}

rm(R2, R2Names, I1Names, reps)

