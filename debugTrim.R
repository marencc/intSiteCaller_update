## this is debugging for:
## run: run20150609
## replicate: GTSP0440-1
##
## problem was trhat
## replicates of GTSP0440 have 165819 reads but only 178
## had all primer, ltrbit, linker.
## see https://microb215.med.upenn.edu/Download/data/share/GTSPReports/bushmanlab/html/run20150609.stat.html
##
## debugging data was located
## microb244:~/run20150609/
##
## The following debugging info seems to suggest that
## the code performed as expected. Sample GTSP0440 might
## have some problem.

codeDir <- get(load("codeDir.RData"))
source(file.path(codeDir, "intSiteLogic.R"))
## here goes the id for GTSP0440-1
sampleID <- 29 

completeMetadata <- get(load("completeMetadata.RData"))[sampleID,]
alias <- completeMetadata$alias
print(t(as.data.frame(completeMetadata)), quote=FALSE)

workingDir <- alias
setwd(workingDir)

stats.bore <- data.frame(sample=alias)



qualityThreshold=completeMetadata$qualityThreshold
badQuality=completeMetadata$badQualityBases
qualityWindow=completeMetadata$qualitySlidingWindow
primer=completeMetadata$primer
ltrbit=completeMetadata$ltrBit
largeLTRFrag=completeMetadata$largeLTRFrag
linker=completeMetadata$linkerSequence
linker_common=completeMetadata$linkerCommon
mingDNA=completeMetadata$mingDNA
read1=completeMetadata$read1
read2=completeMetadata$read2
alias=completeMetadata$alias
vectorSeq=completeMetadata$vectorSeq


read1 <- completeMetadata$read1
read2 <- completeMetadata$read2
read1 <- sub("GTSPRun/", "", read1)
read2 <- sub("GTSPRun/", "", read2)
stopifnot(all(file.exists(read1, read2)))




## inside the function
reads <- lapply(list(read1, read2), sapply, readFastq)

stats.bore$barcoded <- sum(sapply(reads[[1]], length))

r <- lapply(reads, function(x){
    seqs <- x[[1]]
    if(length(seqs) > 0){
        ##remove anything after 5 bases under Q30 in 10bp window
        ##trimmed <- trimTailw(seqs, badQuality, qualityThreshold,
        ##                     round(qualityWindow/2))
        ## this step is not necessary at  all
        ## trim if 5 bases are below '0'(fred score 15) in a window of 10 bases
        ## trimmed <- trimTailw(seqs, 5, '+', 5)
        ## trimmed <- trimTailw(seqs, 5, '#', 5)
        ## this step is necessary because many shortreads functions work on ACGT only
        ##trimmed <- trimmed[width(trimmed) > 65]
        trimmed <- seqs
        trimmed <- trimmed[!grepl('N', sread(trimmed))]
        if(length(trimmed) > 0){
            trimmedSeqs <- sread(trimmed)
            trimmedqSeqs <- quality(quality(trimmed))
            names(trimmedSeqs) <- names(trimmedqSeqs) <- 
            sapply(sub("(.+) .+","\\1",ShortRead::id(trimmed)),
                   function(z){paste0(alias, "%", strsplit(z, "-")[[1]][2])})
        }
    }
    list(trimmedSeqs, trimmedqSeqs)
})

reads <- sapply(r, "[[", 1)
qualities <- sapply(r, "[[", 2)
R1Quality <- qualities[[1]]
rm(r)
gc()

reads.p <- trim_Ltr_side_reads(reads[[2]], primer, ltrbit)
stats.bore$LTRed <- length(reads.p)
## note very few reads left
## debugging trim_Ltr_side_reads()

reads.p = reads[[2]]
maxMisMatch=2
## inside trim_Ltr_side_reads()

stopifnot(class(reads.p) %in% "DNAStringSet")
stopifnot(!any(duplicated(names(reads.p))))
stopifnot(length(primer)==1)
stopifnot(length(ltrbit)==1)

submat1 <- nucleotideSubstitutionMatrix(match=1,
                                        mismatch=0,
                                        baseOnly=TRUE)
## p for primer
## search for primer from the beginning
aln.p <- pairwiseAlignment(pattern=subseq(reads.p, 1, 1+nchar(primer)),
                           subject=primer,
                           substitutionMatrix=submat1,
                           gapOpening = 0,
                           gapExtension = 1,
                           type="overlap")
aln.p.df <- PairwiseAlignmentsSingleSubject2DF(aln.p)

## l for ltrbit
## search for ltrbit fellowing primer
## note, for SCID trial, there are GGG between primer and ltr bit and hence 5
## for extra bases
aln.l <- pairwiseAlignment(pattern=subseq(reads.p, nchar(primer)+1, nchar(primer)+nchar(ltrbit)+1),
                           subject=ltrbit,
                           substitutionMatrix=submat1,
                           gapOpening = 0,
                           gapExtension = 1,
                           type="overlap")
aln.l.df <- PairwiseAlignmentsSingleSubject2DF(aln.l, shift=nchar(primer)-1)


## check primer

print(primer)
## [1] "CCTCGGG"

tdf <- as.data.frame(table(unname(as.character(subseq(reads.p, 1, nchar(primer))))))
head(tdf[order(-tdf$Freq),])

##        Var1  Freq
## 331 CCTCGGG 13494
## 330 CCTCGGC 11596
## 345 CCTCTTG   193
## 332 CCTCGGT    86
## 329 CCTCGGA    57
## 355 CCTGGGC    28
## most reads match primer

## check ltrbit
print(ltrbit)
## [1] "GGTCTTTCA"

tdf <- as.data.frame(table(unname(as.character(subseq(reads.p, nchar(primer)+1, nchar(primer)+nchar(ltrbit))))))
head(tdf[order(-tdf$Freq),])


##           Var1  Freq
## 476  CTCCCAAAG 12698
## 1008 GTGGAGGGT  2908
## 500  CTCCCAGAG  1461
## 551  CTCCGAAAG  1109
## 240  CAGGTTACG  1071
## 663  CTTCCAAAG   748

