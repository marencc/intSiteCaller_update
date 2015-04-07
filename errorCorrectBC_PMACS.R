library("ShortRead")

metadata = read.csv("allRepsBCs.csv")

bushmanJobID = get(load("bushmanJobID.RData"))

I1 = readFastq("Data/Undetermined_S0_L001_I1_001.fastq.gz")
I1 = trimTailw(I1, 2, "0", 12)
I1 = I1[width(I1)==max(width(I1))]

I1 = split(I1, ceiling(seq_along(I1)/500000))

for(chunk in names(I1)){
  writeFasta(I1[[chunk]], file=paste0("Data/trimmedI1-", chunk, ".fasta"))
}

system(paste0('bsub -n1 -q normal -J "BushmanErrorCorrectWorker_', bushmanJobID, '[1-', length(I1),']" -o logs/errorCorrectWorkerOutput%I.txt python ~/EAS/PMACS_scripts/processGolay.py'))

#system(paste0('bsub -n1 -q normal -w "done(BushmanErrorCorrectWorker_', bushmanJobID, ')" -J "BushmanErrorCorrectMerge_', bushmanJobID, '" -o logs/errorCorrectMergeOutput.txt "cat Data/correctedI1-*.fasta > Data/correctedI1.fasta"'))

system(paste0('bsub -n1 -q max_mem64 -w "done(BushmanErrorCorrectWorker_', bushmanJobID, ')" -J "BushmanDemultiplex_', bushmanJobID, '" -o logs/demultiplexOutput.txt Rscript ~/EAS/PMACS_scripts/demultiplex_PMACS.R'))

#trim seqs
system(paste0('bsub -n1 -q plus -w "done(BushmanDemultiplex_', bushmanJobID, ')" -J "BushmanTrimSeqs_', bushmanJobID, '[1-', nrow(metadata), ']" -o logs/trimOutput%I.txt Rscript ~/EAS/PMACS_scripts/trimSeqs_PMACS.R'))

system("sleep 10")

#post-trim processing, also kicks off alignment and int site calling jobs
system(paste0('bsub -n1 -q normal -w "done(BushmanTrimSeqs_', bushmanJobID, ')" -J "BushmanPostTrimProcessing_', bushmanJobID, '" -o logs/postTrimOutput.txt Rscript ~/EAS/PMACS_scripts/postTrimSeqs_PMACS.R'))

