print("started analysis")

workingDir = "~/EAS/primerIDRecoveryOpt2"

dataDir = paste0(workingDir, "/Data")

alignDir = "~/EAS/primerIDRecoveryOpt2/analysis"

maxCores=14

base = list("?", 5, 10, "GAAAATC", "TCTAGCA", "TGCTAGAGATTTTCCACACTGACTAAAAGGGTCT", "CGCCTGCTCGCTAAGCTCTNNNNNNNNNNNNCTCCGCTTAAGGGACT", TRUE, "AGTCCCTTAAGCGGAG", 30, "~/EAS/primerIDRecoveryOpt2/Data/demultiplexedReps/&UNIF_CTRL.1_S0_L001_R1_001.fastq.gz", "~/EAS/primerIDRecoveryOpt2/Data/demultiplexedReps/&UNIF_CTRL.1_S0_L001_R2_001.fastq.gz", TRUE, "~/EAS/primerIDRecoveryOpt2/UNIF_CTRL", "~/EAS/vectorSeqs/MJVector.fasta", 95, 5, 2500, TRUE, FALSE)

parameters = list()

allReps = read.csv("~/EAS/primerIDRecoveryOpt2/allRepsBCs.csv")

vectors = rep("foo", nrow(allReps))
vectors[c(1:nrow(allReps))]="~/EAS/vectorSeqs/MJVector.fasta"

for(i in c(1:nrow(allReps))){
  base[[7]] = as.character(allReps$Linker.Sequence[i])
  base[[11]] = paste0("~/EAS/primerIDRecoveryOpt2/Data/demultiplexedReps/&", allReps$Replicate[i], "_S0_L001_R1_001.fastq.gz")
  base[[12]] = paste0("~/EAS/primerIDRecoveryOpt2/Data/demultiplexedReps/&", allReps$Replicate[i], "_S0_L001_R2_001.fastq.gz")
  base[[14]] = paste0("~/EAS/primerIDRecoveryOpt2/", allReps$Alias[i])
  base[[15]] = vectors[i]
  parameters = append(parameters, list(base))
}

alignMethod = "BLAT"

save(parameters, file="parameters.RData")

#demultiplex seqs
system(paste0('bsub -n1 -q max_mem64 -J "BushmanDemultiplex" -o demultiplexOutput.txt Rscript ~/EAS/PMACS_scripts/demultiplex_PMACS.R'))

#trim seqs
system(paste0('bsub -n1 -q max_mem30 -w "BushmanDemultiplex" -J "BushmanTrimSeqs[1-', length(parameters), ']" -o trimOutput%I.txt Rscript ~/EAS/PMACS_scripts/trimSeqs_PMACS.R'))

system("sleep 30")

#post-trim processing, also kicks off alignment and int site calling jobs
system(paste0('bsub -n1 -q normal -w "BushmanTrimSeqs[1-', length(parameters), ']" -J "BushmanPostTrimProcessing" -o postTrimOutput.txt Rscript ~/EAS/PMACS_scripts/postTrimSeqs_PMACS.R ', length(parameters)))

