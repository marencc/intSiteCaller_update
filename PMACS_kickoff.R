print("started analysis")

workingDir = "~/EAS/intSiteValidation"

dataDir = paste0(workingDir, "/Data")

alignDir = "~/EAS/intSiteValidation/analysis"

base = list("?", 5, 10, "GAAAATC", "TCTAGCA", "TGCTAGAGATTTTCCACACTGACTAAAAGGGTCT", "CGCCTGCTCGCTAAGCTCTNNNNNNNNNNNNCTCCGCTTAAGGGACT", "AGTCCCTTAAGCGGAG", 30, "~/EAS/intSiteValidation/Data/demultiplexedReps/&UNIF_CTRL.1_S0_L001_R1_001.fastq.gz", "~/EAS/intSiteValidation/Data/demultiplexedReps/&UNIF_CTRL.1_S0_L001_R2_001.fastq.gz", "~/EAS/intSiteValidation/UNIF_CTRL", "~/EAS/vectorSeqs/p746vector.fasta", 95, 5, 2500)

parameters = list()

allReps = read.csv("~/EAS/intSiteValidation/allRepsBCs.csv")

vectors = rep("foo", nrow(allReps))
vectors[c(1:nrow(allReps))]="~/EAS/vectorSeqs/p746vector.fasta"

for(i in c(1:nrow(allReps))){
  base[[7]] = as.character(allReps$Linker.Sequence[i])
  base[[10]] = paste0("~/EAS/intSiteValidation/Data/demultiplexedReps/&", allReps$Alias[i], "_S0_L001_R1_001.fastq.gz")
  base[[11]] = paste0("~/EAS/intSiteValidation/Data/demultiplexedReps/&", allReps$Alias[i], "_S0_L001_R2_001.fastq.gz")
  base[[12]] = paste0("~/EAS/intSiteValidation/", allReps$Alias[i])
  base[[13]] = vectors[i]
  parameters = append(parameters, list(base))
}

save(parameters, file="parameters.RData")

suppressWarnings(dir.create("logs"))

bushmanJobID = "intSiteValidation"
blatStartPort = 5560

save(bushmanJobID, file="~/EAS/intSiteValidation/bushmanJobID.RData")
save(blatStartPort, file="~/EAS/intSiteValidation/bushmanBlatStartPort.RData")

#demultiplex seqs
system(paste0('bsub -n1 -q plus -J "BushmanErrorCorrect_', bushmanJobID, '" -o logs/errorCorrectOutput.txt Rscript ~/EAS/PMACS_scripts/errorCorrectBC_PMACS.R'))
