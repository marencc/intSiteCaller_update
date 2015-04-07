from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *
from golay import *
import commands
import sys

seqs = list()

status, fileNumber = commands.getstatusoutput("echo $LSB_JOBINDEX")

fileNumber = str(fileNumber)

out_handle = open("Data/correctedI1-" + fileNumber + ".fasta","a")

for seq_record in SeqIO.parse("Data/trimmedI1-" + fileNumber + ".fasta", "fasta"):  
  res = decode(str(seq_record.seq))
  if res[0] != None:
    SeqIO.write(SeqRecord(Seq(res[0], SingleLetterAlphabet()), id=seq_record.id, description=""), out_handle, "fasta")