from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *
from golay import *
import commands
import sys
import os

seqs = list()

status, fileNumber = commands.getstatusoutput("echo $LSB_JOBINDEX")

fileNumber = str(fileNumber)

#out_handle = open("Data/correctedI1-" + fileNumber + ".fasta","a")
outfile = ( "Data/correctedI1-" + fileNumber + ".fasta" )
os.remove(outfile) if os.path.exists(outfile) else None
out_handle = open(outfile, "a")

for seq_record in SeqIO.parse("Data/trimmedI1-" + fileNumber + ".fasta", "fasta"):  
  res = decode(str(seq_record.seq))
  if res[0] != None:
    SeqIO.write(SeqRecord(Seq(res[0], SingleLetterAlphabet()), id=seq_record.id, description=""), out_handle, "fasta")
