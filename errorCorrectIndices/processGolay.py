from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *
from golay import *
import commands
import sys


out_handle = open("Data/correctedI1-" + str(sys.argv[1]) + ".fasta","w")


for seq_record in SeqIO.parse("Data/trimmedI1-" + str(sys.argv[1]) + ".fasta", "fasta"):  
  res = decode(str(seq_record.seq))
  if res[0] != None:
    SeqIO.write(SeqRecord(Seq(res[0], SingleLetterAlphabet()), id=seq_record.id, description=""), out_handle, "fasta")

done_handle = open("Data/correctedI1-" + str(sys.argv[1]) + ".done","w")
done_handle.close()
