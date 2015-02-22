# PrintTotalAvgQuals.py v1

# Calculates the overall average quality of a sequence
  # For use on the output of FastqInterTag.py

# INPUT:
  # 1) (Gzipped) data piped in to script:
      # zcat data.fq.gz | python PrintTotalAvgQuals.py 
  # 2) Unzipped data
      # python PrintTotalAvgQuals.py data.fq
      
# OUTPUT:
  # Outputs number of reads, total mean quality and standard deviation to standard out

from Bio import SeqIO 
import sys            
import numpy as np
import collections as coll

###################

qual_dict = coll.Counter()
num_seqs = 0


if sys.stdin.isatty() == False:

  for record in SeqIO.parse(sys.stdin, "fastq"):
    
    quals = record.letter_annotations.values()[0]
    
    for b in quals:
      
      qual_dict[b] += 1
  
    num_seqs += 1
    
    if num_seqs % 1000 == 0:
      sys.stderr.write(str(num_seqs) + " | ")
      
elif sys.stdin.isatty() == True: 	# If don't have stdin input
  
  if (len(sys.argv) <> 2):
    print "Please supply a fastq filename (i.e. python script.py file.fastq) or pipe data in (i.e. cat file.fastq | python script.py"
    sys.exit()
  else:
    fastqfile = open(sys.argv[1], "rU")
  
  for record in SeqIO.parse(fastqfile,"fastq"):
    
    quals = record.letter_annotations.values()[0]
    
    for b in quals:
      
      qual_dict[b] += 1
  
    num_seqs += 1
    
    if num_seqs % 1000 == 0:
      sys.stderr.write(str(num_seqs) + " | ")

allquals = []

for q in qual_dict:
  allquals = allquals + ([q] * qual_dict[q])

print "Mean quality:\t" + str(np.mean(allquals))
print "Standard dev:\t" + str(np.std(allquals))
print "Number seqs:\t" + str(num_seqs)
