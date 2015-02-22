# FastqInterTag.fq v1.0

# 6th July 2014

### HISTORY/PURPOSE ###
# Produces fastq reads of just the inter-tag region (i.e. V-tag to J-tag)
# Allows comparison of quality metrics for the actual TCR portion of the sequence

### RUNNING/INPUT ###
# Takes the output of verbose_dcr or cmv_vDCR 
# python FastqInterTag.py FILE.n12

### OUTPUT ###
# A fastq file with the same name as the n12 file

import sys
import os

if (len(sys.argv) <> 2):
  print "Please supply a filename containing vDCR-decombined data, e.g. python CollapseTCRs.py FILENAME.n12"
  sys.exit()
else:
  dcrfilename = str(sys.argv[1])
  dcrfile = open(dcrfilename, "rU")
  
outfilename = dcrfilename.split(".")[0] + "_TCRs.fq"

if outfilename in os.listdir("./"):
  print "Output filename \'" + outfilename + "\' already in use. Delete and re-run"
  sys.exit()
#else:
  #outfile = open(outfilename, "rU")
  #outfile.close()

with open(outfilename, "a") as outfile:
  for line in dcrfile:
    bits = line.rstrip().split(", ")
    readid = "@" + bits[5]
    seq = bits[6]
    qual = bits[7]
    
    if len(seq) <> len(qual):
      print "Error - mismatch between sequence and quality string lengths for read", readid
      
    outread = readid + "\n" + seq + "\n+\n" + qual + "\n"
    
    #print outread
    outfile.write(outread)

#with open(fdf, "a") as a:
  #a.write(result)

  









