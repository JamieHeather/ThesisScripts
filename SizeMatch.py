# SizeMatch.py v1.1

### HISTORY/PURPOSE ###
# Want to randomly select from collapsed decombined data (.freq files) so as to size match to the size of the smallest sample

### RUNNING/INPUT ###
# Takes .freq file (or any comma-space delimited file of format [5-field ID, frequency]) 
# python SizeMatch.py FILE.freq SIZE            
  # e.g.: python SizeMatch.py somefile.freq 10000
  # batch e.g.: for i in *freq; do echo $i; python SizeMatch.py $i 10000; done

### OUTPUT ###
# Outputs same format file, but of DCRs randomly selected from the original population, and their new frequencies
# For now called .matched files, but can be altered by changing "suffix" attribute


from __future__ import division
import collections as coll
import random
import sys
import re


if (len(sys.argv) <> 3):
  print "Please supply a filename for a DCR file (with frequencies) and an integer sample size"
  print "e.g.: python SizeMatch.py FILE.freq 10000"
  sys.exit()
else:
  dcrfilename = str(sys.argv[1])
  dcrfile = open(dcrfilename, "rU")
  samplesize = int(sys.argv[2])
  
  
prefix = str(dcrfilename.split(".")[0])
  
suffix = ".matched"                                     ## output file extension can be specified here

uncollapsed = []

line_count = 0

for line in dcrfile:
  
  line_count += 1
  
  comma = [m.start() for m in re.finditer(',', line)]           #define commas, from which we define all our 

  dcr = str(line[:comma[4]])
  
  freq = int(line[comma[4]+2:-1])

  for n in range(0,freq):
    
    uncollapsed.append(dcr)


sampled = random.sample(uncollapsed, samplesize)

recollapsed = coll.defaultdict(int)

for s in sampled:
  
  recollapsed[s] += 1

outfile = open(prefix+suffix, "w")

for i in recollapsed:
  
  outline = i + ", " + str(recollapsed[i])
  
  print >> outfile, outline

outfile.close()

print '{0:,}'.format(line_count), "records scanned from", str(dcrfilename)
