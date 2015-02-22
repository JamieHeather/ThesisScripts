# Expand.py v1.1

### HISTORY/PURPOSE ###
# Want to randomly select from collapsed decombined data (.freq files) in R 
# For that I need to take a nicely frequency-collapsed data file and re-write it out at one per line
# Built on SizeMatch.py; does the same thing just outputs earlier in the process
### RUNNING/INPUT ###
# Takes .freq file (or any comma-space delimited file of format [5-field ID, frequency]) 
# python Expand.py FILE.freq             

### OUTPUT ###
# Outputs same TCRs as in input file, but one row per collapsed frequency
# For now called .exp files, but can be altered by changing "suffix" attribute


from __future__ import division
import collections as coll
import random
import sys
import re


if (len(sys.argv) <> 2):
  print "Please supply a filename for a DCR file (with frequencies)"
  print "e.g.: python Expand.py FILE.freq "
  sys.exit()
else:
  dcrfilename = str(sys.argv[1])
  dcrfile = open(dcrfilename, "rU")
  
  
prefix = str(dcrfilename.split(".")[0])
  
suffix = ".exp"                                     ## output file extension can be specified here

uncollapsed = []

line_count = 0

for line in dcrfile:
  
  line_count += 1
  
  comma = [m.start() for m in re.finditer(',', line)]           #define commas, from which we define all our 

  dcr = str(line[:comma[4]])
  
  freq = int(line[comma[4]+2:-1])

  for n in range(0,freq):
    
    uncollapsed.append(dcr)


#sampled = random.sample(uncollapsed, samplesize)

#recollapsed = coll.defaultdict(int)

#for s in sampled:
  
  #recollapsed[s] += 1

outfile = open(prefix+suffix, "w")

for i in uncollapsed:
  
  outline = i 
  
  print >> outfile, outline

outfile.close()

print '{0:,}'.format(line_count), "records scanned from", str(dcrfilename)
