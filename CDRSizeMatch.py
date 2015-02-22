  # SizeMatch.py v1.2 CDR3 version

### HISTORY/PURPOSE ###
# Want to randomly select from collapsed decombined data (.freq files) so as to size match to the size of the smallest sample

### RUNNING/INPUT ###
# Takes .freq file (or any comma-space delimited file of format [5-field ID, frequency]) 
# python SizeMatch.py FILE.freq SIZE            e.g. python SizeMatch.py somefile.freq 10000
# for i in *freq; do echo $i; python SizeMatch.py $i 10000; done

### OUTPUT ###
# Outputs same format file, but of DCRs randomly selected from the original population, and their new frequencies
# For now called .matched files, but can be altered by changing "suffix" attribute


from __future__ import division
import collections as coll
import random
import sys
import re

uniq = False
# Provides the option to output a sample from unique CDR3s, which is useful when it comes to size-matching for sharing indexes
  # (due to HIV repertoires' higher Ginis, a sampling will favour the larger clones, thus making them appear less 'shared')


if (len(sys.argv) <> 3):
  print "Please supply a filename for a DCR file (with frequencies) and an integer sample size"
  print "e.g.: python SizeMatch.py FILE.freq 10000"
  sys.exit()
else:
  cdrfilename = str(sys.argv[1])
  cdrfile = open(cdrfilename, "rU")
  samplesize = int(sys.argv[2])
  
print "Unique parameter is set to:", uniq
  
prefix = str(cdrfilename.split(".")[0])
  
suffix = ".cdrmtchd"                                     ## output file extension can be specified here

uncollapsed = []

for line in cdrfile:
  
  comma = [m.start() for m in re.finditer(',', line)]           #define commas, from which we define all our 

  cdr = str(line[:comma[0]])
  
  freq = int(line[comma[0]+2:-1])


  # By setting uniq parameter at top of script one can select between randomly sampling CDR3s from the 'uncollapsed' pool
    # i.e. CDR3s are present multiple times, equal to their frequency
  #... or from the pool of unique CDR3s
  
  if uniq == False:
    
    for n in range(0,freq):
    
      uncollapsed.append(cdr)
  
  elif uniq == True:
    
    uncollapsed.append(cdr)


sampled = random.sample(uncollapsed, samplesize)

recollapsed = coll.defaultdict(int)

for s in sampled:
  
  recollapsed[s] += 1

outfile = open(prefix+suffix, "w")

for i in recollapsed:
  
  outline = i + ", " + str(recollapsed[i])
  
  print >> outfile, outline

outfile.close()





