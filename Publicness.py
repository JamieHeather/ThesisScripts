# Publicness.py
# 1st Nov 2014

# Want to plot the 'publicness' of different rearrangments (i.e. in how many donors a given sequence appears)
## Built on CDR3_Publicness from HIV analysis
# However I want to try and do this not only for ALL CDR3s, but also try and do the same for actual rearrangments

# Tracks the number of three different sequences:
  # DCRs = decombinator assignations; effectively nucleotide sequences
  # CDRs = CDR3s amino acid sequences (conserved C to F)
  # DCDRs = CDR3s + the V and J from the DCR, i.e. the CDR3 region and it's V/J context

from __future__ import division
import os
import collections as coll
import itertools as it
import re
import matplotlib.pyplot as plt
import random
import numpy as np
import math
import time 

plt.rcParams.update({'font.size': 15})

savepath = "/home/jme/TCR/WRITE_UP/THESIS/WorkingPlots/" + time.strftime("%Y %m %d") + " "

def basename(filename):
  # given a .dcrcdr3 filename, returns the base (i.e. the individual)
  return re.split('[_.]', filename)[2]

def getdeets(line):
  # split a line of a .dcrcdr3 file into its components, returned as a list
    # output returned: 0 = DCR, 1 = CDR3, 2 = frequency 
  bits1 = line.split(":")
  dcr = bits1[0]
  bits2 = bits1[1].rstrip().split(", ")
  cdr = bits2[0]
  freq = int(bits2[1])
  return (dcr, cdr, freq)

def subsample(counter,samplesize):
  ############### REQUIRES 'random' package! ################
  # function to subsample a repertoire sample
  # sample must be in a Counter, in format key = sequence (CDR3/DCR etc) : frequency
  if samplesize > sum(counter.values()):
    print "Sample size is larger than dictionary to sample from"
    return
  tmparr = []
  outcount = coll.Counter()
  for seq in counter:
    for x in range(counter[seq]):
      tmparr.append(seq)
  sampled = random.sample(tmparr, samplesize)
  for s in sampled:
    outcount[s] += 1
  return(outcount)


all_a_files = [f for f in os.listdir(os.getcwd()) if f.endswith("a.dcrcdr3") and "HV" in f and "v1" not in f and "CD" not in f]
all_b_files = [f for f in os.listdir(os.getcwd()) if f.endswith("b.dcrcdr3") and "HV" in f and "v1" not in f and "CD" not in f]


all_a_files.sort()
all_b_files.sort()
all_files = all_a_files + all_b_files

all_a_cdr = []
all_a_dcr = []
all_a_dcdr = []
all_b_cdr = []
all_b_dcr = []
all_b_dcdr = []

print "Reading data in to dictionaries..."


# generate a separate donor-specific dictionary for both CDR3s and DCRs 
for f in all_a_files:  
  vars()[basename(f) + "_cdr"] = coll.Counter()
  vars()[basename(f)+"_dcr"] = coll.Counter()  
  vars()[basename(f) + "_dcdr"] = coll.Counter()
  all_a_cdr.append(basename(f) + "_cdr")
  all_a_dcr.append(basename(f) + "_dcr")
  all_a_dcdr.append(basename(f) + "_dcdr")
     
for f in all_b_files:  
  vars()[basename(f) + "_cdr"] = coll.Counter()
  vars()[basename(f)+"_dcr"] = coll.Counter()   
  vars()[basename(f) + "_dcdr"] = coll.Counter()
  all_b_cdr.append(basename(f) + "_cdr")
  all_b_dcr.append(basename(f) + "_dcr")
  all_b_dcdr.append(basename(f) + "_dcdr")

# read all data into those dictionaries
for d in all_files:
  
  thisfile = open(d, "rU")
  cdrout = vars()[basename(d)+"_cdr"]
  dcrout = vars()[basename(d)+"_dcr"]
  dcdrout = vars()[basename(d)+"_dcdr"]
  
  for line in thisfile:
    dcr, cdr, freq = getdeets(line)
    cdrout[cdr] += freq
    dcrout[dcr] += freq
    dcdr = dcr.split(", ")[:2] + [cdr]
    dcdrout[str(dcdr)] += freq
    
    
# for each sequence (DCR or CDR3) we need to know two things:
  # How many donors does that sequence appear in (x axis), and
  # what is the total frequency of that sequence (or average??????????)
  
a_cdr_ppl = coll.Counter()  # how many people (donors) we find these sequences in
b_cdr_ppl = coll.Counter()
a_dcr_ppl = coll.Counter()
b_dcr_ppl = coll.Counter()
a_dcdr_ppl = coll.Counter()  
b_dcdr_ppl = coll.Counter()

#a_cdr_freq = coll.Counter()  # what the cumulative frequency
#b_cdr_freq = coll.Counter()
#a_dcr_freq = coll.Counter()
#b_dcr_freq = coll.Counter()

# go through all donor specific dictionaries

for d in all_a_cdr:
  # each donor
  for s in vars()[d]:
    # each sequence
    a_cdr_ppl[s] += 1

for d in all_a_dcr:
  # each donor
  for s in vars()[d]:
    # each sequence
    a_dcr_ppl[s] += 1

for d in all_a_dcdr:
  # each donor
  for s in vars()[d]:
    # each sequence
    a_dcdr_ppl[s] += 1

for d in all_b_cdr:
  # each donor
  for s in vars()[d]:
    # each sequence
    b_cdr_ppl[s] += 1

for d in all_b_dcr:
  # each donor
  for s in vars()[d]:
    # each sequence
    b_dcr_ppl[s] += 1
    
for d in all_b_dcdr:
  # each donor
  for s in vars()[d]:
    # each sequence
    b_dcdr_ppl[s] += 1

# code to save all public seqs for comparison with BV data
#outf = open("PanAlphaVJCDR.publicness", "w")
#for c in a_dcdr_ppl:
  #print >> outf, c + "|" + str(a_dcdr_ppl[c])
#outf.close()

#sys.exit()    

# convert these dictionaries to y values for plotting (x values = number of donors, i.e. 1 to 10)

acdr = []
adcr = []
adcdr = []
bcdr = []
bdcr = []
bdcdr = []

xs = range(1,11)

for i in xs: 
  acdr.append(a_cdr_ppl.values().count(i))
  adcr.append(a_dcr_ppl.values().count(i))
  adcdr.append(a_dcdr_ppl.values().count(i))
  bcdr.append(b_cdr_ppl.values().count(i))
  bdcr.append(b_dcr_ppl.values().count(i))
  bdcdr.append(b_dcdr_ppl.values().count(i))
  
# Plotting

# alpha = green, beta = purple
# dcr = whole, dcdr = dashed, cdr3 = dotted

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

ax.plot(xs, adcr, color ="green", linewidth=3, alpha=.7, label="Alpha DCR")
ax.plot(xs, adcdr, color ="green", linestyle="dashed", linewidth=3, alpha=.7, label="Alpha VJ-CDR3")
ax.plot(xs, acdr, color ="green", linestyle="dotted", linewidth=3, alpha=.7, label="Alpha CDR3")
ax.plot(xs, bdcr, color ="purple", linewidth=3, alpha=.7, label="Beta DCR")
ax.plot(xs, bdcdr, color ="purple", linestyle="dashed", linewidth=3, alpha=.7, label="Beta VJ-CDR3")
ax.plot(xs, bcdr, color ="purple", linestyle="dotted", linewidth=3, alpha=.7, label="Beta CDR3")

ax.set_yscale('log')
plt.xlabel("Number of donors")
plt.ylabel("Frequency")

plt.legend(loc='upper right', scatterpoints = 1, prop={'size':15})
#plt.show()
#plt.savefig(savepath + " Publicness.svg", bbox_inches='tight')
plt.close()

adc8 = []
avj8 = []

for d in a_dcr_ppl:
  if a_dcr_ppl[d] == 8:
    adc8.append(d)
    
for vj in a_dcdr_ppl:
  if a_dcdr_ppl[vj] == 8:
    avj8.append(vj)
  
##############################

# Calculate overlap (by jaccard) of the various different groups

def jaccard(A, B): # takes 2 sets, outputs the jaccard index of them
    inter = set(A).intersection(B)
    union = set(A).union(B)
    return float(len(inter))/float(len(union))
  

# find out the smallest number of productive rearrangements any donor, so we can sample to that
mina = 1e7
minb = 1e7

for x in all_a_cdr:
  if len(vars()[x].keys()) < mina:
    mina = len(vars()[x].keys())
    
for x in all_a_cdr:
  if len(vars()[x].keys()) < mina:
    mina = len(vars()[x].keys())    

for x in all_b_cdr:
  if len(vars()[x].keys()) < minb:
    minb = len(vars()[x].keys())    

for x in all_b_dcr:
  if len(vars()[x].keys()) < minb:
    minb = len(vars()[x].keys())   
  

mina = int(round(mina, -3))
minb = int(round(minb, -3))

print "Productive Jaccards sampled to:"
print 'Alpha: {0:,}'.format(mina)
print 'Beta: {0:,}'.format(minb)
  
jadcr = []
sjadcr = []
for x in it.combinations(all_a_dcr, 2): 
  jadcr.append(jaccard(vars()[x[0]], vars()[x[1]]))
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), mina)
  #s2 = random.sample(dat2.keys(), mina)
  #sjadcr.append(jaccard(s1, s2))  
  sjadcr.append(jaccard(subsample(dat1, mina), subsample(dat2, mina)))  

jadcdr = []  
sjadcdr = []  
for x in it.combinations(all_a_dcdr, 2): 
  jadcdr.append(jaccard(vars()[x[0]], vars()[x[1]]))
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), mina)
  #s2 = random.sample(dat2.keys(), mina)
  #sjadcdr.append(jaccard(s1, s2))
  sjadcdr.append(jaccard(subsample(dat1, mina), subsample(dat2, mina)))
  
jacdr = []  
sjacdr = []  
for x in it.combinations(all_a_cdr, 2): 
  jacdr.append(jaccard(vars()[x[0]], vars()[x[1]]))
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), mina)
  #s2 = random.sample(dat2.keys(), mina)
  #sjacdr.append(jaccard(s1, s2))
  sjacdr.append(jaccard(subsample(dat1, mina), subsample(dat2, mina)))
  
jbdcr = []    
sjbdcr = []    
for x in it.combinations(all_b_dcr, 2): 
  jbdcr.append(jaccard(vars()[x[0]], vars()[x[1]]))
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), minb)
  #s2 = random.sample(dat2.keys(), minb)
  #sjbdcr.append(jaccard(s1, s2))
  sjbdcr.append(jaccard(subsample(dat1, minb), subsample(dat2, minb)))
  
jbdcdr = []
sjbdcdr = []
for x in it.combinations(all_b_dcdr, 2): 
  jbdcdr.append(jaccard(vars()[x[0]], vars()[x[1]]))
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), minb)
  #s2 = random.sample(dat2.keys(), minb)
  #sjbdcdr.append(jaccard(s1, s2))
  sjbdcdr.append(jaccard(subsample(dat1, minb), subsample(dat2, minb)))
  
jbcdr = []
sjbcdr = []
for x in it.combinations(all_b_cdr, 2): 
  jbcdr.append(jaccard(vars()[x[0]], vars()[x[1]]))
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), minb)
  #s2 = random.sample(dat2.keys(), minb)
  #sjbcdr.append(jaccard(s1, s2))
  sjbcdr.append(jaccard(subsample(dat1, minb), subsample(dat2, minb)))

def jitter(numb):  
  # returns a number with jitter
  val = 0.32  
  return(numb+random.uniform(-val, val))
  
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

ax.scatter([jitter(x) for x in [1]*45], jadcr, color="green", alpha=.5, marker="o", s=50)
ax.scatter([jitter(x) for x in [2]*45], jadcdr, color="green", alpha=.5, marker="*",s=50)
ax.scatter([jitter(x) for x in [3]*45], jacdr, color="green", alpha=.5, marker="x",s=50)
ax.scatter([jitter(x) for x in [4]*45], jbdcr, color="purple", alpha=.5, marker="o",s=50)
ax.scatter([jitter(x) for x in [5]*45], jbdcdr, color="purple", alpha=.5, marker="*",s=50)
ax.scatter([jitter(x) for x in [6]*45], jbcdr, color="purple", alpha=.5, marker="x", s=50)

ax.scatter([jitter(x) for x in [1]*45], sjadcr, color="cyan", alpha=.5, marker="o", s=50)
ax.scatter([jitter(x) for x in [2]*45], sjadcdr, color="cyan", alpha=.5, marker="*", s=50)
ax.scatter([jitter(x) for x in [3]*45], sjacdr, color="cyan", alpha=.5, marker="x", s=50)
ax.scatter([jitter(x) for x in [4]*45], sjbdcr, color="goldenrod", alpha=.5, marker="o", s=50)
ax.scatter([jitter(x) for x in [5]*45], sjbdcdr, color="goldenrod", alpha=.5, marker="*", s=50)
ax.scatter([jitter(x) for x in [6]*45], sjbcdr, color="goldenrod", alpha=.5, marker="x", s=50)

#ax.scatter(-1,1, color="green", alpha=.5, marker="s", label="Whole Alpha", s=50)
#ax.scatter(-1,1, color="cyan", alpha=.5, marker="s", label="Subsampled Alpha", s=50)
#ax.scatter(-1,1, color="purple", alpha=.5, marker="s", label="Whole Beta", s=50)
#ax.scatter(-1,1, color="goldenrod", alpha=.5, marker="s", label="Subsampled Beta", s=50)
ax.scatter(-1,1, color="black", alpha=.5, marker="o", label="DCR", s=50)
ax.scatter(-1,1, color="black", alpha=.5, marker="*", label="VJ-CDR3", s=50)
ax.scatter(-1,1, color="black", alpha=.5, marker="x", label="CDR3", s=50)

ax.set_yscale('log')
ax.set_xticks([2,5])
ax.set_xticklabels(['Alpha', 'Beta'])

cols = ['green', 'cyan', 'purple', 'goldenrod']
gps = ['Whole Alpha', 'Subsampled Alpha', 'Whole Beta', 'Subsampled Beta']
rex = []

for g in range(len(gps)):
  pr = plt.Rectangle((0, 0), 1, 1, fc=cols[g])
  rex.append(pr)
  
legend1 = plt.legend(rex, gps,loc='lower left', prop={'size':14})

plt.ylabel('Jaccard index')
plt.ylim(1e-5,1e-1)
plt.xlim(0,7)
plt.legend(loc='lower right', prop={'size':14})
plt.gca().add_artist(legend1)
#plt.show()
#plt.savefig(savepath + " Productive Jaccards.svg", bbox_inches='tight')
plt.close()

#sys.exit()

# NB mean alpha cdr3 jaccard values are about equivalent to one shared sequence if both had seventeen unique sequences
# eg:
# jaccard([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],[31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,1])
#sys.exit()

# Bit of optional statistics

#import scipy.stats as ss

## all alphas are signif diff
#ss.mannwhitneyu(jadcr,jadcdr)
#(57.0, 6.4625259532266893e-15)

#ss.ttest_ind(jadcr,jadcdr)
#(array(-10.646907112225716), 1.7072236204062232e-17)

#ss.mannwhitneyu(jacdr,jadcdr)
#(46.0, 3.2120027273423209e-15)

#ss.ttest_ind(jacdr,jadcdr)
#(array(11.889845248362867), 5.3010138170968135e-20)

## as are all betas
#ss.mannwhitneyu(jbdcr,jbdcdr)
#(52.0, 4.7076924170404701e-15)

#ss.ttest_ind(jbdcr,jbdcdr)
#(array(-7.15159270880985), 2.4125372627634777e-10)

#ss.mannwhitneyu(jbcdr,jbdcdr)
#(6.0, 2.3673292361000783e-16)

#ss.ttest_ind(jbcdr,jbdcdr)
#(array(22.607392876514435), 2.0397978611303692e-38)

## and most shared beta (cdr3) is signif diff to least shared alpha (dcr)
#ss.mannwhitneyu(jadcr,jbcdr)
#(516.0, 3.1328872785854394e-05)

#ss.ttest_ind(jadcr,jbcdr)
#(array(4.0296262042783075), 0.00011836323599804811)

###########################

# As have been going off the .dcrcdr3 files, everything I've plotted so far is for functional rearrangments
  # As of this week, CDR3ulator also outputs non-productive rearrangements into separate files, so I can compare those too

all_a_npfiles = [f for f in os.listdir(os.getcwd()) if f.endswith("a.np") and "HV" in f and "v1" not in f and "CD" not in f]
all_b_npfiles = [f for f in os.listdir(os.getcwd()) if f.endswith("b.np") and "HV" in f and "v1" not in f and "CD" not in f]


all_a_npfiles.sort()
all_b_npfiles.sort()
all_npfiles = all_a_npfiles + all_b_npfiles

all_a_npdcr = []
all_b_npdcr = []

print "Reading non-productive data in to dictionaries..."

# generate a separate donor-specific dictionary for both CDR3s and DCRs 
for f in all_a_npfiles:  
  vars()[basename(f)+"_npdcr"] = coll.Counter()  
  all_a_npdcr.append(basename(f) + "_npdcr")
     
for f in all_b_files:  
  vars()[basename(f)+"_npdcr"] = coll.Counter()   
  all_b_npdcr.append(basename(f) + "_npdcr")

# read all data into those dictionaries
for d in all_npfiles:
  
  thisfile = open(d, "rU")
  npdcrout = vars()[basename(d)+"_npdcr"]

  for line in thisfile:
    dcr, fail, freq = getdeets(line)
    npdcrout[dcr] += freq
  
    
a_npdcr_ppl = coll.Counter()
b_npdcr_ppl = coll.Counter()

# go through all donor specific dictionaries

for d in all_a_npdcr:
  # each donor
  for s in vars()[d]:
    # each sequence
    a_npdcr_ppl[s] += 1

for d in all_b_npdcr:
  # each donor
  for s in vars()[d]:
    # each sequence
    b_npdcr_ppl[s] += 1
    
# generate plotting points    

anpdcr = []
bnpdcr = []

for i in xs: 
  anpdcr.append(a_npdcr_ppl.values().count(i))
  bnpdcr.append(b_npdcr_ppl.values().count(i))
      

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

ax.plot(xs, adcr, color ="green", linewidth=3, alpha=.7, label="Alpha P DCR")
ax.plot(xs, anpdcr, color ="green", linewidth=3, linestyle='-.', alpha=.7, label="Alpha NP DCR")
ax.plot(xs, bdcr, color ="purple", linewidth=3, alpha=.7, label="Beta P DCR")
ax.plot(xs, bnpdcr, color ="purple", linewidth=3, linestyle='-.', alpha=.7, label="Beta NP DCR")
ax.set_yscale('log')
plt.xlabel("Number of donors")
plt.ylabel("Frequency")

plt.legend(loc='upper right', scatterpoints = 1, prop={'size':15})
#plt.show()
#plt.savefig(savepath + " NP Publicness.svg", bbox_inches='tight')
plt.close()    
    
    
#####



janpdcr = []
jbnpdcr = []

for x in it.combinations(all_a_npdcr, 2): 
  janpdcr.append(jaccard(vars()[x[0]], vars()[x[1]]))
    
for x in it.combinations(all_b_npdcr, 2): 
  jbnpdcr.append(jaccard(vars()[x[0]], vars()[x[1]]))
    

#fig = plt.figure(figsize=(7,5))
#ax = fig.add_subplot(111)

#ax.scatter([jitter(x) for x in [1]*45], jadcr, color="green", alpha=.5, marker="o", label="Alpha P DCR", s=50)
#ax.scatter([jitter(x) for x in [2]*45], janpdcr, color="green", alpha=.5, marker="^", label="Alpha NP DCR", s=50)
#ax.scatter([jitter(x) for x in [3]*45], jbdcr, color="purple", alpha=.5, marker="o", label="Beta P DCR", s=50)
#ax.scatter([jitter(x) for x in [4]*45], jbnpdcr, color="purple", alpha=.5, marker="^", label="Beta NP DCR", s=50)
#ax.set_yscale('log')

#ax.set_xticks([1.5,3.5])
#ax.set_xticklabels(['Alpha', 'Beta'])

#plt.ylabel('Jaccard index')
#plt.ylim(1e-5,1e-1)
#plt.legend(loc='lower left', prop={'size':15})
#plt.show()
###plt.savefig(savepath + " Jaccards.svg", bbox_inches='tight')
##plt.close()
    
    
    
    
# find out the smallest number of non-productive in any donor, so we can sample to that
mina = 1e6
minb = 1e6

for x in all_a_npdcr:
  if len(vars()[x].keys()) < mina:
    mina = len(vars()[x].keys())
    
for x in all_b_npdcr:
  if len(vars()[x].keys()) < minb:
    minb = len(vars()[x].keys()) 
 

mina = int(round(mina, -3))
minb = int(round(minb, -1))

print "Non-productive Jaccards sampled to:"
print 'Alpha: {0:,}'.format(mina)
print 'Beta: {0:,}'.format(minb) 

#sys.exit()
 
sjadcr = []
for x in it.combinations(all_a_dcr, 2): 
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), mina)
  #s2 = random.sample(dat2.keys(), mina)
  #sjadcr.append(jaccard(s1, s2))
  sjadcr.append(jaccard(subsample(dat1, minb), subsample(dat2, minb)))
    
sjbdcr = []    
for x in it.combinations(all_b_dcr, 2):
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  #s1 = random.sample(dat1.keys(), minb)
  #s2 = random.sample(dat2.keys(), minb)
  #sjbdcr.append(jaccard(s1, s2))
  sjbdcr.append(jaccard(subsample(dat1, minb), subsample(dat2, minb)))

sjanpdcr = []
for x in it.combinations(all_a_npdcr, 2): 
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  if len(dat1.keys()) >= mina and len(dat2.keys()) >= mina:
    #s1 = random.sample(dat1.keys(), mina)
    #s2 = random.sample(dat2.keys(), mina)
    #sjanpdcr.append(jaccard(s1, s2))
    sjanpdcr.append(jaccard(subsample(dat1, minb), subsample(dat2, minb)))

sjbnpdcr = []    
for x in it.combinations(all_b_npdcr, 2):
  dat1 = vars()[x[0]]
  dat2 = vars()[x[1]]
  if len(dat1.keys()) >= minb and len(dat2.keys()) >= minb:
    #s1 = random.sample(dat1.keys(), minb)
    #s2 = random.sample(dat2.keys(), minb)
    #sjbnpdcr.append(jaccard(s1, s2))
    sjbnpdcr.append(jaccard(subsample(dat1, minb), subsample(dat2, minb)))
  else:
    print x, len(dat1.keys()), len(dat2.keys())
        

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

ax.scatter([jitter(x) for x in [1]*len(jadcr)], jadcr, color="green", alpha=.5, marker="o", s=50)
ax.scatter([jitter(x) for x in [1]*len(sjadcr)], sjadcr, color="cyan", alpha=.5, marker="o", s=50)

ax.scatter([jitter(x) for x in [2]*len(janpdcr)], janpdcr, color="green", alpha=.5, marker="^", s=50)
ax.scatter([jitter(x) for x in [2]*len(sjanpdcr)], sjanpdcr, color="cyan", alpha=.5, marker="^", s=50)

ax.scatter([jitter(x) for x in [3]*len(jbdcr)], jbdcr, color="purple", alpha=.5, marker="o", s=50)
ax.scatter([jitter(x) for x in [3]*len(sjbdcr)], sjbdcr, color="goldenrod", alpha=.5, marker="o", s=50)

ax.scatter([jitter(x) for x in [4]*len(jbnpdcr)], jbnpdcr, color="purple", alpha=.5, marker="^", s=50)
ax.scatter([jitter(x) for x in [4]*len(sjbnpdcr)], sjbnpdcr, color="goldenrod", alpha=.5, marker="^", s=50)

s1 = ax.scatter(-1,1, color="black", alpha=.5, marker="o", label="P DCR", s=50)
s2 = ax.scatter(-1,1, color="black", alpha=.5, marker="^", label="NP DCR", s=50)

cols = ['green', 'cyan', 'purple', 'goldenrod']
gps = ['Whole Alpha', 'Subsampled Alpha', 'Whole Beta', 'Subsampled Beta']
rex = []

for g in range(len(gps)):
  pr = plt.Rectangle((0, 0), 1, 1, fc=cols[g])
  rex.append(pr)
  
legend1 = plt.legend(rex, gps,loc='lower left', prop={'size':14})

ax.set_yscale('log')
ax.set_xticks([1.5,3.5])
ax.set_xticklabels(['Alpha', 'Beta'])
plt.xlim(0,5)
plt.ylabel('Jaccard index')
plt.ylim(2e-6,1e-1)
plt.legend(loc='lower right', prop={'size':14})
plt.gca().add_artist(legend1)
#plt.show()
#plt.savefig(savepath + " NP Jaccards.svg", bbox_inches='tight')
plt.close()

#ax.scatter(-1,1, color="green", alpha=.5, marker="s", label="Whole Alpha", s=50)
#ax.scatter(-1,1, color="cyan", alpha=.5, marker="s", label="Subsampled Alpha", s=50)
#ax.scatter(-1,1, color="purple", alpha=.5, marker="s", label="Whole Beta", s=50)
#ax.scatter(-1,1, color="goldenrod", alpha=.5, marker="s", label="Subsampled Beta", s=50)
    
# and the stats

#sys.exit()
import scipy.stats as ss

## whole NP > P
  ## alpha
#ss.mannwhitneyu(jadcr,janpdcr)
#(2.0, 1.8136507938454002e-16)

#ss.ttest_ind(jadcr,janpdcr)
#(array(-15.759274516435436), 2.3484398258831157e-27)

  ## beta
#ss.mannwhitneyu(jbdcr,jbnpdcr)
#(783.5, 0.032061428379470965)

#ss.ttest_ind(jbdcr,jbnpdcr)
#(array(-4.228429453150417), 5.7410582221270999e-05)

## subsampled NP > P
  ## alpha
#ss.mannwhitneyu(sjadcr,sjanpdcr)
#(116.0, 2.2179722451897189e-13)

#ss.ttest_ind(sjadcr,sjanpdcr)
#(array(-9.5558568564929), 2.9419702988373262e-15)

  ## beta
#ss.mannwhitneyu(sjbdcr,sjbnpdcr)
#(810.0, 0.0008705444715759209)

#ss.ttest_ind(sjbdcr,sjbnpdcr)
#(array(-2.730237230296352), 0.0076442061743619526)


################

# Want to know whether the more public clones are more frequent in their respective donors

# generate empty lists for lists of frequencies
for i in range(1,11):
  vars()['a_dcr_f' + str(i)] = []
  vars()['a_dcdr_f' + str(i)] = []
  vars()['a_cdr_f' + str(i)] = []
  vars()['b_dcr_f' + str(i)] = []
  vars()['b_dcdr_f' + str(i)] = []
  vars()['b_cdr_f' + str(i)] = []
  vars()['a_npdcr_f' + str(i)] = []
  vars()['b_npdcr_f' + str(i)] = []
  
for seq in a_dcr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_a_dcr:
    if s in vars()[d]:
      vars()['a_dcr_f' + str(ppl)].append(vars()[d][s])

for seq in a_dcdr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_a_dcdr:
    if s in vars()[d]:
      vars()['a_dcdr_f' + str(ppl)].append(vars()[d][s])

for seq in a_cdr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_a_cdr:
    if s in vars()[d]:
      vars()['a_cdr_f' + str(ppl)].append(vars()[d][s])

for seq in b_dcr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_b_dcr:
    if s in vars()[d]:
      vars()['b_dcr_f' + str(ppl)].append(vars()[d][s])

for seq in b_dcdr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_b_dcdr:
    if s in vars()[d]:
      vars()['b_dcdr_f' + str(ppl)].append(vars()[d][s])

for seq in b_cdr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_b_cdr:
    if s in vars()[d]:
      vars()['b_cdr_f' + str(ppl)].append(vars()[d][s])
      
for seq in a_npdcr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_a_npdcr:
    if s in vars()[d]:
      vars()['a_npdcr_f' + str(ppl)].append(vars()[d][s])
      
for seq in b_npdcr_ppl.most_common()[::-1]:
  s = seq[0]
  ppl = seq[1]
  for d in all_b_npdcr:
    if s in vars()[d]:
      vars()['b_npdcr_f' + str(ppl)].append(vars()[d][s])

a_dcr_fmeans = []
a_dcdr_fmeans = []
a_cdr_fmeans = []
b_dcr_fmeans = []
b_dcdr_fmeans = []
b_cdr_fmeans = []

a_dcr_fstd = []
a_dcdr_fstd = []
a_cdr_fstd = []
b_dcr_fstd = []
b_dcdr_fstd = []
b_cdr_fstd = []

a_npdcr_fmeans = []
b_npdcr_fmeans = []
a_npdcr_fstd = []
b_npdcr_fstd = []


for i in range(1,11):
  a_dcr_fmeans.append(np.mean(vars()['a_dcr_f' + str(i)]))
  a_dcdr_fmeans.append(np.mean(vars()['a_dcdr_f' + str(i)]))
  a_cdr_fmeans.append(np.mean(vars()['a_cdr_f' + str(i)]))
  b_dcr_fmeans.append(np.mean(vars()['b_dcr_f' + str(i)]))
  b_dcdr_fmeans.append(np.mean(vars()['b_dcdr_f' + str(i)]))
  b_cdr_fmeans.append(np.mean(vars()['b_cdr_f' + str(i)]))
  a_dcr_fstd.append(np.mean(vars()['a_dcr_f' + str(i)]))
  a_dcdr_fstd.append(np.mean(vars()['a_dcdr_f' + str(i)]))
  a_cdr_fstd.append(np.mean(vars()['a_cdr_f' + str(i)]))
  b_dcr_fstd.append(np.mean(vars()['b_dcr_f' + str(i)]))
  b_dcdr_fstd.append(np.mean(vars()['b_dcdr_f' + str(i)]))
  b_cdr_fstd.append(np.mean(vars()['b_cdr_f' + str(i)]))
  a_npdcr_fmeans.append(np.mean(vars()['a_npdcr_f' + str(i)]))
  b_npdcr_fmeans.append(np.mean(vars()['b_npdcr_f' + str(i)]))
  a_npdcr_fstd.append(np.std(vars()['a_npdcr_f' + str(i)]))
  b_npdcr_fstd.append(np.std(vars()['b_npdcr_f' + str(i)]))


# plot av freqs of productive dcr/vj/cdr seqs

xs = np.arange(1,11)
width=.16

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

ax.set_xticks(xs)
ax.set_xticklabels(xs, rotation=0)

ad = plt.bar(xs-(3*width), a_dcr_fmeans, width, color="yellowgreen", bottom=1e0, label="Alpha DCR", log=True)
adc = plt.bar(xs-(2*width), a_dcdr_fmeans, width, color="green", bottom=1e0, label="Alpha VJ-CDR3", log=True)
ac = plt.bar(xs-(1*width), a_cdr_fmeans, width, color="lime", bottom=1e0, label="Alpha CDR3", log=True)

bd = plt.bar(xs, b_dcr_fmeans, width, color="violet", bottom=1e0, label="Beta DCR", log=True)
bdc = plt.bar(xs+(1*width), b_dcdr_fmeans, width, color="purple", bottom=1e0, label="Beta VJ-CDR3", log=True)
bc = plt.bar(xs+(2*width), b_cdr_fmeans, width, color="indigo", bottom=1e0, label="Beta CDR3", log=True)

ax.errorbar(xs-(3*width)+(width/2), a_dcr_fmeans, yerr=a_dcr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)
ax.errorbar(xs-(2*width)+(width/2), a_dcdr_fmeans, yerr=a_dcdr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)
ax.errorbar(xs-(1*width)+(width/2), a_cdr_fmeans, yerr=a_cdr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)

ax.errorbar(xs+(width/2), b_dcr_fmeans, yerr=b_dcr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)
ax.errorbar(xs+(1*width)+(width/2), b_dcdr_fmeans, yerr=b_dcdr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)
ax.errorbar(xs+(2*width)+(width/2), b_cdr_fmeans, yerr=b_cdr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)

plt.legend(loc='upper left', prop={'size':15})
plt.ylabel("Mean frequency")
plt.xlabel("Number of donors")
#plt.show()
#plt.savefig(savepath + " Avg Freqs.svg", bbox_inches='tight')
plt.close()



# plot av freqs for p/np dcrs


xs = np.arange(1,11)
width=.2

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

ax.set_xticks(xs)
ax.set_xticklabels(xs, rotation=0)

ap = plt.bar(xs-(2*width), a_dcr_fmeans, width, color="green", bottom=1e0, label="Alpha P DCR", log=True)
anp = plt.bar(xs-(1*width), a_npdcr_fmeans, width, color="springgreen", bottom=1e0, label="Alpha NP DCR", log=True)

bp = plt.bar(xs, b_dcr_fmeans, width, color="purple", bottom=1e0, label="Beta P DCR", log=True)
bnp = plt.bar(xs+width, b_npdcr_fmeans, width, color="rosybrown", bottom=1e0, label="Beta NP DCR", log=True)

ax.errorbar(xs-(2*width)+(width/2), a_dcr_fmeans, yerr=a_dcr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)
ax.errorbar(xs-(1*width)+(width/2), a_npdcr_fmeans, yerr=a_npdcr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)
ax.errorbar(xs+(width/2), b_dcr_fmeans, yerr=b_dcr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)
ax.errorbar(xs+(1*width)+(width/2), b_npdcr_fmeans, yerr=b_npdcr_fstd, fmt="none", ecolor="black", lwd=2, capthick=2, alpha=.5)

for ba in range(len(anp)):
  anp[ba].set_hatch('*')


for ba in range(len(anp)):
  bnp[ba].set_hatch('*')


plt.ylim(1e0, 1e3)
#plt.ylim(0, 300)
plt.ylabel("Mean frequency")
plt.xlabel("Number of donors")
plt.legend(loc='upper left', prop={'size':15})
#plt.show()
#plt.savefig(savepath + " NP Avg Freqs.svg", bbox_inches='tight')
plt.close()

#sys.exit()

## stats for frequencies (just p vals)


#def sig(val):
  ## Given a p-value, returns an appropriate asterisk significance indicator
  ##if val <= 0.0001:	# optional additional value, not sure that depth is needed for this analysis
    ##return "****"
  #if val <= 0.001:
    #return "***"
  #elif val <= 0.01:
    #return "**"
  #elif val <= 0.05:
    #return "*"
  #else:
    #return "ns"

## productives

#mwu_a_dcr = []
#mwu_a_dcdr = []
#mwu_a_cdr = []
#mwu_b_dcr = []
#mwu_b_dcdr = []
#mwu_b_cdr = []

#tt_a_dcr = []
#tt_a_dcdr = []
#tt_a_cdr = []
#tt_b_dcr = []
#tt_b_dcdr = []
#tt_b_cdr = []

#for i in range(2,11): 
  #if vars()['a_dcr_f'+str(i)]:
    #mwu_a_dcr.append(ss.mannwhitneyu(vars()['a_dcr_f'+str(i)], vars()['a_dcr_f'+str(i-1)])[1])
    #tt_a_dcr.append(ss.ttest_ind(vars()['a_dcr_f'+str(i)], vars()['a_dcr_f'+str(i-1)])[1])
  #if vars()['a_dcdr_f'+str(i)]:
    #mwu_a_dcdr.append(ss.mannwhitneyu(vars()['a_dcdr_f'+str(i)], vars()['a_dcdr_f'+str(i-1)])[1])
    #tt_a_dcdr.append(ss.ttest_ind(vars()['a_dcdr_f'+str(i)], vars()['a_dcdr_f'+str(i-1)])[1])
  #if vars()['a_cdr_f'+str(i)]:
    #mwu_a_cdr.append(ss.mannwhitneyu(vars()['a_cdr_f'+str(i)], vars()['a_cdr_f'+str(i-1)])[1])
    #tt_a_cdr.append(ss.ttest_ind(vars()['a_cdr_f'+str(i)], vars()['a_cdr_f'+str(i-1)])[1])
  #if vars()['b_dcr_f'+str(i)]:
    #mwu_b_dcr.append(ss.mannwhitneyu(vars()['b_dcr_f'+str(i)], vars()['b_dcr_f'+str(i-1)])[1])
    #tt_b_dcr.append(ss.ttest_ind(vars()['b_dcr_f'+str(i)], vars()['b_dcr_f'+str(i-1)])[1])
  #if vars()['b_dcdr_f'+str(i)]:
    #mwu_b_dcdr.append(ss.mannwhitneyu(vars()['b_dcdr_f'+str(i)], vars()['b_dcdr_f'+str(i-1)])[1])
    #tt_b_dcdr.append(ss.ttest_ind(vars()['b_dcdr_f'+str(i)], vars()['b_dcdr_f'+str(i-1)])[1])
  #if vars()['b_cdr_f'+str(i)]:
    #mwu_b_cdr.append(ss.mannwhitneyu(vars()['b_cdr_f'+str(i)], vars()['b_cdr_f'+str(i-1)])[1])
    #tt_b_cdr.append(ss.ttest_ind(vars()['b_cdr_f'+str(i)], vars()['b_cdr_f'+str(i-1)])[1])
    
#[sig(x) for x in mwu_a_dcr]
#['***', '***', '***', '***', 'ns', '***', '***', '**', '***']
#[3.2404222846199139e-05, 2.0821579010730027e-05, 0.00025791899506505704, 0.00048314613164768911, 0.42991370860839861, 4.1078701231273605e-05, 0.00058690959063189185, 0.0090431280836152702, 3.9431945396284627e-08]

#[sig(x) for x in mwu_a_dcdr]
#['***', '***', '**', '***', 'ns', '**', '*', '***', '**']
#[1.0743650081151436e-13, 1.1230745645625631e-05, 0.008896189129730234, 4.4090763516607324e-06, 0.2370794168161045, 0.0017094028864595521, 0.048686468873711454, 2.1741324446741167e-08, 0.0023127918838069606]

#[sig(x) for x in mwu_a_cdr]
#['***', '***', '***', '***', '***', '***', '**', '***', '***']
#[4.3817618157330213e-17, 9.5740556350326621e-08, 2.0260762351598012e-08, 0.0002122830439529041, 3.7124984976836034e-05, 0.0005397323882328324, 0.0011318119069219945, 4.5544582739162072e-09, 1.7569344821094949e-06]

#[sig(x) for x in mwu_b_dcr]
#['**', '***', 'ns', 'ns']
#[0.003247040394309579, 0.00031281792188566428, 0.37629240786407625, 0.28745131493109255]

#[sig(x) for x in mwu_b_dcdr]
#['**', '***', 'ns', '**']
#[0.00180421793987902, 0.00040685468489148023, 0.43212812721865396, 0.008715686876500997]

#[sig(x) for x in mwu_b_cdr]
#['***', '***', '**', '*', 'ns', 'ns', 'ns', 'ns', '*']
#[2.3847587572789904e-06, 0.00079631496302019018, 0.0038749045789528585, 0.021069902970105098, 0.21417561152209441, 0.32416654179065674, 0.24585066149143819, 0.29640665038371122, 0.036012183039142998]

#[sig(x) for x in tt_a_dcr]
#['***', '*', 'ns', 'ns', 'ns', '***', 'ns', '***', '**']
#[2.7369686409905566e-07, 0.043859178192022161, 0.69331917583253566, 0.31707451447042978, 0.69388248304664357, 0.0006922138280932846, 0.13579200271814582, 5.0416329590372965e-05, 0.0070417332579105112]

#[sig(x) for x in tt_a_dcdr]
#['***', 'ns', 'ns', '***', 'ns', '*', 'ns', '***', '*']
#[4.0074762616052228e-07, 0.098120243351673284, 0.94036689770743076, 5.5184582352429936e-05, 0.43152979699318905, 0.033221983024975005, 0.060831725158149753, 0.00012724882673579267, 0.033331743635195227]

#[sig(x) for x in tt_a_cdr]
#['***', '*', 'ns', '**', '*', '***', '*', '***', '**']
#[2.4885850447570349e-07, 0.028624256779914477, 0.34261487130943136, 0.0039811459802480825, 0.010831040724799042, 0.00065211559043052379, 0.043870159358650032, 1.2198565539836151e-09, 0.0057818183264629732]

#[sig(x) for x in tt_b_dcr]
#['***', 'ns', 'ns', 'ns']
#[1.3517238265695259e-28, 0.75886309845916933, 0.17086328854643754, 0.79878111788554274]

#[sig(x) for x in tt_b_dcdr]
#['***', 'ns', 'ns', 'ns']
#[7.5977442872409216e-14, 0.81230231137570175, 0.10032861750556474, 0.49794572909844625]

#[sig(x) for x in tt_b_cdr]
#['***', 'ns', '**', 'ns', 'ns', 'ns', 'ns', 'ns', '*']
#[1.0129738268202184e-07, 0.35807299810637017, 0.0016253655767205423, 0.11419628290313764, 0.53892441476170472, 0.67185732719046043, 0.78708939594294713, 0.87375183129023926, 0.015073628885003272]


## or by one way anova
#ss.f_oneway(a_dcr_f1,a_dcr_f2,a_dcr_f3,a_dcr_f4,a_dcr_f5,a_dcr_f6,a_dcr_f7,a_dcr_f8,a_dcr_f9,a_dcr_f10)

#ss.f_oneway(a_dcdr_f1,a_dcdr_f2,a_dcdr_f3,a_dcdr_f4,a_dcdr_f5,a_dcdr_f6,a_dcdr_f7,a_dcdr_f8,a_dcdr_f9,a_dcdr_f10)

#ss.f_oneway(a_cdr_f1,a_cdr_f2,a_cdr_f3,a_cdr_f4,a_cdr_f5,a_cdr_f6,a_cdr_f7,a_cdr_f8,a_cdr_f9,a_cdr_f10)


#ss.f_oneway(b_dcr_f1,b_dcr_f2,b_dcr_f3,b_dcr_f4,b_dcr_f5)

#ss.f_oneway(b_dcdr_f1,b_dcdr_f2,b_dcdr_f3,b_dcdr_f4,b_dcdr_f5)

#ss.f_oneway(b_cdr_f1,b_cdr_f2,b_cdr_f3,b_cdr_f4,b_cdr_f5,b_cdr_f6,b_cdr_f7,b_cdr_f8,b_cdr_f9,b_cdr_f10)


# actually probably better to plot them as scatter plot, shows off the features of the distribution better

# productive cdr3s
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

for i in range(1,11):
  adat = vars()['a_cdr_f'+str(i)]
  ax.scatter(np.repeat(i, len(adat))-.1, adat, color='green', alpha=.2)
  ax.scatter(i-.1, np.mean(adat), color="green", marker="_", s=100, linewidth=5, alpha=.7)
  bdat = vars()['b_cdr_f'+str(i)]
  ax.scatter(np.repeat(i, len(bdat))+.1, bdat, color='purple', alpha=.2)
  ax.scatter(i+.1, np.mean(bdat), color="purple", marker="_", s=100, linewidth=5, alpha=.7)

ax.scatter(-1,-1, color='green', alpha=.5, label="Alpha")
ax.scatter(-1,-1, color='purple', alpha=.5, label="Beta")
xs = np.arange(1,11)
ax.set_xticks(xs)
ax.set_xticklabels(xs, rotation=0)
ax.set_yscale('log')
ax.set_xlim( [.5,10.5] )
ax.set_ylim( [9e-1,1e4] )
plt.legend(loc='upper center', prop={'size':15})
plt.xlabel("Number of donors")
plt.ylabel("CDR3 frequency")
#plt.show()
##plt.savefig(savepath + " CDR3 scatter publicness.svg", bbox_inches='tight')
#plt.savefig(savepath + " CDR3 scatter publicness.png", dpi=300, bbox_inches='tight')
plt.close()
 
# non-prod dcrs vs productive

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

for i in range(1,11):
  adat = vars()['a_dcr_f'+str(i)]
  ax.scatter(np.repeat(i, len(adat))-.225, adat, color='green', alpha=.2)
  ax.scatter(i-.225, np.mean(adat), color="green", marker="_", s=100, linewidth=5, alpha=.7)
  bdat = vars()['b_dcr_f'+str(i)]
  ax.scatter(np.repeat(i, len(bdat))-.075, bdat, color='purple', alpha=.2)
  ax.scatter(i-.075, np.mean(bdat), color="purple", marker="_", s=100, linewidth=5, alpha=.7)
  anpdat = vars()['a_npdcr_f'+str(i)]
  ax.scatter(np.repeat(i, len(anpdat))+.075, anpdat, color='springgreen', alpha=.2, marker="*")
  ax.scatter(i+.075, np.mean(anpdat), color="springgreen", marker="_", s=100, linewidth=5, alpha=.7)
  bnpdat = vars()['b_npdcr_f'+str(i)]
  ax.scatter(np.repeat(i, len(bnpdat))+.225, bnpdat, color='rosybrown', alpha=.2, marker="*")
  ax.scatter(i+.225, np.mean(bnpdat), color="rosybrown", marker="_", s=100, linewidth=5, alpha=.7)  

ax.scatter(-1,-1, color='green', alpha=.5, label="Alpha P DCR")
ax.scatter(-1,-1, color='springgreen', alpha=.5, label="Alpha NP DCR", marker="*")
ax.scatter(-1,-1, color='purple', alpha=.5, label="Beta P DCR")
ax.scatter(-1,-1, color='rosybrown', alpha=.5, label="Beta NP DCR", marker="*")

xs = np.arange(1,11)
ax.set_xticks(xs)
ax.set_xticklabels(xs, rotation=0)
ax.set_yscale('log')
ax.set_xlim( [.5,10.5] )
ax.set_ylim( [9e-1,1e4] )
plt.legend(loc='upper center', prop={'size':15})
plt.xlabel("Number of donors")
plt.ylabel("DCR frequency")
#plt.show()
##plt.savefig(savepath + " DCR scatter publicness.svg", bbox_inches='tight')
#plt.savefig(savepath + " DCR scatter NP publicness.png", dpi=300, bbox_inches='tight')

plt.close()
  
 
###################################################

sys.exit()

# looking for invariant sequences

# checking to see whether the potential new invariant is always encoded by the same rearrangement
  # would see same DCR in all donors
  
for x in a_dcr_ppl.most_common():
  b = x[0].split(", ")
  if b[0] == "6" and b[1] == "46" and x[1] >= 5:
    print x, travnam(int(b[0])), trajnam(int(b[1]))

# it's not!

