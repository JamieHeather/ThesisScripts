# Translational_Convergence.py v1.2
# 25th July 2014

# Built on /home/jme/TCR/Sequences/PBTPHIV/freq/CDR3/CDR3_Sharing.py
# Compares the translational convergence of healthy vs HIV v1 samples
  # i.e. how many DCRs are used to encode CDR3s (on average) within a given person
  
# Takes .dcrcdr3 as input (products of CDR3ulator with dcr_output set to True)
  # These files are effectively a combination of .freq and .cdr3 files
  # Each unique DCR has a line, with its CDR3, and frequency
  # NB CDR3s can thus appear on multiple lines
  
# Also note that this analysis (like all others for paper so far) has been run with CDR3ulator.py|include_gxg = False
  # This means that all CDR3 sequences do NOT include the GXG section of the 3' motif
  # This could slightly increase convergence compared to earlier analyses
    # That said, this probably won't overcome the decreased convergence we see due to increased stringency in collapsing
    
# Edited on 20th Aug to tweak the plotting for final figures - just removing inlay plot, make legend bigger
  # v1.2 for final figures = less whitespace, no titles, bigger fonts   
        
from __future__ import division
import os
import collections as coll
import itertools as it
import re
import matplotlib.pyplot as plt
import random
import time
import numpy as np
import scipy.stats as ss

plt.rcParams.update({'font.size': 15})

savepath = "/home/jme//TCR/WRITE_UP/THESIS/WorkingPlots/" + time.strftime("%Y %m %d") + " "

def basename(filename):
  # given a .dcrcdr3 filename, returns the base (i.e. the individual)
  if "HV" in filename:
    return re.split('[_.v]', filename)[3]
  elif "P" in filename:
    return re.split('[_.v]', filename)[3]
  
def getcdr(line):
  # split a line of a .dcrcdr3 file into its components, returned as a list
    # output returned: 0 = dcr, 1 = cdr, 2 = frequence (of that dcr)
  dcr = line.split(":")[0]
  cdrfreq = line.split(":")[1].split(", ")
  cdr = cdrfreq[0]
  freq = int(cdrfreq[1].rstrip())
  return (dcr, cdr, freq)

print "Getting data"

# Get all .dcrcdr3 files (contain (comma-delimited) decombinator index : CDR3 , frequency)

a_hv_files = [f for f in os.listdir(os.getcwd()) if 
                      f.endswith("a.dcrcdr3") and "HV" in f and "CD" not in f and "v1" not in f]
b_hv_files = [f for f in os.listdir(os.getcwd()) if 
                      f.endswith("b.dcrcdr3") and "HV" in f and "CD" not in f and "v1" not in f]
# NB use v2 of early HV files, as look most 'healthy' like, suspect v1 occurred during some seasonal infection

a_hiv_files = [f for f in os.listdir(os.getcwd()) if 
                  f.endswith("v1a.dcrcdr3") and "P0" in f]
b_hiv_files = [f for f in os.listdir(os.getcwd()) if 
                  f.endswith("v1b.dcrcdr3") and "P0" in f]


a_hv_files.sort()
b_hv_files.sort()

a_hiv_files.sort()
b_hiv_files.sort()

all_alpha_files = a_hv_files + a_hiv_files
all_beta_files = b_hv_files + b_hiv_files


# dictionaries to store how many DCRs each CDR3 is encoded byfa
alpha_dcrs = coll.defaultdict(list)
beta_dcrs = coll.defaultdict(list)

alpha_dcrs_c = coll.Counter()
beta_dcrs_c = coll.Counter()

print "Generating donor-specific dictionaries"

# Populate donor-specific dictionaries


#single_bleeds = ["HV07", "HV08", "HV09", "HV10", "HV11", "HV13"]

for d in all_alpha_files:
 
  # if basename(d) not in single_bleeds:
    
  thisfile = open(d, "rU")
  vars()[basename(d)+"a"] = coll.Counter()
  vars()[basename(d)+"a_dcr_c"] = coll.Counter()
  outdict = vars()[basename(d)+"a"]
  indiv_dcr_c_out = basename(d)+"a_dcr_c"      # 8th april indiv convergence counts
  
  for line in thisfile:
    
    cdretc = getcdr(line)
    
    dcr = cdretc[0]
    cdr = cdretc[1]
    freq = cdretc[2]
            
    outdict[cdr] += freq
    
    if dcr not in alpha_dcrs[cdr]:
      alpha_dcrs[cdr].append(dcr)
      alpha_dcrs_c[cdr] += 1
    
    if dcr not in vars()[indiv_dcr_c_out].keys():      # 8th april indiv convergence counts
      vars()[indiv_dcr_c_out][cdr] += 1
      
      
      
for d in all_beta_files:
 
  thisfile = open(d, "rU")
  vars()[basename(d)+"b"] = coll.Counter()
  vars()[basename(d)+"b_dcr_c"] = coll.Counter()
  outdict = vars()[basename(d)+"b"]
  indiv_dcr_c_out = basename(d)+"b_dcr_c"      # 8th april indiv convergence counts
  
  for line in thisfile:
    
    cdretc = getcdr(line)
    
    dcr = cdretc[0]
    cdr = cdretc[1]
    freq = cdretc[2]
            
    outdict[cdr] += freq
    
    if dcr not in beta_dcrs[cdr]:
      beta_dcrs[cdr].append(dcr)
      beta_dcrs_c[cdr] += 1
      
    if dcr not in vars()[indiv_dcr_c_out].keys():      # 8th april indiv convergence counts
      vars()[indiv_dcr_c_out][cdr] += 1   


# Defining groups donors

hva = [x for x in dir() if "HV" in x and "a" in x and not "dcr" in x]
hvb = [x for x in dir() if "HV" in x and "b" in x and not "dcr" in x]

hiva = [x for x in dir() if "P" in x and "a" in x and not "dcr" in x]
hivb = [x for x in dir() if "P" in x and "b" in x and not "dcr" in x]    
    
shiva = random.sample(hiva, len(hva))
shivb = random.sample(hivb, len(hvb))

# Calc convergence


setrange = 1500   # how many CDR3s to sample

hvax = []
hvay = []

rhvaconv = coll.defaultdict(list)       # list of the different convergences for each rank of cdr
rhivaconv = coll.defaultdict(list)       # list of the different convergences for each rank of cdr

for smpl in hva:
  for i in range(setrange):
    cdr = vars()[smpl+"_dcr_c"].most_common(i+1)[i][0]
    #hvax.append(i+1)
    #hvay.append(alpha_dcrs_c[cdr])
    rhvaconv[i+1].append(alpha_dcrs_c[cdr])

hivax = []
hivay = []

for smpl in shiva:
  for i in range(setrange):
    cdr = vars()[smpl+"_dcr_c"].most_common(i+1)[i][0]
    #hivax.append(i+1)
    #hivay.append(alpha_dcrs_c[cdr])
    rhivaconv[i+1].append(alpha_dcrs_c[cdr])

hvameans = []
hivameans = []

for i in range(setrange):
  hvameans.append(np.mean(rhvaconv[i+1]))
  hivameans.append(np.mean(rhivaconv[i+1]))

hvameanmeans = []
hivameanmeans = []
  
for i in range(setrange):
  hvameanmeans.append(np.mean(hvameans[i+1:i+100]))
  hivameanmeans.append(np.mean(hivameans[i+1:i+100]))


############################################################################
############################# PLOTTING #####################################
############################################################################


# NEED THIS TO ADJUST LEGEND FONTSIZE - ONLY RUN ONCE PER SCRIPT (defaults fixes it)
#plt.rc('legend',**{'fontsize':9})

#plt.rcdefaults()

inset_size = 50
####### ALPHA


fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
plt.xlabel("Ranked Frequency")
plt.ylabel("CDR3 translational convergence")
#plt.title("Alpha mean translational convergence")

ax.plot(range(setrange+1)[1:], hvameans, color="blue", alpha=.3, label="Healthy")
ax.plot(range(setrange+1)[1:], hivameans, color="red", alpha=.3, label="HIV")

#fake plots to get right color in legend
#ax.plot(0,0, color="blue", label="Healthy")
#ax.plot(0,0, color="red", label="HIV")

ax.plot(range(setrange+1)[1:], hvameanmeans, color="green", label="Healthy ranked means (100)", linewidth=2)
ax.plot(range(setrange+1)[1:], hivameanmeans, color="purple", label="HIV ranked means (100)", linewidth=2)

#ax2 = plt.axes([.5,.5,.37,.37])
#ax2.plot(range(setrange+1)[1:inset_size+1], hvameans[:inset_size], color="blue", alpha=.5, linewidth=2)
#ax2.plot(range(setrange+1)[1:inset_size+1], hivameans[:inset_size], color="red", alpha=.5, linewidth=2)

#ax2.plot(range(setrange+1)[1:inset_size+1], hvameanmeans[:inset_size], color="green", linewidth=2)
#ax2.plot(range(setrange+1)[1:inset_size+1], hivameanmeans[:inset_size], color="purple", linewidth=2)

ax.axis(xmin=0, xmax=1000, ymin=0, ymax=30)
#ax2.axis(xmin=1, xmax=inset_size, ymin=0, ymax=30)

#leg = ax.legend(loc=[.05,.8])
#for legobj in leg.legendHandles:
  #legobj.set_linewidth(3.0)

plt.legend(loc='upper right', prop={'size':13})
 
plt.tight_layout()
#plt.show()

plt.savefig(savepath + "Alpha translational convergence.svg", bbox_inches='tight')

  
  

######## BETA

rhvbconv = coll.defaultdict(list)       # list of the different convergences for each rank of cdr
rhivbconv = coll.defaultdict(list)       # list of the different convergences for each rank of cdr

for smpl in hvb:
  for i in range(setrange):
    cdr = vars()[smpl+"_dcr_c"].most_common(i+1)[i][0]
    #hvax.append(i+1)
    #hvay.append(beta_dcrs_c[cdr])
    rhvbconv[i+1].append(beta_dcrs_c[cdr])

hivbx = []
hivby = []

for smpl in shivb:
  for i in range(setrange):
    cdr = vars()[smpl+"_dcr_c"].most_common(i+1)[i][0]
    #hivbx.append(i+1)
    #hivby.append(beta_dcrs_c[cdr])
    rhivbconv[i+1].append(beta_dcrs_c[cdr])

hvbmeans = []
hivbmeans = []



for i in range(setrange):
  hvbmeans.append(np.mean(rhvbconv[i+1]))
  hivbmeans.append(np.mean(rhivbconv[i+1]))
  
hvbmeanmeans = []
hivbmeanmeans = []
    
for i in range(setrange):
  hvbmeanmeans.append(np.mean(hvbmeans[i+1:i+100]))
  hivbmeanmeans.append(np.mean(hivbmeans[i+1:i+100]))


  
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
plt.xlabel("Ranked Frequency")
plt.ylabel("CDR3 translational convergence")
#plt.title("Beta mean translational convergence")

ax.plot(range(setrange+1)[1:], hvbmeans, color="blue", alpha=.3, label="Healthy")
ax.plot(range(setrange+1)[1:], hivbmeans, color="red", alpha=.3, label="HIV")



#fake plots to get right color in legend
#ax.plot(0,0, color="blue", label="Healthy")
#ax.plot(0,0, color="red", label="HIV")

ax.plot(range(setrange+1)[1:], hvbmeanmeans, color="green", label="Healthy ranked means (100)", linewidth=2)
ax.plot(range(setrange+1)[1:], hivbmeanmeans, color="purple", label="HIV ranked means (100)", linewidth=2)

# inset plot first x (inset_size) ranks
#ax2 = plt.axes([.5,.5,.37,.37])
#ax2.plot(range(setrange+1)[1:inset_size+1], hvbmeans[:inset_size], color="blue", alpha=.5, linewidth=2)
#ax2.plot(range(setrange+1)[1:inset_size+1], hivbmeans[:inset_size], color="red", alpha=.5, linewidth=2)

#ax2.plot(range(setrange+1)[1:inset_size+1], hvbmeanmeans[:inset_size], color="green", linewidth=2)
#ax2.plot(range(setrange+1)[1:inset_size+1], hivbmeanmeans[:inset_size], color="purple", linewidth=2)

ax.axis(xmin=0, xmax=500, ymin=0, ymax=15)
#ax2.axis(xmin=1, xmax=inset_size, ymin=0, ymax=15)

#leg = ax.legend(loc=[.05,.8])
#for legobj in leg.legendHandles:
  #legobj.set_linewidth(3.0)
  
plt.legend(loc='upper right', prop={'size':13})

plt.tight_layout()  
plt.savefig(savepath + "Beta translational convergence.svg", bbox_inches='tight')
#plt.show()



##################################
####### OUTPUT STATISTICS ########
##################################

outstats = open(savepath + "translational convergence stats.txt", "w")

print >> outstats, "Run on " + time.strftime("%Y %m %d") 
print >> outstats, "Healthy alpha files used: " + str(hva)
print >> outstats, "Healthy beta files used: " + str(hvb)
print >> outstats, "HIV alpha files used: " + str(shiva)
print >> outstats, "HIV beta files used: " + str(shivb)

print >> outstats, "\n----------------Alpha----------------"
print >> outstats, "Healthy"
print >> outstats, "Whole:\tMean TC:\t" + str(np.mean(hvameans)) + "\tMedian TC:\t" + str(np.median(hvameans))
print >> outstats, "top500:\tMean TC:\t" + str(np.mean(hvameans[:500])) + "\tMedian TC:\t" + str(np.median(hvameans[:500]))
print >> outstats, "HIV"
print >> outstats, "Whole:\tMean TC:\t" + str(np.mean(hivameans)) + "\tMedian TC:\t" + str(np.median(hivameans))
print >> outstats, "top500:\tMean TC:\t" + str(np.mean(hivameans[:500])) + "\tMedian TC:\t" + str(np.median(hivameans[:500]))


print >> outstats, "\n----------------Beta----------------"
print >> outstats, "Healthy"
print >> outstats, "Whole:\tMean TC:\t" + str(np.mean(hvbmeans)) + "\tMedian TC:\t" + str(np.median(hvbmeans))
print >> outstats, "top500:\tMean TC:\t" + str(np.mean(hvbmeans[:500])) + "\tMedian TC:\t" + str(np.median(hvbmeans[:500]))
print >> outstats, "HIV"
print >> outstats, "Whole:\tMean TC:\t" + str(np.mean(hivbmeans)) + "\tMedian TC:\t" + str(np.median(hivbmeans))
print >> outstats, "top500:\tMean TC:\t" + str(np.mean(hivbmeans[:500])) + "\tMedian TC:\t" + str(np.median(hivbmeans[:500]))


print >> outstats, "\n----------------Comparison----------------"
print >> outstats, "\nMann-Whitney U test p-values (one sided)"
print >> outstats, "Healthy vs HIV:"
print >> outstats, "Whole alpha:\t" + str(ss.mannwhitneyu(hvameans, hivameans)[1]) + "\ttop500:\t" + str(ss.mannwhitneyu(hvameans[:500], hivameans[:500])[1]) 
print >> outstats, "Whole beta:\t" + str(ss.mannwhitneyu(hvbmeans, hivbmeans)[1]) + "\ttop500:\t" + str(ss.mannwhitneyu(hvbmeans[:500], hivbmeans[:500])[1])
print >> outstats, "Healthy alpha vs beta:"
print >> outstats, "Whole:\t" + str(ss.mannwhitneyu(hvameans, hvbmeans)[1]) + "\ttop500:\t" + str(ss.mannwhitneyu(hvameans[:500], hvbmeans[:500])[1]) 
print >> outstats, "HIV alpha vs beta:"
print >> outstats, "Whole:\t" + str(ss.mannwhitneyu(hvameans, hvbmeans)[1]) + "\ttop500:\t" + str(ss.mannwhitneyu(hvameans[:500], hvbmeans[:500])[1]) 

print >> outstats, "\nT test p-values (indep for inter donor, related for intra)"
print >> outstats, "Healthy vs HIV:"
print >> outstats, "Whole alpha:\t" + str(ss.ttest_ind(hvameans, hivameans)[1]) + "\ttop500:\t" + str(ss.ttest_ind(hvameans[:500], hivameans[:500])[1]) 
print >> outstats, "Whole beta:\t" + str(ss.ttest_ind(hvbmeans, hivbmeans)[1]) + "\ttop500:\t" + str(ss.ttest_ind(hvbmeans[:500], hivbmeans[:500])[1])
print >> outstats, "Healthy alpha vs beta:"
print >> outstats, "Whole:\t" + str(ss.ttest_rel(hvameans, hvbmeans)[1]) + "\ttop500:\t" + str(ss.ttest_rel(hvameans[:500], hvbmeans[:500])[1]) 
print >> outstats, "HIV alpha vs beta:"
print >> outstats, "Whole:\t" + str(ss.ttest_rel(hvameans, hvbmeans)[1]) + "\ttop500:\t" + str(ss.ttest_rel(hvameans[:500], hvbmeans[:500])[1]) 


outstats.close()



