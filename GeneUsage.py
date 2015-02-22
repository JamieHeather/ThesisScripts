# GeneUsage.py v1

  # Built on FindUsedGenes.py v1 (and using results from same document)
  # Designed with a similar strategy to that used in CDR3_JSD.py used in paper (and probably later in thesis)

# 24rd Oct 2014

# Describing healthy repertoires, want to demonstrate V/J gene usage/variance 
  # Only looking at healthy repertoires, v2 for those that have two, ignoring CD4/8

from __future__ import division
import time
import os
import collections as coll 
import matplotlib.pyplot as plt
import numpy as np
from random import uniform
import scipy.stats as ss

plt.rcParams.update({'font.size': 12})

savepath = "/home/jme/TCR/WRITE_UP/THESIS/WorkingPlots/" + time.strftime("%Y %m %d") + " "

all_a_files = [f for f in os.listdir(os.getcwd()) if f.endswith("a.freq") and "HV" in f and "v1" not in f and "CD" not in f]
all_b_files = [f for f in os.listdir(os.getcwd()) if f.endswith("b.freq") and "HV" in f and "v1" not in f and "CD" not in f]

def getVJ(filename, noV, noJ):
  # given a filename (and number of Vs and Js in appr), returns two dictionaries
    # = raw frequency usage of all the possible genes covered in those tag sets, for V and J
  # Initialise dictionaries and ensure have keys for all possible genes
  Vdict = coll.Counter()
  Jdict = coll.Counter()
  for i in range(noV):
    Vdict[i] = 0
  for i in range(noJ):
    Jdict[i] = 0
  # Read through file and add the frequency of each read to the appropriate V and J gene entries in dictionary
  for line in open(filename, "rU"):
    l = line.rstrip().split(", ")        
    Vdict[int(l[0])] += int(l[5])
    Jdict[int(l[1])] += int(l[5])   
  return(Vdict, Jdict)

# Genes used in our data, and their attendent names
trav = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 51, 52, 54, 55]
traj = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 57, 59, 60]
travnam = ['TRAV1-1', 'TRAV1-2', 'TRAV10', 'TRAV12-1', 'TRAV12-2', 'TRAV12-3', 'TRAV13-1', 'TRAV13-2', 'TRAV14/DV4', 'TRAV16', 'TRAV17', 'TRAV18', 'TRAV19', 'TRAV2', 'TRAV20', 'TRAV21', 'TRAV22', 'TRAV23/DV6', 'TRAV24', 'TRAV25', 'TRAV26-1', 'TRAV26-2', 'TRAV27', 'TRAV29/DV5', 'TRAV3', 'TRAV30', 'TRAV34', 'TRAV35', 'TRAV36/DV7', 'TRAV38-1', 'TRAV38-2/DV8', 'TRAV39', 'TRAV4', 'TRAV40', 'TRAV41', 'TRAV5', 'TRAV6', 'TRAV7', 'TRAV8-1', 'TRAV8-2', 'TRAV8-3', 'TRAV8-6', 'TRAV9-1', 'TRAV9-2', 'TRDV1', 'TRDV2', 'TRDV3', 'TRAV8-5', 'TRAV28', 'TRAV31', 'TRAV33', 'TRAV37']
trajnam = ['TRAJ10', 'TRAJ11', 'TRAJ12', 'TRAJ13', 'TRAJ14', 'TRAJ15', 'TRAJ16', 'TRAJ17', 'TRAJ18', 'TRAJ20', 'TRAJ21', 'TRAJ22', 'TRAJ23', 'TRAJ24', 'TRAJ26', 'TRAJ27', 'TRAJ28', 'TRAJ29', 'TRAJ3', 'TRAJ30', 'TRAJ31', 'TRAJ32', 'TRAJ33', 'TRAJ34', 'TRAJ36', 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRAJ4', 'TRAJ40', 'TRAJ41', 'TRAJ42', 'TRAJ43', 'TRAJ44', 'TRAJ45', 'TRAJ46', 'TRAJ47', 'TRAJ48', 'TRAJ49', 'TRAJ5', 'TRAJ50', 'TRAJ52', 'TRAJ53', 'TRAJ54', 'TRAJ56', 'TRAJ57', 'TRAJ6', 'TRAJ7', 'TRAJ8', 'TRAJ9', 'TRAJ1', 'TRAJ2', 'TRAJ25', 'TRAJ35', 'TRAJ58', 'TRAJ60', 'TRAJ61']

trbv = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 57, 58]
trbj = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
trbvnam = ['TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2', 'TRBV20-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1', 'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5', 'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-8', 'TRBV6-9', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8', 'TRBV7-9', 'TRBV9', 'TRBV17', 'TRBV23-1', 'TRBV5-3', 'TRBV5-7', 'TRBV6-7', 'TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV21-1', 'TRBV22-1', 'TRBV5-2', 'TRBV7-5']
trbjnam = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7']


travcol = ['g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'lime', 'lime', 'lime', 'r', 'r', 'r', 'r', 'r']
trajcol = ['g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'yellow', 'yellow', 'yellow', 'yellow', 'yellow', 'r', 'yellow']
trbvcol = ['g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'yellow', 'yellow', 'yellow', 'yellow', 'yellow', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
trbjcol = ['g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g']

numva = len(trav)
numja = len(traj)
numvb = len(trbv)
numjb = len(trbj)


# loop through various genes

chains = ['a', 'b']

# variables used to store standard deviations of various groups, for last plot
fhigh = []
flow = []
ophigh = []
oplow = []
allmean = []
allstd = []

for c in chains:
  
  # read V and J data from files

  vars()['all' + c + 'v'] = []
  vars()['all' + c + 'j'] = []
  vars()['all' + c] = []

  for f in vars()['all_' + c + '_files']:
    nam = f.split("_")[2].split(".")[0]
    vnam = nam + "v"
    jnam = nam + "j"
    vars()[vnam], vars()[jnam] = getVJ(f, vars()["numv" + c], vars()["numj" + c])
    vars()['all' + c + 'v'].append(vnam)
    vars()['all' + c + 'j'].append(jnam)
    vars()['all' + c].append(nam)

  # loop through all genes and record the proportions with which each donor uses that gene
  vusage = coll.defaultdict(list)
  jusage = coll.defaultdict(list)  
  
  for g in vars()['tr' + c + 'v']:	# loop through genes
    for s in vars()['all' + c + 'v']:	# loop through samples
      dat = vars()[s]
      vusage[g].append(dat[g]/sum(dat.values()))
      
  for g in vars()['tr' + c + 'j']:	# loop through genes
    for s in vars()['all' + c + 'j']:	# loop through samples
      dat = vars()[s]
      jusage[g].append(dat[g]/sum(dat.values()))

  # loop through all genes again, and this time calculate mean proportional usage, and standard deviation
  vmeans = []
  jmeans = []

  vstd = []
  jstd = []

  for g in vars()['tr' + c + 'v']:	# loop through genes
    vmeans.append(np.mean(vusage[g]))
    vstd.append(np.std(vusage[g]))

  for g in vars()['tr' + c + 'j']:	# loop through genes
    jmeans.append(np.mean(jusage[g]))
    jstd.append(np.std(jusage[g]))
    
  # need to account for the fact that some rarely used genes have high std dev, causing error bars to go off axis
    # just need to substract the lower limit of the plots, i.e. 1e-6
  vstdchk = [x-y for x,y in zip(vmeans,vstd)]
  jstdchk = [x-y for x,y in zip(jmeans,jstd)]
  
  vlow = []
  jlow = []
  
  for i in range(len(vstdchk)):
    if vstdchk[i] == 0:
      vlow.append(0)
    elif vstdchk[i] > 0:
      vlow.append(vstd[i])
    else:
      vlow.append(vmeans[i]-1e-6)
      
  for i in range(len(jstdchk)):
    if jstdchk[i] == 0:
      jlow.append(0)
    elif jstdchk[i] > 0:
      jlow.append(jstd[i])
    else:
      jlow.append(jmeans[i]-1e-6)
      
  
  # also want to demonstrate that the lower used genes (i.e. below median) have greater variance 
    # furthermore want to distinguish between F and ORF or P genes (f vs op), showing whether they follow this rule
  
  for i in range(len(vmeans)):
    if vmeans[i] > np.median(vmeans):
      if vars()['tr'+c+'vcol'][i] == 'g' or vars()['tr'+c+'vcol'][i] == 'lime': 
	# infer functionality from the color tuples
	fhigh.append(vstd[i])
      else:
	ophigh.append(vstd[i])
    else:
      if vars()['tr'+c+'vcol'][i] == 'g' or vars()['tr'+c+'vcol'][i] == 'lime': 
	# infer functionality from the color tuples
	flow.append(vstd[i])
      else:
	oplow.append(vstd[i])      
  
  for i in range(len(jmeans)):
    if jmeans[i] > np.median(jmeans):
      if vars()['tr'+c+'jcol'][i] == 'g': 
	# infer functionality from the color tuples
	fhigh.append(jstd[i])
      else:
	ophigh.append(jstd[i])
    else:
      if vars()['tr'+c+'jcol'][i] == 'g': 
	# infer functionality from the color tuples
	flow.append(jstd[i])
      else:
	oplow.append(jstd[i])      


  # ... having done this and noted the opposite of what's expected (i.e. the more expressed the higher the variance)
    # now want to plot this, i.e. mean usage vs standard deviation
  
  allmean = allmean + vmeans + jmeans
  allstd = allstd + vstd + jstd       

  # plot 

  fig = plt.figure(figsize=(10,3.5))
  ax = fig.add_subplot(111)
  ax.bar(range(vars()['numv'+c]), vmeans, log=True, color=vars()['tr' + c + 'vcol'], bottom=1e-6)
  ax.errorbar([x+.45 for x in range(vars()['numv'+c])], vmeans, yerr=[vlow,vstd], xerr=[0]*len(vlow), fmt="none", log=True, color='black')
  plt.ylabel("Proportional gene usage")
  ax.set_xticks([x+.45 for x in range(vars()['numv'+c])])
  ax.set_xticklabels(vars()['tr'+c+'vnam'], rotation=90)
  plt.xlim(0,vars()['numv'+c])
  #if c=="a":	# legend looked awful on logged plot - will just crop it out and paste it over in inkscape
    #ax.bar(100,1, color="g", label="F")
    #ax.bar(100,1, color="lime", label="F (TRDV)")
    #ax.bar(100,1, color="yellow", label="ORF")
    #ax.bar(100,1, color="r", label="P")
    #plt.legend(loc='lower center', prop={'size':14})	
  #plt.show()
  plt.savefig(savepath + "TR" + c.upper() + "V HV Prop Usage.svg", bbox_inches='tight')
  plt.close()
  
  if c == "a":
    fig = plt.figure(figsize=(10,3.5))
  elif c == "b":
    # can make the beta J plot narrower, as it has fewer columns, leaving space for the variance analysis
    fig = plt.figure(figsize=(4,3.5))
  ax = fig.add_subplot(111)
  ax.bar(range(vars()['numj'+c]), jmeans, log=True, color=vars()['tr' + c + 'jcol'], bottom=1e-6)
  ax.errorbar([x+.45 for x in range(vars()['numj'+c])], jmeans, yerr=[jlow,jstd], xerr=[0]*len(jlow), fmt="none", log=True, color='black')
  plt.ylabel("Proportional gene usage")
  ax.set_xticks([x+.45 for x in range(vars()['numj'+c])])
  ax.set_xticklabels(vars()['tr'+c+'jnam'], rotation=90)
  plt.xlim(0,vars()['numj'+c])
  #plt.show()
    
  plt.savefig(savepath + "TR" + c.upper() + "J HV Prop Usage.svg", bbox_inches='tight')
    
  plt.close()

  #sys.exit()
  
  
# plot standard deviation of various genes
  # NB - ALL P/ORF genes average proportional usage less than the median, which is probably to be expected
  
  
jitter = 0.1 
  # add random noise to allow production of nice R-style stripcharts via scatterplots

fhighx = [1+uniform(-jitter, jitter) for x in range(len(fhigh))]
flowx = [2+uniform(-jitter, jitter) for x in range(len(flow))]
oplowx = [3+uniform(-jitter, jitter) for x in range(len(oplow))]


fig = plt.figure(figsize=(4,3.5))
ax = fig.add_subplot(111)
ax.scatter(fhighx, fhigh, color='g', alpha=.5, label="F")
ax.scatter(flowx, flow, color='g', alpha=.5)
ax.scatter(oplowx, oplow, color='orange', alpha=.5, label="ORF/P")
ax.set_yscale('log')
ax.set_xticks(range(1,4,1))
ax.set_xticklabels(['F high', 'F low', 'ORF/P low'])
# add median lines
ax.plot([0.75, 1.25], [np.median(fhigh)]*2, color="g", linewidth=3, alpha=.5)
ax.plot([1.75, 2.25], [np.median(flow)]*2, color="g", linewidth=3, alpha=.5)
ax.plot([2.75, 3.25], [np.median(oplow)]*2, color="orange", linewidth=3, alpha=.5)
plt.ylim(9e-6,1e-1)
plt.ylabel('Standard deviation')
plt.legend(loc="upper right", prop={'size':14})

plt.savefig(savepath + "std distributions.svg", bbox_inches='tight')
#plt.show()
plt.close()

# statistics
ss.mannwhitneyu(fhigh, flow)
# (915.0, 1.0566879686569566e-13)
# standard deviation of F high signif greater than F low

ss.mannwhitneyu(flow, oplow)
# (208.0, 4.0799875936820635e-08)
# standard deviation of P/ORF signif lower than those of F low

# So basically this means the more expressed a gene is, the more variable the usage of that gene 
# (confirmed with a quick plt.scatter(vmeans,vstd))

# plotting that properly

allcol = travcol + trajcol + trbvcol + trbjcol

coefficients = np.polyfit(allmean, allstd, 1)
polynomial = np.poly1d(coefficients)
ys = polynomial(allmean)

fig = plt.figure(figsize=(4,3.5))
ax = fig.add_subplot(111)
ax.scatter(allmean, allstd, color=allcol, alpha = .7)

#ax.set_xscale('log')
##ax.set_yscale('log')
#plt.xlim(1e-6,1e0)
#plt.ylim(1e-6,1e0)

plt.xlim(-.01,.16)
plt.ylim(-.002,.04)

plt.plot(allmean, ys)
plt.xlabel('Mean proportional usage')
plt.ylabel('Standard deviation')
plt.savefig(savepath + "Mean Prop Usage vs std.svg", bbox_inches='tight')
#plt.show()
plt.close()


ss.linregress(allmean,allstd)
#slope, intercept, r_value, p_value, std_err
#(0.19947942005866842, -3.0833303562486855e-05, 0.90303236375583307, 7.3883657851553767e-67, 0.0071325469665176177)














