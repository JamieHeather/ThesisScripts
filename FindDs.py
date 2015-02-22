# FindDs_AllFiles.py
  # = FindDs.py v2, or rather a version of that applied to all files in a directory.

# This file just designed to track TRBD (as have no delta data)

from __future__ import division
import collections as coll
import difflib as dl
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import time 
import scipy.stats as ss

fontsize = 12
plt.rcParams.update({'font.size': fontsize})

savepath = "/home/jme/TCR/WRITE_UP/THESIS/WorkingPlots/" + time.strftime("%Y %m %d") + " "

all_b_files = [f for f in os.listdir(os.getcwd()) if f.endswith("b.freq") and "HV" in f and "v1" not in f and "CD" not in f]

if len(sys.argv) <> 2:
  print "Please supply threshold length of D region allowed for assignation"
  print "e.g. python CollapseTCRs.py FILENAME.freq b 8 (default is 7)"
  sys.exit()
else:
  len_threshold = int(sys.argv[1])


def sig(val):
  # Given a p-value, returns an appropriate asterisk significance indicator
  #if val <= 0.0001:	# optional additional value, not sure that depth is needed for this analysis
    #return "****"
  if val <= 0.001:
    return "***"
  elif val <= 0.01:
    return "**"
  elif val <= 0.05:
    return "*"
  else:
    return ""



chain = 'b'

bdnames = ['TRBD1','TRBD2','TRBD2*2'] 

bd = ['GGGACAGGGGGC', 'GGGACTAGCGGGGGG', 'GGGACTAGCGGGAGGG']

ddnames = ['TRDD1', 'TRDD2', 'TRDD3']
dd = ['GAAATAGT', 'CCTTCCTAC', 'ACTGGGGGATACG']


bd1vdicts = []
bd2vdicts = []
bd1jdicts = []
bd2jdicts = []

uniddicts = []
d1dicts = []
d2dicts = []

d1c = 0
d2c = 0

############ 

# Used genes, i.e. those that occur in our data and thus are the only ones we need both plotting

trbv = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 57, 58]
trbvnam = ['TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2', 'TRBV20-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1', 'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5', 'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-8', 'TRBV6-9', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8', 'TRBV7-9', 'TRBV9', 'TRBV17', 'TRBV23-1', 'TRBV5-3', 'TRBV5-7', 'TRBV6-7', 'TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV21-1', 'TRBV22-1', 'TRBV5-2', 'TRBV7-5']
trbjnam = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7']

# TRBJ is easier as it's just 0-12

###########

# Need to loop through all files

d1 = 0
d2 = 0
  


for f in all_b_files:
  nam = f.split("_")[2].split(".")[0]

  uniq_dcr_count = 0
  tot_dcr_count = 0 
  uniq_with_d = 0
  uniq_with_dd = 0
  tot_with_d = 0
  tot_with_dd = 0

  # dictionaries to store the frequency of all the associated V genes 
    # also initialise with empty values for all possible genes (i.e. number of extended tags)
  vars()[nam+"_1vs"] = coll.Counter()
  vars()[nam+"_1js"] = coll.Counter()
  vars()[nam+"_2vs"] = coll.Counter()
  vars()[nam+"_2js"] = coll.Counter()
  
  bd1vdicts.append(nam+"_1vs")
  bd2vdicts.append(nam+"_2vs")
  bd1jdicts.append(nam+"_1js")
  bd2jdicts.append(nam+"_2js")

  vars()[nam+"_unid"] = coll.Counter()
  vars()[nam+"_d1"] = coll.Counter()
  vars()[nam+"_d2"] = coll.Counter()
  
  uniddicts.append(nam+"_unid")
  d1dicts.append(nam+"_d1")
  d2dicts.append(nam+"_d2")

  for g in trbv:
    vars()[nam+"_1vs"][g] = 0
    vars()[nam+"_2vs"][g] = 0

  for i in range(13):
    vars()[nam+"_1js"][str(i)] = 0
    vars()[nam+"_2js"][str(i)] = 0

  # keep each unit in temp dicts first though, as will be recording proportional gene usage
  tmp1vs = coll.Counter()
  tmp2vs = coll.Counter()
  tmp1js = coll.Counter()
  tmp2js = coll.Counter()
  
  # dictionaries to store length distribution (only for NON TANDEM genes)
  tmpunidlens = coll.Counter()	# those we can't ID
  tmpd1lens = coll.Counter()
  tmpd2lens = coll.Counter()
    
  unid = 0 	# counter for number of unidentifiable D genes
  canid = 0	# and those we can ID
  
  for dcr in open(f, "rU"):
    # loop through all records in Decombined file
        
    bits = dcr.rstrip().split(", ")
    insert = bits[4]
    
    uniq_dcr_count += 1
    tot_dcr_count += int(bits[5])
    
    possds = coll.Counter()
    
    found = False 	# check that only count each matching dcr once
    double_chk = False # need to check haven't already counted a DCR (otherwise tandems get counted twice)

        
    for g in range(len(vars()[chain + "d"])):
      # look for all known D genes of that chain
      s = dl.SequenceMatcher(None, insert, vars()[chain + "d"][g])
      match = s.find_longest_match(0, len(insert), 0, len(vars()[chain + "d"][g]))
      
      if match[2] >= len_threshold:
	possds[vars()[chain + "dnames"][g]] = match[2]
	
	#print insert, vars()[chain + "d"][g], vars()[chain + "dnames"][g], insert[match[0]:match[0]+match[2]], bits[0]
      
      if possds and found==False: 
	uniq_with_d += 1
	tot_with_d += int(bits[5])
	found = True
      # check for tandem Ds, i.e. VDDJ rearrangements
      
      if chain == "b" and double_chk == False:
	if ('TRBD1' in possds.keys() and 'TRBD2' in possds.keys()) or ('TRBD1' in possds.keys() and 'TRBD2*2' in possds.keys()): 
	  uniq_with_dd += 1
	  tot_with_dd += int(bits[5])
	  double_chk = True
	  #print dcr.rstrip()
	  #sys.exit()
      elif chain == "d" and double_chk == False:
	if len(possds.keys()) > 1:
	  uniq_with_dd += 1
	  tot_with_dd += int(bits[5])
	  double_chk = True
      #if len(possds.keys()) >1: sys.exit()
	##sys.exit()
	
    # want to record the length of the different Ds (longest D match available)
    
    # count lengths of Ds    
    if len(possds.keys()) == 1 and double_chk == False:
      if possds.most_common()[0][0] == 'TRBD1':
	tmpd1lens[possds.most_common()[0][1]] += int(bits[5])
	d1c += int(bits[5])
      elif possds.most_common()[0][0] in ['TRBD2', 'TRBD2*2']:
	tmpd2lens[possds.most_common()[0][1]] += int(bits[5])
	d2c += int(bits[5])
    elif len(possds.keys()) > 1:
      if sorted([possds.most_common()[0][0], possds.most_common()[1][0]]) == ['TRBD2', 'TRBD2*2']:
	# many times TRDB2 will get two equal hits
	tmpd2lens[possds.most_common()[0][1]] += int(bits[5])     
	d2c += int(bits[5])
      elif possds.most_common()[0][1] == possds.most_common()[1][1]:      
	tmpunidlens[possds.most_common()[0][1]] += int(bits[5])
      else:
	if possds.most_common()[0][0] == 'TRBD1':
	  tmpd1lens[possds.most_common()[0][1]] += int(bits[5])
	  d1c += int(bits[5])
	elif possds.most_common()[0][0] in ['TRBD2', 'TRBD2*2']:
	  tmpd2lens[possds.most_common()[0][1]] += int(bits[5])     
	  d2c += int(bits[5])
	
    # want to record which Js different Ds associate with (only for singles)
    
    if found == True and double_chk == False:
	
      if len(possds.keys()) == 1:
	canid += 1
	tempgene = possds.most_common()[0][0]
      
      else:
	firstmatchlen = possds.most_common()[0][1]
	secondmatchlen = possds.most_common()[1][1]
      
	if firstmatchlen == secondmatchlen:
	  unid += 1
	  continue

	else:
	  canid += 1
	  tempgene = possds.most_common()[0][0]
	
      if tempgene == 'TRBD1':
	tmp1vs[bits[0]] += int(bits[5])
	tmp1js[bits[1]] += int(bits[5])
	d1 += int(bits[5])
      elif tempgene == 'TRBD2' or tempgene == 'TRBD2*2':
	tmp2vs[bits[0]] += int(bits[5])
	tmp2js[bits[1]] += int(bits[5])
	d2 += int(bits[5])

  tot_dcr = tot_dcr_count
  # load proportional values of gene usage into donor specific dictionaries
  for g in trbv:
    vars()[nam+"_1vs"][g] += tmp1vs[str(g)]/tot_dcr
    vars()[nam+"_2vs"][g] += tmp2vs[str(g)]/tot_dcr

  for i in range(13):
    vars()[nam+"_1js"][i] += tmp1js[str(i)]/tot_dcr
    vars()[nam+"_2js"][i] += tmp2js[str(i)]/tot_dcr
  
  longest = max(max(tmpd1lens.keys()),max(tmpd2lens.keys()))
  
  for i in range(len_threshold, longest):
    vars()[nam+"_unid"][i] += tmpunidlens[i]/tot_dcr
    vars()[nam+"_d1"][i] += tmpd1lens[i]/tot_dcr
    vars()[nam+"_d2"][i] += tmpd2lens[i]/tot_dcr
    
########## plotting D len distributions

unidmeans = []
d1means = []
d2means = []

unidstds = []
d1stds = []
d2stds = []

up = max([len(x) for x in bd])

for i in range(len_threshold, up):
  tmpun = []
  tmp1 = []
  tmp2 = []
  for d in uniddicts:
    dat = vars()[d]
    tmpun.append(dat[i])
  for d in d1dicts:
    dat = vars()[d]
    tmp1.append(dat[i])
  for d in d2dicts:
    dat = vars()[d]
    tmp2.append(dat[i])
  unidmeans.append(np.mean(tmpun))
  d1means.append(np.mean(tmp1))
  d2means.append(np.mean(tmp2))
  unidstds.append(np.std(tmpun))
  d1stds.append(np.std(tmp1))
  d2stds.append(np.std(tmp2))
  

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

xs = np.arange(up - len_threshold)

#ax.bar(xs, d2means, bottom=0, color="lightgray")	# old version - bar plot
#ax.bar(xs, d1means, bottom=d2means, color="black")
#ax.bar(xs, unidmeans, bottom=[x+y for x,y in zip(d2means, d1means)], color="white")
#ax.errorbar([x+.45 for x in xs], d2means, yerr=d2stds, fmt="none", ecolor='grey', lw=1.5, capthick=1.5)
#ax.errorbar([x+.45 for x in xs], [x+y for x,y in zip(d1means, d2means)], \
  #yerr=d1stds, fmt="none", ecolor='grey', lw=1.5, capthick=1.5)
#ax.errorbar([x+.45 for x in xs], [x+y+z for x,y,z in zip(d1means, d2means, unidmeans)], \
  #yerr=unidstds, fmt="none", ecolor='grey', lw=1.5, capthick=1.5)

ax.plot([x+.45 for x in xs], d1means, linewidth=2, alpha=.9, color="black", label="TRBD1")
ax.plot([x+.45 for x in xs], d2means, linewidth=2, alpha=.9, color="darkgray", label="TRBD2")		# new version - line plot, remove confusion
ax.plot([x+.45 for x in xs], unidmeans, linewidth=2, alpha=.9, linestyle="dashed", color="black", label="Unidentifiable")

ax.errorbar([x+.45 for x in xs], d1means, yerr=d1stds, fmt="none", ecolor='black', lw=1.5, capthick=1.5)
ax.errorbar([x+.45 for x in xs], d2means, yerr=d2stds, fmt="none", ecolor='darkgrey', lw=1.5, capthick=1.5)
ax.errorbar([x+.45 for x in xs], unidmeans, yerr=unidstds, fmt="none", ecolor='black', lw=1.5, capthick=1.5)

## fake plots for legend in nice order
#ax.bar(100,0, color="white", label="Unidentifiable")
#ax.bar(100,00, color="black", label="TRBD1")
#ax.bar(100,00, color="lightgray", label="TRBD2")
plt.xlim(0,up-len_threshold)
ax.set_xticks([x+.45 for x in range(up-len_threshold)])
ax.set_xticklabels(range(len_threshold, up))
plt.ylabel('Proportion of rearrangements')
plt.xlabel('Length of longest TRBD match')
plt.legend(loc="upper right", prop={'size':fontsize})
plt.ylim(-.001,.085)

#plt.show()

plt.savefig(savepath + str(len_threshold) + "_" + "TRBD_Lengths.svg", bbox_inches='tight')
plt.close()



  
#sys.exit()

######### plotting gene usage

# first need to generate the average values for mean gene usage, and std dev

v1means = []
v2means = []
v1std = []
v2std = []
vsig = []

for g in trbv:
  tmp1 = []
  tmp2 = []
  for d in bd1vdicts:
    dat = vars()[d]
    tmp1.append(dat[g])
  for d in bd2vdicts:
    dat = vars()[d]
    tmp2.append(dat[g])
  v1means.append(np.mean(tmp1))
  v1std.append(np.std(tmp1))
  v2means.append(np.mean(tmp2))
  v2std.append(np.std(tmp2))  
  if tmp1 <> tmp2:
    vsig.append(ss.mannwhitneyu(tmp1,tmp2)[1])
  else:
    vsig.append(1)
    
j1means = []
j2means = []
j1std = []
j2std = []
jsig = []

for i in range(13):
  tmp1 = []
  tmp2 = []
  for d in bd1jdicts:
    dat = vars()[d]
    tmp1.append(dat[i])
  for d in bd2jdicts:
    dat = vars()[d]
    tmp2.append(dat[i])
  j1means.append(np.mean(tmp1))
  j1std.append(np.std(tmp1))
  j2means.append(np.mean(tmp2))
  j2std.append(np.std(tmp2))  
  if tmp1 <> tmp2:  
    jsig.append(ss.mannwhitneyu(tmp1,tmp2)[1])
  else:
    jsig.append(1)
    
# plot J
#sys.exit()

xs = np.arange(13)
width=.4
trbj = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7']
sigmx = [sig(jsig[x]) for x in range(13)]

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

ax.set_xticks([x+.45 for x in range(13)])
ax.set_xticklabels(trbj, rotation=90)

plt.bar(xs, j1means, width, color="black", label="TRBD1", bottom=0)
plt.bar(xs+width, j2means, width, color="lightgray", label="TRBD2", bottom=0)

ax.errorbar([x+.22 for x in xs], j1means, yerr=[j1std, [0 for x in j1std]], fmt="none", ecolor='white', lw=2, capthick=2)
ax.errorbar([x+.22 for x in xs], j1means, yerr=[[0 for x in j1std],j1std], fmt="none", ecolor='black', lw=2, capthick=2)
ax.errorbar([x+width+.22 for x in xs], j2means, yerr=j2std, fmt="none", ecolor='black', lw=2, capthick=2)

for s in range(len(sigmx)):
  ax.text(xs[s]+width+.19, .057, sigmx[s], ha="center", rotation="vertical", size=12)
  
plt.xlim(0,13)
plt.ylim(0,.06)

#plt.legend(loc="upper left", prop={'size':fontsize})
plt.ylabel("Proportion of rearrangements")

plt.savefig(savepath + str(len_threshold) + "_" + "TRBDJ_Pairing.svg", bbox_inches='tight')
#plt.show()
plt.close()

#sys.exit()

###############

# TRBV pairing

xs = np.arange(len(trbv))
width=.4
sigmx = [sig(vsig[x]) for x in range(len(trbv))]

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)

ax.set_xticks([x+.45 for x in range(len(trbv))])
ax.set_xticklabels(trbvnam, rotation=90)

plt.bar(xs, v1means, width, color="black", label="TRBD1", bottom=0)
plt.bar(xs+width, v2means, width, color="lightgray", label="TRBD2", bottom=0)

for s in range(len(sigmx)):
  ax.text(xs[s]+width+.25, .035, sigmx[s], ha="center", rotation="vertical", size=12)

ax.errorbar([x+.21 for x in xs], v1means, yerr=[v1std, [0 for x in v1std]], fmt="none", ecolor='white', capsize=2)
ax.errorbar([x+.21 for x in xs], v1means, yerr=[[0 for x in v1std],v1std], fmt="none", ecolor='black', capsize=2)
ax.errorbar([x+width+.21 for x in xs], v2means, yerr=v2std, fmt="none", ecolor='black', capsize=2)

plt.ylabel("Proportion of rearrangements")
plt.xlim(0,len(trbv))
plt.ylim(0,.037)
#plt.legend(loc="upper left", prop={'size':fontsize})

plt.savefig(savepath + str(len_threshold) + "_TRBDV_Pairing.svg", bbox_inches='tight')
#plt.show()
plt.close()











