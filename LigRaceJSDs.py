# LigRaceJSDs.py v1
  # JH, UCL Oct 2014
  
# Want to by how much our Lig RACE protocol and other published sets differ in their V gene usage
  # For a figure for thesis - whether our protocol introduces bias
  # Calculate by comparing the Jensen Shannons Divergences of the V distributions of the various groups

from __future__ import division
import os  
import re
import collections as coll
import matplotlib.pyplot as plt
import numpy as np
import time
import itertools as it
import random

# key 

# _f = file
# _i = index (used to keep track of which sample is which in the big props tuples)
# _js = JSDs 
# _p = plot points
# _s = source (i.e. which file that instance belongs to)

# a_ = alpha 
# b_ = beta 

# xVy = JSD comparison of an x distribution versus a y distribution

###########################
#### DEFINE GENE LISTS ####
###########################

#savepath = "/media/jme/SAMSUNG/ThesisAnalysis/others_data/uniqs/" + time.strftime("%Y %m %d") + " "
savepath = "/media/jme/KINGSTON/uniqs/" + time.strftime("%Y %m %d") + " "

outstats = open(savepath + time.strftime("%Y %m %d") + " JSD outstats.txt", "w")

# Number of genes per gene type

numb_trav = 56
numb_trbv = 63



##########################
#### DEFINE FUNCTIONS ####
##########################

  
def basename(filename):
  # given a .freq filename, returns the base (i.e. the individual)
  if "HV" in filename:
    temp = re.split('[_.v]', filename)[3]
    return re.split('[_.ab]', temp)[0]
  else:
    return re.split('[_.v]', filename)[3] + re.split('[_.v]', filename)[4]   

def jsd(x,y): #Jensen-shannon divergence
  # from Jonathan Friedman
    # http://stats.stackexchange.com/questions/29578/jensen-shannon-divergence-calculation-for-3-prob-distributions-is-this-ok
    # author of pysurvey package, which uses JS
  # see 4th equation on http://www.mtm.ufsc.br/~taneja/book/node21.html
    import warnings
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    x = np.array(x)
    y = np.array(y)
    d1 = x*np.log2(2*x/(x+y))
    d2 = y*np.log2(2*y/(x+y))
    d1[np.isnan(d1)] = 0
    d2[np.isnan(d2)] = 0
    d = 0.5*np.sum(d1+d2)    
    return d  

def flatten(lst): # Flatten a list 
  return([thing for sublist in lst for thing in sublist])



##################
#### ANALYSIS ####
##################

print "Getting data"


# Get all .freq files (contain (comma-delimited) decombinator index : CDR3 , frequency)
# ONLY deal with "v2" healthy files, i.e.  only one bleed per donor

a_us_f = ['vDCRe_alpha_HV07a.freq', 'vDCRe_alpha_HV08a.freq', 'vDCRe_alpha_HV09a.freq']

a_wang_f = ['vDCRe_alpha_wang_1_SRR030416.uniq', 'vDCRe_alpha_wang_2_SRR030417.uniq', 'vDCRe_alpha_wang_3_SRR030709.uniq']


b_us_f = ['vDCRe_beta_HV07b.freq', 'vDCRe_beta_HV08b.freq', 'vDCRe_beta_HV09b.freq']

b_wang_f = ['vDCRe_beta_wang_1_SRR030416.uniq', 'vDCRe_beta_wang_2_SRR030417.uniq', 'vDCRe_beta_wang_3_SRR030709.uniq']

b_warren_f = ['vDCRe_beta_warren_1_SRR060699.uniq', 'vDCRe_beta_warren_2_SRR060700.uniq', 'vDCRe_beta_warren_3_SRR060701.uniq']

b_bolotin_f = ['vDCRe_beta_bol_1_SRR494403.uniq', 'vDCRe_beta_bol_2_SRR494404.uniq', 'vDCRe_beta_bol_3_SRR494405.uniq']


#####################
####### ALPHA ####### #############################################################################
#####################

#all_a_files = a_us_f + a_wang_f

#all_a = [basename(x) for x in all_a_files]

#print "TRAV"

#gene_numb = numb_trav

## get indexes of which files (in order) relate to different groupings

#a_us_i = []
#a_wang_i = []

#for i in range(len(all_a_files)):
  #if "HV" in all_a_files[i]: 
    #a_us_i.append(i)
  #else: 
    #a_wang_i.append(i)


## calculate the proportional of each V gene within a repertoire
#props = []

#for f in all_a_files:
  #infile = open(f, "rU")
  #outdict = "vb" + basename(f)
  #vars()[outdict] = coll.Counter()
  #for i in range(gene_numb):
    #vars()[outdict][i] = 0 	# populate counter with empty values for each V region  
  #for line in infile:
    ## read in frequencies of each V region usage
    #split1 = line.rstrip().split(", ")
    #v = int(split1[0])
    #freq = int(split1[5])
    #vars()[outdict][v] += freq
  #temp_props = []
  #for i in range(gene_numb):
    #temp_props.append(vars()[outdict][i]/sum(vars()[outdict].values()))
  #props.append(temp_props)





#usVus_js = []
#usVwang_js = []
#wangVwang_js = []

#same_chk = 0

#for i1 in range(len(props)):
  #temp_usVus_js = []
  #temp_usVwang_js = []
  #temp_wangVwang_js = []

  #for i2 in range(len(props)):
    #if i1 == i2:
      #same_chk += 1
    #elif i1 in a_us_i and i2 in a_us_i:
      #temp_usVus_js.append(jsd(props[i1], props[i2]))
    #elif ((i1 in a_us_i and i2 in a_wang_i) or (i1 in a_wang_i and i2 in a_us_i)):
      #temp_usVwang_js.append(jsd(props[i1], props[i2]))
    #elif i1 in a_wang_i and i2 in a_wang_i:
      #temp_wangVwang_js.append(jsd(props[i1], props[i2]))
    #else:
      ## just check I haven't missed any eventuality
      #print "Error, something\'s missing"
      
  #usVus_js.append(temp_usVus_js)
  #usVwang_js.append(temp_usVwang_js)
  #wangVwang_js.append(temp_wangVwang_js)

#if same_chk <> len(props): print "Error, same check failed"


#fig = plt.figure(figsize=(6,5))

#ax = fig.add_subplot(111)

#uVu_p = [[x]*2 for x in range(3)] + [[]] * 3
#ax.scatter(uVu_p[:3], usVus_js[:3], color="blue")

#uVw_p = [[x]*3 for x in range(6)] 
#ax.scatter(uVw_p, usVwang_js, color="cyan")

#wVw_p = [[]] * 3 + [[x]*2 for x in range(3,6,1)] 
#ax.scatter(wVw_p[3:], wangVwang_js[3:], color="green")

#plt.show()

###################
###### BETA ####### ###############################################################################
###################



all_b_files = b_us_f + b_wang_f + b_warren_f + b_bolotin_f

all_b = [basename(x) for x in all_b_files]

print "TRBV"

gene_numb = numb_trbv

# get indexes of which files (in order) relate to different groupings

b_us_i = []
b_wang_i = []
b_warren_i = []
b_bolotin_i = []

index = coll.defaultdict(list)

for i in range(len(all_b_files)):
  if "HV" in all_b_files[i]: 
    b_us_i.append(i)
    index[i].append("us")
  elif "wang" in all_b_files[i]: 
    b_wang_i.append(i)
    index[i].append("wang")
  elif "warren" in all_b_files[i]: 
    b_warren_i.append(i)
    index[i].append("warren")
  elif "bol" in all_b_files[i]: 
    b_bolotin_i.append(i)
    index[i].append("bol")

# calculate the proportional of each V gene within a repertoire
props = []

for f in all_b_files:
  infile = open(f, "rU")
  outdict = "vb" + basename(f)
  vars()[outdict] = coll.Counter()
  for i in range(gene_numb):
    vars()[outdict][i] = 0 	# populate counter with empty values for each V region  
  for line in infile:
    # read in frequencies of each V region usage
    split1 = line.rstrip().split(", ")
    v = int(split1[0])
    freq = int(split1[5])
    vars()[outdict][v] += freq
  temp_props = []
  for i in range(gene_numb):
    temp_props.append(vars()[outdict][i]/sum(vars()[outdict].values()))
  props.append(temp_props)

sources = ['us', 'wang', 'warren', 'bol']

#combos = it.combinations(sources, 2)
#for x in combos: print x[0] + "V" + x[1] + "_js"

combos = it.product(sources, repeat=2)
for x in combos: vars()[x[0] + "V" + x[1] + "_js"] = []


combos = it.product(sources, repeat=2)
jsd_lists = [x[0] + "V" + x[1] + "_js" for x in combos]


#usVus_js = []
#wangVwang_js = []
#warrenVwarren_js = []
#bolotinVbolotin_js = []

same_chk = 0

for i1 in range(len(props)):
  #temp_usVus_js = []
  #temp_wangVwang_js = []
  #temp_warrenVwarren_js = []
  #temp_bolotinVbolotin_js = []
  combos = it.product(sources, repeat=2)
  for x in combos: vars()["temp_" + x[0] + "V" + x[1] + "_js"] = []

  for i2 in range(len(props)):
    if i1 == i2:
      same_chk += 1    
      
    else: 
      i1_s = index[i1][0]
      i2_s = index[i2][0]
      
      vars()["temp_" + i1_s + "V" + i2_s + "_js"].append(jsd(props[i1], props[i2]))

  #for i2 in range(len(props)):
    #if i1 == i2:
      #same_chk += 1
    
    #elif i1 in a_us_i and i2 in a_us_i:
      #temp_usVus_js.append(jsd(props[i1], props[i2]))
    #elif ((i1 in a_us_i and i2 in a_wang_i) or (i1 in a_wang_i and i2 in a_us_i)):
      #temp_usVwang_js.append(jsd(props[i1], props[i2]))
    #elif i1 in a_wang_i and i2 in a_wang_i:
      #temp_wangVwang_js.append(jsd(props[i1], props[i2]))
    #else:
      ## just check I haven't missed any eventuality
      #print "Error, something\'s missing"
  
  for j in jsd_lists:
    vars()[j].append(vars()["temp_" + j])
  #usVus_js.append(temp_usVus_js)
  #usVwang_js.append(temp_usVwang_js)
  #wangVwang_js.append(temp_wangVwang_js)

if same_chk <> len(props): print "Error, same check failed"

names = ['HV07', 'HV08', 'HV09', 'Wang1', 'Wang2', 'Wang3', 'Warren1', 'Warren2', 'Warren3', 'Bolotin1', 'Bolotin2', 'Bolotin3']


fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
c = 0
jitter = 0.2

for js in jsd_lists:
  #print js, c
  dat = vars()[js]
  for x in range(len(dat)):
    #print x
    if dat[x] <> []:
      #print dat[x], x, js
      if c == 0: 
	col = "cyan"		# us V us
	nam = "BLR vs BLR"
      elif c == 1 or c == 4: 
	col = "purple"		# us V wang
	nam = "BLR vs Wang"
      elif c == 2 or c == 8: 
	col = "green"		# us V warren
	nam = "BLR vs Warren"
      elif c == 3 or c == 12: 
	col = "darkblue"	# us V bol
	nam = "BLR vs Bolotin"
      elif c == 5: 
	col = "red"		# wang V wang
	nam = "Wang vs Wang"
      elif c == 6 or c == 9: 
	col = "orange"		# wang V warren
	nam = "Wang vs Warren"
      elif c == 7 or c == 13: 
	col = "darkred"		# wang V bol
	nam = "Wang vs Bolotin"
      elif c == 10: 
	col = "yellow"		# warren V warren
	nam = "Warren vs Warren"
      elif c == 11 or c == 14: 
	col = "darkgoldenrod"	# warren V bol
	nam = "Warren vs Bolotin"
      elif c == 15: 
	col = "black"		# bol V bol
	nam = "Bolotin vs Bolotin"
	
      rawpts = [x]*len(dat[x])
      jitted = [y+random.uniform(-jitter, jitter) for y in rawpts]
      
      if c in [0,1,2,3,5,6,7,10,11,15] and x%3 == 0:
	ax.scatter(jitted, dat[x], color=col, s=50, alpha=.7, label=nam)
      else:
	ax.scatter(jitted, dat[x], color=col, s=50, alpha=.7)
  c += 1     

### add average lines to each group, across relevant sections
# and add legend 

avg_usVus = np.mean(flatten(usVus_js))
avg_usVwang = np.mean([flatten(usVwang_js),flatten(wangVus_js)])
avg_usVwar = np.mean([flatten(usVwarren_js),flatten(warrenVus_js)])
avg_usVbol = np.mean([flatten(usVbol_js),flatten(bolVus_js)])
avg_wangVwang = np.mean(flatten(wangVwang_js))
avg_wangVwar = np.mean([flatten(wangVwarren_js),flatten(warrenVwang_js)])
avg_wangVbol = np.mean([flatten(wangVbol_js),flatten(bolVwang_js)])
avg_warVwar = np.mean(flatten(warrenVwarren_js))
avg_warVbol = np.mean([flatten(warrenVbol_js),flatten(bolVwarren_js)])
avg_bolVbol =np.mean(flatten(bolVbol_js))

plt.legend(loc="right")
ax.plot([-.6,2.4],[avg_usVus]*2, linewidth=6, color="cyan", alpha=.3)
ax.plot([2.6,5.4],[avg_wangVwang]*2, linewidth=6, color="red", alpha=.3)
ax.plot([5.6,8.4],[avg_warVwar]*2, linewidth=6, color="yellow", alpha=.3)
ax.plot([8.6,11.4],[avg_bolVbol]*2, linewidth=6, color="black", alpha=.3)
ax.plot([-.6,2.4],[avg_usVbol]*2, linewidth=6, color="darkblue", alpha=.3)
ax.plot([8.6,11.4],[avg_usVbol]*2, linewidth=6, color="darkblue", alpha=.3)
ax.plot([-.6,5.4],[avg_usVwang]*2, linewidth=6, color="purple", alpha=.3)
ax.plot([-.6,2.4],[avg_usVwar]*2, linewidth=6, color="green", alpha=.3)
ax.plot([5.6,8.4],[avg_usVwar]*2, linewidth=6, color="green", alpha=.3)
ax.plot([5.6,11.4],[avg_warVbol]*2, linewidth=6, color="darkgoldenrod", alpha=.3)
ax.plot([2.6,8.4],[avg_wangVwar]*2, linewidth=6, color="orange", alpha=.3)
ax.plot([2.6,5.4],[avg_wangVbol]*2, linewidth=6, color="darkred", alpha=.3)
ax.plot([8.6,11.4],[avg_wangVbol]*2, linewidth=6, color="darkred", alpha=.3)

plt.tight_layout()
#.subplots_adjust(bottom=0.25)
plt.gca().set_ylim(ymin=-.02)
#plt.xticks(rotation=70)  
#plt.xticks(range(len(jsd_lists)), names)
ax.set_xticks(range(12))
ax.set_xticklabels(names, rotation=65)
plt.xlim(-1,17)
plt.ylabel("Jensen-Shannon Divergence")
#plt.ylim(0, .3)
#plt.show()


plt.savefig(savepath + "V JSD.svg", bbox_inches='tight')


sys.exit()






fig = plt.figure(figsize=(6,5))

ax = fig.add_subplot(111)

uVu_p = [[x]*2 for x in range(3)] + [[]] * 3
ax.scatter(uVu_p[:3], usVus_js[:3], color="blue")

uVw_p = [[x]*3 for x in range(6)] 
ax.scatter(uVw_p, usVwang_js, color="cyan")

wVw_p = [[]] * 3 + [[x]*2 for x in range(3,6,1)] 
ax.scatter(wVw_p[3:], wangVwang_js[3:], color="green")

plt.show()






