# ConsensusMerged.py
# 28th Nov 2014

# Takes merged reads from flash (that were first demultipled with DID, so have indexes and barcodes at each end)

from __future__ import division
import collections as coll
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys 
import re
import numpy as np
import Levenshtein as lev
from acora import AcoraBuilder
from string import Template

##from time import time, clock
#import decimal as dec
#import string
#import operator as op


inputfile = sys.argv[1]
chain = sys.argv[2]



def ldist(one, two):
  # Return the Levenshtein distance of two strings
  return lev.distance(one, two)

###########################################################################################################
###########################################################################################################
                                                                             ##############################
# Reading data in                                                            ##############################
                                                                             ##############################
###########################################################################################################
###########################################################################################################


barcodes = coll.defaultdict(list)
ids = coll.defaultdict(list)

#print "Reading in barcoded data..."

counts = coll.Counter()

for read in SeqIO.parse(open(inputfile, "rU"), "fastq"):
  
  counts['total_seq_in'] += 1
  
  # check both barcodes are the same (i.e. demultiplexing worked) 
  bc1 = str(read.seq)[:12]
  bc2 = str(Seq(str(read.seq)).reverse_complement()[:12])
  
  if bc1 != bc2:
    counts['bc_diff'] += 1
    continue
  
  # check barcode is of sufficient quality
    
  bc_q_thresh = 25 # specify a phred score that each base call in the barcode must beat
  
  if all(b > bc_q_thresh for b in read.letter_annotations.values()[0][:12]): 
    transcript = str(read.seq)[24:-24]
  else:
    counts['bc_fail'] += 1
    continue
  
  # check read length is sensible

  min_length = 200
 
  if len(transcript) < min_length:
    counts['len_fail'] += 1
    continue
  
  # check no ambiguous base calls 
  
  if "N" in transcript:
    counts['N'] += 1
    continue
  
  # Do I need to institute a quality score on actual reads? 
    # Currently don't think so, as I'm going to have a minimum number of reads
    
  counts['seqs_passed_to_barcodes'] += 1
  
  # Put successfully checked reads into dict based on their barcodes
  barcodes[bc1].append(transcript)
  ids[bc1].append(str(read.id))



###########################################################################################################
###########################################################################################################
                                                                             ##############################
# Generate consensus sequences within barcode families, then collapse        ##############################
                                                                             ##############################
###########################################################################################################
###########################################################################################################

#print "Generating consensus sequences..."

read_threshold = 3 	# number of reads a barcode has to have to be collapsed

cons_seqs = coll.Counter()
cons_ids = coll.defaultdict(list)
cons_bc = coll.defaultdict(list)

for bc in barcodes:
  
  # check that this barcode has the threshold number of reads
  if len(barcodes[bc]) < read_threshold: 
    counts['bc_too_few_reads'] += 1
    continue
  
  # read all sequences into a counter
    # allows us to see which is the biggest
    # also allows us to see if all reads are the same and skip collapsing if that's the case
  
  seqs = coll.Counter()
  for x in barcodes[bc]: seqs[x] += 1
  
  if len(seqs) == 1:
    counts['bc_all_seq_same'] += 1
    cons_seqs[x] += len(barcodes[bc])
    cons_bc[x].append(bc)
    for t in ids[bc]:
      cons_ids[x].append(t)
    continue

  # if the most common sequence is a certain amount bigger than the next biggest clone, take that as genuine
  
  amount_bigger = 3
  
  top2 = seqs.most_common(2)

  if top2[0][1] > ( top2[1][1] * amount_bigger ):
    counts['bc_from_proto'] += 1
    cons_seqs[top2[0][0]] += len(barcodes[bc])
    cons_bc[top2[0][0]].append(bc)
    for t in ids[bc]:
      cons_ids[top2[0][0]].append(t)    
    continue
  
  # otherwise will have to generate a consensus the old fashioned way
    # generate an array of all reads and take the most common base at each position
  
    # first see whether all reads are the same length
  
  lens = coll.Counter()
  
  for s in seqs:
    lens[len(s)] += seqs[s]
    
    # make sure that there are enough of same length to pass filter
  
  if lens.most_common(1)[0][1]:
    counts['bc_too_few_same_len'] += 1
    continue

  # get consensus seq from any left
    
  arr_len = lens.most_common(1)[0][0]	# set the length of array as most frequently used length
  arr_height = lens.most_common(1)[0][1]	# set the length of array as most frequently used length  
  
  consensus = ""
  
    # make array: rows = # reads of same length, columns = # bases in reads
  read_arr = np.chararray(( arr_height, arr_len ))
  
  # enter sequence characters into rows
  row = 0
  for sq in barcodes[bc]:
    if len(sq) == arr_len:
      read_arr[row:row+1] = sq
      row += 1
      
  # take consensus character and add to overall string
  for c in range(arr_len):
    base_count = coll.Counter(read_arr[:,c])
    # check that most frequent base at a position is uniquely common (i.e. no draw)
    if len(base_count) > 1:
      if base_count.most_common(2)[0][1] == base_count.most_common(2)[1][1]:
	base_count['N'] += base_count.most_common(1)[0][1] + 1
    # otherwise, having passed the tests, add the most common base at that position to the consensus!
    consensus = consensus + base_count.most_common(1)[0][0]
  
  if 'N' in consensus:
    counts['bc_with_N'] += 1
  
  if len(consensus) == arr_len:
    cons_seqs[consensus] += len(barcodes[bc])
    cons_bc[consensus].append(bc)
    for t in ids[bc]:
      cons_ids[consensus].append(t)
    counts['bc_passed'] += 1
  else:
    counts['bc_wrong_len_cons'] += 1
  
  #print "Erroroororororororr00r0r0rrr0r0R"

# Then collapse on barcodes

dist_thresh = 3
col_reads = coll.Counter()
  
for cons in cons_bc.keys():
  tracking = coll.Counter()
  clustered = coll.Counter()
  cbc = cons_bc[cons]         # find all barcodes associated with a given consensus      
  for b1 in cbc:
    if tracking[b1] == 0:
      ref_bc = b1
      tracking[b1] = 1
      clustered[b1] = 1
    for b2 in cbc:
      if tracking[b2] == 0:
	if ldist(ref_bc, b2) < dist_thresh:
	  tracking[b2] = 1
	  clustered[ref_bc] += 1
  col_reads[cons] = len(clustered)


###########################################################################################################
###########################################################################################################
                                                                             ##############################
# Decombine collapsed reads                                                  ############################## 
                                                                             ##############################
###########################################################################################################
###########################################################################################################

# Can operate a very trimmed down version of Decombinator
  # Don't need to work with fastqs, as just have sequences
  # Don't need to output anything additional, as sequences are pre-collapsed
  # Don't need to worry about reverse complements, as have merged all reads UTR first, so are in correct orientation
  # Don't need to worry about choice of tags, as just using original (as are internal and therefore less likely to be spliced)

###################
# Decombinator FUNCTIONS #
###################

def get_v_deletions( rc, v_match, temp_end_v, v_regions_cut ):
    # This function determines the number of V deletions in sequence rc
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in rc.
    function_temp_end_v = temp_end_v
    pos = -1
    is_v_match = 0
    while is_v_match == 0 and 0 <= function_temp_end_v < len(rc):
        if str(v_regions_cut[v_match])[pos] == str(rc)[function_temp_end_v] and str(v_regions_cut[v_match])[pos-1] == str(rc)[function_temp_end_v-1] and str(v_regions_cut[v_match])[pos-2] == str(rc)[function_temp_end_v-2]:
            is_v_match = 1
            deletions_v = -pos - 1
            end_v = function_temp_end_v
        else:
            pos -= 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        return []

def get_j_deletions( rc, j_match, temp_start_j, j_regions_cut ):
    # This function determines the number of J deletions in sequence rc
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in rc.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= function_temp_start_j+2 < len(str(rc)):
        if str(j_regions_cut[j_match])[pos] == str(rc)[function_temp_start_j] and str(j_regions_cut[j_match])[pos+1] == str(rc)[function_temp_start_j+1] and str(j_regions_cut[j_match])[pos+2] == str(rc)[function_temp_start_j+2]:
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1
            
    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        return []

def get_v_tags(file_v, half_split):
    import string
    
    v_seqs = [] # Holds all V tags
    jump_to_end_v = [] # Holds the number of jumps to make to look for deletions for each V region once the corresponding tag has been found
    for line in file_v:
        elements = line.rstrip("\n") # Turns every element in a text file line separated by a space into elements in a list
        v_seqs.append(string.split(elements)[0]) # Adds elements in first column iteratively
        jump_to_end_v.append(int(string.split(elements)[1])) # Adds elements in second column iteratively

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])
    
    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]

def get_j_tags(file_j, half_split):
    import string
    
    j_seqs = [] # Holds all J tags
    jump_to_start_j = [] # Holds the number of jumps to make to look for deletions for each J region once the corresponding tag has been found

    for line in file_j:
        elements = line.rstrip("\n")
        j_seqs.append(string.split(elements)[0])
        jump_to_start_j.append(int(string.split(elements)[1]))

    half1_j_seqs = []
    half2_j_seqs = []

    for j in range(len(j_seqs)):
        half1_j_seqs.append(j_seqs[j][0:half_split])
        half2_j_seqs.append(j_seqs[j][half_split:])

    return [j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j]


# make a dictionary that contains all the details about a given transcript, including DCRs

def nest():
  # Allows the creation of nested defaultdicts
  return coll.defaultdict(list)

transcripts = coll.defaultdict(nest)

# also have a list of rearrangements so we can see if we get multiples of the same
dcrs = coll.defaultdict(list)
dcrfreqs = coll.Counter()

dcrcounts = coll.Counter()
 
#print "Assigning recombination features using Decombinator functions"

# Do not change - V tags are split at position 10, J at position 10, to look for half tags if no full tag is found.
if chain == "a":
  v_half_split, j_half_split = [10,10] 
elif chain == "b":
  v_half_split, j_half_split = [10,6]

# Look for tag and V/J fasta files: if cannot find locally, sources from GitHub repositories
handle = open("human_TR"+ chain.upper() + "V_region.fasta", "rU")
v_genes = list(SeqIO.parse(handle, "fasta"))
handle.close()

handle = open("human_TR"+ chain.upper() + "J_region.fasta", "rU")
j_genes = list(SeqIO.parse(handle, "fasta"))
handle.close()

v_regions = []
v_nams = []
for v in range(0, len(v_genes)):
  v_regions.append(str(v_genes[v].seq).upper())
  v_nams.append(v_genes[v].id.split("|")[1])

j_regions = []
j_nams = [] 
for j in range(0, len(j_genes)):
  j_regions.append(str(j_genes[j].seq).upper())
  j_nams.append(v_genes[v].id.split("|")[1])

## Build keyword tries of V and J tags for fast assignment
v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("tags_tr"+ chain.lower() + "v.txt", "rU"), v_half_split)
j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("tags_tr"+ chain.lower() + "j.txt", "rU"), j_half_split)

v_builder = AcoraBuilder()
for i in range(0,len(v_seqs)):
    v_builder.add(str(v_seqs[i])) # Add all V tags to keyword trie

v_key = v_builder.build()

j_builder = AcoraBuilder()
for i in range(0,len(j_seqs)):
    j_builder.add(str(j_seqs[i])) # Add all J tags to keyword trie

j_key = j_builder.build()

## Build keyword tries for first and second halves of both V and J tags
v_half1_builder = AcoraBuilder()
for i in range(0,len(half1_v_seqs)):
    v_half1_builder.add(str(half1_v_seqs[i]))
half1_v_key = v_half1_builder.build()

v_half2_builder = AcoraBuilder()
for i in range(0,len(half2_v_seqs)):
    v_half2_builder.add(str(half2_v_seqs[i]))
half2_v_key = v_half2_builder.build()

j_half1_builder = AcoraBuilder()
for i in range(0,len(half1_j_seqs)):
    j_half1_builder.add(str(half1_j_seqs[i]))
half1_j_key = j_half1_builder.build()

j_half2_builder = AcoraBuilder()
for i in range(0,len(half2_j_seqs)):
    j_half2_builder.add(str(half2_j_seqs[i]))
half2_j_key = j_half2_builder.build()


## Begin analysing sequences

for cseq in col_reads:
  
  transcripts[cseq]['freq'] = col_reads[cseq]

  found_v_match = 0
  found_j_match = 0
  dcrcounts['dcr_seqs'] += 1

  v_seq_start = 0     
  j_seq_end = 0                          
      
  hold_v = v_key.findall(cseq)
  hold_j = j_key.findall(cseq)

  if hold_v:
    v_match = v_seqs.index(hold_v[0][0]) # Assigns V
    temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
    v_seq_start = hold_v[0][1]      
    
    if get_v_deletions( cseq, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
	[ end_v, deletions_v] = get_v_deletions( cseq, v_match, temp_end_v, v_regions )
  else:
    found_v_match = 0
    hold_v1 = half1_v_key.findall(cseq)
    hold_v2 = half2_v_key.findall(cseq)
    
    for i in range(len(hold_v1)):
      indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
      for k in indices:
	if len(v_seqs[k]) == len(cseq[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
	  if lev.hamming( v_seqs[k], cseq[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
	    v_match = k
	    temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
	    found_v_match += 1
	    v_seq_start = hold_v1[i][1]      

    for i in range(len(hold_v2)):
      indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
      for k in indices:
	if len(v_seqs[k]) == len(cseq[hold_v2[i][1]-v_half_split:hold_v2[i][1]-v_half_split+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
	  if lev.hamming( v_seqs[k], cseq[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
	    v_match = k
	    temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
	    found_v_match += 1
	    v_seq_start = hold_v2[i][1] - v_half_split      

  if hold_j:
    j_match = j_seqs.index(hold_j[0][0]) # Assigns J
    temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
    
    j_seq_end = hold_j[0][1] + len(hold_j[0][0])      
    
    if get_j_deletions( cseq, j_match, temp_start_j, j_regions ):
      [ start_j, deletions_j] = get_j_deletions( cseq, j_match, temp_start_j, j_regions )
      
  else:
    found_j_match = 0
    hold_j1 = half1_j_key.findall(cseq)
    hold_j2 = half2_j_key.findall(cseq)
    for i in range(len(hold_j1)):
      indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
      for k in indices:
	if len(j_seqs[k]) == len(cseq[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
	  if lev.hamming( j_seqs[k], cseq[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
	    j_match = k
	    temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
	    found_j_match += 1
	    j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split                                              

    for i in range(len(hold_j2)):
      indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
      for k in indices:
	if len(j_seqs[k]) == len(cseq[hold_j2[i][1]-j_half_split:hold_j2[i][1]-j_half_split+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
	  if lev.hamming( j_seqs[k], cseq[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
	    j_match = k
	    temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
	    found_j_match += 1
	    j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])                                                

  TCRseq = str(cseq[v_seq_start:j_seq_end])       

  if hold_v and hold_j:
    if get_v_deletions( cseq, v_match, temp_end_v, v_regions ) and get_j_deletions( cseq, j_match, temp_start_j, j_regions ):
  
      if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
	dcrcounts['imposs_del'] += 1                    
      elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
	  dcrcounts['tag_overlap'] += 1                     
      else: 
	dcr = tuple([ v_match, j_match, deletions_v, deletions_j, cseq[end_v+1:start_j] ])
	dcrs[dcr].append(cseq)
	dcrfreqs[dcr] += 1
	transcripts[cseq]['dcr'] = dcr
	transcripts[cseq]['tcrseq'] = TCRseq
	dcrcounts['assigned'] += 1                   
	    
  elif hold_v and found_j_match == 1:
    if get_v_deletions( cseq, v_match, temp_end_v, v_regions ) and get_j_deletions( cseq, j_match, temp_start_j, j_regions ):
      [ start_j, deletions_j] = get_j_deletions( cseq, j_match, temp_start_j, j_regions )
	  
      if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
	dcrcounts['imposs_del'] += 1                    
      elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
	dcrcounts['tag_overlap'] += 1                     
      else:        
	dcr = tuple([ v_match, j_match, deletions_v, deletions_j, cseq[end_v+1:start_j] ])
	dcrs[dcr].append(cseq)
	dcrfreqs[dcr] += 1
	transcripts[cseq]['dcr'] = dcr
	transcripts[cseq]['tcrseq'] = TCRseq
	dcrcounts['assigned'] += 1                   
	  
  elif found_v_match == 1 and hold_j:
    if get_v_deletions( cseq, v_match, temp_end_v, v_regions ) and get_j_deletions( cseq, j_match, temp_start_j, j_regions ):
      [ end_v, deletions_v] = get_v_deletions( cseq, v_match, temp_end_v, v_regions )
      
      if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
	dcrcounts['imposs_del'] += 1                    
      elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
	dcrcounts['tag_overlap'] += 1                     
      else:        
	dcr = tuple([ v_match, j_match, deletions_v, deletions_j, cseq[end_v+1:start_j] ])
	dcrs[dcr].append(cseq)
	dcrfreqs[dcr] += 1
	transcripts[cseq]['dcr'] = dcr
	transcripts[cseq]['tcrseq'] = TCRseq       
	dcrcounts['assigned'] += 1                   

  elif found_v_match == 1 and found_j_match == 1:
    if get_v_deletions( cseq, v_match, temp_end_v, v_regions ) and get_j_deletions( cseq, j_match, temp_start_j, j_regions ):
      [ end_v, deletions_v] = get_v_deletions( cseq, v_match, temp_end_v, v_regions )
      [ start_j, deletions_j] = get_j_deletions( cseq, j_match, temp_start_j, j_regions )
	  
      if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
	dcrcounts['imposs_del'] += 1                    
      elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
	dcrcounts['tag_overlap'] += 1                     
      else:        
	dcr = tuple([ v_match, j_match, deletions_v, deletions_j, cseq[end_v+1:start_j] ])
	dcrs[dcr].append(cseq)
	dcrfreqs[dcr] += 1
	transcripts[cseq]['dcr'] = dcr
	transcripts[cseq]['tcrseq'] = TCRseq
	dcrcounts['assigned'] += 1                   


###########################################################################################################
###########################################################################################################
                                                                             ##############################
# See whether rearrangment is functional or not (CDR3ulating)                ############################## 
                                                                             ##############################
###########################################################################################################
###########################################################################################################

#print "Checking rearrangment functionality..."


# Get conserved V regions (genes imported earlier)

with open("TR" + chain.upper() + "V_ConservedC.txt", 'r') as f:
  vconservedc = [int(line.rstrip('\n')) for line in f]
  
cdr3counts = coll.Counter()

for tr in transcripts:
  
  # Checks the productivity of a given DCR-assigned rearrangement 
  # NB: A productively rearranged receptor does not necessarily mean that it is the working receptor used in a cell
    # It could be a silenced chain that isn't used, or could have inactivating mutations upstream of the sequenced region
  
  # only need to bother if that transcript has been decombined
  
  if 'dcr' not in transcripts[tr].keys():
    cdr3counts['no_dcr'] += 1
    continue
  else:
    dcr = transcripts[tr]['dcr']
  
  # 0.5 Set up check variables
  
  # Boolean productivity checks that CDR3s must pass 
  in_frame = 0
  no_stop = 0
  found_c = 0
  found_fgxg = 0
  
  # CDR3-defining positions 
  start_cdr3 = 0
  end_cdr3 = 0 
  
  # 1. Rebuild whole nucleotide sequence from Decombinator assignment
  v = int(dcr[0])
  j = int(dcr[1])
  vdel = int(dcr[2])
  jdel = int(dcr[3])
  ins_nt = dcr[4]

  if vdel == 0:
    v_used = v_regions[v]
  else:
    v_used = v_regions[v][:-vdel]
  j_used = j_regions[j][jdel:]

  nt = ''.join([v_used, ins_nt, j_used])
  
  # 2. Translate
  aa = str(Seq(nt, generic_dna).translate())

  # 3. Check whether whole rearrangement is in frame
  if (len(nt)-1) % 3 == 0:
    in_frame = 1
  else:
    if '*' in aa:
      cdr3counts["OOF_with_stop"] += 1
      transcripts[tr]['cdr3fail'] = "OOF_with_stop"
      continue
    else:
      cdr3counts["OOF_without_stop"] += 1
      transcripts[tr]['cdr3fail'] = "OOF_without_stop"
      continue


  # 4. Check for stop codons in the in-frame rearrangements
  if '*' not in aa:
    no_stop = 1
  else:
    cdr3counts["IF_with_stop"] += 1
    transcripts[tr]['cdr3fail'] = "IF_with_stop"
    continue

  # 5. Check for conserved cysteine in the V gene
  if aa[vconservedc[v]-1] == 'C':
    found_c = 1
    start_cdr3 = vconservedc[v]-1
  elif chain == "b" and v in [45, 56]: # Allow for TRBV17 and TRBV26, which use Y instead of C to start CDR3s
    if aa[vconservedc[v]-1] == 'Y':
      found_c = 1
      start_cdr3 = vconservedc[v]-1
    
  else:
    cdr3counts["No_conserved_cysteine"] +=1
    transcripts[tr]['cdr3fail'] = "No_conserved_cysteine"
    continue
  
  # 5.5 Having found conserved cysteine, only need look downstream to find other end of CDR3
  downstream_c = aa[start_cdr3:]

  # 6. Check for presence of FGXG motif (or equivalent)
  if chain == 'a':      

    # TRAJs get a little more interesting with their conserved sequences than the other genes
      # All TRAJ FGXGs are at the -11 position apart from one
        # TRAJ59 (DCR 58) is truncated in its 3', and thus its FGXG is at -9:-5
      # All TRAJS use FGXG as their CDR3 ending motif, apart from the following:
        #TRAJ16 (DCR 6) = FARG
        #TRAJ33 (DCR 22) = WGAG
        #TRAJ38 (DCR 26) = WGLG
        #TRAJ35 (DCR 54) = CGSG
        #TRAJ51 (DCR 55) = FGKE
        #TRAJ55 (DCR 56) = WGKG
        #TRAJ61 (DCR 60) = FGAN      
    
    if j <> 58:
      site = downstream_c[-11:-7]
      
      if re.findall('FG.G', site):
	end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 1
      
      else: # Allow for the non-canonical FGXGs in TRAJ      

	awkward_ajs = [6, 22, 26, 54, 55, 56, 60]
	alt_aj_motifs = ['FARG', 'WGAG', 'WGLG', 'CGSG', 'FGKE', 'WGKG', 'FGAN']
	
	if j in awkward_ajs and site == alt_aj_motifs[awkward_ajs.index(j)]:
	  end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 1
   
	else:
	  cdr3counts["No_conserved_FGXG"] += 1
	  transcripts[tr]['cdr3fail'] = "No_conserved_FGXG"
	  continue
	  
    else: # TRAJ59
      site = downstream_c[-9:-5]     
      
      if re.findall('FG.G', site):
	end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 1 
    
  elif chain == 'b':
    # All F TRBJ FGXGs being at -10 position
      # TRBJ2-2P is at -8:-4 (although I am unaware of any evidence in our own or others' data that it gets recombined)
    
    site = downstream_c[-10:-6]
    
    if re.findall('FG.G', site):
      end_cdr3 = len(downstream_c) - 10 + start_cdr3 + 1

    elif j == 13: # TRBJ2-2Ps
      site = downstream_c[-8:-4]
      
      if re.findall('LGGG', site):
	end_cdr3 = len(downstream_c) - 8 + start_cdr3 + 1      

    else:
      cdr3counts["No_conserved_FGXG"] += 1  
      transcripts[tr]['cdr3fail'] = "No_conserved_FGXG"
      continue
  
  transcripts[tr]['cdr3'] =  aa[start_cdr3:end_cdr3]    
  transcripts[tr]['cdr3fullaa'] = aa


###########################################################################################################
###########################################################################################################
                                                                             ##############################
# See whether transcript is functional, and whether is spliced               ############################## 
                                                                             ##############################
###########################################################################################################
###########################################################################################################

#print "Checking transcript functionality..."

# Functions to find open reading frames

def findsss(aaseq):
  # Once given an amino acid sequence, returns the indexes of any start and stop codons amino acids it finds (M/*)
  starts = [i for i, l in enumerate(aaseq) if l == "M"]
  stops = [i for i, l in enumerate(aaseq) if l == "*"]
  return([starts, stops])

def orfscan(sss):
  # Takes the output of findsss function, i.e. two part list detailing start and stop codon matches 
  # Outputs any possible complete uORFs based on those criteria
  # NB some previous papers (e.g. Calvo et al PNAS 09 paper) disallow overlapping uORFs
    # Given the principle of leaky scanning, that doesn't seem sensible to me, so I'm currently allowing them
  ORF_length_filter = 10
  temp_orfs = []
  for i in range(len(sss[0])):
    if [x for x in sss[1] if x > sss[0][i]]:
      temp_orfs.append([[sss[0][i], x] for x in sss[1] if x > sss[0][i]])
  flattened = [orf for lst in temp_orfs for orf in lst]  
  out_orfs = []
  for y in flattened:
    if y[1] - y[0] >= ORF_length_filter:
      out_orfs.append(y)   
  return(out_orfs)	# ensure output flattened dictionary in order for UTRs to be correctly counted


translcounts = coll.Counter()

# sequences for the end of constant regions, to look for correct frame

trac = "IQNPDPA"
trbc1 = "DLKNVF"
trbc2 = "DLNKVF"

# genes that have shown to be spliceable
if chain == "a":
  spgenes = [44, 8, 34, 0, 16, 38, 7, 40]
elif chain == "b":
  spgenes = [43, 4, 1, 42, 12, 38, 37, 14, 40]



# look through each transcript, translating (taking frame with constant region sequence as right), noting splice sites

for tr in transcripts:  
  
  # find the correct reading frame (with respect to in-frame constant region)
  
  found = -1
  for f in range(3):
    aa = str(Seq(tr[f:], generic_dna).translate())    
    if chain == "a":    
      if trac in aa:
	# be wary might find constant region sequences cropping up unexpectly elsewhere
	if found == -1:
	  found = f
	else:
	  found = [found, f]
    elif chain == "b":
      if trbc1 in aa:
	if found == -1:
	  found = f
	else:
	  found = [found, f]
      if trbc2 in aa:
	if found == -1:
	  found = f
	else:
	  found = [found, f]
	  
  if found == -1:
    transcripts[tr]['trXc'] = "not_found"
    translcounts['c_not_found'] += 1
    
  elif isinstance(found,list):
    transcripts[tr]['trXc'] = "not_unambig_found"
    translcounts['c_not_unambig_found'] += 1
  
  else:
    aa = str(Seq(tr[found:], generic_dna).translate())   
    transcripts[tr]['frame'] = found
    transcripts[tr]['aa'] = aa
    translcounts['c_found'] += 1
    
    # first define the orf that contains the constant region (if there is one)
    
    if chain == "a":
      
      # move back from edge of constant region look for last M (ie start) or stop
      c_lim = aa.find(trac)
      
      start = -1
      for a in range(c_lim, -1, -1):
	if aa[a] == "M":
	  start = a
	#elif aa[a] == "*":
    
  # find start and stop codons and ORFs
  
  ss = findsss(aa)
  
  #orfs = orfscan(ss)
  
  if ss[0] and ss[1]:
    # both start and stop sites
    # potentially productive
    transcripts[tr]['pot_prod'] = "start+stop"
    
  elif ss[0] and not ss[1]:
    # starts no stops
    # potenially productive
    transcripts[tr]['prob_prod'] = "startnostop"
    
  elif not ss[0] and ss[1]:
    # no starts but stops
    # not productive
    transcripts[tr]['not_prod'] = "nostart+stop"
    
  elif not ss[0] and not ss[1]:
    # no starts no stops    
    # not productive (unless uses non-canonical start site)
    transcripts[tr]['unlikely_prod'] = "nostartnostop"
  
  #print "---------------------------------------------------------------"
  #for x in orfs:
    #print aa[x[0]:x[1]]
  #if len(orfs) > 4:
    #sys.exit()
  
  # Check for known splice sites
  
  if chain == "a":
    splicenams = ['TRAV1-1|D1A2', 'TRAV1-1|D2A2', 'TRAV1-1|D3A1.5', 'TRAV1-1|D3A2', 'TRAV1-1|D4A2', 'TRAV13-2|D2A2', 'TRAV14/DV4|D2A2', 'TRAV14/DV4|D2A2', 'TRAV22|D1.5A2', 'TRAV41|D2A2', 'TRAV8-1|D2A2', 'TRAV8-3|D2A2', 'TRDV1|D2A2', 'TRDV1|D2A3'] 
    splices = ['CCATGAAGATGGGAGGAGCTCCAGATGAAA', 'TTGTCCAGATAAACTGAGCTCCAGATGAAA', 'TTATGGGCTGTCCTGGGAGCACCCACATTT', 'TTATGGGCTGTCCTGGAGCTCCAGATGAAA', 'TTCTTTCTTACAATGGAGCTCCAGATGAAA', 'ACAAGCAAGAATCTGGAAAAGCTACTCAACCTGGAGACTCA', 'TCAAAGTTATGGTCTATTCTGAAGGCAAGAAAATCCGCCAA', 'TCCAAGTTATGGTCTATTCTGAAGGCAAGAAAATCCGCCAA', 'TGTGAGTTGGGGCTGTCCAAGGGACAAAACAGAATGGAAGA', 'ACAATCAACTGCAGTTACTCGGAAAAGCACAGCTCCCTGCA', 'TGGAACTGTTAATCTCTTCTGGAAACCCTCTGTGCAGTGG', 'GGCAACACCTTATCTCTTCTGGAAACCCTCTGTGCATTGGA', 'AGCAGTCACCCTGAACTGCCTGGTTCTGATGAACAGAATGC', 'AGCAGTCACCCTGAACTGCCTAATGCAAAAAGTGGTCGCTA']
  elif chain == "b":
    splicenams = ['TRBV10-2|D2A2', 'TRBV10-2|D3A3', 'TRBV11-2|D1A2', 'TRBV11-2|D2A3', 'TRBV18|D2A2', 'TRBV2|D2A2', 'TRBV7-2|D2A2', 'TRBV7-2|D2A3', 'TRBV7-2|D2A4', 'TRBV7-3|D2A2', 'TRBV7-3|D2A3', 'TRBV7-6|D-1A2', 'TRBV7-6|D1A2', 'TRBV7-6|D2A2', 'TRBV7-6|D3A2', 'TRBV7-8|D1A3', 'TRBV7-8|D2A2', 'TRBV7-8|D2A3', 'TRBV7-9|D2A2', 'TRBV7-9|D2A3']
    splices = ['CAAGATCACAGAGACAGGAAGATAAAGGAGAAGTCCCCGAT', 'GAGCCACAGCTATATGTTCTGATCCAAGACAGAGAATTTCC', 'CCCTCTGTCTCCTGGGAGCAGAGAGGCTCAAAGGAGTAGAC', 'TGGCCATGCTACCCTTTACTGGCTCAAAGGAGTAGACTCCA', 'ATCGGCAGCTCCCAGAGGAAGAGGGCCCCAGCATCCTGAGG', 'TAATCACTTATACTTCTATTGTTGAAAGGCCTGATGGATCA', 'GGGAAAGGATGTAGAGCTCAGGCAACAGTGCACCAGACAAA', 'GGGAAAGGATGTAGAGCTCAGTGCACCAGACAAATCAGGGC', 'GGGAAAGGATGTAGAGCTCAGAGAGGACTGGGGGATCCGTC', 'GGGAAAATATGTAGAGCTCAGGCACGGGTGCGGCAGATGAC', 'GGGAAAATATGTAGAGCTCAGTCAGGCCTGAGGGATCCGTC', 'GCTCACAGTGACACTGATCTGAGAGGCCTGAGGGATCCATC', 'TCCTGGGTTTCCTAGGGACAGAGAGGCCTGAGGGATCCATC', 'TGGAGTCTCCCAGTCTCCCAGAGAGGCCTGAGGGATCCATC', 'GGGACAGGATGTAGCTCTCAGAGAGGCCTGAGGGATCCATC', 'TCCTGGGTTTCCTAGGGACAGAAAGGCCTGAGGGATCCGTC', 'AGGACAGGATGTAGCTCTCAGAATGAAGCTCAACTAGACAA', 'AGGACAGGATGTAGCTCTCAGAAAGGCCTGAGGGATCCGTC', 'GGGACAGAATGTAACTTTCAGAATGAAGCTCAACTAGAAAA', 'GGGACAGAATGTAACTTTCAGAGAGGCCTAAGGGATCTTTC']
    
  for i in range(len(splices)):
    if splices[i] in tr:
      transcripts[tr]['splice'].append(splicenams[i])
      #if ss[0] and not ss[1]:
	#print "yes"
	#print "----------------------\n", tr, transcripts[tr], aa[ss[0][0]:]

  ## need to make this do soemthing with translated seqs, i.e check there's a functional one, no stops, right L1, etc
  ## then finally need to collate all data
  
## TOTAL NUMBER TRANSCRIPTS
#tot = sum([transcripts[x]['freq'] for x in transcripts])
  
## TOTAL NUMBER TRANSCRIPTS THAT DECOMBINED (unique * their frequencies)
#tot_dcr = sum([transcripts[x]['freq'] for x in transcripts if "dcr" in transcripts[x].keys() ])

## TOTAL THAT WERE RECOMBINATORIALLY PRODUCTIVE
#tot_cdr = sum([transcripts[x]['freq'] for x in transcripts if "cdr3" in transcripts[x].keys() ])

## TOTAL THAT WEREN'T RECOMBINATORIALLY PRODUCTIVE
#tot_cdrfail = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() ])

## TOTAL THAT WERE RECOMBINATORIALLY PRODUCTIVE AND SPLICED
#tot_sp_cdr = sum([transcripts[x]['freq'] for x in transcripts if "cdr3" in transcripts[x].keys() and "splice" in  transcripts[x].keys()])

## TOTAL THAT WEREN'T RECOMBINATORIALLY PRODUCTIVE AND SPLICED
#tot_sp_cdrfail = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in  transcripts[x].keys()])


#print "\n---------------------------", inputfile, "----------------------------------------\n"

##print "Number collapsed transcripts\t\t\t\t", str(tot)
##print "Number decombined transcripts\t\t\t\t", str(tot_dcr)

##print "Number decombined transcripts with cdr3s\t\t", str(tot_cdr)
##print "Number decombined transcripts without cdr3s\t\t", str(tot_cdrfail)

##print "Number decombined+spliced transcripts with cdr3s\t", str(tot_sp_cdr)
##print "Number decombined+spliced transcripts without cdr3s\t", str(tot_sp_cdrfail)



## TOTAL THAT WEREN'T RECOMBINATORIALLY PRODUCTIVE AND SPLICED BY REASON
#tot_sp_cdrfail_nof = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "No_conserved_FGXG"])
#tot_sp_cdrfail_noc = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "No_conserved_cysteine"])
#tot_sp_cdrfail_ifw = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "IF_with_stop"])
#tot_sp_cdrfail_oofwo = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "OOF_without_stop"])
#tot_sp_cdrfail_oofw = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "OOF_with_stop"])

#tot_cdrfail_nof = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" not in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "No_conserved_FGXG"])
#tot_cdrfail_noc = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" not in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "No_conserved_cysteine"])
#tot_cdrfail_ifw = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" not in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "IF_with_stop"])
#tot_cdrfail_oofwo = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" not in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "OOF_without_stop"])
#tot_cdrfail_oofw = sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" not in  transcripts[x].keys() and transcripts[x]['cdr3fail'] == "OOF_with_stop"])

##print "Numbers of specific fails (unspliced/spliced)"
##print "No FGXG\t", str(tot_cdrfail_nof), "\t", str(tot_sp_cdrfail_nof)
##print "No C\t", str(tot_cdrfail_noc), "\t", str(tot_sp_cdrfail_noc)
##print "IF+stop[\t", str(tot_cdrfail_ifw), "\t", str(tot_sp_cdrfail_ifw)
##print "OOF no stop\t", str(tot_cdrfail_oofwo), "\t", str(tot_sp_cdrfail_oofwo)
##print "OOF+stop\t", str(tot_cdrfail_oofw), "\t", str(tot_sp_cdrfail_oofw)

## Breakdown of probably transcriptional functionality of productive CDR3s spliced vs unspliced

#unspliced = [x for x in transcripts if "cdr3" in transcripts[x].keys() and "splice" not in transcripts[x].keys() and "trXc" not in transcripts[x].keys() ]
#spliced = [x for x in transcripts if "cdr3" in transcripts[x].keys() and "splice" in transcripts[x].keys() and "trXc" not in transcripts[x].keys() ]

#tot_unspliced = sum( [transcripts[x]['freq'] for x in unspliced] )
#tot_spliced = sum( [transcripts[x]['freq'] for x in spliced] )

#un_prob = sum( [transcripts[x]['freq'] for x in unspliced if "prob_prod" in transcripts[x].keys()] )
#sp_prob = sum( [transcripts[x]['freq'] for x in spliced if "prob_prod" in transcripts[x].keys()] )

## little function to assess whether a potentially transcriptionally productive transcript (ie both * and Ms)
#def potprob(foundsss):
  #if max(foundsss[0]) > max(foundsss[1]):
    #return True
  #else:
    #return False

#un_potprob = sum( [transcripts[x]['freq'] for x in unspliced if "pot_prod" in transcripts[x].keys() and potprob(findsss(transcripts[x]['aa'])) == True ] )
#sp_potprob = sum( [transcripts[x]['freq'] for x in spliced if "pot_prod" in transcripts[x].keys() and potprob(findsss(transcripts[x]['aa'])) == True ] ) 

#un_potprobnot = sum( [transcripts[x]['freq'] for x in unspliced if "pot_prod" in transcripts[x].keys() and potprob(findsss(transcripts[x]['aa'])) == False ] )
#sp_potprobnot = sum( [transcripts[x]['freq'] for x in spliced if "pot_prod" in transcripts[x].keys() and potprob(findsss(transcripts[x]['aa'])) == False ] ) 

#un_nostart = sum( [transcripts[x]['freq'] for x in unspliced if "unlikely_prod" in transcripts[x].keys() or "not_prod" in transcripts[x].keys()] )
#sp_nostart = sum( [transcripts[x]['freq'] for x in spliced if "unlikely_prod" in transcripts[x].keys() or "not_prod" in transcripts[x].keys()] )

#print "\nNumber P unspliced:\t", str(tot_unspliced)
#print "Probable translated (M no *):\t", str(un_prob)
#print "Probable translated (M + *):\t", str(un_potprob)
#print "Improbable translated (M + *):\t", str(un_potprobnot)
#print "Not translated  (no M):\t", str(un_nostart)

#print "\nNumber P spliced:\t", str(tot_spliced)
#print "Probable translated (M no *):\t", str(sp_prob)
#print "Probable translated (M + *):\t", str(sp_potprob)
#print "Improbable translated (M + *):\t", str(sp_potprobnot)
#print "Not translated  (no M):\t", str(sp_nostart)

## TOTAL THAT WERE SPLICED AND REMAIN PROBABLY TRANSCRIPTIONALLY PRODUCTIVE
##tot_sp_prob = sum([transcripts[x]['freq'] for x in transcripts if "prob_prod" in transcripts[x].keys() and "splice" in  transcripts[x].keys()])
##tot_sp_pot = sum([transcripts[x]['freq'] for x in transcripts if "pot_prod" in transcripts[x].keys() and "splice" in  transcripts[x].keys()])
##tot_sp_unl = sum([transcripts[x]['freq'] for x in transcripts if "unlikely_prod" in transcripts[x].keys() and "splice" in  transcripts[x].keys()])
##tot_sp_not = sum([transcripts[x]['freq'] for x in transcripts if "not_prod" in transcripts[x].keys() and "splice" in  transcripts[x].keys()])




#sys.exit()


## output some examples of each splice type for plotting
#outf = open("exampleBsplices.fa", "a")

##t = [x for x in transcripts if "splice" in transcripts[x].keyTTATGGGCTGTCCTGGGAGCACCCACATTTs() and transcripts[x]['splice'][0] == splicenams[2]]
##for l in t:
  ##outl = ">" + transcripts[l]['splice'][0] + "_" + str(transcripts[l]['freq']) + "\n" + l
  ##print >> outf, outl

#for i in range(len(splicenams)):
  ##q = []
  #count = 0
  #for t in transcripts:
    #if count < 10:
      #if 'splice' in transcripts[t].keys():
	#if transcripts[t]['splice'] == [splicenams[i]]:
	  #outl = ">" + str(transcripts[t]['splice']) + "_" + str(transcripts[t]['freq']) + "\n" + t
	  #print >> outf, outl 
	#count += 1
      #else:
	#continue
   
#for i in range(0):
  ##q = []
  #count = 0
  #for t in transcripts:
    #if count < 10:
      #if 'splice' in transcripts[t].keys():
        #if transcripts[t]['splice'] == [splicenams[3]]:
	  #outl = ">" + str(transcripts[t]['splice']) + "_" + str(transcripts[t]['freq']) + "\n" + t
	  #print >> outf, outl 
	#count += 1
      #else:
	#continue 
  ##t = [x for x in transcripts if "splice" in transcripts[x].keys() and transcripts[x]['splice'][0] == splicenams[i]]
  ##count=0
  ##for l in t:
    ##if count < 6:
      ##outl = ">" + transcripts[l]['splice'] + "_" + str(transcripts[l]['freq']) + "\n" + l
      ##print >> outf, outl
    ##else:
      ##continue
    ##count += 1


    
#outf.close()






#sys.exit()  

#[transcripts[x]['freq'] for x in transcripts if "prob_prod" in transcripts[x].keys() and "splice" in  transcripts[x].keys() and "cdr3" in transcripts[x].keys()]


#[x for x in transcripts if "cdr3" in transcripts[x].keys() and transcripts[x]['cdr3'] == "CADPGGYNKLIF"]


##len([x for x in transcripts if "dcr" in transcripts[x].keys()])

##len([x for x in transcripts if "cdr3" in transcripts[x].keys() and "splice" in transcripts[x].keys()])

##len([x for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in transcripts[x].keys()])

##for y in ([x for x in transcripts if "cdr3" in transcripts[x].keys() and "splice" in transcripts[x].keys()]): print "-----------------\n", transcripts[y]['aa']  




##c = coll.Counter()


##for y in [x for x in transcripts if "splice" not in transcripts[x].keys()]: 
  ##c['tot'] += 1
  ##if "*" in transcripts[y]['aa']:
    ##c['stop'] += 1
  ##else:
    ##print y, transcripts[y]
    
  


#print sum([transcripts[x]['freq'] for x in transcripts if "cdr3" in transcripts[x].keys() and "splice" in transcripts[x].keys()])

#print sum([transcripts[x]['freq'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in transcripts[x].keys()])

#sys.exit()


#[transcripts[x]['dcr'] for x in transcripts if "dcr" in transcripts[x].keys() and "splice" in transcripts[x].keys()]



#q = coll.Counter()

#for e in [transcripts[x]['cdr3fail'] for x in transcripts if "cdr3fail" in transcripts[x].keys() and "splice" in transcripts[x].keys()]:
  #q[e] += 1

##for d in dcrs:
  ##if len(dcrs[d]) > 1:
    ##both = coll.Counter()
    ##for t in dcrs[d]:
      
## want to know which clones (dcrs) exist as both spliced and unspliced in same person
    


## find big clones that don't have DCRs 
#bigdcrless = [x for x in transcripts if "dcr" not in transcripts[x].keys() and transcripts[x]['freq'] > 3]

#for b in bigdcrless:
  #print b, transcripts[b]
  

##oasis reversible shopper

#gens = coll.Counter()
#for n in splicenams:
  #gens[n.split("|")[0]] += 1
  
#spliceable = [travnam(g) for g in gens.keys()]


#spt = []
#notspt=[]
#for t in transcripts:
  #if 'dcr' in transcripts[t].keys():
    #if transcripts[t]['dcr'][0] in spgenes:
      #spt.append(t)
    #else:
      #notspt.append(t)

## READS THAT ARE RECOMBINATORIALLY PRODUCTIVE, BUT NOT TRANSCRIPTIONALLY
  ## THESE CAN BE MINED FOR NEW SPLICE SITES


#[transcripts[x]['freq'] for x in transcripts if "cdr3" in transcripts[x].keys() and "not_prod" in transcripts[x].keys()]


#print sum([transcripts[x]['freq'] for x in spt if "cdr3" in transcripts[x].keys() and "splice" in transcripts[x].keys()])

#print sum([transcripts[x]['freq'] for x in spt if "cdr3fail" in transcripts[x].keys() and "splice" in transcripts[x].keys()])

#[x for x in transcripts if "splice" in transcripts[x].keys() and transcripts[x]['splice'] == ['TRBV10-2|D2A2']]

##[transcripts[x]['splice'][0] for x in transcripts if "splice"]# in transcripts[x].keys() and transcripts[x]['splice'][0] == 'TRBV10-2|D2A2']

##[transcripts[x]['splice'][0] for x in transcripts if "splice" in transcripts[x].keys()] # in transcripts[x].keys() and transcripts[x]['splice'][0] == 'TRBV10-2|D2A2']




### see how many DCRs have both spliced and unspliced versions present

numb_diff = []

for d in dcrs:
  if len(dcrs[d]) > 1:
    b = coll.Counter()
    for t in dcrs[d]:
      for k in transcripts[t].keys():
	b[k] += 1
    if b['splice'] > 0 and b['splice'] == b['aa'] and b['cdr3'] > 0:
      numb_diff.append(d)
      #print d, b
  else:
    continue
  
print len(numb_diff)/len(dcrs.keys()), "\n---------------"

