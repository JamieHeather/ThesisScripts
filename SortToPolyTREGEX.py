## Simple script to loop through a fastq file and remove any records that contain an ambiguous N nucleotide in their sequence
# Jamie Heather, 10th March 2013
# modified from earlier versions to find polyT starting reads in 454 data

from Bio import SeqIO # module needed for sequence input
import sys            # module needed for command line argument to get fastq name
import re	     # module for regex
import collections as coll
###################

countseqs = 0
countNs = 0
countNT = 0
countNnonT = 0
polyT = 0
nonpolyT = 0

# pass command line argument to fastqfile var

if (len(sys.argv) > 1):
     fastqfile1 = open(sys.argv[1], "rU")

 

elif (len(sys.argv) == 1):
     print "Please supply a fastq filename (i.e. python RemoveNs.py file.fastq)"
     sys.exit()


pToutfilename = (sys.argv[1])[:15] + "..._pT.fastq"
pToutput_file = open(pToutfilename, "w")
pToutput_file.close()
RCoutfilename = (sys.argv[1])[:15] + "..._RC.fastq"
RCoutput_file = open(RCoutfilename, "w")
RCoutput_file.close()


def len_hp(seq):
  # function to find the length of a homopolymer tract
  # only really applies for the MiSeq data
    # 454 data is too prone to substitutions in the homopolymer itself
  
  hp_len = 8
  hp_base = 'T'
  hp = hp_base * hp_len
  
  start = seq.index(hp)
  end = start
  for pos in range(start + hp_len, len(seq), 1):
    
    base = seq[pos]
    if base == hp_base:
      end += 1
    else:
      return(end-start + hp_len)
  

####################

ptlens = coll.Counter()

ptc = 0
pt1c = 0
pt2c = 0

record_iterator = SeqIO.parse(fastqfile1, "fastq")

for record in SeqIO.parse(fastqfile1,"fastq"):

    countseqs += 1
    sequence = str(record.seq)[:70]
   
    pT = re.search('TTTTTTTT', sequence)
    pT1 = re.search('TT+[ACGN]T+', sequence)
    pT2 = re.search('TT+[ACGN]T*[ACGN]T+', sequence)
#    pA = re.search('AAAA+[TCGN]AAAA+', sequence)

    if pT:
      ptlens[len_hp(str(record.seq))] += 1
      ptc += 1
    elif pT1:
      if len(str((pT1.group(0)))) >= 8: 
	#ptlens[len(str((pT1.group(0))))] += 1
	pt1c += 1
    elif pT2:
      if len(str((pT2.group(0)))) >= 8: 
	#ptlens[len(str((pT2.group(0))))] += 1
	pt2c += 1

#    if pT1:#or pT2:# or pA:
    if pT or pT1 and len(str((pT1.group(0)))) >= 8 or pT2 and len(str((pT2.group(0)))) >= 8:
        polyT += 1
        output_handle = open(str(pToutfilename), "a")
        SeqIO.write(record, output_handle, "fastq-sanger")
        output_handle.close() 
        
    
    else:
        nonpolyT += 1
        output_handle = open(str(RCoutfilename), "a")
        SeqIO.write(record, output_handle, "fastq-sanger")
        output_handle.close()      


print "Reading from", fastqfile1
print countseqs, "records analysed."
print polyT, "poly-T sequences"
print nonpolyT, "non-poly-T sequences"
    
sys.exit()    
################### picture of poly t lengths (only really applies to MiSeq)

# for miseq, count the length of homopols
# looks in first 70 nt for the 8nt 'seed' hp
# then looks downstream of that for more instances of same base
# keeps counting until it finds a base that isn't the same, i.e. the end of the hp



import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
suffix = ".svg"

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)

  
ax.plot(ptlens.keys(), ptlens.values())
ax.vlines(24, 0, 100000, linestyles='dashed')
ax.plot(test.keys(), test.values())
plt.ylabel("Frequency")
plt.xlabel("Homopolymer length (bp)")
plt.ylim(-500, 100000)
plt.show()


#plt.legend(loc='upper left')
plt.savefig('d2polyTlengths' + suffix, bbox_inches='tight')

test = coll.Counter()
test[25] = 80000


  
