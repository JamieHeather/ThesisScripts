# PrintErrorProbabilities.py v1

# Similar to PrintTotalAvgQuals.py v1
  # Instead of printing average Q score, prints the calculated probability of base-call errors, and their resultant accuracies
    # (As averaging Q-scores does not accurately represent the data - see http://drive5.com/usearch/manual/avgq.html)


from Bio import SeqIO 
import sys            
import numpy as np


#Q = -10 log10 P /

#therefore

#P = 10 ^ (-Q/10)

def prob(Q):
  # Returns the probability of the base being incorrect, based on quality score
  return(10**(-Q/10))
  # NB 1 - prob therefore gives the estimated base call accuracy

num_seqs = 0

probs = []
test = []

if sys.stdin.isatty() == False:

  for record in SeqIO.parse(sys.stdin, "fastq"):
    
    quals = record.letter_annotations.values()[0]
    
    probs.append([prob(q) for q in quals])
        
    num_seqs += 1
    
    if num_seqs % 1000 == 0:
      sys.stderr.write(str(num_seqs) + " | ")
      
elif sys.stdin.isatty() == True: 	# If don't have stdin input
  
  if (len(sys.argv) <> 2):
    print "Please supply a fastq filename (i.e. python script.py file.fastq) or pipe data in (i.e. cat file.fastq | python script.py"
    sys.exit()
  else:
    fastqfile = open(sys.argv[1], "rU")
  
  for record in SeqIO.parse(fastqfile,"fastq"):
        
    quals = record.letter_annotations.values()[0]
    
    probs.append([prob(q) for q in quals])  

    num_seqs += 1
    
    if num_seqs % 1000 == 0:
      sys.stderr.write(str(num_seqs) + " | ")

flattened_probs = [item for sublist in probs for item in sublist]

mean_prob = np.mean(flattened_probs)
pred_acc = (1 - mean_prob) * 100

if len(sys.argv) > 1:
  print "Running PrintErrorProbabilities, reading from\n", sys.argv[1] + ",", '{0:,}'.format(num_seqs), "sequences processed"
  
print "Mean prob:\t\tPredicted accuracy"
print str(mean_prob) + "\t" + str(pred_acc) + "%"
