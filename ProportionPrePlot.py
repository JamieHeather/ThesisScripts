# ProportionPrePlot.py v1

# Built on NaiveMemory_PrePlot.py
# Takes two given files (DCR or CDR3 say) and outputs two columns detailing the frequency with which
  # each identifier appears in both files

import sys


if (len(sys.argv) <> 3):
  print "Please supply two files, e.g. python ProportionPrePlot.py file1.txt file2.txt"
  sys.exit()
else:
  file1 = list(open(sys.argv[1], "rU"))		
  file2 = list(open(sys.argv[2], "rU"))		        

file2 = file2[:]

sum_file1 = 0		# Record the total number of clones
sum_file2 = 0

n_array = []		# Arrays that frequencies will be stored in
m_array = []

m_counter = 0		# Counter for file2-clone iteration number

shared = 0		# Boolean, if clone has been found shared between file1/mem yet

for n in file1:				# First loops through all file1

  n_freq = float(n.split(',')[5])		# Takes number before first comma as frequency, and all after as clone

  n_clone = n.split(',')[:5]

  sum_file1 += n_freq
  
  m_counter = 0

  for m in file2:			# Then scrolls through file2 looking for that same clone

    m_freq = float(m.split(',')[5])

    m_clone = m.split(',')[:5]
    
    if n_clone == m_clone:		# If there's a match, the frequencies both go into the right lists
    
      n_array.append(n_freq)
      m_array.append(m_freq)
    
      sum_file2 += m_freq
      
      shared = 1 

      del file2[m_counter]	# Memory clone gets deleted 
      
    m_counter += 1      

  if shared == 0:			# If not file1 goes into file1, and file2 gets a zero

    n_array.append(n_freq)
    m_array.append(0)

  shared = 0
      


for mt in file2:			# All leftover file2 clones (i.e. the unique ones) get added to array, with zeros for their file1

  m_freq = float(mt.split(',')[5])

  n_array.append(0)

  m_array.append(m_freq)

  sum_file2 += m_freq

print "File1, File2"			# Headers for csv

for i in range(len(n_array)):

  print str(n_array[i]/sum_file1) + "," + str(m_array[i]/sum_file2)		# Outputs x,y coordinates of n/m usage



