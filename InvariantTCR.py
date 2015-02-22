# InvariantTCR.py
# Wed 5th Nov

# Built from /mnt/chain/TcRSeq_Processed/JH_HIV_Paper_Analysis/freq/Published_CDR3s_VJ.py

# Want to simply enumerate the number of detectable TCR sequences from invariant (and semi-invariant) TCRs 
  # Or rather the TCRs of non-MHC restricted T-cell subpopulations
  # Invariably these sequences are alpha chain, but sometimes there are associated beta chains (typically not invariant)

# CDR3s are stored in lists of lists in format: [[V1, J1, 'CDR3-1'], [V2, J2, 'CDR3-2']...]
  # This allows us to play with either VJ-linked CDR3s, or individual CDR3 sequences as we please
  # Will also try to allow redudant TCRs (i.e. from multiple papers) and just process it better than in Published_CDR3s_VJ

###############################################

## FIX! need to copy over the top of publicness and tweak so that I read all the right data in




###############################################################################################################################
###    ##   ###  #  ######  ###  ######      ###    ######  ######   ###  ##        ###########################################
####  ###    ##  ##  ####  ###    #####  ###  ###  ######    #####    ##  #####  ##############################################
####  ###  #  #  ###  ##  ###  ##  ####     #####  #####  ##  ####  #  #  #####  ##############################################
####  ###  ##    ####    ###        ###  ##  ####  ####        ###  ##    #####  ##############################################
###    ##  ###   #####  ###  ######  ##  ###  ##    ##  ######  ##  ###   #####  ##############################################
###############################################################################################################################


#############################
########## GEM ##############
#############################

# GEM are TRAV1-2+ TRAJ9+

gem_alpha_vanrhijn= [[1, 49, 'CAVLNTGGFKTIFGXG'], [1, 49, 'CAVRRTGGFKTIFGXG'], [1, 49, 'CAVRNTGGFKTIFGXG'], \
	  [1, 49, 'CAVRVTGGFKTIFGXG'], [1, 49, 'CAVRGTGGFKTIFGXG']]
  # NB CASMYTGGFKTIFGXG removed, as while is CD1b restricted, isn't GEM

gem_beta_vanrhijn= [[23, 8, 'CASSQPIGGGEQFFGXG'], [31, 7, 'CASRPPLTARGLKHTGELFFGXG'], [22, 8, 'CAWAKTGFGGDTQYFGXG'], [31, 6, 'CASSPRLAGDEQFFGXG'], [22, 8, 'CAWSPTSLRLASTDTQYFGXG']]


all_apub_lists = [x for x in dir() if "gem_alpha_" in x] 
all_bpub_lists = [x for x in dir() if "gem_beta_" in x] 

all_apub_raw = []
all_bpub_raw = []

for l in all_apub_lists:
  for c in vars()[l]:
    all_apub_raw.append(c)
    
for l in all_bpub_lists:
  for c in vars()[l]:
    all_bpub_raw.append(c)

# Prevent duplicate CDR3s making it through

apub = coll.Counter()
bpub = coll.Counter()

for x in all_apub_raw:
  apub[x] += 1
  
for x in all_bpub_raw:
  bpub[x] += 1

uniqa = apub.keys()
uniqb = bpub.keys()

all_apub = ['uniqa']
all_bpub = ['uniqb']

# NB - invariant GEM alpha chains = TRAV1-2 + TRAJ9
  # = 1 + 49 (DCR)

sys.exit()

elif pathogen == "m":


#############################
########## MAIT #############
#############################

### Predominantly TRAV1-2 + TRAJ33 ###

# Alpha DCR indexes
  # V 1 = TRAV1-2
  # J 22 = TRAJ33
  # J 2 = TRAJ12
  
# Beta DCR indexes
  # V 15 = TRBV20-1
  # J 6 = TRBJ2-1


# NB Gold et al also provides some paired beta chain phenotyping
  # In order, their betas used: TRBV30, TRBV6-5, TRBV6-5, TRBV6-2, TRBV19 and TRBV20-1

mait_alpha_gold = [[1, 22, 'CAVLDSNYQLIWGAG'], [1, 22, 'CAVRDSNYQLIQWGAG'], [1, 22, 'CAVRDSNYQLIQWGAG'], \
		   [1, 22, 'CARSDSNYQLIWGAG'], [1, 22, 'CASMDSNYQLIWGAG'], [1, 9, 'CAVNGDDYKLSFGAG']]

# NB - invariant MAIT alpha chains = TRAV1-2 + TRAJ33 
    # = decombinator indexes 1 and 22 respectively

mait_alpha_vanrhijn = [[1, 22, 'CAFMDSNYQLIWGAG']]

mait_beta_vanrhijn = [[15, 6, 'CSARTSGDFGEQFFGXG']]



mait_alpha_greenaway = [[1, 22, 'CAVRDSNYQLIWGAG'], [1, 22, 'CAVMDSNYQLIWGAG'], [1, 22, 'CAVLDSNYQLIWGAG'], \
	    [1, 22, 'CAVTDSNYQLIWGAG'], [1, 22, 'CAVKDSNYQLIWGAG'], [1, 22, 'CAVIDSNYQLIWGAG'], [1, 22, 'CAVGDSNYQLIWGAG'], \
	    [1, 22, 'CAGLDSNYQLIWGAG'], [1, 22, 'CAALDSNYQLIWGAG'], [1, 22, 'CAVQDSNYQLIWGAG'], [1, 22, 'CALLDSNYQLIWGAG'], \
	    [1, 22, 'CAAVDSNYQLIWGAG'], [1, 22, 'CAAMDSNYQLIWGAG'], [1, 22, 'CAGMDSNYQLIWGAG'], [1, 22, 'CAPLDSNYQLIWGAG'], \
	    [1, 22, 'CASMDSNYQLIWGAG'], [1, 22, 'CATMDSNYQLIWGAG'], [1, 22, 'CARSDSNYQLIWGAG'], [1, 22, 'CAIMDSNYQLIWGAG'], \
	    [1, 22, 'CAAEDSNYQLIWGAG'], [1, 22, 'CAGWDSNYQLIWGAG'], [1, 22, 'CAAIDSNYQLIWGAG'], [1, 22, 'CAFMDSNYQLIWGAG'], \
	    [1, 22, 'CAVRDRDYQLIWGAG'], [1, 22, 'CASIDSNYQLIWGAG'], [1, 22, 'CAVMDDNYQLIWGAG'], [1, 22, 'CAVVDSNYQLIWGAG']]

# Lepore beta chains = mostly sequences from TRAV1-2+TRAJ33+ MAIT cells

mait_beta_lepore = [[33, 12, 'CASRLMSGSSYEQYFGXG'],  [32, 8, 'CASSAASGGADTQYFGXG'], [15, 0, 'CSARDRRETEAFFGXG'], \
	    [15, 6, 'CSARGDREAYNEQFFGXG'], [15, 21, 'CSARGIDRVTNEQFFGXG'], [25, 10, 'CASSQERGSQETQYFGXG'], \
	    [9, 2, 'CASSLGSSGNTIYFGXG'], \ # end of donor A TRAV1-2+TRAJ33+ MAIT
	    [33, 10, 'CASSPSGGGAQETQYFGXG'], [15, 7, 'CSATGTGDTGELFFGXG'], [15, 3, 'CSAPWAGVNEKLFFGXG'], \
	    [16, 12, 'CATSREGRASNEQYFGXG'], [31, 8, 'CASSYSTSGADTQYFGXG'], [31, 12, 'CASSVGQENYEQYFGXG'], \
	    [24, 8, 'CASSQEGQGAPTDTQYFGXG'], [25, 8, 'CASSQEGQGAPTDTQYFGXG'], [34, 8, 'CASSNRVTSTDTQYFGXG'], \ 
	      # end of 2B TRAV1-2+TRAJ33+ MAIT
	    [31, 10, 'CASSDGPAEETQYFGXG'], [32, 8, 'CASSPGTASTDTQYFGXG']] # end of 2B TRAV1-2+TRAJ12+

# Lepore alpha chains sequences however = from non-TRAJ33 using MAIT (TRAJ12 and TRAJ20)

mait_alpha_lepore = [[1, 2, 'CVWMDSSYKLIFGXG'], [1, 2, 'CAVRDSSYKLIFGXG'], [1, 2, 'CAVLDSSYKLIFGXG'], \
	    [1, 2, 'CAVMDSSYKLIFGXG'], [1, 2, 'CAVVDSSYKLIFGXG'], [1, 2, 'CAVLDSSYKLIFGXG'], [1, 2, 'CAVMDSSYKLIFGXG'], \
	    [1, 2, 'CAVMDSSYKLIFGXG'], [1, 2, 'CAVMDSSYKLIFGXG'], [1, 2, 'CAVMDSSYKLIFGXG'], [1, 2, 'CAVMDSSYKLIFGXG'], \
	    [1, 2, 'CAVMDSSYKLIFGXG'], [1, 9, 'CAVRDFDYKLSFGXG'], [1, 9, 'CAVRDGDYKLSFGXG'], [1, 9, 'CAVRDNDYLSFGXG'], \
	    [1, 9, 'CAVRDNDYKLSFGXG'], [1, 9, 'CAVRDNDYKLSFGXG'], [1, 9, 'CAVRLNDYKLSFGXG'], [1, 9, 'CAVRLNDYKLSFGXG'], \
	    [1, 9, 'CAVRDADYKLSFGXG'], [1, 9, 'CAVNSNDYKLSFGXG']]

mait_alpha_reantragoon = [[1, 22, 'CAVMDSNYQLIWGXG'], [1, 2, 'CAVMDSSYKLIFGXG'], [1, 9, 'CAVRDGDYKLSFGXG']]


all_apub_lists = [x for x in dir() if "mait_alpha_" in x]

all_bpub_lists = [x for x in dir() if "mait_beta_" in x]

all_apub_raw = []
all_bpub_raw = []

for l in all_apub_lists:
  for c in vars()[l]:
    all_apub_raw.append(c)
    
for l in all_bpub_lists:
  for c in vars()[l]:
    all_bpub_raw.append(c)

# Prevent duplicate CDR3s making it through

apub = coll.Counter()
bpub = coll.Counter()

for x in all_apub_raw:
  apub[x] += 1
  
for x in all_bpub_raw:
  bpub[x] += 1

preuniqa = apub.keys()
uniqb = bpub.keys()

uniqa = []

for x in preuniqa:
  if x[:-3].endswith("QLIW"):
    uniqa.append("1-22-" + x)
  elif x[:-3].endswith("LSF"):
    uniqa.append("1-9-" + x)
  elif x[:-3].endswith("KLIF"):
    uniqa.append("1-2-" + x)
    
all_apub = ['uniqa']
all_bpub = ['uniqb']








#############################
########## NKT ##############
#############################

elif pathogen == "n":

nkt_alpha_vanrhijn = [[2, 8, 'CVVSDRGSTLGRLYFGXG']]

nkt_beta_vanrhijn = [[17, 6, 'CASSESQYGRAAYNEQFFGXG']]


nkt_alpha_greenaway = [[2, 8, 'CVVSDRGSTLGRLYFGXG'], [2, 8, 'CVVSRGSTLGRLYFGXG'], [2, 8, 'CVVTDRGSTLGRLYFGXG'], \
	[2, 8, 'CVVNDRGSTLGRLYFGXG'], [2, 8, 'CVVGDRGSTLGRLYFGXG'], [2, 8, 'CVVRADRGSTLGRLYFGXG'], \
	[2, 8, 'CVVFDRGSTLGRLYFGXG'], [2, 8, 'CVASDRGSTLGRLYFGXG'], [2, 8, 'CVVARDRGSTLGRLYFGXG'], \
	[2, 8, 'CVSVDRGSTLGRLYFGXG'], [2, 8, 'CVVSDKGSTLGRLYFGXG'], [2, 8, 'CVVSATNRGSTLGRLYFGXG'], \
	[2, 8, 'CVVPGRLYFGXG']]

nkt_beta_exley = [[17, 7, 'CASREGAMGTGELFFGEG'], [17, 8, 'CASSATRALTGSDTQYFGPG'], [17, 6, 'CASSFLDRDYSYNEQFFGPG'], \
	[17, 12, 'CASSENRQGAGYEQYFGPG'], [17, 7, 'CASSERTTNTGELFFGEG'], [17, 6, 'CASSVRPGGNEQFFGPG'], \
	[17, 0, 'CASSDGEQANTEAFFGQG'], [17, 1, 'CASSATIRDRASGYTFGSG'], [17, 7, 'CASSDTRVGGELFFGEG'], \
	[17, 4, 'CASSLGESNQPQHFGDG'], [17, 12, 'CASSVPGPAYEQYFGPG'], [17, 7, 'CASSDTRVGGELFFGEG'], \
	[17, 0, 'CASGTQGNTEAFFGQG'], [17, 1, 'CASEYGGPSYGYTFGSG']]

nkt_alpha_han = [[20, 17, 'CILRDVGNTPLVFGXG'], [34, 37, 'CAVKNFGNEKLTFGXG'], [2, 8, 'CVVSDRGSTLGRLYFGXG']]
  # removed CAVMDSNYQLIFGXG as it actually a MAIT sequence ( TRAV1-2 TRAJ33 )
    # isnt even CAVMDSNYQLIFGXG, should be CAVMDSNYQLIWGXG
  
nkt_alpha_demoulins = [[2, 13, 'CVVTLPMTTDSWGKFQFGXG'], [2, 8, 'CVVSDRGSTLGRLYFGXG'], [2, 8, 'CVVTDRGSTLGRLYFGXG'], \
	[2, 24, 'CVVSGRQTGANNLFFGXG'], [2, 48, 'CVVSIAGFQKLVFGXG'], [2, 49, 'CVVINTGGFKTIFGXG'], [2, 24, 'CVVIATGANNLFFGXG'], \
	[2, 8, 'CVVSNDRGSTLGRLYFGXG'], [36, 0, 'CALDMVAGGGNKLTFGXG'], [36, 46, 'CALDMSGGSYIPTFGXG'], [36, 0, 'CALEITGGGNKLTFGXG'],\
	[36, 3, 'CALSNSGGYQKVTFGXG'], [36, 22, 'CALDREGSNYQLIWGXG'], [36, 34, 'CAPDSGGGADGLTFGXG'], [36, 16, 'CALEVAGSYQLTFGXG'], \
	[36, 39, 'CALDTGRRALTFGXG'], [36, 20, 'CALDQSNARLMFGXG'], [36, 27, 'CALDNAGNMLTFGXG'], [36, 29, 'CALKAGTYKYIFGXG'], \
	[36, 57, 'CATETSGSRLTFGXG'], [2, 8, 'CVVRADRGSTLGRLYFGXG'], [2, 20, 'CVVSRENNARLMFGXG'], [2, 46, 'CVVSSGGSYIPTFGXG'], \
	[2, 28, 'CVVSPISGGYNKLIFGXG'], [2, 20, 'CVVSRENNARLMFGXG'], [2, 33, 'CVVSARAGTASKLTFGXG'], [2, 20, 'CVVSRENNARLMFGXG'], \
	[2, 9, 'CVVSVDYKLSFGXG'], [2, 20, 'CVVSRENNARLMFGXG'], [2, 0, 'CVVQTGGGNKLTFGXG'], [2, 31, 'CVSLKGSQGNLIFGXG'], \
	[2, 30, 'CVVSRENNARLMFGXG'], [36, 11, 'CALGHLRSARQLTFGXG'], [36, 16, 'CALYSGAGSYQLTFGXG'], [36, 25, 'CALHGSNNTGKLIFGXG'],\
	[36, 25, 'CALHGSSNTGKLIFGXG'], [36, 6, 'CALEGPGDQKLLFGXG'], [36, 20, 'CALGKLNARLMFGXG'], [36, 15, 'CALITNAGKSTFGXG'], \
	[36, 0, 'CALHGVTGGGNKLTFGXG'], [36, 25, 'CALDDYSGNTGKLIFGXG'], [36, 46, 'CALDLSGGSYIPTFGXG'], [36, 16, 'CALYSGAGSYQLTFGXG'],\
	[36, 25, 'CALHGSGNTGKLIFGXG'], [36, 25, 'CALHGSSNTGKLIFGXG'], [36, 33, 'CALGYTGTASKLTFGXG'], [36, 12, 'CALVNQGGKLIFGXG'],\
	[2, 8, 'CVVRADRGSTLGRLYFGXG'], [2, 8, 'CVVSDRGSTLGRLYFGXG'], [2, 8, 'CVVSDKGSTLGRLYFGXG'], [2, 8, 'CVVSATNRGSTLGRLYFGXG']]

################ add http://onlinelibrary.wiley.com/doi/10.1002/eji.201141956/full ##################

all_apub_lists = [x for x in dir() if "nkt_alpha_" in x] 

all_bpub_lists = [x for x in dir() if "nkt_beta_" in x] 

all_apub_raw = []
all_bpub_raw = []

for l in all_apub_lists:
  for c in vars()[l]:
    all_apub_raw.append(c)
    
for l in all_bpub_lists:
  for c in vars()[l]:
    all_bpub_raw.append(c)

# Prevent duplicate CDR3s making it through

apub = coll.Counter()
bpub = coll.Counter()

for x in all_apub_raw:
  apub[x] += 1
  
for x in all_bpub_raw:
  bpub[x] += 1

uniqa = apub.keys()
uniqb = bpub.keys()

all_apub = ['uniqa']
all_bpub = ['uniqb']




# FUCKED UP!
# Accidentally overwrote this script with an old version, basically wiping the majority of the code.
# Trying to recover the code off my home drive, not optimistic

# Still have the figures, but might have to regenerate the code

# Basically did same as in HIV analysis, just on 1 bleed of all healthies

# had some more sequences though
  # had more nkt from sandserson
  # also had new putative from van schaik (and some more mait)
  
  # then also output the most common from publicness, subtracted those that are present in any of the other groups, then look for those as well

  # see plotting for details












