#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

import sys, random, re, math
from decimal import *



# TODO: for PCR, identification of primers on the edge of a circular sequence



dna_alphabet = {'A':'A', 'C':'C', 'G':'G', 'T':'T',
                'R':'AG', 'Y':'CT', 'W':'AT', 'S':'CG', 'M':'AC', 'K':'GT',
                'H':'ACT', 'B':'CGT', 'V':'ACG', 'D':'AGT',
                'N':'ACGT',
                'a': 'a', 'c': 'c', 'g': 'g', 't': 't',
                'r':'ag', 'y':'ct', 'w':'at', 's':'cg', 'm':'ac', 'k':'gt',
                'h':'act', 'b':'cgt', 'v':'acg', 'd':'agt',
                'n':'acgt'}

complement_alphabet = {'A':'T', 'T':'A', 'C':'G', 'G':'C','R':'Y', 'Y':'R',
                       'W':'W', 'S':'S', 'M':'K', 'K':'M', 'H':'D', 'D':'H',
                       'B':'V', 'V':'B', 'N':'N','a':'t', 'c':'g', 'g':'c',
                       't':'a', 'r':'y', 'y':'r', 'w':'w', 's':'s','m':'k',
                       'k':'m', 'h':'d', 'd':'h', 'b':'v', 'v':'b', 'n':'n'}
gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def translate( sequence ):
    """Return the translated protein from 'sequence' assuming +1 reading frame"""
    return ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])



# Read in all enzymes:
def EnzymeDictionary():
	EnzymeDictionary = {}
	fh = open('REases.tsv', 'rU')
	for line in fh:
		card = line.rstrip().split('\t')
		card[0] = re.sub(r'\-','_',card[0])
		EnzymeDictionary[card[0]] = restrictionEnzyme(card[0],card[1],card[2],card[3],card[4],
													card[5],card[6],card[7],card[8],card[9])
	return EnzymeDictionary

# Suffix Tree implementation from: http://chipsndips.livejournal.com/2005/12/07/
inf = 1000000

# Define a  for a node in the suffix tree
class SuffixNode(dict):
	def __init__(self):
		self.suffixLink = None # Suffix link as defined by Ukkonen

class LCS:
	def __init__(self,str1,str2):
		str = str1 + str2
		inf = len(str)
		self.str = str   #Keep a reference to str to ensure the string is not garbage collected
		self.seed = SuffixNode() #Seed is a dummy node. Suffix link of root points to seed. For any char,there is a link from seed to root
		self.root = SuffixNode() # Root of the suffix tree
		self.root.suffixLink = self.seed
		self.root.depth = 0
		self.deepest = 0,0

		# For each character of str[i], create suffixtree for str[0:i]
		s = self.root; k=0
		for i in range(len(str)):
			self.seed[str[i]] = -2,-2,self.root
			oldr = self.seed
			t = str[i]
			#Traverse the boundary path of the suffix tree for str[0:i-1]
			while True:
				# Descend the suffixtree until state s has a transition for the stringstr[k:i-1]
				while i>k:
					kk,pp,ss = s[str[k]]
					if pp-kk < i-k:
						k = k + pp-kk+1
						s = ss
					else:
						break
			   # Exit this loop if s has a transition for the string str[k:i] (itmeans str[k:i] is repeated);
			   # Otherwise, split the state if necessary
				if i>k:
					tk = str[k]
					kp,pp,sp = s[tk]
					if t.lower() == str[kp+i-k].lower():
						break
					else: # Split the node
						r = SuffixNode()
						j = kp+i-k
						tj = str[j]
						r[tj] = j, pp, sp
						s[str[kp]] = kp,j-1, r
						r.depth = s.depth + (i-k)
						sp.depth = r.depth + pp - j + 1
						if j<len(str1)<i and r.depth>self.deepest[0]:
							self.deepest = r.depth,j-1
				elif s.has_key(t):
					break
				else:
					r = s
			   # Add a transition from r that starts with the letter str[i]
				tmp = SuffixNode() 
				r[t] = i,inf,tmp
				# Prepare for next iteration
				oldr.suffixLink = r
				oldr = r
				s = s.suffixLink
			# Last remaining endcase
			oldr.suffixLink = s

	def LongestCommonSubstring(self):
		start = self.deepest[1]-self.deepest[0]+1
		end = self.deepest[1]+1
		return (self.str[start:end],start,end)

class PrimerError(Exception):
    """Exception raised for errors in the primer(s) input.

    Attributes:
        primer -- sequence for one (or both, in tuple form) of the given input primers
        template -- sequence for the given PCR template
        msg  -- explanation of the error
    """

    def __init__(self, primer, template, msg):
        self.primer = primer
        self.template = template
        self.msg = msg

#Note: PCR() product is not case preserving
def PCR(primer1DNA, primer2DNA, templateDNA):
	template = templateDNA.sequence + '$'
	primer_1 = primer1DNA.sequence + '$'
	primer_2 = primer2DNA.sequence + '$'
	try:
		indices = [0,0,0,0,0,0]		# List for storing primer annealing region start/stop indices and strand association
		counter = 0
		rightStub = ''				# "stub" = non-annealing regions of the input primers. Want to append these 
		leftStub = ''					# to the right (5') and left (3') ends of the output PCR product.
		nextOrientation = 0
		for currentPrimer in (primer_1, primer_2):		# not making any assumptions about the directionality of the input primers
			fwdMatch = LCS(template.upper(),currentPrimer.upper())
			fwdTuple = fwdMatch.LongestCommonSubstring()
			first = re.compile(fwdTuple[0], re.IGNORECASE)
			fwd_stub = currentPrimer[0:len(currentPrimer)-len(fwdTuple[0])-1]
			fList = first.findall(template)
			matchCount = 0
			matchedAlready = 0
			start = 0
			stop = 0
			for match in fList:
				if primerTm(match) >= 45:						# switched this to Tm >= 45 C for matches
					matchCount = matchCount + 1
			if matchCount == 0 & len(fList) > 0:				# no matches in forward direction
				tooShort1 = True
			else:
				tooShort1 = False
			if matchCount > 1:									# if matches in forward direction more than once
				if nextOrientation == 2: 							# ... but was supposed to match in reverse direction
					raise PrimerError(currentPrimer,template,'Primers have same forward (5\'->3\') orientation AND primer anneals to multiple sites in template:')
				raise PrimerError(currentPrimer,template,'Primer anneals to multiple sites in template:')
			elif matchCount == 1:								# if matches in the forward direction exactly once
				if nextOrientation == 2: 							# ... but was supposed to match in reverse direction
					raise PrimerError(currentPrimer,template,'Primers have same forward (5\'->3\') orientation')
				matchedAlready = 1
			revcomp = reverseComplement(currentPrimer)
			revMatch = LCS(template.upper(),revcomp.upper()+'$')
			revTuple = revMatch.LongestCommonSubstring()
			last = re.compile(revTuple[0], re.IGNORECASE)
			rev_stub = currentPrimer[len(revTuple[0]):len(currentPrimer)]
			lList = last.findall(template)
			matchCount = 0
			for match in lList:
				if primerTm(match) >= 45:						# switched this to Tm >= 45 C for matches
					matchCount = matchCount + 1
			if matchCount == 0 & len(lList) > 0:				# no matches in forward direction
				tooShort2 = True
			else:
				tooShort2 = False
			if matchCount > 1:									# if matches in reverse direction more than once
				if matchedAlready == 1:								# ... and already matched in forward direction
					if nextOrientation == 1: 							# ... but was supposed to match in forward direction
						raise PrimerError(currentPrimer,template,'Primers have same reverse (3\'->5\') orientation AND primer anneals to multiple sites in template AND it primes in both orientations:')
					raise PrimerError(currentPrimer,template,'Primer anneals to multiple sites in template AND it primes in both orientations:')
				if nextOrientation == 1: 
					raise PrimerError(currentPrimer,template,'Primers have same reverse (3\'->5\') orientation AND primer anneals to multiple sites in template:')
				raise PrimerError(currentPrimer,template,'Primer anneals to multiple sites in template:')
			elif matchCount == 1: 								# if matches in the reverse direction exactly once
				if matchedAlready == 1:								# ... and already matched in forward direction
					if nextOrientation == 1: 							# ... but was supposed to match in forward direction
						raise PrimerError(currentPrimer,template,'Primers have same reverse (3\'->5\') orientation AND primer primes in both orientations:')
					raise PrimerError(currentPrimer,template,'Primer primes in both orientations:')
				else:
					matchedAlready = 2
			if matchedAlready == 0:								# if no matches
				if tooShort1 and tooShort2:							# ... it may be because the annealing region has a tm < 45 C
					raise PrimerError(currentPrimer, template,'Primer is too short and does not anneal')
				raise PrimerError(currentPrimer,template,'Primer does not prime in either orientation:') 	# ... or not.
			if matchedAlready == 1:
				indices[counter] = fwdTuple[1]
				counter = counter + 1
				indices[counter] = fwdTuple[2]
				counter = counter + 1
				indices[counter] = 'fwd'
				counter = counter + 1
				nextOrientation = 2
				leftStub = fwd_stub
			if matchedAlready == 2:
				indices[counter] = revTuple[1]
				counter = counter + 1
				indices[counter] = revTuple[2]
				counter = counter + 1
				indices[counter] = 'rev'
				counter = counter + 1
				nextOrientation = 1
				rightStub = reverseComplement(rev_stub)
		if indices[2] == 'fwd':
			fwdStart = indices[0]
			fwdEnd = indices[1]
			revStart = indices[3]
			revEnd = indices[4]
			if fwdStart < revStart and fwdEnd < revEnd:
				return DNA(leftStub+template[fwdStart:revEnd]+rightStub,'PCR product')
			else:
				if templateDNA.topology == 'circular':	# circular template is exception to the fwdStart < revStart and fwdEnd < revEnd rule
					return DNA(leftStub+template[fwdStart:len(template)-1]+template[:revStart]+rightStub,'PCR product')
				else:
					raise PrimerError((primer1DNA.sequence, primer2DNA.sequence),template,'Forward primer beginning and ending indices must be before those of the reverse:')
		elif indices[2] == 'rev':
			fwdStart = indices[3]
			fwdEnd = indices[4]
			revStart = indices[0]
			revEnd = indices[1]
			if fwdStart < revStart and fwdEnd < revEnd:
				return DNA(leftStub+template[fwdStart:revEnd]+rightStub,'PCR product')
			else:
				if templateDNA.topology == 'circular':
					return DNA(leftStub+template[fwdStart:len(template)-1]+template[:revStart]+rightStub,'PCR product')
				else:
					raise PrimerError((primer1DNA.sequence, primer2DNA.sequence),template,'Forward primer beginning and ending indices must be before those of the reverse:')
	except PrimerError as error:
		print 'EXCEPTION: '+ error.msg
		print 'primer: ' 
		print error.primer
		#print 'template: '
		#print error.template
		sys.exit()

# Note: reverseComplement() is case preserving
def reverseComplement(sequence):
	basecomplement = {'G':'C', 'A':'T', 'T':'A', 'C':'G', 'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 'S':'S', 'W':'W', 'H':'D', 'B':'V', 'V':'B', 'D':'H', 'N':'N','g':'c', 'a':'t', 't':'a', 'c':'g', 'r':'y', 'y':'r', 'm':'k', 'k':'m', 's':'s', 'w':'w', 'h':'d', 'b':'v', 'v':'b', 'd':'h', 'n':'n'}
  	return "".join([basecomplement.get(nucleotide.lower(), '') for nucleotide in sequence[::-1]])

# Note: reverseComplement() is case preserving
def Complement(sequence):
	basecomplement = {'G':'C', 'A':'T', 'T':'A', 'C':'G', 'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 'S':'S', 'W':'W', 'H':'D', 'B':'V', 'V':'B', 'D':'H', 'N':'N','g':'c', 'a':'t', 't':'a', 'c':'g', 'r':'y', 'y':'r', 'm':'k', 'k':'m', 's':'s', 'w':'w', 'h':'d', 'b':'v', 'v':'b', 'd':'h', 'n':'n'}
  	return "".join([basecomplement.get(nucleotide.lower(), '') for nucleotide in sequence[0:]])

def primerTm(sequence):
	milliMolarSalt = 50
	milliMolarMagnesium = 1.5
	nanoMolarPrimerTotal = 200
	molarSalt = milliMolarSalt/1000
	molarMagnesium = milliMolarMagnesium/1000
	molarPrimerTotal = Decimal(nanoMolarPrimerTotal)/Decimal(1000000000)
	re.sub(r'\s','', sequence)
 	return nearestNeighborTmNonDegen(sequence, molarSalt, molarPrimerTotal, molarMagnesium)

def primerTmsimple(sequence):
  	return 64.9+41*(GCcontent(sequence)*len(sequence) - 16.4)/len(sequence)

# phusion notes on Tm
# https://www.finnzymes.fi/optimizing_tm_and_annealing.html

# get substring from the beginning of input that is 55C Tm
def get_55_primer(sequence):
	lastChar = 17
	myPrimer = sequence.substring(0,lastChar)
	while( primerTmsimple(myPrimer) < 54.5 or lastChar > 60):
		lastChar = lastChar + 1
		myPrimer = sequence[0:lastChar]
	return myPrimer

def nearestNeighborTmNonDegen (sequence, molarSalt, molarPrimerTotal, molarMagnesium):
	# The most sophisticated Tm calculations take into account the exact sequence and base stacking parameters, not just the base composition.
	# m = ((1000* dh)/(ds+(R * Math.log(primer concentration))))-273.15;
	# Borer P.N. et al. (1974)  J. Mol. Biol. 86, 843.
	# SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
	# Allawi, H.T. and SantaLucia, J. Jr. (1997) Biochemistry 36, 10581.
	# von Ahsen N. et al. (1999) Clin. Chem. 45, 2094.

	sequence = sequence.lower()

	R = 1.987 # universal gas constant in Cal/degrees C * mol
	ds = 0 	  # cal/Kelvin/mol
	dh = 0    # kcal/mol

	# perform salt correction
	correctedSalt = molarSalt + molarMagnesium * 140 # adjust for greater stabilizing effects of Mg compared to Na or K. See von Ahsen et al 1999
	ds = ds + 0.368 * (len(sequence) - 1) * math.log(correctedSalt) # from von Ahsen et al 1999

	#  perform terminal corrections
	termDsCorr = getTerminalCorrectionsDsHash()
	ds = ds + termDsCorr[sequence[0]]
	ds = ds + termDsCorr[sequence[len(sequence) - 1]]

	termDhCorr = getTerminalCorrectionsDhHash()
	dh = dh + termDhCorr[sequence[0]]
	dh = dh + termDhCorr[sequence[len(sequence) - 1]]

	dsValues = getDsHash()
	dhValues = getDhHash()

	for i in range(len(sequence)-1):
		ds = ds + dsValues[sequence[i] + sequence[i + 1]]
		dh = dh + dhValues[sequence[i] + sequence[i + 1]]
	return (((1000 * dh) / (ds + (R * math.log(molarPrimerTotal / 2)))) - 273.15)

def getTerminalCorrectionsDsHash():
	# SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
	dictionary = {'g' : -2.8,'a': 4.1,'t' : 4.1,'c' : -2.8}
	return dictionary

def getTerminalCorrectionsDhHash():
	# SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
	dictionary = {'g':0.1,'a' : 2.3,'t' : 2.3,'c' : 0.1}
	return dictionary

def getDsHash():
	# SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
	dictionary = {
	'gg' : -19.9,
	'ga' : -22.2,
	'gt' : -22.4,
	'gc' : -27.2,
	'ag' : -21.0,
	'aa' : -22.2,
	'at' : -20.4,
	'ac' : -22.4,
	'tg' : -22.7,
	'ta' : -21.3,
	'tt' : -22.2,
	'tc' : -22.2,
	'cg' : -27.2,
	'ca' : -22.7,
	'ct' : -21.0,
	'cc' : -19.9}
	return dictionary

def getDhHash():
	# SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
	dictionary = {'gg': -8.0,
	'ga' : -8.2,
	'gt' : -8.4,
	'gc' : -10.6,
	'ag' : -7.8,
	'aa' : -7.9,
	'at' : -7.2,
	'ac' : -8.4,
	'tg' : -8.5,
	'ta' : -7.2,
	'tt' : -7.9,
	'tc' : -8.2,
	'cg' : -10.6,
	'ca' : -8.5,
	'ct' : -7.8,
	'cc' : -8.0}
	return dictionary

def Digest(InputDNA, Enzymes):
	# TODO: Error/Exception handling
		# fixed: DNA('TGGGGACTGCCGTTCATGGTGAGATGAGTGAAGGCGAGCTGGTGGATGCATTCCGCCATGTGAGTGATGCGTTTGAGCAAACCAGCGAAACCATCGGCGTGCGCGCCAATAACGCGATCAACGACATGGTGCGTCAACGTCTGCTGAACCGCTTTACCAGCGAGCAGGCGGAAGGGAACGCAATTTACCGTCTGACGCCGCTCGGCATCGGCATTACTGACTACNNNATCCGTCAGCGCGAGTTTTCTACGCTGCGTCTTTCTATGCAGTTGTCGATTGTGGCGGGTGAGCTCAAACGCGCAGCAGATGCCGCCGAAGAGGGCGGTGATGAATTTCACTGGCACCGTAATGTCTATGCGCCACTGAAATATTCGGTAGCAGAAATTTTCGACAGTATCGACCTGACGCAACGTCTGATGGACGAACAGCAGCAGCAGGTGAAGGACGATATCGCCCAGTTGCTGAACAAAGACTGGCGGGCGGCGATTTCCAGCTGGATCCTGAATTGTTGCTTTCGGAAACTTCCGGAACGCTGCGTGAATTGCAGGATACGCTGGAAGCGGCAGGCGACAAATTGCAGGCTAATCTGTTGCGCATTCAGGATGCGACGATGACCCATGACGATCTGCATTTCGTCGATCGTCTGGTGTTCGATCTGCAGAGCAAACTCGATCGTATTATCAGTTGGGGCCAGCAATCCATCGACTTGTGGATTGGCTACGACCGCCACGTACACAAATTTATTCGTACCGCGATCGATATGGATAAAAACCGCGTCTTTGCTCAGCGGTTACGTCAGTCGGTACAAACCTATTTTGATGAACGGCGGGCGCTAACTTATGCCAATGCCGATCGTCTGCTGGATATGCGTGACGAAGAGATGGCACTGCGCGATGAAGAAGTGACTGGGGAACTTCCTGAGGATCTGGAATACGAAGAGTTTAACGAGATCCGCGAACAGCTGGCGGCGATCATCGAAGAACAACTTGCCGTGTACAAAACCAGACAAGTGCCGCTGGATCTTGGTCTGGTGGTACGCGAATATCTGTCACAGTATCCGCGTGCACGTCACTTTGACGTTGCGCGTATTGTTATTGATACCTGACGCAACGTCTGCGAATTCCTGCAGTA','plasmid')
		# fixed: DNA('TGACTGCCGTTCATGGTGAGATGAGTGAAGGCGAGCTGGTGGATGCATTCCGCCATGTGAGTGATGCGTTTGAGCAAACCAGCGAAACCATCGGCGTGCGCGCCAATAACGCGATCAACGACATGGTGCGTCAACGTCTGCTGAACCGCTTTACCAGCGAGCAGGCGGAAGGGAACGCAATTTACCGTCTGACGCCGCTCGGCATCGGCATTACTGACTACNNNATCCGTCAGCGCGAGTTTTCTACGCTGCGTCTTTCTATGCAGTTGTCGATTGTGGCGGGTGAGCTCAAACGCGCAGCAGATGCCGCCGAAGAGGGCGGTGATGAATTTCACTGGCACCGTAATGTCTATGCGCCACTGAAATATTCGGTAGCAGAAATTTTCGACAGTATCGACCTGACGCAACGTCTGATGGACGAACAGCAGCAGCAGGTGAAGGACGATATCGCCCAGTTGCTGAACAAAGACTGGCGGGCGGCGATTTCCAGCTGGATCCTGAATTGTTGCTTTCGGAAACTTCCGGAACGCTGCGTGAATTGCAGGATACGCTGGAAGCGGCAGGCGACAAATTGCAGGCTAATCTGTTGCGCATTCAGGATGCGACGATGACCCATGACGATCTGCATTTCGTCGATCGTCTGGTGTTCGATCTGCAGAGCAAACTCGATCGTATTATCAGTTGGGGCCAGCAATCCATCGACTTGTGGATTGGCTACGACCGCCACGTACACAAATTTATTCGTACCGCGATCGATATGGATAAAAACCGCGTCTTTGCTCAGCGGTTACGTCAGTCGGTACAAACCTATTTTGATGAACGGCGGGCGCTAACTTATGCCAATGCCGATCGTCTGCTGGATATGCGTGACGAAGAGATGGCACTGCGCGATGAAGAAGTGACTGGGGAACTTCCTGAGGATCTGGAATACGAAGAGTTTAACGAGATCCGCGAACAGCTGGCGGCGATCATCGAAGAACAACTTGCCGTGTACAAAACCAGACAAGTGCCGCTGGATCTTGGTCTGGTGGTACGCGAATATCTGTCACAGTATCCGCGTGCACGTCACTTTGACGTTGCGCGTATTGTTATTGATACCTGACGCAACGTCTGCGAATTCCTGCAGTA','PCR product')
		# fixed: DNA('GATGACTGAATTCTCATGGTGAGATGAGTGAAGGCGAGCTGGTGGATGCATTCCGCCATGTGAGTGATGCGTTTGAGCAAACCAGCGAAACCATCGGCGTGCGCGCCAATAACGCGATCAACGACATGGTGCGTCAACGTCTGCTGAACCGCTTTACCAGCGAGCAGGCGGAAGGGAACGCAATTTACCGTCTGACGCCGCTCGGCATCGGCATTACTGACTACNNNATCCGTCAGCGCGAGTTTTCTACGCTGCGTCTTTCTATGCAGTTGTCGATTGTGGCGGGTGAGCTCAAACGCGCAGCAGATGCCGCCGAAGAGGGCGGTGATGAATTTCACTGGCACCGTAATGTCTATGCGCCACTGAAATATTCGGTAGCAGAAATTTTCGACAGTATCGACCTGACGCAACGTCTGATGGACGAACAGCAGCAGCAGGTGAAGGACGATATCGCCCAGTTGCTGAACAAAGACTGGCGGGCGGCGATTTCCAGCTGGATCCTGAATTGTTGCTTTCGGAAACTTCCGGAACGCTGCGTGAATTGCAGGATACGCTGGAAGCGGCAGGCGACAAATTGCAGGCTAATCTGTTGCGCATTCAGGATGCGACGATGACCCATGACGATCTGCATTTCGTCGATCGTCTGGTGTTCGATCTGCAGAGCAAACTCGATCGTATTATCAGTTGGGGCCAGCAATCCATCGACTTGTGGATTGGCTACGACCGCCACGTACACAAATTTATTCGTACCGCGATCGATATGGATAAAAACCGCGTCTTTGCTCAGCGGTTACGTCAGTCGGTACAAACCTATTTTGATGAACGGCGGGCGCTAACTTATGCCAATGCCGATCGTCTGCTGGATATGCGTGACGAAGAGATGGCACTGCGCGATGAAGAAGTGACTGGGGAACTTCCTGAGGATCTGGAATACGAAGAGTTTAACGAGATCCGCGAACAGCTGGCGGCGATCATCGAAGAACAACTTGCCGTGTACAAAACCAGACAAGTGCCGCTGGATCTTGGTCTGGTGGTACGCGAATATCTGTCACAGTATCCGCGTGCACGTCACTTTGACGTTGCGCGTATTGTTATTGATACCTGACGCAACGTCTGCGAATTCCTGCAGTA', 'plasmid')
	print 'Input DNA ('+InputDNA.topology+'): '+InputDNA.sequence
	(indices, frags, sites, totalLength) = ([], [], "", len(InputDNA.sequence)) # Initialization
	if InputDNA.topology == "linear":	
		# Initialize indices array with start and end indices of the linear fragment
			# Add dummy REase to avoid null pointers
		dummy = restrictionEnzyme("dummy", "", "", "", "", "", 0, 0, "(0/0)","")
		indices = [(0,0,'',dummy), (totalLength,0,'',dummy)]
	# Identify restriction sites, fill in indices array
	for enzyme in Enzymes:
		sites = enzyme.find_sites(InputDNA)
		for site in sites:
			# WARNING: end proximity for linear fragments exception
			if InputDNA.topology == 'linear' and int(site[0]) - int(enzyme.endDistance) < 0 or int(site[1]) + int(enzyme.endDistance) > totalLength:
				print 'WARNING: end proximity for '+enzyme.name+' restriction site at indices '+str(site[0]%totalLength)+','+str(site[1]%totalLength)+' for input with length '+str(totalLength)
				if InputDNA.topology == 'linear' and site[2] == 'antisense' and site[1] - max(enzyme.bottom_strand_offset,enzyme.top_strand_offset) < 0:
					print 'WARNING: restriction cut site for '+enzyme.name+' with recognition indices '+str(site[0]%totalLength)+','+str(site[1]%totalLength)+' out of bounds for input with length '+str(totalLength)
				else:
					pass
			# WARNING: restriction index out of bounds exception
			elif InputDNA.topology == 'linear' and site[2] == 'antisense' and site[1] - max(enzyme.bottom_strand_offset,enzyme.top_strand_offset) < 0:
				print 'WARNING: restriction cut site for '+enzyme.name+' with recognition indices '+str(site[0]%totalLength)+','+str(site[1]%totalLength)+' out of bounds for input with length '+str(totalLength)
				pass
			else: 
				site = site + (enzyme, )
				indices.append(site)
		indices.sort()
	# If you have overlapping restriction sites, choose the first one and discard they
		# second (TODO: there may be a better, non-greedy way to do this... not sure)
	filtered = []
	n = 0
	while n < len(indices):
		try:
			(currentTuple, nextTuple) = (indices[n], indices[n+1])
			(currentStart, nextStart) = (currentTuple[0], nextTuple[0])
			filtered.append(indices[n])
			if currentStart + len(enzyme.alpha_only_site) >= nextStart:
				currentIndex = indices[n+1]
				if currentIndex[0] == len(InputDNA.sequence):
					pass
				else:
					print 'WARNING: overlapping restriction sites '+currentTuple[3].name+' (indices '+str(currentTuple[0])+','+str(currentTuple[1])+') and '+nextTuple[3].name+' (indices '+str(nextTuple[0])+','+str(nextTuple[1])+')'
					n = n + 1
			n = n + 1
		except:
			# got to end of list,
			filtered.append(indices[n])
			n = n + 1
	indices = filtered
	# If it's linear, only act on the first n - 1 fragments until you hit the blunt ending
		# If it's circular, then the 'last' segment is adjacent to the 'first' one, so you
		# need to consider the adjacency relationships among the full n fragments
	if InputDNA.topology == "linear":
		lastIt = len(indices) - 1
	else:
		lastIt = len(indices)
	# Consider enzyme for the current restriction site as well as the next restriction
		# site, so that you can generate overhangs for both sides of the current fragment
	for n in range(lastIt):
		currentTuple = indices[n]
		if n+1 > len(indices) - 1:
			n = -1
		nextTuple = indices[n+1]
		(currentStart, currentEnd, direction, currentEnzyme) = currentTuple
		(nextStart, nextEnd, nextDirection, nextEnzyme) = nextTuple
		# If it's on the sense strand, then overhang is positive
			# If it's on the antisense strand, then you have to go back towards the 5'
			# to generate the overhang (so multiply by -1)
		# CT(B)O = current top (bottom) overhang, AL(R)L = add left (right) length, NT(B)O = next top (bottom) overhang
		(ALL, ARL) = (0,0)
		if direction == "sense":
			(CTO, CBO) = (currentEnzyme.top_strand_offset, currentEnzyme.bottom_strand_offset)
			ALL = max(CTO,CBO)
		else:
			(CTO, CBO) = (-1 * currentEnzyme.top_strand_offset, -1 * currentEnzyme.bottom_strand_offset)
			ALL = max(CTO,CBO)
		if nextDirection == "sense":
			(NTO, NBO) = (nextEnzyme.top_strand_offset, nextEnzyme.bottom_strand_offset)
			ARL = min(NTO,NBO)
		else:
			(NTO, NBO) = (-1 * nextEnzyme.top_strand_offset + 1, -1 * nextEnzyme.bottom_strand_offset + 1)
			ARL = min(NTO,NBO)-1
		# Update start value currentStart and apply ( mod length ) to deal with edge cases
			# Also, update end value digEnd for fragment indices
		currentStart = currentStart+ALL
		currentStart = currentStart % totalLength
		digEnd = nextStart + ARL
		if currentEnzyme.reach and direction == "sense":
			currentStart = currentStart + len(currentEnzyme.alpha_only_site)
		if nextEnzyme.reach and nextDirection == "sense":
			digEnd = digEnd + len(nextEnzyme.alpha_only_site)
		# Loop around fragment case for circular InputDNA's
		if digEnd > 0 and currentStart > 0 and digEnd < currentStart and InputDNA.topology == 'circular':
			digested = DNA(InputDNA.sequence[currentStart:]+InputDNA.sequence[:digEnd],'digest')
		else:
			digested = DNA(InputDNA.sequence[currentStart:digEnd],'digest')
		# Adjust top and bottom overhang values based on the orientation of the restriction site
		if direction == "sense":
			(TO, BO) = (CTO, CBO)
		else:
			(TO, BO) = (CBO, CTO)
		difference = abs(abs(BO) - abs(TO))
		# Generate TLO and BLO fragment overhangs
		if abs(TO) < abs(BO) and direction == "sense" or abs(TO) > abs(BO) and direction == "antisense":
			if currentStart - len(currentEnzyme.alpha_only_site) < 0:
				digested.topLeftOverhang = Overhang(InputDNA.sequence[currentStart-difference:]+InputDNA.sequence[:currentStart])
			else:
				digested.topLeftOverhang = Overhang(InputDNA.sequence[currentStart-difference:currentStart])
			digested.bottomLeftOverhang = Overhang('')
		else:
			digested.topLeftOverhang = Overhang('') 
			# Edge case statement
			if currentStart - len(currentEnzyme.alpha_only_site) < 0:
				digested.bottomLeftOverhang = Overhang(Complement(InputDNA.sequence[currentStart-difference:]+InputDNA.sequence[:currentStart]))
			else:
				digested.bottomLeftOverhang = Overhang(Complement(InputDNA.sequence[currentStart-difference:currentStart]))
		print 'Fragment: '+digested.sequence
		print 'Fragment.TLO: '+ digested.topLeftOverhang.sequence + " ("+str(currentStart-difference)+","+str(currentStart)+")"
		print 'Fragment.BLO: '+ digested.bottomLeftOverhang.sequence + " ("+str(currentStart-difference)+","+str(currentStart)+")"
		# Adjust top and bottom overhang values based on the orientation of the restriction site
		if direction == "sense":
			(TO, BO) = (NTO, NBO)
		else:
			(TO, BO) = (NBO, NTO)
		difference = abs(abs(BO) - abs(TO))
		# Apply ( mod length ) operator to end index value digDiff to deal with edge cases
		digDiff = digEnd + difference
		digDiff = digDiff % totalLength
		# Generate TRO and BRO fragment overhangs
		if abs(TO) < abs(BO) and direction == "sense" or abs(TO) > abs(BO) and direction == "antisense":
			digested.topRightOverhang = Overhang('')
			# Edge case statement
			if digDiff - len(nextEnzyme.alpha_only_site) < 0:
				digested.bottomRightOverhang = Overhang(Complement(InputDNA.sequence[digEnd:]+InputDNA.sequence[:digDiff]))
			else:
				digested.bottomRightOverhang = Overhang(Complement(InputDNA.sequence[digEnd:digDiff])) 
		else:
			# Edge case statement
			if digDiff - len(nextEnzyme.alpha_only_site) < 0:
				digested.topRightOverhang = Overhang(InputDNA.sequence[digEnd:]+InputDNA.sequence[:digDiff])
			else:
				digested.topRightOverhang = Overhang(InputDNA.sequence[digEnd:digDiff])
			digested.bottomRightOverhang = Overhang('')
		print 'Fragment.TRO: '+ digested.topRightOverhang.sequence + " ("+str(digEnd)+","+str(digDiff)+")"
		print 'Fragment.BRO: '+ digested.bottomRightOverhang.sequence + " ("+str(digEnd)+","+str(digDiff)+")"
		# Discard small fragments
		# TODO: what is the right length for this? For a zymo we will discard all small frags, but without Zymo not clear what is best
		if len(digested.sequence) < 4:
			pass
		else:
			# frags.append((currentStart,digested))
			frags.append(digested)
		# frags.sort()
	return frags

def revcomp(string):
       letters = list(string)
       letters = [complement_alphabet[base] for base in letters]
       rcomp = ''.join(letters)
       return rcomp[::-1] #reverses string

class Overhang(object):
	def __init__(self, seq=""):
		self.sequence = seq

class DNA(object):
	#for linear DNAs, this string should include the entire sequence (5' and 3' overhangs included
	def __init__(self, seq="",DNAclass=""):
		self.sequence = seq
		self.length = len(seq)
		notDNA = re.compile('([^gatcrymkswhbvdn])')
		isnotDNA = False
		exceptionText = "" 
		for m in notDNA.finditer(self.sequence.lower()):
			exceptionText = exceptionText + m.group()+ " at position "+ str( m.start()) + " is not valid IUPAC DNA; "
			isnotDNA = True
		if(isnotDNA):
			raise Exception(exceptionText)
		self.name = "pbca1256" #would be pbca1256 for vectors or pbca1256-Bth8199 for plasmids
		self.description = "SpecR pUC" #this is for humans to read
		self.dam_methylated = True
		self.overhang = "circular" #blunt, 3', 5', circular... should be a class in itself?
		self.topLeftOverhang = ""
		self.bottomLeftOverhang = ""
		self.topRightOverhang = ""
		self.bottomRightOverhang = ""
		#PCR product, miniprep, genomic DNA
		self.provenance = ""
		#Here is the linked list references for building up action-chains
		# an action chain would be something like do PCR on day 1, do transformation on day 2, etc
		self.head = None
		self.tail = None
		if DNAclass == "primer" or DNAclass == "genomic" or DNAclass == "PCR product" or DNAclass == "digest":
			self.topology = "linear"
		elif DNAclass == 'plasmid':
			self.topology = "circular" #circular or linear, genomic should be considered linear
		else:
			raise Exception("Invalid molecule class. Acceptable classes are 'digest', genomic', 'PCR product', 'plasmid' and 'primer'.")
	def reversecomp(self):
		return revcomp(self.sequence) #reverses string
		#code to handle the overhangs & other object attributes
	def find(self, string):
		return 0
	def prettyPrint(self):
		#prints out top and bottom strands, truncates middle so length is ~100bp
		#example:
		# TTATCG...[1034bp]...GGAA
		#   ||||              ||||
		#   TAGC..............CCTTAA
		return 0
	
#taken from BioPython
def BaseExpand(base):
    """BaseExpand(base) -> string.

    given a degenerated base, returns its meaning in IUPAC alphabet.

    i.e:
        b= 'A' -> 'A'
        b= 'N' -> 'ACGT'
        etc..."""
    base = base.upper()
    return dna_alphabet[base]

#function to convert recog site into regex, from Biopython
def regex(site):
    """regex(site) -> string.

    Construct a regular expression from a DNA sequence.
    i.e.:
        site = 'ABCGN'   -> 'A[CGT]CG.'"""
    reg_ex = site
    for base in reg_ex:
        if base in ('A', 'T', 'C', 'G', 'a', 'c', 'g', 't'):
            pass
        if base in ('N', 'n'):
            reg_ex = '.'.join(reg_ex.split('N'))
            reg_ex = '.'.join(reg_ex.split('n'))
        if base in ('R', 'Y', 'W', 'M', 'S', 'K', 'H', 'D', 'B', 'V'):
            expand = '['+ str(BaseExpand(base))+']'
            reg_ex = expand.join(reg_ex.split(base))
    return reg_ex

def ToRegex(site, name):
	sense = ''.join(['(?P<', name, '>', regex(site.upper()), ')'])
	antisense = ''.join(['(?P<', name, '_as>', regex(revcomp( site.upper() )), ')'])
	rg = sense + '|' + antisense
	return rg	

class restrictionEnzyme(object):
	def __init__(self,name="", buffer1="", buffer2="", buffer3="", buffer4="", bufferecori="", heatinact="", incubatetemp="", recognitionsite="",distance=""):
		self.name = name
		self.buffer_activity =[buffer1, buffer2, buffer3, buffer4, bufferecori]
		self.inactivate_temp = heatinact
		self.incubate_temp = incubatetemp
		#human-readable recognition site
		self.recognition_site = recognitionsite
		self.endDistance = distance
		#function to convert recog site into regex
		alpha_only_site = re.sub('[^a-zA-Z]+', '', recognitionsite)
		self.alpha_only_site = alpha_only_site
		# print ToRegex(alpha_only_site, name)
		self.compsite = ToRegex(alpha_only_site, name)
		self.reach = False
		#convert information about where the restriction happens to an offset on the top and bottom strand
		#for example, BamHI -> 1/5 with respect to the start of the site match
		hasNum = re.compile('(-?\d+/-?\d+)')
		not_completed = 1
		for m in hasNum.finditer(recognitionsite):
			(top, bottom) = m.group().split('/')
		  	self.top_strand_offset = int(top)
		  	self.bottom_strand_offset = int(bottom)
		  	self.reach = True
		  	not_completed = 0
		p = re.compile("/")
		for m in p.finditer(recognitionsite):
			if not_completed:
				self.top_strand_offset = int(m.start())
				self.bottom_strand_offset =  len(recognitionsite) - 1 - self.top_strand_offset	

	def prettyPrint(self):
		print "Name: ", self.name, "Recognition Site: ", self.recognition_site
	def find_sites(self, DNA):
		seq = DNA.sequence
		(fwd, rev) = self.compsite.split('|')
		fwd_rease_re = re.compile(fwd)
		rev_rease_re = re.compile(rev)
		indices = []
		seen = {}
		if DNA.topology == "circular":
			searchSequence = seq.upper() + seq[0:len(self.recognition_site)-2]
		else:
			searchSequence = seq.upper()
		for m in fwd_rease_re.finditer(searchSequence):
			span = m.span()
			span = (span[0] % len(seq), span[1] % len(seq))
			seen[span[0]] = 1
			span = span + ('sense',)
			indices.append(span)
		for m in rev_rease_re.finditer(searchSequence):
			span = m.span()
			try:
				seen[span[0]]
			except:
				span = span + ('antisense',)
				indices.append(span)	
		return indices

#accepts two primers and list of input template DNAs
def SOE(primer1, primer2, templates):
	return 0

def Distinguish2DNABands(a, b):
        #case of 2
	#for a standard 1-2% agarose gel,
        #we can distinguish a and b if
        return ( abs(a.length - b.length) > (0.208*a.length+42))  & (min(a.length, b.length) > 250 )

def DistinguishDNABands(list_of_dnas):
	ret_val = True
	for i in range(len(list_of_dnas)-1):
		ret_val = ret_val & Distinguish2DNABands(list_of_dnas[i], list_of_dnas[i+1])
	return ret_val

def FindDistinguishingEnzyme(list_of_dnas):
	#find the REase that can distinguish between the input DNAs
	#DistinguishDNABands(a, b) returns true if we can
	# tell apart bands a, b on a gel and a and b are both > 300bp, < 7kb
	#Let n be the number of DNAs in the list.  Let E be the enzyme under question
	# Then we construct a n-dimensional matrix
	# where the dimensions have max value defined by the number of fragments generated by E
	# E can be used to distinguish between the DNAs if there is a complete row or column 
	# that is distinguishable (all True by DistinguishDNABands)
	
	#iterate over all enzymes (enzyme list should be prioritized by availability and "goodness")
		#execute find good enz
	#iterate over all combinations of 2 enzymes
		#execute find good enz

	##find good enz
	#for each enzyme/combo in the list
		#calculate fragments for each input DNA
		#skip if any DNA has # fragments > 6
		#n-length list, each character represents the DNA fragment currently under investigation
		#iterate to fill in the hypermatrix values

		#find if the hypermatrix has a column/row that has all True
	return 0
def FindDistEnz():
	return FindDistinguishingEnzyme(list_of_dnas)

#accepts list of DNA, outputs list of DNA
def Ligate(inputDNAs):
	products = []
	# self ligation
	for fragment in inputDNAs:
		if fragment.topology == 'circular':
			print 'Invalid input molecule removed -- only linear digest fragments accepted.'
			inputDNAs.remove(fragment)
		if fragment.topLeftOverhang.sequence != '':
			if fragment.topLeftOverhang.sequence.lower() == Complement(fragment.bottomRightOverhang.sequence.lower()):
				products.append(DNA(fragment.topLeftOverhang.sequence+fragment.sequence,'plasmid'))
		elif fragment.bottomLeftOverhang.sequence != '':
			if fragment.topLeftOverhang.sequence.lower() == Complement(fragment.topRightOverhang.sequence.lower()):
				products.append(DNA(fragment.sequence+fragment.topRightOverhang.sequence,'plasmid'))
	# TODO: pairwise ligation
	# for fragOne in inputDNAs:
	# 	for fragTwo in inputDNAs:
	# There is a better way to this for sure...
	# 		if fragOne.topLeftOverhang.sequence != '' and fragOne.bottomRightOverhang.sequence  != '' and fragTwo.topLeftOverhang.sequence  != '' and fragTwo.bottomRightOverhang.sequence != '':
	# 			if fragOne.topLeftOverhang.sequence.lower() == Complement(fragTwo.bottomRightOverhang.sequence.lower()) and fragTwo.topLeftOverhang.sequence.lower() == Complement(fragOne.bottomRightOverhang.sequence.lower()):
	# 				products.append(DNA(fragOne.topLeftOverhang.sequence+fragOne.sequence+fragTwo.topLeftOverhang.sequence+fragTwo.sequence,'plasmid'))
	return products


###checks for presence of regex-encoded feature in seq
def HasFeature(regex, seq):
	#Regex must be lower case!
	return bool( re.search(regex, seq.lower()) ) | bool( re.search(regex, reverseComplement(seq.lower()) ) )

#####Origins#####
def HasColE2(seq):
	#has ColE2 origin, data from PMID 16428404
	regexp = '....tga[gt]ac[ct]agataagcc[tgc]tatcagataacagcgcccttttggcgtctttttgagcacc' 
	return HasFeature(regexp, seq)
	#necessary and sufficient element for ColE2 replication, however a longer sequence is needed for stable replication
	# 'AGCGCCTCAGCGCGCCGTAGCGTCGATAAAAATTACGGGCTGGGGCGAAACTACCATCTGTTCGAAAAGGTCCGTAAATGGGCCTACAGAGCGATTCGTCAGGGCTGGCCTGTATTCTCACAATGGCTTGATGCCGTTATCCAGCGTGTCGAAATGTACAACGCTTCGCTTCCCGTTCCGCTTTCTCCGGCTGAATGTCGGGCTATTGGCAAGAGCATTGCGAAATATACACACAGGAAATTCTCACCAGAGGGATTTTCCGCTGTACAGGCCGCTCGCGGTCGCAAGGGCGGAACTAAATCTAAGCGCGCAGCAGTTCCTACATCAGCACGTTCGCTGAAACCGTGGGAGGCATTAGGCATCAGTCGAGCGACGTACTACCGAAAATTAAAATGTGACCCAGACCTCGCnnnntga'
	#longer element shown in the Anderson lab that stably replicates


def HasR6K(seq):
	#has R6k, data from Anderson lab observations
	regexp = 'gcagttcaacctgttgatagtacgtactaagctctcatgtttcacgtactaagctctcatgtttaacgtactaagctctcatgtttaacgaactaaaccctcatggctaacgtactaagctctcatggctaacgtactaagctctcatgtttcacgtactaagctctcatgtttgaacaataaaattaatataaatcagcaacttaaatagcctctaaggttttaagttttataagaaaaaaaagaatatataaggcttttaaagcttttaaggtttaacggttgtggacaacaagccagggatgtaacgcactgagaagcccttagagcctctcaaagcaattttgagtgacacaggaacacttaacggctgacatggg'.lower()
	return HasFeature(regexp, seq)

def HasP15A(seq):
	regex = 'aatattttatctgattaataagatgatcttcttgagatcgttttggtctgcgcgtaatctcttgctctgaaaacgaaaaaaccgccttgcagggcggtttttcgaaggttctctgagctaccaactctttgaaccgaggtaactggcttggaggagcgcagtcaccaaaacttgtcctttcagtttagccttaaccggcgcatgacttcaagactaactcctctaaatcaattaccagtggctgctgccagtggtgcttttgcatgtctttccgggttggactcaagacgatagttaccggataaggcgcagcggtcggactgaacggggggttcgtgcatacagtccagcttggagcgaactgcctacccggaactgagtgtcaggcgtggaatgagacaaacgcggccataacagcggaatgacaccggtaaaccgaaaggcaggaacaggagagcgcacgagggagccgccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccaccactgatttgagcgtcagatttcgtgatgcttgtcaggggggcggagcctatggaaaaacggctttgccgcggccctctcacttccctgttaagtatcttcctggcatcttccaggaaatctccgccccgttcgtaagccatttccgctcgccgcagtcgaacgaccgagcgtagcgagtcagtgagcgaggaagcggaatatatcctgtatcacatattctgctgacgcaccggtgcagccttttttctcctgccacatgaagcacttcactgacaccctcatcagtgccaacatagtaag'
	return HasFeature(regex, seq)

def HaspUC(seq):
	regex = 'cccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagcattgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacat'
	return HasFeature(regex, seq)

def HasAAFeature(regex, DNAseq):
	#must be uppercase, checks all six possibilities, fwd, rev x 3 frames
	seq = DNAseq
        retval = bool( re.search(regex, translate(seq.upper() )) ) | bool( re.search(regex,translate(seq[1:].upper() ) ) ) |  bool( re.search(regex,translate(seq[2:].upper() ) ) )
	seq = reverseComplement(seq)
	retval = retval | bool( re.search(regex, translate(seq.upper() )) ) | bool( re.search(regex,translate(seq[1:].upper() ) ) ) |  bool( re.search(regex,translate(seq[2:].upper() ) ) )
	return retval
 
def HasSpecR(seq):
	regex='MRSRNWSRTLTERSGGNGAVAVFMACYDCFFGVQSMPRASKQQARYAVGRCLMLWSSNDVTQQGSRPKTKLNIMREAVIAEVSTQLSEVVGVIERHLEPTLLAVHLYGSAVDGGLKPHSDIDLLVTVTVRLDETTRRALINDLLETSASPGESEILRAVEVTIVVHDDIIPWRYPAKRELQFGEWQRNDILAGIFEPATIDIDLAILLTKAREHSVALVGPAAEELFDPVPEQDLFEALNETLTLWNSPPDWAGDERNVVLTLSRIWYSAVTGKIAPKDVAADWAMERLPAQYQPVILEARQAYLGQEEDRLASRADQLEEFVHYVKGEITKVVGK'
	return HasAAFeature(regex, seq)
def HasAmpR(seq):
	regex='MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW'
	return HasAAFeature(regex, seq)
def HasKanR(seq):
	regex='MSHIQRETSCSRPRLNSNMDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDVTDEMVRLNWLTEFMPLPTIKHFIRTPDDAWLLTTAIPGKTAFQVLEEYPDSGENIVDALAVFLRRLHSIPVCNCPFNSDRVFRLAQAQSRMNNGLVDASDFDDERNGWPVEQVWKEMHKLLPFSPDSVVTHGDFSLDNLIFDEGKLIGCIDVGRVGIADRYQDLAILWNCLGEFSPSLQKRLFQKYGIDNPDMNKLQFHLMLDEFF'
	return HasAAFeature(regex, seq)
def HasCmR(seq):
	regex='MEKKITGYTTVDISQWHRKEHFEAFQSVAQCTYNQTVQLDITAFLKTVKKNKHKFYPAFIHILARLMNAHPEFRMAMKDGELVIWDSVHPCYTVFHEQTETFSSLWSEYHDDFRQFLHIYSQDVACYGENLAYFPKGFIENMFFVSANPWVSFTSFDLNVANMDNFFAPVFTMGKYYTQGDKVLMPLAIQVHHAVCDGFHVGRMLNELQQYCDEWQGGA'
	return HasAAFeature(regex, seq)
def HasResistance(seq):
	retval = []
	if HasCmR(seq):
		retval.append( 'CmR' )
	if HasKanR(seq):
		retval.append('KanR')
	if HasAmpR(seq):
		retval.append('AmpR')
	if HasSpecR(seq):
		retval.append('SpecR')
	return retval
def HasReplicon(seq):
	retval = []
	if HasColE2(seq):
		retval.append('ColE2')
	if HasR6K(seq):
		retval.append('R6K')
	if HasP15A(seq):
		retval.append('P15A')
	if HaspUC(seq):
		retval.append('pUC')
	return retval
class Strain(object):
	def __init__(self, name="", replication="", resistance="", plasmid=DNA("",'plasmid')):
		#pass everything in as a comma separated list
		self.name = name
		self.replication = replication.split(",")
		self.resistance = resistance.split(",") #should include the plasmid resistance!
		self.plasmid = plasmid #DNA object
	

#accepts list of dnas and a strain, it should output a list of DNAs that survive the transformation
def TransformPlate(DNAs, strain, selection_antibiotic):
	#strain is an object
	transformed = []
	for dna in DNAs:
		#check if circular, confers new resistance on strain, and doesn't compete with existing plasmid in strain
		if dna.topology == 'circular':
			newR = False
			replicon_ok = False
			no_existing_plasmid = False
			err_msg = ""
			success_msg = ""
			resistances = HasResistance(dna.sequence)
			replicons = HasReplicon(dna.sequence)
			#just need one resistance not already in strain
			for resistance in resistances:
				if not(resistance in strain.resistance):
					newR = True
					success_msg += "Use "+resistance
			for replicon in replicons:
				#has the pir/repA necessary for ColE2/R6K?
				if replicon in strain.replication:
					replicon_ok = True
			for replicon in replicons:
				#check if existing plasmid would compete
				if not(replicon in HasReplicon(strain.plasmid.sequence) ):
					no_existing_plasmid = True
			if(newR & replicon_ok & no_existing_plasmid):
				transformed.append(dna)	
				print success_msg
			else:
				if not(newR):
					print "Plasmid either doesn't have an antibiotic resistance or doesn't confer a new one on this strain"
				if not(replicon_ok):
					print "Plasmid replicon won't function in this strain"
				if not(no_existing_plasmid):
					print "Transformed plasmid replicon competes with existing plasmid in strain"
	if len(transformed)<1:
		print "No DNAs successfully transformed.  DNAs may be linear."
        return 0

yes = DNA('GATCCtaaCTCGAcgtgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatCAATTCGACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAATAATTGCTCATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTGCCATCCAGCTGATATCCCCTATAGTGAGTCGTATTACATGGTCATAGCTGTTTCCTGGCAGCTCTGGCCCGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGCCTCCTCTAGACCAGCCAGGACAGAAATGCCTCGACTTCGCTGCTGCCCAAGGTTGCCGGGTGACGCACACCGTGGAAACGGATGAAGGCACGAACCCAGTGGACATAAGCCTGTTCGGTTCGTAAGCTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGGGGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAAACATCATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATCGAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGCGGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAAACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGCGAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGTTATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGTATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAACATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAGGATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCTGGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGCAAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTATCAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCCTCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTAGTCGGCAAATAACCCTCGAGCCACCCATGACCAAAATCCCTTAACGTGAGTTACGCGTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCATTGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTcctaggTGTAAAACGACGGCCAGTCTTAAGCTCGGGCCCCAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCTCCGAATTGgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctgGAATTCATGA', 'plasmid')
pth1601kc = DNA('GATCCtaaCTCGAcgtgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatCAATTCGACCCAGCTTTCTTGTACAAAGTGGTTGATCcttacAGATCCcggttatccacagaatcaggggAGGCCTtagaaatattttatctgattaataagatgatcttcttgagatcgttttggtctgcgcgtaatctcttgctctgaaaacgaaaaaaccgccttgcagggcggtttttcgaaggttctctgagctaccaactctttgaaccgaggtaactggcttggaggagcgcagtcaccaaaacttgtcctttcagtttagccttaaccggcgcatgacttcaagactaactcctctaaatcaattaccagtggctgctgccagtggtgcttttgcatgtctttccgggttggactcaagacgatagttaccggataaggcgcagcggtcggactgaacggggggttcgtgcatacagtccagcttggagcgaactgcctacccggaactgagtgtcaggcgtggaatgagacaaacgcggccataacagcggaatgacaccggtaaaccgaaaggcaggaacaggagagcgcacgagggagccgccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccaccactgatttgagcgtcagatttcgtgatgcttgtcaggggggcggagcctatggaaaaacggctttgccgcggccctctcacttccctgttaagtatcttcctggcatcttccaggaaatctccgccccgttcgtaagccatttccgctcgccgcagtcgaacgaccgagcgtagcgagtcagtgagcgaggaagcggaatatatcctgtatcacatattctgctgacgcaccggtgcagccttttttctcctgccacatgaagcacttcactgacaccctcatcagtgccaacatagtaagccagtatacactccgctagcgctgatgtccggcggtgcGCATGCcgttaagggattttggtcatgaACTAGCttgatcgggcacgtaagaggttccaactttcaccataatgaaataagatcactaccgggcgtattttttgagttatcgagattttcaggagctaaggaagctaaaatggagaaaaaaatcactggatataccaccgttgatatatcccaatggcatcgtaaagaacattttgaggcatttcagtcagttgctcaatgtacctataaccagaccgttcagctggatattacggcctttttaaagaccgtaaagaaaaataagcacaagttttatccggcctttattcacattcttgcccgcctgatgaatgctcatccggaatttcgtatggcaatgaaagacggtgagctggtgatatgggatagtgttcacccttgttacaccgttttccatgagcaaactgaaacgttttcatcgctctggagtgaataccacgacgatttccggcagtttctacacatatattcgcaagatgtggcgtgttacggtgaaaacctggcctatttccctaaagggtttattgagaatatgtttttcgtctcagccaatccctgggtgagtttcaccagttttgatttaaacgtggccaatatggacaacttcttcgcccccgttttcaccatgggcaaatattatacgcaaggcgacaaggtgctgatgccgctggcgattcaggttcatcatgccgtttgtgatggcttccatgtcggcagaatgcttaatgaattacaacagtactgcgatgagtggcagggcggggcgtaatttgatatcgagctcgcttggactcctgttgatagatccagtaatgacctcagaactccatctggatttgttcagaacgctcggttgccgccgggcgttttttattggtgagaatccaagcCTGCAGataacttcgtatagcatacattatacgaagttatctcgagctgatccttcaactcagcaaaagttcgatttattcaacaaagccacgttgtgtctcaaaatctctgatgttacattgcacaagataaaaatatatcatcatgaacaataaaactgtctgcttacataaacagtaatacaaggggtgttatgagccatattcaacgggaaacgtcttgctccaggccgcgattaaattccaacatggatgctgatttatatgggtataaatgggctcgcgataatgtcgggcaatcaggtgcgacaatctatcgattgtatgggaagcccgatgcgccagagttgtttctgaaacatggcaaaggtagcgttgccaatgatgttacagatgagatggtcagactaaactggctgacggaatttatgcctcttccgaccatcaagcattttatccgtactcctgatgatgcatggttactcaccactgcgatccccgggaaaacagcattccaggtattagaagaatatcctgattcaggtgaaaatattgttgatgcgctggcagtgttcctgcgccggttgcattcgattcctgtttgtaattgtccttttaacagcgatcgcgtatttcgtctcgctcaggcgcaatcacgaatgaataacggtttggttgatgcgagtgattttgatgacgagcgtaatggctggcctgttgaacaagtctggaaagaaatgcataagcttttgccattctcaccggattcagtcgtcactcatggtgatttctcacttgataaccttatttttgacgaggggaaattaataggttgtattgatgttggacgagtcggaatcgcagaccgataccaggatcttgccatcctatggaactgcctcggtgagttttctccttcattacagaaacggctttttcaaaaatatggtattgataatcctgatatgaataaattgcagtttcatttgatgctcgatgagtttttctaatcagaattggttaattggttgtaacactggcagagcattacgctgacttgacggGaattgCCATTATCAACAAGTTTGTACAAAAAAGCAGGCTCCGAATTGgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctgGAATTCATGA', 'plasmid')
lol = Strain("name", "pUC,R6K,ColE2", "KanR,CmR", yes)
TransformPlate([yes], lol, "SpecR")

