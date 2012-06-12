#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

import sys, random, re, math
from DNA import DNA, restrictionEnzyme, Overhang
from decimal import *

# TODO: for PCR, identification of primers on the edge of a circular sequence
# TODO: Digest function Error/Exception handling, e.g. proximity to terminal sequence on a linear fragment
# TODO: Comment digest function

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

# Define a class for a node in the suffix tree
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

def digest(InputDNA, Enzymes):
	# TODO: Error/Exception handling, e.g. proximity to terminal sequence on a linear fragment
	print 'Input DNA ('+InputDNA.topology+'): '+InputDNA.sequence
	indices = []			# for restriction sites
	frags = []				# producted DNA fragments
	sites = ""
	totalLength = len(InputDNA.sequence)
	if InputDNA.topology == "linear":	
		# Initialize indices array with start and end indices of the linear fragment
			# Add dummy REase for consistency
		dummy = restrictionEnzyme("dummy", "", "", "", "", "", 0, 0, "(0/0)","")
		indices = [(0,0,'',dummy), (len(InputDNA.sequence),0,'',dummy)]
	# Identify restriction sites, fill in indices array
	for enzyme in Enzymes:
		sites = enzyme.find_sites(InputDNA)
		for site in sites:
			site = site + (enzyme, )
			indices.append(site)
		indices.sort()
	# If you have overlapping restriction sites, choose the first one and discard the
		# second (TODO: there may be a better, non-greedy way to do this... not sure)
	for n in range(len(indices)-1):
		(currentTuple, nextTuple) = (indices[n], indices[n+1])
		(currentStart, nextStart) = (currentTuple[0], nextTuple[0])
		if currentStart + len(enzyme.recognition_site) >= nextStart:
			indices.pop(n+1)
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
			CTO = currentEnzyme.top_strand_offset
			CBO = currentEnzyme.bottom_strand_offset
			ALL = max(CTO,CBO)
		else:
			CTO = len(currentEnzyme.alpha_only_site) -1 * currentEnzyme.top_strand_offset
			CBO  = len(currentEnzyme.alpha_only_site) -1 * currentEnzyme.bottom_strand_offset
			ALL = max(CTO,CBO)-1
		if nextDirection == "sense":
			NTO = nextEnzyme.top_strand_offset
			NBO = nextEnzyme.bottom_strand_offset
			ARL = min(NTO,NBO)
		else:
			NTO = len(nextEnzyme.alpha_only_site) -1 * nextEnzyme.top_strand_offset
			NBO  = len(nextEnzyme.alpha_only_site) -1 * nextEnzyme.bottom_strand_offset
			ARL = min(NTO,NBO)-1
			# TODO: off by errors
			# from test case plasmid = DNA('TTCATGGTGAGATGAGTGAAGGCGAGCTGGTGGATGCATTCCGCCATGTGAGTGATGCGTTTGAGCAAACCAGCGAAACCATCGGCGTGCGCGCCAATAACGCGATCAACGACATGGTGCGTCAACGTCTGCTGAACCGCTTTACCAGCGAGCAGGCGGAAGGGAACGCAATTTACCGTCTGACGCCGCTCGGCATCGGCATTACTGACTACNNNATCCGTCAGCGCGAGTTTTCTACGCTGCGTCTTTCTATGCAGTTGTCGATTGTGGCGGGTGAGCTCAAACGCGCAGCAGATGCCGCCGAAGAGGGCGGTGATGAATTTCACTGGCACCGTAATGTCTATGCGCCACTGAAATATTCGGTAGCAGAAATTTTCGACAGTATCGACCTGACGCAACGTCTGATGGACGAACAGCAGCAGCAGGTGAAGGACGATATCGCCCAGTTGCTGAACAAAGACTGGCGGGCGGCGATTTCCAGCTGGATCCTGAATTGTTGCTTTCGGAAACTTCCGGAACGCTGCGTGAATTGCAGGATACGCTGGAAGCGGCAGGCGACAAATTGCAGGCTAATCTGTTGCGCATTCAGGATGCGACGATGACCCATGACGATCTGCATTTCGTCGATCGTCTGGTGTTCGATCTGCAGAGCAAACTCGATCGTATTATCAGTTGGGGCCAGCAATCCATCGACTTGTGGATTGGCTACGACCGCCACGTACACAAATTTATTCGTACCGCGATCGATATGGATAAAAACCGCGTCTTTGCTCAGCGGTTACGTCAGTCGGTACAAACCTATTTTGATGAACGGCGGGCGCTAACTTATGCCAATGCCGATCGTCTGCTGGATATGCGTGACGAAGAGATGGCACTGCGCGATGAAGAAGTGACTGGGGAACTTCCTGAGGATCTGGAATACGAAGAGTTTAACGAGATCCGCGAACAGCTGGCGGCGATCATCGAAGAACAACTTGCCGTGTACAAAACCAGACAAGTGCCGCTGGATCTTGGTCTGGTGGTACGCGAATATCTGTCACAGTATCCGCGTGCACGTCACTTTGACGTTGCGCGTATTGTTATTGATACCTGACGCAACGTCTGCGAA','plasmid')
				# left side of frag with currentEnzyme == 'sense' strand BceAI
				# right side of frag with nextEnzyme == 'sense' strand BceAI
		# Update start value currentStart and apply ( mod length ) to deal with edge cases
			# Also, update end value digEnd for fragment indices
		currentStart = currentStart+ALL
		currentStart = currentStart % totalLength
		digEnd = nextStart + ARL
		digested = DNA(InputDNA.sequence[currentStart:digEnd],'digest')
		# Adjust top and bottom overhang values based on the orientation of the restriction site
		if direction == "sense":
			TO = CTO
			BO = CBO
		else:
			TO = CBO
			BO = CTO
		difference = abs(abs(BO) - abs(TO))
		if abs(TO) < abs(BO) and direction == "sense" or abs(TO) > abs(BO) and direction == "antisense":
			# Edge case statement
			print 'got in here'
			print currentEnzyme.compsite
			print nextEnzyme.compsite
			print currentStart
			print difference
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
		# Adjust top and bottom overhang values based on the orientation of the restriction site
		if direction == "sense":
			TO = NTO
			BO = NBO
		else:
			TO = NBO
			BO = NTO
		difference = abs(abs(BO) - abs(TO))
		# Apply ( mod length ) operator to end index value digDiff to deal with edge cases
		digDiff = digEnd + difference
		digDiff = digDiff % totalLength
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
		frags.append((currentStart,digested))
		frags.sort()
	for frag in frags:
		print 'Fragment: '+frag[1].sequence
		print 'Fragment.TRO: '+frag[1].topRightOverhang.sequence
		print 'Fragment.BRO: '+frag[1].bottomRightOverhang.sequence
		print 'Fragment.TLO: '+frag[1].topLeftOverhang.sequence
		print 'Fragment.BLO: '+frag[1].bottomLeftOverhang.sequence
	return sites