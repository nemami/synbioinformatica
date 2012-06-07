#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

import sys, random, re

# Suffix Tree implementation from: http://chipsndips.livejournal.com/2005/12/07/
inf = 1000000

# Define a class for a node in the suffix tree
class SuffixNode(dict):
	def __init__(self):
		self.suffixLink = None # Suffix link as defined by Ukkonen

	def Print(self,str,ws=""):
		for t in self:
			k,p,s = self[t]
			if p == inf:
				print "%s%s" % (ws, str[k:])
			else:
				print "%s%s" % (ws, str[k:p+1])
				s.Print(str,ws+"|"*(p-k+1))
		
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
				# Decend the suffixtree until state s has a transition for the stringstr[k:i-1]
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
					if t == str[kp+i-k]:
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

	def Print(self):
		self.root.Print(self.str)

def PCR(template, fwd_primer, rev_primer):
	template = template + '$'
	fwd_primer = fwd_primer + '$'
	rev_primer = rev_primer + '$'
	fwdMatch = LCS(template,fwd_primer)
	fwdTuple = fwdMatch.LongestCommonSubstring()
	first = re.compile(fwdTuple[0])
	fList = first.findall(template)
	if len(fList) > 1:
		return ('',0,0)
	emails = re.findall(r'', fwdTuple[0])
	revcomp = reverseComplement(rev_primer)
	revMatch = LCS(template,revcomp+'$')
	revTuple = revMatch.LongestCommonSubstring()
	last = re.compile(revTuple[0])
	lList = last.findall(template)
	if len(lList) > 1:
		return ('',0,0)
	if fwdTuple[1] < revTuple[1] and fwdTuple[2] < revTuple[2]:
		return (template[fwdTuple[1]:revTuple[2]], fwdTuple[1], revTuple[2])
	else:
		return ('',0,0)

def reverseComplement(sequence):
  complement = {'a':'t','c':'g','g':'c','t':'a','n':'n'}
  return "".join([complement.get(nucleotide.lower(), '') for nucleotide in sequence[::-1]])

# def main():
# 	print PCR('aatgctactactattagtagaattgatgccaccttttcagctcgcgccccaaatgaaaatatagctaaacaggttattgaccatttgcgaaatgtatctaatggtcaaactaaatctactcgttcgcagaattgggaatcaactgttatatggaatgaaacttccagacaccgtactttagttgcatatttaaaacatgttgagctacagcattatattcagcaattaagctctaagccatccgcaaaaatgacctcttatcaaaaggagcaattaaaggtactctctaatcctgacctgttggagtttgcttccggtctggttcgctttgaagctcgaattaaaacgcgatatttgaagtctttcgggcttcctcttaatctttttgatgcaatccgctttgcttctgactataatagtcagggtaaagacctgatttttgatttatggtcattctcgttttctgaactgtttaaagcatttgagggggattcaatgaatatttatgacgattccgcagtattggacgctatccagtctaaacattttactattaccccctctggcaaaacttcttttgcaaaagcctctcgctattttggtttttatcgtcgtctggtaaacgagggttatgatagtgttgctcttactatgcctcgtaattccttttggcgttatgtatctgcattagttgaatgtggtattcctaaatctcaactgatgaatctttctacctgtaataatgttgttccgttagttcgttttattaacgtagatttttcttcccaacgtcctgactggtataatgagccagttcttaaaatcgcataaggtaattcacaatgattaaagttgaaattaaaccatctcaagcccaatttactactcgttctggtgtttctcgtcagggcaagccttattcactgaatgagcagctttgttacgttgatttgggtaatgaatatccggttcttgtcaagattactcttgatgaaggtcagccagcctatgcgcctggtctgtacaccgttcatctgtcctctttcaaagttggtcagttcggttcccttatgattgaccgtctgcgcctcgttccggctaagtaacatggagcaggtcgcggatttcgacacaatttatcaggcgatgatacaaatctccgttgtactttgtttcgcgcttggtataatcgctgggggtcaaagatgagtgttttagtgtattcttttgcctctttcgttttaggttggtgccttcgtagtggcattacgtattttacccgtttaatggaaacttcctcatgaaaaagtctttagtcctcaaagcctctgtagccgttgctaccctcgttccgatgctgtctttcgctgctgagggtgacgatcccgcaaaagcggcctttaactccctgcaagcctcagcgaccgaatatatcggttatgcgtgggcgatggttgttgtcattgtcggcgcaactatcggtatcaagctgtttaagaaattcacctcgaaagcaagctgataaaccgatacaattaaaggctccttttggagccttttttttggagattttcaacgtgaaaaaattattattcgcaattcctttagtggtacctttctattctcactcggccgaaactgttgaaagttgtttagcaaaatcccatacagaaaattcatttactaacgtctggaaagacgacaaaactttagatcgttacgctaactatgagggctgtctgtggaatgctacaggcgttgtagtttgtactggtgacgaaactcagtgttacggtacatgggttcctattgggcttgctatccctgaaaatgagggtggtggctctgagggtggcggttctgagggtggcggttctgagggtggcggtactaaacctcctgagtacggtgatacacctattccgggctatacttatatcaaccctctcgacggcacttatccgcctggtactgagcaaaaccccgctaatcctaatccttctcttgaggagtctcagcctcttaatactttcatgtttcagaataataggttccgaaataggcagggggcattaactgtttatacgggcactgttactcaaggcactgaccccgttaaaacttattaccagtacactcctgtatcatcaaaagccatgtatgacgcttactggaacggtaaattcagagactgcgctttccattctggctttaatgaggatttatttgtttgtgaatatcaaggccaatcgtctgacctgcctcaacctcctgtcaatgctggcggcggctctggtggtggttctggtggcggctctgagggtggtggctctgagggtggcggttctgagggtggcggctctgagggaggcggttccggtggtggctctggttccggtgattttgattatgaaaagatggcaaacgctaataagggggctatgaccgaaaatgccgatgaaaacgcgctacagtctgacgctaaaggcaaacttgattctgtcgctactgattacggtgctgctatcgatggtttcattggtgacgtttccggccttgctaatggtaatggtgctactggtgattttgctggctctaattcccaaatggctcaagtcggtgacggtgataattcacctttaatgaataatttccgtcaatatttaccttccctccctcaatcggttgaatgtcgcccttttgtctttggcgctggtaaaccatatgaattttctattgattgtgacaaaataaacttattccgtggtgtctttgcgtttcttttatatgttgccacctttatgtatgtattttctacgtttgctaacatactgcgtaataaggagtcttaatcatgccagttcttttgggtattccgttattattgcgtttcctcggtttccttctggtaactttgttcggctatctgcttacttttcttaaaaagggcttcggtaagatagctattgctatttcattgtttcttgctcttattattgggcttaactcaattcttgtgggttatctctctgatattagcgctcaattaccctctgactttgttcagggtgttcagttaattctcccgtctaatgcgcttccctgtttttatgttattctctctgtaaaggctgctattttcatttttgacgttaaacaaaaaatcgtttcttatttggattgggataaataatatggctgtttattttgtaactggcaaattaggctctggaaagacgctcgttagcgttggtaagattcaggataaaattgtagctgggtgcaaaatagcaactaatcttgatttaaggcttcaaaacctcccgcaagtcgggaggttcgctaaaacgcctcgcgttcttagaataccggataagccttctatatctgatttgcttgctattgggcgcggtaatgattcctacgatgaaaataaaaacggcttgcttgttctcgatgagtgcggtacttggtttaatacccgttcttggaatgataaggaaagacagccgattattgattggtttctacatgctcgtaaattaggatgggatattatttttcttgttcaggacttatctattgttgataaacaggcgcgttctgcattagctgaacatgttgtttattgtcgtcgtctggacagaattactttaccttttgtcggtactttatattctcttattactggctcgaaaatgcctctgcctaaattacatgttggcgttgttaaatatggcgattctcaattaagccctactgttgagcgttggctttatactggtaagaatttgtataacgcatatgatactaaacaggctttttctagtaattatgattccggtgtttattcttatttaacgccttatttatcacacggtcggtatttcaaaccattaaatttaggtcagaagatgaaattaactaaaatatatttgaaaaagttttctcgcgttctttgtcttgcgattggatttgcatcagcatttacatatagttatataacccaacctaagccggaggttaaaaaggtagtctctcagacctatgattttgataaattcactattgactcttctcagcgtcttaatctaagctatcgctatgttttcaaggattctaagggaaaattaattaatagcgacgatttacagaagcaaggttattcactcacatatattgatttatgtactgtttccattaaaaaaggtaattcaaatgaaattgttaaatgtaattaattttgttttcttgatgtttgtttcatcatcttcttttgctcaggtaattgaaatgaataattcgcctctgcgcgattttgtaacttggtattcaaagcaatcaggcgaatccgttattgtttctcccgatgtaaaaggtactgttactgtatattcatctgacgttaaacctgaaaatctacgcaatttctttatttctgttttacgtgcaaataattttgatatggtaggttctaacccttccattattcagaagtataatccaaacaatcaggattatattgatgaattgccatcatctgataatcaggaatatgatgataattccgctccttctggtggtttctttgttccgcaaaatgataatgttactcaaacttttaaaattaataacgttcgggcaaaggatttaatacgagttgtcgaattgtttgtaaagtctaatacttctaaatcctcaaatgtattatctattgacggctctaatctattagttgttagtgctcctaaagatattttagataaccttcctcaattcctttcaactgttgatttgccaactgaccagatattgattgagggtttgatatttgaggttcagcaaggtgatgctttagatttttcatttgctgctggctctcagcgtggcactgttgcaggcggtgttaatactgaccgcctcacctctgttttatcttctgctggtggttcgttcggtatttttaatggcgatgttttagggctatcagttcgcgcattaaagactaatagccattcaaaaatattgtctgtgccacgtattcttacgctttcaggtcagaagggttctatctctgttggccagaatgttccttttattactggtcgtgtgactggtgaatctgccaatgtaaataatccatttcagacgattgagcgtcaaaatgtaggtatttccatgagcgtttttcctgttgcaatggctggcggtaatattgttctggatattaccagcaaggccgatagtttgagttcttctactcaggcaagtgatgttattactaatcaaagaagtattgctacaacggttaatttgcgtgatggacagactcttttactcggtggcctcactgattataaaaacacttctcaggattctggcgtaccgttcctgtctaaaatccctttaatcggcctcctgtttagctcccgctctgattctaacgaggaaagcacgttatacgtgctcgtcaaagcaaccatagtacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgatttgggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcgggctattcttttgatttataagggattttgccgatttcggaaccaccatcaaacaggattttcgcctgctggggcaaaccagcgtggaccgcttgctgcaactctctcagggccaggcggtgaagggcaatcagctgttgcccgtctcactggtgaaaagaaaaaccaccctggcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagcttgcatgcctgcaggtcctcgaattcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgctttgcctggtttccggcaccagaagcggtgccggaaagctggctggagtgcgatcttcctgaggccgatactgtcgtcgtcccctcaaactggcagatgcacggttacgatgcgcccatctacaccaacgtgacctatcccattacggtcaatccgccgtttgttcccacggagaatccgacgggttgttactcgctcacatttaatgttgatgaaagctggctacaggaaggccagacgcgaattatttttgatggcgttcctattggttaaaaaatgagctgatttaacaaaaatttaatgcgaattttaacaaaatattaacgtttacaatttaaatatttgcttatacaatcttcctgtttttggggcttttctgattatcaaccggggtacatatgattgacatgctagttttacgattaccgttcatcgattctcttgtttgctccagactctcaggcaatgacctgatagcctttgtagatctctcaaaaatagctaccctctccggcattaatttatcagctagaacggttgaatatcatattgatggtgatttgactgtctccggcctttctcacccttttgaatctttacctacacattactcaggcattgcatttaaaatatatgagggttctaaaaatttttatccttgcgttgaaataaaggcttctcccgcaaaagtattacagggtcataatgtttttggtacaaccgatttagctttatgctctgaggctttattgcttaattttgctaattctttgccttgcctgtatgatttattggatgtt','attagtagaattgatgccacctttt','actattatagtcagaagcaaag')

# if __name__ == "__main__":
# 	main()