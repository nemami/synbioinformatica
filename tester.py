#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

# PCR function tester script
# Takes 3 invocation arguments: 
# template, forward primer, reverse primer (in that order)
# returns tuple: (PCR product, start index, end index)

import synbioinformatica, sys
from DNA import DNA, restrictionEnzyme

def main(template, primer_1, primer_2):
	print PCR(template, primer_1, primer_2)

if __name__ == "__main__":
	InitializeEnzymes = getattr(synbioinformatica,'EnzymeDictionary')
	PCR = getattr(synbioinformatica,'PCR')
	digest = getattr(synbioinformatica,'digest')
	# main(sys.argv[1],sys.argv[2],sys.argv[3])
	EnzymeDictionary = InitializeEnzymes()
	plasmid = DNA('TTCATGGTGAGATGAGTGAAGGCGAGCTGGTGGATGCATTCCGCCATGTGAGTGATGCGTTTGAGCAAACCAGCGAAACCATCGGCGTGCGCGCCAATAACGCGATCAACGACATGGTGCGTCAACGTCTGCTGAACCGCTTTACCAGCGAGCAGGCGGAAGGGAACGCAATTTACCGTCTGACGCCGCTCGGCATCGGCATTACTGACTACNNNATCCGTCAGCGCGAGTTTTCTACGCTGCGTCTTTCTATGCAGTTGTCGATTGTGGCGGGTGAGCTCAAACGCGCAGCAGATGCCGCCGAAGAGGGCGGTGATGAATTTCACTGGCACCGTAATGTCTATGCGCCACTGAAATATTCGGTAGCAGAAATTTTCGACAGTATCGACCTGACGCAACGTCTGATGGACGAACAGCAGCAGCAGGTGAAGGACGATATCGCCCAGTTGCTGAACAAAGACTGGCGGGCGGCGATTTCCAGCTGGATCCTGAATTGTTGCTTTCGGAAACTTCCGGAACGCTGCGTGAATTGCAGGATACGCTGGAAGCGGCAGGCGACAAATTGCAGGCTAATCTGTTGCGCATTCAGGATGCGACGATGACCCATGACGATCTGCATTTCGTCGATCGTCTGGTGTTCGATCTGCAGAGCAAACTCGATCGTATTATCAGTTGGGGCCAGCAATCCATCGACTTGTGGATTGGCTACGACCGCCACGTACACAAATTTATTCGTACCGCGATCGATATGGATAAAAACCGCGTCTTTGCTCAGCGGTTACGTCAGTCGGTACAAACCTATTTTGATGAGCCGTGGGCGCTAACTTATGCCAATGCCGATCGTCTGCTGGATATGCGTGACGAAGAGATGGCACTGCGCGATGAAGAAGTGACTGGGGAACTTCCTGAGGATCTGGAATACGAAGAGTTTAACGAGATCCGCGAACAGCTGGCGGCGATCATCGAAGAACAACTTGCCGTGTACAAAACCAGACAAGTGCCGCTGGATCTTGGTCTGGTGGTACGCGAATATCTGTCACAGTATCCGCGTGCACGTCACTTTGACGTTGCGCGTATTGTTATTGATACCTGACGCAACGTCTGCGAA','plasmid')
	# primer_1 = DNA("CCAGCTCGCCTTCACTCATCTCAwwwwwwwwwwwwww",'primer')
	# primer_2 = DNA("wwwwwwwwwwwwwwwwwATAACGCGATCAACGACATGGTGCGTCA",'primer')
	# pcrProduct = PCR(plasmid, primer_1, primer_2)
	# print pcrProduct.sequence
	product = digest(plasmid,[EnzymeDictionary["EcoRI"], EnzymeDictionary["BceAI"]])
	# print product