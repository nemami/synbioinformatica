#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

import synbioinformatica, sys
from synbioinformatica import DNA, restrictionEnzyme, Ligate

def main(template, primer_1, primer_2):
	print PCR(template, primer_1, primer_2)

if __name__ == "__main__":
	InitializeEnzymes = getattr(synbioinformatica,'EnzymeDictionary')
	PCR = getattr(synbioinformatica,'PCR')
	Digest = getattr(synbioinformatica,'Digest')
	Ligate = getattr(synbioinformatica,'Ligate')
	GelPurify = getattr(synbioinformatica,'GelPurify')
	ZymoPurify = getattr(synbioinformatica,'ZymoPurify')
	EnzymeDictionary = InitializeEnzymes()
	plasmid = DNA('ATGGTGAGATGAGTGAAGGCGAGCTGGTGGATGCATTCCGCCATGTGAGTGATGCGTTTGAGCAAACCAGCGAAACCATCGGCGTGCGCGCCAATAACGCGATCAACGACATGGTGCGTCAACGTCTGCTGAACCGCTTTACCAGCGAGCAGGCGGAAGGGAACGCAATTTACCGTCTGACGCCGCTCGGCATCGGCATTACTGACTACNNNATCCGTCAGCGCGAGTTTTCTACGCTGCGTCTTTCTATGCAGTTGTCGATTGTGGCGGGTGAGCTCAAACGCGCAGCAGATGCCGCCGAAGAGGGCGGTGATGAATTTCACTGGCACCGTAATGTCTATGCGCCACTGAAATATTCGGTAGCAGAAATTTTCGACAGTATCGACCTGACGCAACGTCTGATGGACGAACAGCAGCAGCAGGTGAAGGACGATATCGCCCAGTTGCTGAACAAAGACTGGCGGGCGGCGATTTCCAGCTGTGAATTGTTGCTTTCGGAAACTTCCGGAACGCTGCGTGAATTGCAGGATACGCTGGAAGCGGCAGGCGACAAATTGCAGGCTAATCTGTTGCGCATTCAGGATGCGACGATGACCCATGACGATCTGCATTTCGTCGATCGTCTGGTGTTCGATCTGCAGAGCAAACTCGATCGTATTATCAGTTGGGGCCAGCAATCCATCGACTTGTGGATTGGCTACGACCGCCACGTACACAAATTTATTCGTACCGCGATCGATATGGATAAAAACCGCGTCTTTGCTCAGCGGTTACGTCAGTCGGTACAAACCTATTTTGATGAGCCGTGGGCGCTAACTTATGCCAATGCCGATCGTCTGCTGGATATGCGTGACGAAGAGATGGCACTGCGCGATGAAGAAGTGACTGGGGAACTTCCTGAGGATCTGGAATACGAAGAGTTTAACGAGATCCGCGAACAGCTGGCGGCGATCATCGAAGAACAACTTGCCGTGTACAAAACCAGACAAGTGCCGCTGGATCTTGGTCTGGTGGTACGCGAATATCTGTCACAGTATCCGCGTGCACGTCACTTTGACGTTGCGCGTATTGTTATTGATACCTGACGCAACGTCTGC','plasmid','bNE015')
	primer_1 = DNA("GTGAAGGCGAGCTGGTGGATGC",'primer','ne001')
	primer_2 = DNA("CTTCTTCATCGCGCAGTGCCATCT",'primer','ne002')
	pcrProduct = PCR(primer_1, primer_2, plasmid)
	pcrProduct.prettyPrint()
	print "\n\tPCR product instructions: "+pcrProduct.instructions+"\n"
	print "\t"+str(pcrProduct.parents)+"\n"
	print "\t"+str(pcrProduct.children)+"\n"
	# products = Digest(plasmid,[EnzymeDictionary["BsaAI"]])
	# products = Digest(plasmid,[EnzymeDictionary["EcoRI"], EnzymeDictionary["BceAI"]])
	plasmid = DNA('aatgctactactattagtagaattgatgccaccttttcagctcgcgccccaaatgaaaatatagctaaacaggttattgaccatttgcgaaatgtatctaatggtcaaactaaatctactcgttcgcagaattgggaatcaactgttatatggaatgaaacttccagacaccgtactttagttgcatatttaaaacatgttgagctacagcattatattcagcaattaagctctaagccatccgcaaaaatgacctcttatcaaaaggagcaattaaaggtactctctaatcctgacctgttggagtttgcttccggtctggttcgctttgaagctcgaattaaaacgcgatatttgaagtctttcgggcttcctcttaatctttttgatgcaatccgctttgcttctgactataatagtcagggtaaagacctgatttttgatttatggtcattctcgttttctgaactgtttGGATCCaaagcatttgagggggattcaatgaatatttatgacgattccgcagtattggacgctatccagtctaaacattttactattaccccctctggcaaaacttcttttgcaaaagcctctcgctattttggtttttatcgtcgtctggtaaacgagggttatgatagtgttgctcttactatgcctcgtaattccttttggcgttatgtatctgcattagttgaatgtggtattcctaaatctcaactgatgaatctttctacctgtaataatgttgttccgttagttcgttttattaacgtagatttttcttcccaacgtcctgactggtataatgagccagttcttaaaatcgcataaggtaattcacaatgattaaagttgaaattaaaccatctcaagcccaatttactactcgttctggtgtttctcgtcagggcaagccttattcactgaatgagcagctttgttacgttgatttgggtaatgaatatccggttcttgtcaagattactcttgatgaaggtcagccagcctatgcgcctggtctgtacaccgttcatctgtcctctttcaaagttggtcagttcggttcccttatgattgaccgtctgcgccAGATCTcttcagaccAGATCTagacaccgtactttagttgcatatttaaaacatgttgagctacagcattatattcagcaattaagctctaagccatccgcaaaaatgacctcttatcaaaaggagcaattaaaggtactctctaatcctgacctgttggagtttgcttccggtctggttcgctttgaagctcgaattaaaacgcgatatttgaagtctttcgggcttcctcttaatctttttgatgcaatccgctttgcttctgactataatagtcagggtaaagacctgatttttgatttatggtcattctcgttttctgaactgtttGGATCCaaagcatttgagggggattcaatgaatatttatgacgattccgcagtattggacgctatccagtctaaacattttactattaccccctctggcaaaacttcttttgcaaaagcctctcgctattttggtttttatcgtcgtctggtaaacgagggttatgatagtgttgctcttactatgcctcgtaattccttttggcgttatgtatctgcattagttgaatgtggtattcctaaatctcaactgatgaatctttctacctgtaataatgttgttccgttagttcgttttattaacgtagatttttcttcccaacgtcctgactggtataatgagccagttcttaaaatcgcataaggtaattcacaatgattaaagttgaaattaaaccatctcaagcccaatttactactcgttctggtgtttctcgtcagggcaagccttattcactgaatgagcagctttgttacgttgatttgggtaatgaatatccggttcttgtcaagattactcttgatgaaggtcagccagcctatgcgcctggtctgtacaccgttcatctgtcctctttcaaagttggtcagttcggttcccttatgattgaccgtctgcgccAGATCTaatgctactactattagtagaattgatgccaccttttcagctcgcgccccaaatgaaaatatagctaaacaggttattgaccatttgcgaaatgtatctaatggtcaaactaaatctactcgttcgcagaattgggaatcaactgttatatggaatgaaacttccagacaccgtactttagttgcatatttaaaacatgttgagctacagcattatattcagcaattaagctctaagccatccg','plasmid','bNE100')
	products = Digest(plasmid,[EnzymeDictionary["BamHI"], EnzymeDictionary["BglII"]])
	for product in products:
		print "\tDigest Fragment instructions: "+product.instructions
	bands = GelPurify(products,'L')
	zymoed = ZymoPurify(bands)
	ligated = Ligate(zymoed)
	for lig in ligated:
		print "\tLigation product instructions: "+lig.instructions+"\n\n"
	print '\nInput DNA:\n'
	plasmid.prettyPrint()
	print '\nDigestion product(s):\n'
	n = 1
	for product in products:
		print 'digest #'+str(n)+':'
		product.prettyPrint()
		n = n + 1
	print '\nGel Purification product(s):\n'
	n = 1
	for product in zymoed:
		print 'purified #'+str(n)+':'
		product.prettyPrint()
		n = n + 1
	print "\nLigation product(s):\n"
	n = 1
	for product in ligated:
		print 'ligation product #'+str(n)+':'
		product.prettyPrint()
		n = n + 1
