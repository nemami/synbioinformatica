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
	DigestBuffer = getattr(synbioinformatica,'DigestBuffer')
	GelAndZymoPurify = getattr(synbioinformatica,'GelAndZymoPurify')
	EnzymeDictionary = InitializeEnzymes()
	print 'buffers:'
	# print DigestBuffer("BglII", "BamHI")
	plasmid = DNA('ATGGTGAGATGAGTGAAGGCGAGCTGGTGGATGCATTCCGCCATGTGAGTGATGCGTTTGAGCAAACCAGCGAAACCATCGGCGTGCGCGCCAATAACGCGATCAACGACATGGTGCGTCAACGTCTGCTGAACCGCTTTACCAGCGAGCAGGCGGAAGGGAACGCAATTTACCGTCTGACGCCGCTCGGCATCGGCATTACTGACTACNNNATCCGTCAGCGCGAGTTTTCTACGCTGCGTCTTTCTATGCAGTTGTCGATTGTGGCGGGTGAGCTCAAACGCGCAGCAGATGCCGCCGAAGAGGGCGGTGATGAATTTCACTGGCACCGTAATGTCTATGCGCCACTGAAATATTCGGTAGCAGAAATTTTCGACAGTATCGACCTGACGCAACGTCTGATGGACGAACAGCAGCAGCAGGTGAAGGACGATATCGCCCAGTTGCTGAACAAAGACTGGCGGGCGGCGATTTCCAGCTGTGAATTGTTGCTTTCGGAAACTTCCGGAACGCTGCGTGAATTGCAGGATACGCTGGAAGCGGCAGGCGACAAATTGCAGGCTAATCTGTTGCGCATTCAGGATGCGACGATGACCCATGACGATCTGCATTTCGTCGATCGTCTGGTGTTCGATCTGCAGAGCAAACTCGATCGTATTATCAGTTGGGGCCAGCAATCCATCGACTTGTGGATTGGCTACGACCGCCACGTACACAAATTTATTCGTACCGCGATCGATATGGATAAAAACCGCGTCTTTGCTCAGCGGTTACGTCAGTCGGTACAAACCTATTTTGATGAGCCGTGGGCGCTAACTTATGCCAATGCCGATCGTCTGCTGGATATGCGTGACGAAGAGATGGCACTGCGCGATGAAGAAGTGACTGGGGAACTTCCTGAGGATCTGGAATACGAAGAGTTTAACGAGATCCGCGAACAGCTGGCGGCGATCATCGAAGAACAACTTGCCGTGTACAAAACCAGACAAGTGCCGCTGGATCTTGGTCTGGTGGTACGCGAATATCTGTCACAGTATCCGCGTGCACGTCACTTTGACGTTGCGCGTATTGTTATTGATACCTGACGCAACGTCTGC','plasmid','bNE015')
	primer_1 = DNA("NNNNNNNNNCTTCTTCATCGCGCAGTGCCATCT",'primer','ne001')
	primer_2 = DNA("NNNNNNNNNGTGAAGGCGAGCTGGTGGATGC",'primer','ne002')
	pcrProduct = PCR(primer_1, primer_2, plasmid)
	pcrProduct.prettyPrint()
	print pcrProduct.sequence
	print "\n\tPCR product instructions: "+pcrProduct.instructions+"\n"
	print "\t"+str(pcrProduct.parents)+"\n"
	print "\t"+str(pcrProduct.children)+"\n"
	# # products = Digest(plasmid,[EnzymeDictionary["BsaAI"]])
	# # products = Digest(plasmid,[EnzymeDictionary["EcoRI"], EnzymeDictionary["BceAI"]])
	plasmid = DNA('aatgctactactattagtagaattgatgccaccttttcagctcgcgccccaaatgaaaatatagctaaacaggttattgaccatttgcgaaatgtatctaatggtcaaactaaatctactcgttcgcagaattgggaatcaactgttatatggaatgaaacttccagacaccgtactttagttgcatatttaaaacatgttgagctacagcattatattcagcaattaagctctaagccatccgcaaaaatgacctcttatcaaaaggagcaattaaaggtactctctaatcctgacctgttggagtttgcttccggtctggttcgctttgaagctcgaattaaaacgcgatatttgaagtctttcgggcttcctcttaatctttttgatgcaatccgctttgcttctgactataatagtcagggtaaagacctgatttttgatttatggtcattctcgttttctgaactgtttGGATCCaaagcatttgagggggattcaatgaatatttatgacgattccgcagtattggacgctatccagtctaaacattttactattaccccctctggcaaaacttcttttgcaaaagcctctcgctattttggtttttatcgtcgtctggtaaacgagggttatgatagtgttgctcttactatgcctcgtaattccttttggcgttatgtatctgcattagttgaatgtggtattcctaaatctcaactgatgaatctttctacctgtaataatgttgttccgttagttcgttttattaacgtagatttttcttcccaacgtcctgactggtataatgagccagttcttaaaatcgcataaggtaattcacaatgattaaagttgaaattaaaccatctcaagcccaatttactactcgttctggtgtttctcgtcagggcaagccttattcactgaatgagcagctttgttacgttgatttgggtaatgaatatccggttcttgtcaagattactcttgatgaaggtcagccagcctatgcgcctggtctgtacaccgttcatctgtcctctttcaaagttggtcagttcggttcccttatgattgaccgtctgcgccAGATCTcttcagaccAGATCTagacaccgtactttagttgcatatttaaaacatgttgagctacagcattatattcagcaattaagctctaagccatccgcaaaaatgacctcttatcaaaaggagcaattaaaggtactctctaatcctgacctgttggagtttgcttccggtctggttcgctttgaagctcgaattaaaacgcgatatttgaagtctttcgggcttcctcttaatctttttgatgcaatccgctttgcttctgactataatagtcagggtaaagacctgatttttgatttatggtcattctcgttttctgaactgtttGGATCCaaagcatttgagggggattcaatgaatatttatgacgattccgcagtattggacgctatccagtctaaacattttactattaccccctctggcaaaacttcttttgcaaaagcctctcgctattttggtttttatcgtcgtctggtaaacgagggttatgatagtgttgctcttactatgcctcgtaattccttttggcgttatgtatctgcattagttgaatgtggtattcctaaatctcaactgatgaatctttctacctgtaataatgttgttccgttagttcgttttattaacgtagatttttcttcccaacgtcctgactggtataatgagccagttcttaaaatcgcataaggtaattcacaatgattaaagttgaaattaaaccatctcaagcccaatttactactcgttctggtgtttctcgtcagggcaagccttattcactgaatgagcagctttgttacgttgatttgggtaatgaatatccggttcttgtcaagattactcttgatgaaggtcagccagcctatgcgcctggtctgtacaccgttcatctgtcctctttcaaagttggtcagttcggttcccttatgattgaccgtctgcgccAGATCTaatgctactactattagtagaattgatgccaccttttcagctcgcgccccaaatgaaaatatagctaaacaggttattgaccatttgcgaaatgtatctaatggtcaaactaaatctactcgttcgcagaattgggaatcaactgttatatggaatgaaacttccagacaccgtactttagttgcatatttaaaacatgttgagctacagcattatattcagcaattaagctctaagccatccg','plasmid','bNE100')
	products = Digest(plasmid,[EnzymeDictionary["BamHI"], EnzymeDictionary["BglII"]])
	for product in products:
		print "\tDigest Fragment instructions: "+product.instructions+"\n\n"
	bands = GelAndZymoPurify(products,'L')
	for band in bands:
		print "\tGelAndZymo Fragment instructions: "+band.instructions+"\n\n"
	ligated = Ligate(bands)
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
	for product in bands:
		print 'purified #'+str(n)+':'
		product.prettyPrint()
		n = n + 1
	print "\nLigation product(s):\n"
	n = 1
	for product in ligated:
		print 'ligation product #'+str(n)+':'
		product.prettyPrint()
		n = n + 1
